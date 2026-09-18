"""CRYSTAL electronic DOS format, spin normalization, and validation regressions."""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from quantum_macaroni import CrystalDOSSParser, DOSResult, available_parsers
from quantum_macaroni.core.constants import HTR_TO_EV

_ALPHA_ROWS = """ -1.0000D-01 1.0000D+00 1.0000D+01
 0.0000E+00 2.0000E+00 2.0000E+01
 1.0000E-01 3.0000E+00 3.0000E+01
"""
_BETA_ROWS = """ -1.0000E-01 -4.0000E+00 -4.0000E+01
 0.0000E+00 -5.0000E+00 -5.0000E+01
 1.0000E-01 -6.0000E+00 -6.0000E+01
"""
_FOOTER = "# EFERMI (HARTREE) -0.20000D+00\n"


def _dat(jspins: int) -> str:
    """Create DOSS.DAT tables using CRYSTAL's NEPTS/NPROJ/NSPIN header."""
    content = f"# NEPTS 3 NPROJ 2 NSPIN {jspins}\n#\n"
    content += '@ XAXIS LABEL "E-EFERMI (HARTREE)"\n'
    content += '@ YAXIS LABEL "DENSITY OF STATES (STATES/HARTREE/CELL)"\n'
    content += _ALPHA_ROWS
    if jspins == 2:  # noqa: PLR2004
        content += "#\n# BETA\n#\n" + _BETA_ROWS
    return content + _FOOTER


def _block(
    flag: int,
    values: list[float],
    minimum: float = -0.3,
    step: float = 0.1,
    fermi: float = -0.2,
) -> str:
    """Write a fort.25 projection with adjoining 12-character numeric fields."""
    content = f"-%-{flag}DOSS{1:5d}{len(values):5d}{0.0:12.5E}{step:12.5E}{fermi:12.5E}\n"
    content += f"{0.0:12.5E}{minimum:12.5E}\n"
    content += "    1   36    0    0    0    0\n"
    # A negative field follows a positive one without whitespace.
    content += "".join(f"{value:12.5E}" for value in values[:2]) + "\n"
    content += "".join(f"{value:12.5E}" for value in values[2:]) + "\n"
    return content


def _fort25(flag: int) -> str:
    """Create alpha projections followed by beta projections when spin-polarized."""
    content = _block(flag, [1.0, 2.0, 3.0]) + _block(flag, [10.0, 20.0, 30.0])
    if flag % 2:
        content += _block(flag, [-4.0, -5.0, -6.0]) + _block(flag, [-40.0, -50.0, -60.0])
    return content


def _parse(tmp_path: Path, content: str) -> DOSResult:
    """Write a fixture with a neutral filename to verify content-based detection."""
    filepath = tmp_path / "dos.txt"
    filepath.write_text(content)
    return CrystalDOSSParser().parse(filepath)


@pytest.mark.parametrize("jspins", [1, 2])
def test_dat_spins_and_units(tmp_path: Path, jspins: int) -> None:
    """Normalize each spin's DOS without shifting the already relative energy grid."""
    parsed = _parse(tmp_path, _dat(jspins))
    assert isinstance(parsed, DOSResult)
    assert (parsed.jspins, parsed.nenergy, parsed.nprojections) == (jspins, 3, 2)
    assert parsed.dos.shape == (jspins, 3, 2)
    np.testing.assert_allclose(parsed.energies, np.array([-0.1, 0.0, 0.1]) * HTR_TO_EV)
    np.testing.assert_allclose(parsed.absolute_energies, np.array([-0.3, -0.2, -0.1]) * HTR_TO_EV)
    assert parsed.fermi_energy == pytest.approx(-0.2 * HTR_TO_EV)
    np.testing.assert_allclose(parsed.dos[0] * HTR_TO_EV, [[1, 10], [2, 20], [3, 30]])
    # Converting DOS inversely to energy preserves the integrated state count.
    np.testing.assert_allclose(trapezoid(parsed.dos[0], parsed.energies, axis=0), [0.4, 4.0])
    if jspins == 2:  # noqa: PLR2004
        np.testing.assert_allclose(parsed.dos[1] * HTR_TO_EV, [[4, 40], [5, 50], [6, 60]])
    assert parsed.energies.flags.c_contiguous
    assert parsed.dos.flags.c_contiguous
    assert parsed.energies.dtype == np.float64
    assert parsed.dos.dtype == np.float64


@pytest.mark.parametrize("flag", [0, 1, 2, 3])
def test_fort25_spins_and_units(tmp_path: Path, flag: int) -> None:
    """Read flag parity and normalize absolute fort.25 grids to the DAT convention."""
    parsed = _parse(tmp_path, _fort25(flag))
    expected = _parse(tmp_path, _dat(flag % 2 + 1))
    assert parsed.jspins == flag % 2 + 1
    np.testing.assert_allclose(parsed.energies, expected.energies, atol=1e-14)
    np.testing.assert_allclose(parsed.dos, expected.dos)
    assert parsed.fermi_energy == expected.fermi_energy


def test_combined_fort25(tmp_path: Path) -> None:
    """Extract DOS when BAND records precede and follow the DOSS projections."""
    band = "-%-0BAND    2    3\nnot DOS metadata\nnot DOS values\n"
    parsed = _parse(tmp_path, band + _fort25(1) + band)
    assert parsed.dos.shape == (2, 3, 2)
    np.testing.assert_allclose(parsed.dos[:, 0, 0] * HTR_TO_EV, [1, 4])


def test_negative_projected_dos(tmp_path: Path) -> None:
    """Preserve signed projected densities instead of applying absolute values."""
    parsed = _parse(tmp_path, _dat(2).replace("2.0000E+01", "-2.0000E+01").replace("-5.0000E+01", "5.0000E+01"))
    np.testing.assert_allclose(parsed.dos[:, 1, 1] * HTR_TO_EV, [-20, -50])
    parsed_fixed = _parse(tmp_path, _block(0, [1.0, -2.0, 3.0]))
    np.testing.assert_allclose(parsed_fixed.dos[0, :, 0] * HTR_TO_EV, [1, -2, 3])


def test_dat_whitespace_and_grace_separators(tmp_path: Path) -> None:
    """Allow formatting lines without depending on hard-coded beta row offsets."""
    content = "\n\n" + _dat(2).replace("# BETA\n", "# BETA\n\n&\n@ legend off\n")
    parsed = _parse(tmp_path, content)
    assert parsed.dos.shape == (2, 3, 2)


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("", "Empty"),
        ("not DOS data", "Unrecognized"),
        (_dat(1).replace("NPROJ", "NBND"), "header"),
        (_dat(1).replace("NSPIN 1", "NSPIN 0"), "spin count"),
        (_dat(1).replace("NEPTS 3", "NEPTS 0"), "dimensions"),
        (_dat(1).replace("NPROJ 2", "NPROJ 0"), "dimensions"),
        (_dat(1).replace(_FOOTER, ""), "Missing EFERMI"),
        (_dat(1) + "# EFERMI = -0.3\n", "Inconsistent EFERMI"),
        (_dat(1).replace("3.0000E+01", ""), "column count"),
        (_dat(2).replace(_BETA_ROWS, ""), "spin-channel rows"),
        (_dat(1) + _ALPHA_ROWS, "spin-channel rows"),
        (_dat(2).replace("-1.0000E-01 -4", "-2.0000E-01 -4"), "grids differ"),
        (_dat(1).replace("0.0000E+00 2", "-1.0000E-01 2"), "strictly increasing"),
        (_dat(1).replace("3.0000E+01", "NaN"), "Invalid numeric"),
        (_dat(1).replace("3.0000E+01", "1.0E+999"), "Non-finite"),
        (_dat(1).replace("-0.20000D+00", "1.0E+308"), "Non-finite converted"),
        ("-%-0PDOS    1    3\n", "No electronic DOSS"),
        ("-%-0DOSS    1\n", "DOSS header"),
        (_fort25(1) + _block(1, [-1.0, -2.0, -3.0]), "Unpaired"),
        (_block(0, [1.0, 2.0, 3.0]) + _block(1, [1.0, 2.0, 3.0]), "spin flags"),
        (_fort25(0) + _block(0, [1.0, 2.0, 3.0], fermi=-0.3), "Fermi energies"),
        (_fort25(0) + _block(0, [1.0, 2.0, 3.0], step=0.2), "energy grids"),
        (_fort25(0) + _block(0, [1.0, 2.0]), "energy grids"),
        (_block(0, [1.0, 2.0, 3.0], step=0.0), "step must be positive"),
        (_block(0, [1.0, 2.0, 3.0]).replace(" 3.00000E+00", ""), "values"),
        (_block(0, [1.0, 2.0, 3.0]) + " 4.00000E+00\n", "values"),
        (_block(0, [1.0, 2.0, 3.0]) + "-%-0DOSS broken\n", "DOSS header"),
    ],
)
def test_invalid_dos(tmp_path: Path, content: str, message: str) -> None:
    """Reject missing spin data and inconsistent metadata instead of guessing."""
    with pytest.raises(ValueError, match=message):
        _parse(tmp_path, content)


def test_dos_is_separate_from_band_parsers() -> None:
    """Keep DOS files out of the k-point transport parser selector."""
    assert CrystalDOSSParser.name not in available_parsers()
