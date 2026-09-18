"""Regression tests for CRYSTAL properties output and transport integration."""

from pathlib import Path

import numpy as np
import pytest

from quantum_macaroni import CrystalOutpParser, available_parsers, calculate_spin_polarized_transport, get_parser
from quantum_macaroni.core.constants import HTR_TO_EV

_HEADER = """DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)
 3.0 0.0 0.0
 0.0 4.0 0.0
 0.0 0.0 5.0
 N. OF SCF CYCLES 18 FERMI ENERGY -0.15D+00
 NUMBER OF K POINTS IN THE IBZ 2
 ***** 2 SYMMOPS - TRANSLATORS IN ANGSTROM
 ***** MATRICES AND TRANSLATORS IN THE CARTESIAN REFERENCE FRAME
 NO. 1 INVERSE 1 ORDER 1       NO. 2 INVERSE 2 ORDER 2
 1.000 0.000 0.000 0.000      -1.000 0.000 0.000 0.000
 0.000 1.000 0.000 0.000       0.000 -1.000 0.000 0.000
 0.000 0.000 1.000 0.000       0.000 0.000 -1.000 0.000
 TTTTTTTTTTTTTTTTTTTT SYMMOPS TELAPSE 0.44 TCPU 0.18
 *** K POINTS COORDINATES (OBLIQUE COORDINATES IN UNITS OF IS = 4)
 1-R( 0 0 0) 2-C( -1 1 0)

"""
_ALPHA = """ EIGENVALUES - K= 1 ( 0 0 0)
 -2.0D-01(Ag ) -1.0E-01(B1u) +2.0E-01(Ag )

 EIGENVALUES - K= 2 ( -1 1 0)
 -1.8E-01 -0.8E-01 2.2E-01

"""
_BETA = """ EIGENVALUES - K= 2 ( -1 1 0)
 -1.7E-01 -0.7E-01 2.3E-01

 EIGENVALUES - K= 1 ( 0 0 0)
 -1.9E-01(Ag ) -0.9E-01(B1u) 2.1E-01(Ag )

"""


def _write_outp(tmp_path: Path, content: str) -> Path:
    """Write a small CRYSTAL output fixture."""
    path = tmp_path / "test.outp"
    path.write_text(content)
    return path


def test_afm_example() -> None:
    """Read the supplied AFM example without collapsing its spin channels."""
    path = Path(__file__).resolve().parents[1] / "examples" / "outp" / "mno2afm.outp"
    parsed = CrystalOutpParser().parse(path)
    assert parsed.jspins == 2  # noqa: PLR2004
    assert parsed.eigenvalues.shape == (2, 301, 152)
    assert (parsed.nk, parsed.nbands) == (301, 152)
    assert parsed.symops.shape == (8, 3, 3)
    np.testing.assert_allclose(parsed.kpoints[[0, 1, -1]], [[0, 0, 0], [1 / 12, 0, 0], [0.5, 0.5, 0.5]])
    np.testing.assert_allclose(parsed.lattice, np.diag([4.39371, 4.39371, 2.88061]))
    np.testing.assert_allclose(parsed.eigenvalues[:, 0, 0], -235.8 * HTR_TO_EV)
    assert parsed.fermi_energy == pytest.approx(-HTR_TO_EV)
    assert not np.array_equal(parsed.eigenvalues[0], parsed.eigenvalues[1])


@pytest.mark.parametrize("polarized", [False, True])
def test_channels_and_units(tmp_path: Path, polarized: bool) -> None:
    """Convert coordinates and energies and align out-of-order beta points."""
    blocks = "ALPHA ELECTRONS\n" + _ALPHA + "BETA ELECTRONS\n" + _BETA if polarized else _ALPHA
    path = _write_outp(tmp_path, _HEADER + blocks + "FERMI ENERGY AND DENSITY MATRIX CALCULATION\n 9.0E+00\n")
    parsed = get_parser("crystal-outp").parse(path)
    assert "crystal-outp" in available_parsers()
    assert parsed.jspins == (2 if polarized else 1)
    assert parsed.eigenvalues.shape == (parsed.jspins, 2, 3)
    np.testing.assert_allclose(parsed.kpoints, [[0, 0, 0], [-0.25, 0.25, 0]])
    np.testing.assert_allclose(parsed.eigenvalues[0] / HTR_TO_EV, [[-0.2, -0.1, 0.2], [-0.18, -0.08, 0.22]])
    assert parsed.fermi_energy == pytest.approx(-0.15 * HTR_TO_EV)
    np.testing.assert_array_equal(parsed.symops, [np.eye(3, dtype=int), -np.eye(3, dtype=int)])
    if polarized:
        np.testing.assert_allclose(parsed.eigenvalues[1] / HTR_TO_EV, [[-0.19, -0.09, 0.21], [-0.17, -0.07, 0.23]])


def test_cartesian_symmetry_conversion(tmp_path: Path) -> None:
    """Convert rounded Cartesian rotations correctly for a nonorthogonal cell."""
    header = _HEADER.replace("3.0 0.0 0.0", "1.0 0.0 0.0").replace("0.0 4.0 0.0", "-0.5 0.8660254 0.0")
    header = header.replace("-1.000 0.000 0.000 0.000", "-0.500 -0.866 0.000 0.000")
    header = header.replace("0.000 -1.000 0.000 0.000", "0.866 -0.500 0.000 0.000")
    header = header.replace("0.000 0.000 -1.000 0.000", "0.000 0.000 1.000 0.000")
    parsed = CrystalOutpParser().parse(_write_outp(tmp_path, header + _ALPHA))
    np.testing.assert_array_equal(parsed.symops[1], [[0, -1, 0], [1, -1, 0], [0, 0, 1]])


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (_HEADER + "ALPHA ELECTRONS\n" + _ALPHA, "Missing or mixed spin channels"),
        (_HEADER + "SPIN POLARIZED SYSTEM\n" + _ALPHA, "requires ALPHA and BETA"),
        (_HEADER + _ALPHA.replace("-1.8E-01 -0.8E-01 2.2E-01", "-1.8E-01"), "inconsistent band counts"),
        (_HEADER + _ALPHA.replace("K= 2 ( -1 1 0)", "K= 2 ( 1 1 0)"), "coordinates do not match"),
        (_HEADER + _ALPHA + _ALPHA, "Duplicate eigenvalue block"),
        (_HEADER + _ALPHA.split(" EIGENVALUES - K= 2", maxsplit=1)[0], "Incomplete k-point eigenvalues"),
        (_HEADER.replace("IS = 4", "IS = 0") + _ALPHA, "coordinate divisor"),
        (_HEADER.replace("1-R( 0 0 0) 2-C( -1 1 0)", "1-R( 0 0 0)") + _ALPHA, "Incomplete k-point coordinate"),
        (_HEADER.replace("FERMI ENERGY", "REFERENCE ENERGY") + _ALPHA, "Missing numeric FERMI"),
        (_HEADER.replace("SYMMOPS - TRANSLATORS", "OPERATORS - TRANSLATORS") + _ALPHA, "enable SYMMOPS"),
        (_HEADER.replace("***** 2 SYMMOPS", "***** 3 SYMMOPS") + _ALPHA, "Incomplete SYMMOPS"),
        (_HEADER.replace("DIRECT LATTICE VECTOR COMPONENTS", "CELL") + _ALPHA, "enable COORPRT"),
    ],
)
def test_invalid_output(tmp_path: Path, content: str, message: str) -> None:
    """Reject malformed or partial output before passing it to interpolation."""
    with pytest.raises(ValueError, match=message):
        CrystalOutpParser().parse(_write_outp(tmp_path, content))


def test_iteration_selector(tmp_path: Path) -> None:
    """Accept the single dataset selector and reject unsupported history indices."""
    path = _write_outp(tmp_path, _HEADER + _ALPHA)
    parsed = CrystalOutpParser().parse(path, iteration="1")
    assert parsed.jspins == 1
    with pytest.raises(ValueError, match="single dataset"):
        CrystalOutpParser().parse(path, iteration="2")


@pytest.mark.parametrize("polarized", [False, True])
def test_transport_workflow(tmp_path: Path, polarized: bool) -> None:
    """Run real interpolation and transport with both CRYSTAL spin layouts."""
    blocks = "ALPHA ELECTRONS\n" + _ALPHA + "BETA ELECTRONS\n" + _BETA if polarized else _ALPHA
    path = _write_outp(tmp_path, _HEADER + blocks)
    result = calculate_spin_polarized_transport(path, parser="crystal-outp", kpoint_mesh=(3, 3, 3), lr_ratio=2)
    assert result["parser"] == "crystal-outp"
    assert result["jspins"] == (2 if polarized else 1)
    assert np.isfinite(result["sigma"]).all()
    assert np.isfinite(result["seebeck"]).all()
    assert np.isfinite(result["kappa"]).all()
