"""Independent CRYSTAL electronic DOSS.DAT and fort.25 density-of-states reader."""

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

from quantum_macaroni.core.constants import HTR_TO_EV
from quantum_macaroni.parsers.base import DOSResult

_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_DAT_HEADER = re.compile(r"^#\s*\S+\s+(\d+)\s+NPROJ\s+(\d+)\s+NSPIN\s+(\d+)\s*$", re.IGNORECASE)
_FERMI = re.compile(rf"\bEFERMI\s*(?:\([^)]*\)\s*)?(?:[=:]\s*)?({_NUMBER})", re.IGNORECASE)
_FORT_HEADER = re.compile(r"^-\%-([0-9])DOSS\s+(\d+)\s+(\d+)")
_FIELD_WIDTH = 12
_FORT_METADATA_ROWS = 3
_FORT_STEP_START = 30
_FORT_FERMI_START = 42
_SPIN_CHANNELS = {1, 2}


def _float(token: str) -> float:
    """Read one finite number, including Fortran D exponents."""
    token = token.strip()
    if re.fullmatch(_NUMBER, token) is None:
        raise ValueError(f"Invalid numeric value {token!r} in CRYSTAL DOS")
    value = float(token.replace("D", "E").replace("d", "e"))
    if not np.isfinite(value):
        raise ValueError("Non-finite numeric value in CRYSTAL DOS")
    return value


def _result(energies: npt.NDArray[np.float64], dos: npt.NDArray[np.float64], fermi_energy: float) -> DOSResult:
    """Validate the energy grid and normalize Hartree units and beta signs."""
    if len(energies) == 0 or not np.isfinite(energies).all() or np.any(np.diff(energies) <= 0):
        raise ValueError("CRYSTAL DOS energy grid must be finite and strictly increasing")
    normalized = dos.copy()
    if normalized.shape[0] == 2:  # noqa: PLR2004
        normalized[1] *= -1
    with np.errstate(over="ignore"):
        result = DOSResult(
            energies=np.ascontiguousarray(energies * HTR_TO_EV),
            dos=np.ascontiguousarray(normalized / HTR_TO_EV),
            fermi_energy=fermi_energy * HTR_TO_EV,
        )
    if (
        not np.isfinite(result.energies).all()
        or not np.isfinite(result.dos).all()
        or not np.isfinite(result.fermi_energy)
    ):
        raise ValueError("Non-finite converted values in CRYSTAL DOS")
    return result


def _parse_dat(lines: list[str]) -> DOSResult:
    """Read named dimensions, numeric spin tables, and the electronic Fermi footer."""
    header = _DAT_HEADER.fullmatch(lines[0].strip())
    if header is None:
        raise ValueError("Invalid CRYSTAL DOSS.DAT header; expected NPROJ and NSPIN")
    nenergy, nprojections, jspins = (int(value) for value in header.groups())
    if nenergy <= 0 or nprojections <= 0 or jspins not in _SPIN_CHANNELS:
        raise ValueError("Invalid CRYSTAL DOS dimensions or spin count")
    fermi_matches = _FERMI.findall("\n".join(lines))
    if not fermi_matches:
        raise ValueError("Missing EFERMI metadata in electronic DOSS.DAT; phonon DOS is unsupported")
    fermi_energy = _float(fermi_matches[0])
    if any(not np.isclose(_float(value), fermi_energy, rtol=1e-8, atol=1e-10) for value in fermi_matches[1:]):
        raise ValueError("Inconsistent EFERMI metadata in CRYSTAL DOSS.DAT")
    rows = []
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")) or stripped == "&":
            continue
        tokens = stripped.split()
        if len(tokens) != nprojections + 1:
            raise ValueError("Wrong projection column count in CRYSTAL DOSS.DAT")
        rows.append([_float(token) for token in tokens])
    if len(rows) != jspins * nenergy:
        raise ValueError("Incomplete or extra spin-channel rows in CRYSTAL DOSS.DAT")
    tables = np.array(rows, dtype=np.float64).reshape(jspins, nenergy, nprojections + 1)
    energies = tables[0, :, 0]
    if any(not np.allclose(table[:, 0], energies, rtol=1e-8, atol=1e-10) for table in tables[1:]):
        raise ValueError("Alpha and beta energy grids differ in CRYSTAL DOSS.DAT")
    # DOSS.DAT already stores E - E_F in Hartree.
    return _result(energies, tables[:, :, 1:], fermi_energy)


@dataclass(slots=True)
class _DOSBlock:
    """One fort.25 projection in the file's original Hartree units."""

    energies: npt.NDArray[np.float64]
    values: npt.NDArray[np.float64]
    fermi_energy: float
    jspins: int


def _parse_block(lines: list[str]) -> _DOSBlock:
    """Read fixed-width metadata and values for one DOSS projection."""
    header = _FORT_HEADER.match(lines[0])
    if header is None or len(lines) < _FORT_METADATA_ROWS:
        raise ValueError("Incomplete CRYSTAL fort.25 DOSS header")
    jspins = int(header.group(1)) % 2 + 1
    nenergy = int(header.group(3))
    if nenergy <= 0:
        raise ValueError("Invalid point count in CRYSTAL fort.25 DOSS")
    step = _float(lines[0][_FORT_STEP_START : _FORT_STEP_START + _FIELD_WIDTH])
    fermi_energy = _float(lines[0][_FORT_FERMI_START : _FORT_FERMI_START + _FIELD_WIDTH])
    minimum = _float(lines[1][_FIELD_WIDTH : 2 * _FIELD_WIDTH])
    if step <= 0:
        raise ValueError("CRYSTAL fort.25 DOSS energy step must be positive")
    values = []
    for line in lines[_FORT_METADATA_ROWS:]:
        if not line.strip():
            continue
        record = line.rstrip()
        for start in range(0, len(record), _FIELD_WIDTH):
            values.append(_float(record[start : start + _FIELD_WIDTH]))
    if len(values) != nenergy:
        raise ValueError("Incomplete or extra values in CRYSTAL fort.25 DOSS projection")
    # Unlike DOSS.DAT, fort.25 stores an absolute starting energy.
    energies = minimum + step * np.arange(nenergy, dtype=np.float64) - fermi_energy
    return _DOSBlock(energies, np.array(values, dtype=np.float64), fermi_energy, jspins)


def _parse_fort25(lines: list[str]) -> DOSResult:
    """Gather DOS projections, skipping other properties in a combined fort.25 file."""
    headers = [index for index, line in enumerate(lines) if line.startswith("-%-")]
    blocks = []
    for start, end in zip(headers, [*headers[1:], len(lines)], strict=True):
        if "DOSS" in lines[start].split()[0]:
            blocks.append(_parse_block(lines[start:end]))
    if not blocks:
        raise ValueError("No electronic DOSS blocks found in CRYSTAL fort.25")
    first = blocks[0]
    for block in blocks[1:]:
        if block.jspins != first.jspins:
            raise ValueError("Inconsistent spin flags in CRYSTAL fort.25 DOSS")
        if not np.isclose(block.fermi_energy, first.fermi_energy, rtol=1e-8, atol=1e-10):
            raise ValueError("Inconsistent Fermi energies in CRYSTAL fort.25 DOSS")
        if block.energies.shape != first.energies.shape or not np.allclose(
            block.energies, first.energies, rtol=1e-8, atol=1e-10
        ):
            raise ValueError("Inconsistent projection energy grids in CRYSTAL fort.25 DOSS")
    if len(blocks) % first.jspins:
        raise ValueError("Unpaired alpha/beta projections in CRYSTAL fort.25 DOSS")
    nprojections = len(blocks) // first.jspins
    values = np.stack([block.values for block in blocks]).reshape(first.jspins, nprojections, len(first.energies))
    return _result(first.energies, values.transpose(0, 2, 1), first.fermi_energy)


class CrystalDOSSParser:
    """Read CRYSTAL electronic DOS directly from formatted output files."""

    name = "crystal-doss"

    def parse(self, filepath: str | Path) -> DOSResult:  # noqa: PLR6301
        """Parse DOSS.DAT or electronic DOSS blocks in a fort.25 file.

        Args:
            filepath: Path to a CRYSTAL DOSS.DAT, .DOSS, or fort.25 file.
                The format is detected from its contents, including combined
                fort.25 files containing band data before the DOS blocks.

        Returns:
            DOS with one restricted channel or separate alpha/beta channels.
            Energies are in eV relative to the Fermi level, densities in states/eV/cell,
            and the absolute Fermi energy is retained in eV. Projections keep
            their file order; no extra spin or orbital multiplicity is applied.

        Raises:
            ValueError: If the format, dimensions, values, or spin grids are invalid.
                Phonon DOS files and missing electronic Fermi metadata are rejected.

        """
        lines = Path(filepath).read_text().splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        if not lines:
            raise ValueError("Empty CRYSTAL DOS file")
        if lines[0].startswith("-%-"):
            return _parse_fort25(lines)
        if lines[0].strip().startswith("#"):
            return _parse_dat(lines)
        raise ValueError("Unrecognized CRYSTAL DOS format; expected DOSS.DAT or fort.25")
