"""CRYSTAL properties ``outp`` parser for one or two spin channels."""

import re
from pathlib import Path

import numpy as np
import numpy.typing as npt

from quantum_macaroni.core.constants import HTR_TO_EV
from quantum_macaroni.parsers.base import ParserResult

_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_COORDS = r"\(\s*([+-]?\d+)\s+([+-]?\d+)\s+([+-]?\d+)\s*\)"
_KPOINT = re.compile(r"(\d+)-[RC]\s*" + _COORDS)
_EIGENVALUES = re.compile(r"EIGENVALUES\s*-\s*K\s*=\s*(\d+)\s*" + _COORDS)
_SPIN = re.compile(r"^\s*(ALPHA|BETA)\s+ELECTRONS\s*$")
_MATRIX_SIZE = 3
_SYMMETRY_ROW_SIZE = 4
_MIN_LATTICE_VOLUME = 1e-12


def _numbers(line: str) -> list[float]:
    """Read a numeric row, allowing Fortran exponents and band symmetry labels."""
    cleaned = re.sub(r"\([^)]*\)", " ", line).strip()
    if not cleaned or re.fullmatch(rf"{_NUMBER}(?:\s+{_NUMBER})*", cleaned) is None:
        return []
    return [float(token.replace("D", "E").replace("d", "e")) for token in cleaned.split()]


def _lattice(lines: list[str]) -> npt.NDArray[np.float64]:
    """Read the direct lattice in angstrom, with vectors stored as rows."""
    for index, line in enumerate(lines):
        if "DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)" in line:
            rows = [_numbers(row) for row in lines[index + 1 : index + 1 + _MATRIX_SIZE]]
            if any(len(row) != _MATRIX_SIZE for row in rows):
                raise ValueError("Invalid direct lattice vectors in CRYSTAL outp")
            lattice = np.array(rows, dtype=np.float64)
            if not np.isfinite(lattice).all() or abs(np.linalg.det(lattice)) < _MIN_LATTICE_VOLUME:
                raise ValueError("Singular or non-finite direct lattice in CRYSTAL outp")
            return lattice
    raise ValueError("Missing DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM) in CRYSTAL outp; enable COORPRT")


def _fractional_rotation(cartesian: npt.NDArray[np.float64], lattice: npt.NDArray[np.float64]) -> npt.NDArray[np.int_]:
    """Convert a rounded Cartesian matrix to a valid fractional rotation."""
    fractional = np.linalg.solve(lattice.T, cartesian @ lattice.T)
    integer = np.rint(fractional).astype(int)
    # CRYSTAL prints Cartesian rotations to only three decimal places.
    if not np.allclose(lattice.T @ integer, cartesian @ lattice.T, atol=2e-3, rtol=2e-3):
        raise ValueError("Symmetry matrix is incompatible with the direct lattice in CRYSTAL outp")
    if not np.isclose(abs(np.linalg.det(integer)), 1.0):
        raise ValueError("Invalid symmetry rotation in CRYSTAL outp")
    return integer


def _symops(lines: list[str], lattice: npt.NDArray[np.float64]) -> npt.NDArray[np.int_]:
    """Convert printed Cartesian rotations to integer fractional rotations."""
    header = next((i for i, line in enumerate(lines) if "SYMMOPS - TRANSLATORS IN ANGSTROM" in line), None)
    if header is None:
        raise ValueError("Missing symmetry matrices in CRYSTAL outp; enable SYMMOPS")
    count = re.search(r"(\d+)\s+SYMMOPS", lines[header])
    if count is None:
        raise ValueError("Invalid SYMMOPS header in CRYSTAL outp")
    rotations = []
    for index in range(header + 1, len(lines)):
        line = lines[index]
        if "TTTT" in line and "SYMMOPS" in line:
            break
        operators = re.findall(r"NO\.\s+(\d+)\s+INVERSE", line)
        if not operators:
            continue
        rows = [_numbers(row) for row in lines[index + 1 : index + 1 + _MATRIX_SIZE]]
        if any(len(row) != len(operators) * _SYMMETRY_ROW_SIZE for row in rows):
            raise ValueError("Incomplete symmetry matrix in CRYSTAL outp")
        for column in range(len(operators)):
            start = column * _SYMMETRY_ROW_SIZE
            cartesian = np.array([row[start : start + _MATRIX_SIZE] for row in rows])
            rotations.append(_fractional_rotation(cartesian, lattice))
    if len(rotations) != int(count.group(1)):
        raise ValueError("Incomplete SYMMOPS section in CRYSTAL outp")
    return np.ascontiguousarray(rotations, dtype=int)


def _kpoints(lines: list[str]) -> tuple[int, dict[int, tuple[int, int, int]]]:
    """Read the integer k-point table and its fractional coordinate divisor."""
    headers = [i for i, line in enumerate(lines) if "K POINTS COORDINATES" in line]
    if len(headers) != 1:
        raise ValueError("Expected one K POINTS COORDINATES table in CRYSTAL outp")
    header = headers[0]
    factor_match = re.search(r"IS\s*=\s*(\d+)", lines[header])
    if factor_match is None or int(factor_match.group(1)) <= 0:
        raise ValueError("Missing or invalid k-point coordinate divisor IS in CRYSTAL outp")
    points: dict[int, tuple[int, int, int]] = {}
    for line in lines[header + 1 :]:
        matches = list(_KPOINT.finditer(line))
        if not matches:
            if points:
                break
            continue
        for match in matches:
            point_id = int(match.group(1))
            if point_id in points:
                raise ValueError(f"Duplicate k-point {point_id} in CRYSTAL outp coordinate table")
            points[point_id] = (int(match.group(2)), int(match.group(3)), int(match.group(4)))
    counts = re.findall(r"POINTS IN THE IBZ\s+(\d+)", "\n".join(lines[:header]))
    if not points or (counts and len(points) != int(counts[-1])):
        raise ValueError("Incomplete k-point coordinate table in CRYSTAL outp")
    return int(factor_match.group(1)), points


def _band_blocks(lines: list[str], points: dict[int, tuple[int, int, int]]) -> dict[str, dict[int, list[float]]]:
    """Read eigenvalue blocks keyed by spin label and k-point ID."""
    channels: dict[str, dict[int, list[float]]] = {}
    spin = "restricted"
    values: list[float] | None = None
    for line in lines:
        spin_match = _SPIN.match(line)
        point_match = _EIGENVALUES.search(line)
        if spin_match:
            spin = spin_match.group(1)
            values = None
            channels.setdefault(spin, {})
        elif point_match:
            point_id = int(point_match.group(1))
            coords = (int(point_match.group(2)), int(point_match.group(3)), int(point_match.group(4)))
            if point_id not in points or coords != points[point_id]:
                raise ValueError(f"Eigenvalue coordinates do not match k-point {point_id} in CRYSTAL outp")
            channel = channels.setdefault(spin, {})
            if point_id in channel:
                raise ValueError(f"Duplicate eigenvalue block for {spin} k-point {point_id} in CRYSTAL outp")
            values = []
            channel[point_id] = values
        elif values is not None:
            row = _numbers(line)
            if row:
                values.extend(row)
            else:
                values = None
    return channels


def _bands(lines: list[str], points: dict[int, tuple[int, int, int]]) -> npt.NDArray[np.float64]:
    """Align each spin channel by k-point ID, validating completeness."""
    channels = _band_blocks(lines, points)
    expected = ["ALPHA", "BETA"] if "ALPHA" in channels or "BETA" in channels else ["restricted"]
    if set(channels) != set(expected):
        raise ValueError("Missing or mixed spin channels in CRYSTAL outp")
    if expected == ["restricted"] and any("SPIN POLARIZED" in line for line in lines):
        raise ValueError("Spin-polarized CRYSTAL outp requires ALPHA and BETA eigenvalue blocks")
    for channel in channels.values():
        if set(channel) != set(points):
            raise ValueError("Incomplete k-point eigenvalues in CRYSTAL outp spin channel")
    sizes = {len(values) for channel in channels.values() for values in channel.values()}
    if len(sizes) != 1 or 0 in sizes:
        raise ValueError("Empty or inconsistent band counts in CRYSTAL outp")
    energies = np.array([[channels[spin][key] for key in sorted(points)] for spin in expected], dtype=np.float64)
    if not np.isfinite(energies).all():
        raise ValueError("Non-finite eigenvalues in CRYSTAL outp")
    return np.ascontiguousarray(energies * HTR_TO_EV)


def _fermi_energy(lines: list[str]) -> float:
    """Read the last numeric Fermi energy in atomic units and convert to eV."""
    matches = re.findall(rf"FERMI ENERGY\s*(?:[=:]\s*)?({_NUMBER})", "\n".join(lines), flags=re.IGNORECASE)
    if not matches:
        raise ValueError("Missing numeric FERMI ENERGY in CRYSTAL outp")
    energy = float(matches[-1].replace("D", "E").replace("d", "e")) * HTR_TO_EV
    if not np.isfinite(energy):
        raise ValueError("Non-finite Fermi energy in CRYSTAL outp")
    return energy


class CrystalOutpParser:
    """Parser plugin for CRYSTAL properties output, including AFM systems."""

    name = "crystal-outp"

    def parse(self, filepath: str | Path, iteration: str = "last") -> ParserResult:  # noqa: PLR6301
        """Parse one CRYSTAL properties dataset into normalized transport input.

        Args:
            filepath: Path to a CRYSTAL ``outp`` file with COORPRT and SYMMOPS output.
            iteration: ``"last"`` or ``"1"``; properties output contains one dataset.

        Returns:
            Fractional k-points, energies in eV, lattice in angstrom, and symmetries.
            AFM alpha and beta blocks remain separate spin channels. Restricted
            output has one channel; transport supplies its spin degeneracy.

        Raises:
            ValueError: If the dataset is incomplete, inconsistent, or unsupported.

        """
        if iteration not in {"last", "1"}:
            raise ValueError("CRYSTAL outp supports a single dataset: iteration must be 'last' or '1'")
        lines = Path(filepath).read_text().splitlines()
        lattice = _lattice(lines)
        factor, points = _kpoints(lines)
        energies = _bands(lines, points)
        nspin, nk, nbands = energies.shape
        return ParserResult(
            kpoints=np.ascontiguousarray([points[key] for key in sorted(points)], dtype=np.float64) / factor,
            eigenvalues=energies,
            fermi_energy=_fermi_energy(lines),
            jspins=nspin,
            nbands=nbands,
            nk=nk,
            lattice=lattice,
            symops=_symops(lines, lattice),
        )
