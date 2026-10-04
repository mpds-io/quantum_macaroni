"""Base parser protocol, data container, and parser registry."""

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
import numpy.typing as npt


@dataclass(slots=True)
class ParserResult:
    """Normalized electronic-structure payload returned by parsers.

    Attributes:
        kpoints: Fractional k-point coordinates with shape ``(nk, 3)``.
        eigenvalues: Band energies with shape ``(nspin, nk, nbands)`` in eV.
        fermi_energy: Fermi energy in eV.
        jspins: Number of spin channels.
        nbands: Number of bands.
        nk: Number of k-points.
        lattice: Real-space lattice matrix in angstrom.
        symops: Integer symmetry operations with shape ``(nsym, 3, 3)``.

    """

    kpoints: npt.NDArray[np.float64]
    eigenvalues: npt.NDArray[np.float64]
    fermi_energy: float
    jspins: int
    nbands: int
    nk: int
    lattice: npt.NDArray[np.float64]
    symops: npt.NDArray[np.int_]


@dataclass(slots=True)
class DOSResult:
    """Electronic density of states with a separate axis for each spin channel.

    Attributes:
        energies: Energy grid relative to the Fermi energy, shape ``(nenergy,)``, in eV.
        dos: DOS values with shape ``(jspins, nenergy, nprojections)`` in states/eV/cell.
            Projection order is preserved. CRYSTAL's plotted beta sign is reversed
            so alpha and beta use the same physical sign convention. Restricted
            values are retained without applying an additional spin factor.
        fermi_energy: Absolute Fermi energy in eV.

    """

    energies: npt.NDArray[np.float64]
    dos: npt.NDArray[np.float64]
    fermi_energy: float

    @property
    def jspins(self) -> int:
        """Return the number of spin channels (one or two)."""
        return self.dos.shape[0]

    @property
    def nenergy(self) -> int:
        """Return the number of energy sampling points."""
        return len(self.energies)

    @property
    def nprojections(self) -> int:
        """Return the number of DOS columns per spin channel."""
        return self.dos.shape[2]

    @property
    def absolute_energies(self) -> npt.NDArray[np.float64]:
        """Return the energy grid before alignment to the Fermi energy, in eV."""
        return self.energies + self.fermi_energy


class ElectronicStructureParser(Protocol):
    """Protocol implemented by parser plugins."""

    name: str

    def parse(self, filepath: str | Path, iteration: str = "last") -> ParserResult:
        """Parse an input file and return normalized data.

        Args:
            filepath: Path to parser input file.
            iteration: Iteration selector used by parser implementation.

        Returns:
            Normalized parser output.

        """
        ...


_PARSERS: dict[str, ElectronicStructureParser] = {}


def register_parser(parser: ElectronicStructureParser) -> None:
    """Register a parser instance by its ``name`` attribute.

    Args:
        parser: Parser instance to register.

    """
    _PARSERS[parser.name] = parser


def get_parser(name: str) -> ElectronicStructureParser:
    """Return a parser instance registered under ``name``.

    Args:
        name: Registered parser name.

    Returns:
        Parser instance.

    Raises:
        ValueError: If parser name is unknown.

    """
    try:
        return _PARSERS[name]
    except KeyError as exc:
        available = ", ".join(sorted(_PARSERS))
        raise ValueError(f"Unknown parser '{name}'. Available: {available}") from exc


def available_parsers() -> tuple[str, ...]:
    """Return sorted names of all registered parsers.

    Returns:
        Tuple with parser names.

    """
    return tuple(sorted(_PARSERS))
