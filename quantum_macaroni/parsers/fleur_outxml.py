"""FLEUR XML parser and structure extraction helpers."""

import importlib
from pathlib import Path
from typing import TypedDict

import numpy as np
import numpy.typing as npt
from lxml import etree  # ty:ignore[unresolved-import]
from ase.data import chemical_symbols as ase_symbols

from quantum_macaroni.core.constants import BOHR_TO_ANG, HTR_TO_EV
from quantum_macaroni.parsers.base import ParserResult


class FleurRawData(TypedDict):
    """Raw electronic data parsed from XML iteration."""

    kpoints: npt.NDArray[np.float64]
    eigenvalues: npt.NDArray[np.float64]
    fermi_energy: float
    jspins: int
    nbands: int
    nk: int


def _parse_coord(token: str) -> float:
    """Parse numeric token that may be a fraction.

    Args:
        token: Number token in decimal or ``num/den`` form.

    Returns:
        Parsed floating-point value.

    """
    if "/" in token:
        num, den = token.split("/")
        return float(num) / float(den)
    return float(token)


def _sym2z() -> dict[str, int]:
    """Return map from chemical symbol to atomic number."""
    return {value: idx + 1 for idx, value in enumerate(ase_symbols[1:])}


def structure_from_outxml(
    filepath: str | Path,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.int_]]:
    """Extract structure information from FLEUR XML.

    Args:
        filepath: Path to FLEUR ``out.xml`` file.

    Returns:
        Tuple with lattice matrix, fractional positions, and atomic numbers.

    """
    try:
        # Prefer masci-tools when available because it already knows FLEUR schema quirks,
        # so we avoid re-encoding fragile XML semantics in our own parser.
        fleur_xml = importlib.import_module("masci_tools.io.fleur_xml")
        xml_getters = importlib.import_module("masci_tools.util.xml.xml_getters")
        load_outxml = fleur_xml.load_outxml
        get_structure_data = xml_getters.get_structure_data

        xmltree, schema = load_outxml(str(filepath))
        atoms, cell, _ = get_structure_data(xmltree, schema)
        cell = np.array(cell)
        frac_pos = np.array([atom.position for atom in atoms]) @ np.linalg.inv(cell)
        mapping = _sym2z()
        return cell, frac_pos, np.array([mapping.get(atom.symbol, 0) for atom in atoms])
    except ImportError:
        # Keep a pure-lxml fallback so the package remains usable in lightweight setups
        # where the full FLEUR tooling stack is intentionally not installed.
        pass

    tree = etree.parse(str(filepath))
    root = tree.getroot()
    bravais = root.find(".//bravaisMatrix")
    if bravais is None:
        raise ValueError("No <bravaisMatrix> found")

    rows = [[float(x) for x in bravais.find(name).text.split()] for name in ["row-1", "row-2", "row-3"]]
    lattice = np.array(rows, dtype=np.float64) * BOHR_TO_ANG
    inv_lat = np.linalg.inv(lattice)

    species_to_z = {}
    for species in root.findall(".//species"):
        atomic_number = species.get("atomicNumber", "")
        if atomic_number:
            species_to_z[species.get("name", "")] = int(atomic_number)

    frac_pos = []
    numbers = []
    for atom_group in root.findall(".//atomGroup"):
        atomic_number = species_to_z.get(atom_group.get("species", ""), 0)
        for rel_pos in atom_group.findall("relPos"):
            frac_pos.append([_parse_coord(token) for token in rel_pos.text.split()])
            numbers.append(atomic_number)
        for abs_pos in atom_group.findall("absPos"):
            coords = np.array([_parse_coord(token) for token in abs_pos.text.split()], dtype=np.float64) * BOHR_TO_ANG
            frac_pos.append((coords @ inv_lat).tolist())
            numbers.append(atomic_number)
        for film_pos in atom_group.findall("filmPos"):
            frac_pos.append([_parse_coord(token) for token in film_pos.text.split()])
            numbers.append(atomic_number)

    return lattice, np.array(frac_pos, dtype=np.float64), np.array(numbers, dtype=int)


def parse_fleur_outxml(filepath: str | Path, iteration: str = "last") -> FleurRawData:
    """Parse electronic data for one iteration from FLEUR XML.

    Args:
        filepath: Path to FLEUR ``out.xml`` file.
        iteration: Iteration index as string or ``"last"``.

    Returns:
        Dictionary with k-points, eigenvalues, and metadata.

    """
    tree = etree.parse(str(filepath))
    root = tree.getroot()

    magnetism = root.find(".//magnetism")
    jspins = int(magnetism.get("jspins", "1"))

    iterations = root.findall(".//iteration")
    # FLEUR keeps the self-consistent history in the same file; selecting by iteration lets
    # callers inspect intermediate states without needing separate postprocessing files.
    selected_iteration = iterations[-1] if iteration == "last" else iterations[int(iteration) - 1]

    fermi_element = selected_iteration.find(".//FermiEnergy")
    fermi_htr = float(fermi_element.get("value"))
    fermi_energy = fermi_htr * HTR_TO_EV

    eigenvalue_nodes = selected_iteration.findall(".//eigenvaluesAt")
    data_by_spin = {}
    for node in eigenvalue_nodes:
        spin = int(node.get("spin"))
        ikpt = int(node.get("ikpt"))
        kx = float(node.get("k_x"))
        ky = float(node.get("k_y"))
        kz = float(node.get("k_z"))
        values = np.array([float(x) * HTR_TO_EV for x in node.text.split()], dtype=np.float64)

        if spin not in data_by_spin:
            data_by_spin[spin] = {}
        data_by_spin[spin][ikpt] = {"kpoint": np.array([kx, ky, kz], dtype=np.float64), "eigenvalues": values}

    nspin = len(data_by_spin)
    nk = len(data_by_spin[1])
    nbands = len(data_by_spin[1][1]["eigenvalues"])

    kpoints = np.empty((nk, 3), dtype=np.float64)
    eigenvalues = np.empty((nspin, nk, nbands), dtype=np.float64)

    sorted_k = sorted(data_by_spin[1].keys())
    for idx, ikpt in enumerate(sorted_k):
        kpoints[idx] = data_by_spin[1][ikpt]["kpoint"]

    for spin in range(nspin):
        spin_key = spin + 1
        for idx, ikpt in enumerate(sorted(data_by_spin[spin_key].keys())):
            eigenvalues[spin, idx] = data_by_spin[spin_key][ikpt]["eigenvalues"]

    return {
        "kpoints": kpoints,
        "eigenvalues": eigenvalues,
        "fermi_energy": fermi_energy,
        "jspins": jspins,
        "nbands": nbands,
        "nk": nk,
    }


def read_symops_from_outxml(filepath: str | Path) -> npt.NDArray[np.int_]:
    """Read integer symmetry operations from FLEUR XML.

    Args:
        filepath: Path to FLEUR ``out.xml`` file.

    Returns:
        Symmetry operations with shape ``(nsym, 3, 3)``.

    """
    tree = etree.parse(str(filepath))
    root = tree.getroot()
    symops = []
    for symop in root.findall(".//symOp"):
        rows = []
        for name in ["row-1", "row-2", "row-3"]:
            row = symop.find(name)
            parts = row.text.split()
            rows.append([int(float(x)) for x in parts[:3]])
        symops.append(rows)
    return np.array(symops, dtype=int)


class FleurOutxmlParser:
    """Parser plugin for FLEUR ``out.xml`` files."""

    name = "fleur-outxml"

    def parse(self, filepath: str | Path, iteration: str = "last") -> ParserResult:  # noqa: PLR6301
        """Parse a FLEUR output file into :class:`ParserResult`.

        Args:
            filepath: Path to FLEUR ``out.xml`` file.
            iteration: Iteration index as string or ``"last"``.

        Returns:
            Normalized parser output used by calculators.

        """
        raw = parse_fleur_outxml(filepath, iteration=iteration)
        cell = structure_from_outxml(filepath)
        symops = read_symops_from_outxml(filepath)
        return ParserResult(
            kpoints=raw["kpoints"],
            eigenvalues=raw["eigenvalues"],
            fermi_energy=raw["fermi_energy"],
            jspins=raw["jspins"],
            nbands=raw["nbands"],
            nk=raw["nk"],
            lattice=cell[0],
            symops=symops,
        )
