"""Project command-line entry point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from quantum_macaroni import (
    CrystalDOSSParser,
    DOSResult,
    available_calculators,
    available_parsers,
    calculate_spin_polarized_transport,
    get_parser,
)

GRID_SPEC_LEN = 3
DEFAULT_CHECKPOINT_PATH = "transport_state.npz"
INPUT_DETECTION_BYTES = 16384


def _parse_scalar_or_grid(values: list[float], name: str) -> float | np.ndarray:
    """Parse either a scalar value or a [start, stop, npoints] specification."""
    if len(values) == 1:
        return float(values[0])
    if len(values) != GRID_SPEC_LEN:
        raise ValueError(f"{name} must have 1 or 3 numbers")

    start, stop, npoints_raw = values
    npoints = int(round(npoints_raw))
    if not np.isclose(npoints_raw, npoints):
        raise ValueError(f"{name} third value must be an integer number of points")
    if npoints <= 0:
        raise ValueError(f"{name} number of points must be positive")
    if npoints == 1:
        return float(start)
    return np.linspace(float(start), float(stop), npoints, dtype=np.float64)


def _build_parser() -> argparse.ArgumentParser:
    """Create CLI parser for transport runs."""
    parser = argparse.ArgumentParser(description="Boltzmann transport calculator and CRYSTAL DOS reader")
    parser.add_argument("filepath", help="Path to band input (out.xml/outp) or CRYSTAL DOSS file")
    parser.add_argument("dos_filepath", nargs="?", help="Optional CRYSTAL DOSS file accompanying an outp file")
    parser.add_argument("--dos-file", "--doss", help="CRYSTAL DOSS file accompanying an outp file")
    parser.add_argument(
        "--fermi-source",
        choices=("outp", "doss"),
        default="doss",
        help="Fermi-energy reference for a combined outp/DOSS run (default: doss)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        nargs="+",
        default=[300.0],
        metavar="T",
        help="Temperature input: one value (T) or three values (T_START T_END N)",
    )
    parser.add_argument(
        "--chemical-potential",
        type=float,
        nargs="+",
        default=[0.0],
        metavar="MU",
        help="Chemical potential shift(s) in eV: one value (MU) or three values (MU_START MU_END N)",
    )
    parser.add_argument("--tau", type=float, default=1e-14, help="Relaxation time in seconds")
    parser.add_argument(
        "--kmesh",
        type=int,
        nargs=3,
        default=[80, 80, 80],
        metavar=("NX", "NY", "NZ"),
        help="k-point mesh dimensions",
    )
    parser.add_argument("--lr-ratio", type=int, default=20, help="Interpolator star-vector ratio")
    parser.add_argument(
        "--band-window",
        type=float,
        nargs=2,
        default=[-3.0, 3.0],
        metavar=("EMIN", "EMAX"),
        help="Band window relative to Fermi level in eV",
    )
    parser.add_argument("--chunk-size", type=int, default=4096, help="Chunk size for batched evaluations")
    parser.add_argument(
        "--parser",
        choices=(*available_parsers(), CrystalDOSSParser.name),
        default=None,
        help="Input parser (default: detect the file format)",
    )
    parser.add_argument(
        "--calculator",
        choices=available_calculators(),
        default="boltzmann",
        help="Transport calculator",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output JSON path (default: transport_results.json or dos_results.json)",
    )
    parser.add_argument(
        "--checkpoint",
        default=DEFAULT_CHECKPOINT_PATH,
        help=f"Restartable .npz checkpoint path (default: {DEFAULT_CHECKPOINT_PATH})",
    )
    parser.add_argument(
        "--no-checkpoint",
        action="store_true",
        help="Disable checkpoint writing and resume for this run",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Ignore an existing checkpoint and overwrite it as the run progresses",
    )
    return parser


def _to_jsonable(value: Any) -> Any:
    """Convert nested numpy-rich results to JSON-serializable structure."""
    converted: Any = value
    if isinstance(value, dict):
        converted = {str(k): _to_jsonable(v) for k, v in value.items()}
    elif isinstance(value, (list, tuple)):
        converted = [_to_jsonable(v) for v in value]
    elif isinstance(value, np.ndarray):
        converted = value.tolist()
    elif isinstance(value, np.floating):
        converted = float(value)
    elif isinstance(value, np.integer):
        converted = int(value)
    elif isinstance(value, np.bool_):
        converted = bool(value)
    return converted


def _input_parser(filepath: str) -> str:
    """Detect band or DOS input from its header and the outp extension."""
    with Path(filepath).open() as file_obj:
        prefix = file_obj.read(INPUT_DETECTION_BYTES).lstrip()
    if prefix.startswith("-%-") or (prefix.startswith("#") and "NPROJ" in prefix.splitlines()[0].upper()):
        return CrystalDOSSParser.name
    if Path(filepath).suffix.lower() == ".outp" or (not prefix.startswith("<") and "CRYSTAL" in prefix):
        return "crystal-outp"
    return "fleur-outxml"


def _dos_payload(parsed: DOSResult, filepath: str) -> dict[str, Any]:
    """Return a JSON-compatible schema for parsed electronic DOS."""
    return {
        "parser": CrystalDOSSParser.name,
        "filepath": str(Path(filepath).resolve()),
        "energies": parsed.energies,
        "absolute_energies": parsed.absolute_energies,
        "dos": parsed.dos,
        "fermi_energy": parsed.fermi_energy,
        "jspins": parsed.jspins,
        "nenergy": parsed.nenergy,
        "nprojections": parsed.nprojections,
        "units": {"energies": "eV relative to Fermi energy", "fermi_energy": "eV", "dos": "states/eV/cell"},
    }


def _run_transport(args: argparse.Namespace, parser_name: str, dos: DOSResult | None) -> dict[Any, Any]:
    """Run band transport, optionally using the accompanying DOS Fermi reference."""
    temperature = _parse_scalar_or_grid(args.temperature, "temperature")
    chemical_potential = _parse_scalar_or_grid(args.chemical_potential, "chemical_potential")
    fermi_energy = None
    if dos is not None:
        bands = get_parser(parser_name).parse(args.filepath)
        if bands.jspins != dos.jspins:
            raise ValueError(f"Spin channels differ between outp ({bands.jspins}) and DOSS ({dos.jspins})")
        if args.fermi_source == "doss":
            fermi_energy = dos.fermi_energy
        print(f"Combined CRYSTAL input: transport Fermi reference from {args.fermi_source}")
    return calculate_spin_polarized_transport(
        args.filepath,
        temperature=temperature,
        chemical_potential=chemical_potential,
        tau=args.tau,
        kpoint_mesh=tuple(args.kmesh),
        lr_ratio=args.lr_ratio,
        band_window=tuple(args.band_window),
        chunk_size=args.chunk_size,
        parser=parser_name,
        calculator=args.calculator,
        checkpoint_path=None if args.no_checkpoint else args.checkpoint,
        resume_checkpoint=not args.no_resume,
        fermi_energy=fermi_energy,
    )


def _run_input(args: argparse.Namespace) -> tuple[dict[Any, Any], str]:
    """Dispatch DOS export or transport with an optional companion DOS file."""
    parser_name = args.parser or _input_parser(args.filepath)
    if args.dos_filepath and args.dos_file:
        raise ValueError("Specify the accompanying DOSS file either positionally or with --dos-file")
    dos_filepath = args.dos_filepath or args.dos_file
    if dos_filepath and parser_name != "crystal-outp":
        raise ValueError("An accompanying DOSS file requires a CRYSTAL outp band input")
    if parser_name == CrystalDOSSParser.name:
        parsed_dos = CrystalDOSSParser().parse(args.filepath)
        return _dos_payload(parsed_dos, args.filepath), "dos_results.json"
    dos = CrystalDOSSParser().parse(dos_filepath) if dos_filepath else None
    result = _run_transport(args, parser_name, dos)
    if dos is not None:
        result["dos"] = _dos_payload(dos, dos_filepath)
        result["meta"]["fermi_source"] = args.fermi_source
    return result, "transport_results.json"


def main() -> None:
    """Run band transport or export CRYSTAL DOS from CLI input."""
    parser = _build_parser()
    args = parser.parse_args()
    try:
        result, default_output = _run_input(args)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    output = args.output or default_output
    json_payload = _to_jsonable(result)
    with open(output, "w", encoding="utf-8") as fobj:
        json.dump(json_payload, fobj, indent=2, ensure_ascii=False)
    print(f"\nSaved JSON results to {output}")


if __name__ == "__main__":
    main()
