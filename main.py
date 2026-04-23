"""Project command-line entry point."""

from __future__ import annotations

import argparse
import json
from typing import Any

import numpy as np

from quantum_macaroni import available_calculators, available_parsers, calculate_spin_polarized_transport

GRID_SPEC_LEN = 3


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
    parser = argparse.ArgumentParser(description="Boltzman semiclassical transport calculator")
    parser.add_argument("filepath", help="Path to input electronic-structure file")
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
        choices=available_parsers(),
        default="fleur-outxml",
        help="Electronic-structure parser",
    )
    parser.add_argument(
        "--calculator",
        choices=available_calculators(),
        default="boltzmann",
        help="Transport calculator",
    )
    parser.add_argument(
        "--output",
        default="transport_results.json",
        help="Output JSON file path",
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


def main() -> None:
    """Run transport calculation from CLI."""
    parser = _build_parser()
    args = parser.parse_args()

    try:
        temperature = _parse_scalar_or_grid(args.temperature, "temperature")
        chemical_potential = _parse_scalar_or_grid(args.chemical_potential, "chemical_potential")
    except ValueError as exc:
        parser.error(str(exc))

    result = calculate_spin_polarized_transport(
        args.filepath,
        temperature=temperature,
        chemical_potential=chemical_potential,
        tau=args.tau,
        kpoint_mesh=tuple(args.kmesh),
        lr_ratio=args.lr_ratio,
        band_window=tuple(args.band_window),
        chunk_size=args.chunk_size,
        parser=args.parser,
        calculator=args.calculator,
    )

    json_payload = _to_jsonable(result)
    with open(args.output, "w", encoding="utf-8") as fobj:
        json.dump(json_payload, fobj, indent=2, ensure_ascii=False)
    print(f"\nSaved JSON results to {args.output}")


if __name__ == "__main__":
    main()
