"""Backward-compatible facade exposing imports from modular internals."""

import sys
from pathlib import Path
from typing import Any, cast

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from quantum_macaroni.calculators.transport import (  # noqa: E402
    BoltzmannTransportCalculator,
    calculate_spin_polarized_transport,
)
from quantum_macaroni.interpolation import SKWInterpolator  # noqa: E402
from quantum_macaroni.mesh import TetrahedronMesh  # noqa: E402
from quantum_macaroni.parsers.fleur_outxml import (  # noqa: E402
    FleurOutxmlParser,
    parse_fleur_outxml,
    structure_from_outxml,
)
from quantum_macaroni.parsers.fleur_outxml import (  # noqa: E402
    read_symops_from_outxml as _read_symops_from_outxml,
)

__all__ = [
    "SKWInterpolator",
    "TetrahedronMesh",
    "BoltzmannTransportCalculator",
    "parse_fleur_outxml",
    "structure_from_outxml",
    "_read_symops_from_outxml",
    "calculate_spin_polarized_transport",
    "FleurOutxmlParser",
]


if __name__ == "__main__":
    example_dir = Path(__file__).resolve().parent
    default_file = example_dir / "out-nospin.xml"
    checkpoint_path = example_dir / "pbte_transport_state.npz"
    temperature = [300.0]
    chemical_potential = np.linspace(-0.5, 0.5, 3)
    tau = 1e-14
    mesh = (96, 96, 96)
    lr_ratio = 25
    band_window = (-3, 3)
    chunk_size = 4096

    print("=" * 60)
    print("TEST: Non-spin-polarized PbTe (checkpointed chemical potential scan)")
    print("=" * 60)
    print(f"  Input: {default_file}")
    print(f"  Checkpoint: {checkpoint_path}")
    print("  Re-run with the same settings to resume from the latest compatible checkpoint.")
    result = calculate_spin_polarized_transport(
        default_file,
        temperature=temperature,
        chemical_potential=chemical_potential,
        tau=tau,
        kpoint_mesh=mesh,
        lr_ratio=lr_ratio,
        band_window=band_window,
        chunk_size=chunk_size,
        checkpoint_path=checkpoint_path,
        resume_checkpoint=True,
    )
    result_by_mu = cast(dict[float | str, Any], result)

    meta = result_by_mu["meta"]
    print(f"  E_Fermi = {meta['fermi_energy']:.4f} eV")

    for mu in chemical_potential:
        mu_key = float(mu)
        print(f"\n  mu = E_F + {mu_key:+.4f} eV")
        for temp in temperature:
            data = result_by_mu[mu_key][float(temp)]
            print(
                f"    T={temp:6.1f} K: "
                f"sigma={data['sigma_avg']:.4e} S/m, "
                f"seebeck={data['seebeck_avg'] * 1e6:.2f} uV/K, "
                f"kappa={data['kappa_avg']:.4f} W/(m*K)"
            )
