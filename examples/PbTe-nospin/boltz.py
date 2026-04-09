"""Backward-compatible facade exposing historical imports from modular internals."""

import numpy as np

from boltzpy.calculators.transport import (
    BoltzmannTransportCalculator,
    calculate_spin_polarized_transport,
)
from boltzpy.interpolation import SKWInterpolator
from boltzpy.mesh import TetrahedronMesh
from boltzpy.parsers.fleur_outxml import (
    FleurOutxmlParser,
    parse_fleur_outxml,
    structure_from_outxml,
)
from boltzpy.parsers.fleur_outxml import (
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
    default_file = "out-nospin.xml"
    temperature = [300.0, 600.0, 900.0]
    chemical_potential = np.linspace(-0.5, 0.5, 11)
    tau = 1e-14
    mesh = (80, 80, 80)
    lr_ratio = 20
    band_window = (-3, 3)
    chunk_size = 4096

    print("=" * 60)
    print("TEST: Non-spin-polarized PbTe (chemical potential scan)")
    print("=" * 60)
    result = calculate_spin_polarized_transport(
        default_file,
        temperature=temperature,
        chemical_potential=chemical_potential,
        tau=tau,
        kpoint_mesh=mesh,
        lr_ratio=lr_ratio,
        band_window=band_window,
        chunk_size=chunk_size,
    )

    meta = result["meta"]
    print(f"  E_Fermi = {meta['fermi_energy']:.4f} eV")

    for mu in chemical_potential:
        mu_key = float(mu)
        print(f"\n  mu = E_F + {mu_key:+.4f} eV")
        for temp in temperature:
            data = result[mu_key][float(temp)]
            print(
                f"    T={temp:6.1f} K: "
                f"sigma={data['sigma_avg']:.4e} S/m, "
                f"seebeck={data['seebeck_avg'] * 1e6:.2f} uV/K, "
                f"kappa={data['kappa_avg']:.4f} W/(m*K)"
            )
