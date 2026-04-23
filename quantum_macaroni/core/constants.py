"""Physical constants and conversion factors used by quantum_macaroni."""

import math

from ase.units import Bohr, Hartree
from scipy.constants import physical_constants

BOHR_TO_ANG = Bohr
HBAR = physical_constants["Planck constant over 2 pi in eV s"][0]
E_CHARGE = physical_constants["elementary charge"][0]
KB_EV = physical_constants["Boltzmann constant in eV/K"][0]
HTR_TO_EV = Hartree
TWO_PI = 2.0 * math.pi
