"""Calculator protocol and calculator registry helpers."""

from typing import Protocol

import numpy as np
import numpy.typing as npt

TransportTuple = tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]


class TransportCalculator(Protocol):
    """Protocol for transport calculator implementations."""

    name: str

    def calculate_transport(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> TransportTuple:
        """Return transport tensors for the given state point.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional evaluation chunk size.

        Returns:
            Tuple ``(sigma, seebeck, kappa, l0, l1, l2)``.

        """
        ...


_CALCULATORS: dict[str, type[TransportCalculator]] = {}


def register_calculator(name: str, calculator_cls: type[TransportCalculator]) -> None:
    """Register a calculator class under ``name``."""
    _CALCULATORS[name] = calculator_cls


def get_calculator(name: str) -> type[TransportCalculator]:
    """Return calculator class registered under ``name``."""
    try:
        return _CALCULATORS[name]
    except KeyError as exc:
        available = ", ".join(sorted(_CALCULATORS))
        raise ValueError(f"Unknown calculator '{name}'. Available: {available}") from exc


def available_calculators() -> tuple[str, ...]:
    """Return sorted names of all registered calculators."""
    return tuple(sorted(_CALCULATORS))


def tensor_average(tensor: npt.NDArray[np.float64]) -> float:
    """Return isotropic average as one-third of tensor trace."""
    return float(np.trace(tensor) / 3.0)
