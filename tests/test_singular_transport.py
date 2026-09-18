"""Regression tests for transport windows with inactive conduction directions."""

import numpy as np

from quantum_macaroni.calculators.transport import _onsager_to_transport
from quantum_macaroni.core.constants import E_CHARGE


def test_empty_transport_window() -> None:
    """An empty thermal energy window produces zero tensors without inversion errors."""
    zero = np.zeros((3, 3))
    for tensor in _onsager_to_transport(zero, zero, zero, 300.0):
        np.testing.assert_array_equal(tensor, zero)


def test_inactive_direction() -> None:
    """Retain conducting directions while projecting a singular null direction out."""
    l0 = np.diag([2.0, 4.0, 0.0])
    l1 = np.diag([1.0, 2.0, 0.0])
    l2 = np.diag([3.0, 6.0, 0.0])
    sigma, seebeck, kappa = _onsager_to_transport(l0, l1, l2, 300.0)
    np.testing.assert_allclose(sigma, E_CHARGE * l0, atol=0)
    np.testing.assert_allclose(seebeck, np.diag([-0.5, -0.5, 0.0]) / 300.0, atol=0)
    np.testing.assert_allclose(kappa, np.diag([2.5, 5.0, 0.0]) * E_CHARGE / 300.0, atol=0)
