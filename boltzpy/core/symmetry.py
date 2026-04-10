"""Symmetry utilities for point-group construction."""

import numpy as np
import numpy.typing as npt


def point_group(symops: npt.ArrayLike, time_reversal: bool) -> npt.NDArray[np.int_]:
    """Return unique point-group operations with optional time-reversal completion.

    Args:
        symops: Symmetry operations convertible to ``(nsym, 3, 3)`` integer array.
        time_reversal: Whether to add ``-R`` operations if inversion is missing.

    Returns:
        Unique point-group operations.

    """
    symops = np.reshape(np.asarray(symops, dtype=int), (-1, 3, 3))
    unique = [symops[0]]
    for op in symops[1:]:
        if not any(np.array_equal(op, ref) for ref in unique):
            unique.append(op)
    has_inv = any(np.array_equal(op, -np.eye(3, dtype=int)) for op in unique)
    # When inversion is absent, adding -R mimics the symmetry completion that time reversal
    # would otherwise contribute, which keeps the interpolation basis physically consistent.
    if not has_inv and time_reversal:
        return np.array(list(unique) + [-op for op in unique], dtype=int)
    return np.array(unique, dtype=int)
