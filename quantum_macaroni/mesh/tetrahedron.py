"""Regular k-space mesh and tetrahedral decomposition utilities."""

import numpy as np
import numpy.typing as npt

from quantum_macaroni.core.constants import TWO_PI


class TetrahedronMesh:
    """Generate full k-point mesh and tetrahedra for Brillouin-zone integration."""

    # For each possible body-diagonal choice, list the 6 tetrahedra as local cube vertex indices
    # (numbering follows the `vertices` construction in `_generate_tetrahedra`).
    _TETRAHEDRA_CONFIGS = np.array(
        [
            [[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]],
            [[1, 0, 2, 6], [1, 0, 4, 6], [1, 3, 2, 6], [1, 3, 7, 6], [1, 5, 4, 6], [1, 5, 7, 6]],
            [[2, 0, 1, 5], [2, 0, 4, 5], [2, 3, 1, 5], [2, 3, 7, 5], [2, 6, 4, 5], [2, 6, 7, 5]],
            [[3, 0, 1, 4], [3, 0, 2, 4], [3, 1, 5, 4], [3, 1, 7, 4], [3, 2, 6, 4], [3, 2, 7, 4]],
        ],
        dtype=np.int32,
    )

    def __init__(self, lattice: npt.ArrayLike, mesh: tuple[int, int, int]) -> None:
        """Build mesh objects for a given real-space lattice and mesh dimensions.

        Args:
            lattice: Real-space lattice matrix.
            mesh: Number of divisions along ``(a, b, c)`` reciprocal axes.

        """
        if any(m <= 0 for m in mesh):
            raise ValueError(f"mesh dimensions must be positive, got {mesh}")

        self.lattice = np.ascontiguousarray(np.asarray(lattice, dtype=np.float64))
        self.mesh = np.array(mesh, dtype=np.int32)
        self.recip_lattice = TWO_PI * np.linalg.inv(self.lattice).T
        self.vol_bz = abs(np.linalg.det(self.recip_lattice))
        self._generate_kpoints()
        self._generate_tetrahedra()

    def _generate_kpoints(self) -> None:
        """Generate full regular fractional k-point grid."""
        nkx, nky, nkz = [int(x) for x in self.mesh]
        ix, iy, iz = np.meshgrid(
            np.arange(nkx, dtype=np.float64),
            np.arange(nky, dtype=np.float64),
            np.arange(nkz, dtype=np.float64),
            indexing="ij",
        )
        self.nk_full = nkx * nky * nkz
        self.full_kpoints = np.empty((self.nk_full, 3), dtype=np.float64)
        self.full_kpoints[:, 0] = ix.reshape(-1) / nkx
        self.full_kpoints[:, 1] = iy.reshape(-1) / nky
        self.full_kpoints[:, 2] = iz.reshape(-1) / nkz
        self.full_kpoints = np.ascontiguousarray(self.full_kpoints)

    def _generate_tetrahedra(self) -> None:
        """Generate tetrahedral decomposition and tetrahedron volume."""
        nkx, nky, nkz = [int(x) for x in self.mesh]
        n_cubes = nkx * nky * nkz

        diagonals = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.float64)
        diag_lengths = np.linalg.norm((self.recip_lattice @ diagonals.T).T, axis=1)
        # Picking the shortest reciprocal-space body diagonal reduces skinny tetrahedra, which
        # improves the numerical behavior of the tetrahedron integration near sharp band features.
        config_idx = int(np.argmin(diag_lengths))

        ix, iy, iz = np.meshgrid(
            np.arange(nkx, dtype=np.int32),
            np.arange(nky, dtype=np.int32),
            np.arange(nkz, dtype=np.int32),
            indexing="ij",
        )
        ix = ix.reshape(-1)
        iy = iy.reshape(-1)
        iz = iz.reshape(-1)
        jx = (ix + 1) % nkx
        jy = (iy + 1) % nky
        jz = (iz + 1) % nkz

        def flat_index(
            x: npt.NDArray[np.int32], y: npt.NDArray[np.int32], z: npt.NDArray[np.int32]
        ) -> npt.NDArray[np.int32]:
            return x * nky * nkz + y * nkz + z

        vertices = np.stack(
            [
                flat_index(ix, iy, iz),
                flat_index(jx, iy, iz),
                flat_index(ix, jy, iz),
                flat_index(jx, jy, iz),
                flat_index(ix, iy, jz),
                flat_index(jx, iy, jz),
                flat_index(ix, jy, jz),
                flat_index(jx, jy, jz),
            ],
            axis=1,
        )

        tetra_config = self._TETRAHEDRA_CONFIGS[config_idx]
        self.tetrahedra = np.ascontiguousarray(vertices[:, tetra_config].reshape(-1, 4), dtype=np.int32)
        self.n_tetrahedra = self.tetrahedra.shape[0]
        self.tetra_vol = self.vol_bz / (6.0 * n_cubes)
