"""SKW interpolation of band energies and velocities on arbitrary k-points."""

import itertools
from collections import deque

import numpy as np
import numpy.typing as npt
import scipy.linalg
from scipy.special import erfc

from boltzpy.core.constants import TWO_PI
from boltzpy.core.numerics import (
    nb_eval_energy_from_star,
    nb_eval_energy_velocity_from_star,
    nb_star_and_grad_batch,
    nb_star_batch,
)
from boltzpy.core.symmetry import point_group


class SKWInterpolator:
    """Star-function interpolator for periodic band structures."""

    def __init__(
        self,
        kpoints: npt.ArrayLike,
        eigenvalues: npt.ArrayLike,
        cell: npt.ArrayLike,
        symops: npt.ArrayLike,
        time_reversal: bool = True,
        lr_ratio: int = 5,
        filter_params: tuple[float, float] | None = None,
    ) -> None:
        """Fit interpolation coefficients from sampled k-point eigenvalues.

        Args:
            kpoints: Sampled fractional k-points with shape ``(nk, 3)``.
            eigenvalues: Sampled band energies with shape ``(nspin, nk, nbands)``.
            cell: Lattice matrix or tuple-like object whose first item is lattice.
            symops: Integer symmetry operations.
            time_reversal: Whether to include time-reversal completion for symmetries.
            lr_ratio: Number of star generators relative to input k-points.
            filter_params: Optional ``(rcut_scale, sigma)`` filter for high-R coefficients.

        Raises:
            ValueError: If k-point array has unexpected shape.

        """
        eigenvalues = np.ascontiguousarray(np.atleast_3d(eigenvalues), dtype=np.float64)
        self.nspin, self.nk, self.nbands = eigenvalues.shape
        kpoints = np.ascontiguousarray(np.asarray(kpoints, dtype=np.float64))
        # We require one 3D fractional coordinate per sampled eigenvalue and at least two
        # k-points so the reference-subtracted linear system below is well-defined.
        # Fool protection.
        if kpoints.shape != (self.nk, 3) or self.nk < 2:  # noqa: PLR2004
            raise ValueError("kpoints must have shape (nk, 3) with nk >= 2")

        cell_arr = np.asarray(cell, dtype=np.float64)
        lat = np.ascontiguousarray(cell_arr if cell_arr.shape == (3, 3) else np.asarray(cell[0], dtype=np.float64))
        self._lat = lat
        self._frac_to_cart_t = np.ascontiguousarray((lat / TWO_PI).T, dtype=np.float64)
        self._recip_metric = (TWO_PI * TWO_PI) * np.linalg.inv(lat @ lat.T)

        self._pg = np.ascontiguousarray(point_group(symops, time_reversal), dtype=np.int64)
        self._npg = len(self._pg)

        n_need = lr_ratio * self.nk
        has_inv = any(np.array_equal(rot, -np.eye(3, dtype=int)) for rot in self._pg)
        vf = 0.5 if has_inv else 1.0
        # The search radius is estimated from the target number of star generators so we start
        # near a workable basis size instead of repeatedly exploring obviously too-small shells.
        rm = int((1.0 + (lr_ratio * self.nk * self._npg * vf) / 2.0) ** (1.0 / 3.0)) + 1
        rmax = np.array([rm, rm, rm], dtype=int)

        while True:
            self._rpts, r2, ok = self._find_stars(n_need, rmax)
            self.nr = len(self._rpts)
            if ok:
                break
            rmax *= 2

        self._rot_rpts = np.ascontiguousarray(np.array([rot @ self._rpts.T for rot in self._pg], dtype=np.float64))
        self._rot_rpts_flat = np.ascontiguousarray(self._rot_rpts.reshape(-1), dtype=np.float64)

        c1 = 0.25
        c2 = 0.25
        if self.nr > 1 and r2[1] > 0:
            ratio = r2 / r2[1]
            # Damping long-range stars regularizes the linear system and reduces the tendency
            # of SKW fits to oscillate when the basis grows faster than the input data quality.
            inv_rho = 1.0 / ((1.0 - c1 * ratio) ** 2 + c2 * ratio**3)
            inv_rho[0] = 1.0
        else:
            inv_rho = np.ones(self.nr, dtype=np.float64)

        sk = nb_star_batch(kpoints, self._rot_rpts_flat, self._npg, self.nr)

        nm = self.nk - 1
        # Building the system relative to the last k-point removes the null-mode associated with
        # a constant energy shift and makes the Hermitian solve better conditioned.
        dsk = sk[:nm, 1:] - sk[nm, 1:]
        hmat = (dsk * inv_rho[1:]) @ dsk.conj().T
        np.fill_diagonal(hmat, hmat.diagonal().real)

        de = np.empty((nm, self.nbands, self.nspin), dtype=np.complex128)
        for spin in range(self.nspin):
            for band in range(self.nbands):
                de[:, band, spin] = eigenvalues[spin, :nm, band] - eigenvalues[spin, nm, band]

        lam = scipy.linalg.solve(
            hmat,
            de.reshape(nm, -1),
            assume_a="her",
            check_finite=False,
        ).reshape(nm, self.nbands, self.nspin)

        self._coeffs = np.empty((self.nspin, self.nbands, self.nr), dtype=np.complex128)
        inner = np.einsum("mr,mbs->rbs", dsk.conj(), lam)
        self._coeffs[:, :, 1:] = np.einsum("rbs->sbr", inner * inv_rho[1:, None, None])
        self._coeffs[:, :, 0] = eigenvalues[:, nm, :] - np.einsum("sbr,r->sb", self._coeffs[:, :, 1:], sk[nm, 1:])
        self._coeffs = np.ascontiguousarray(self._coeffs)

        if filter_params is not None:
            rcut = filter_params[0] * np.sqrt(r2[-1])
            sigma = filter_params[1]

            for ir in range(1, self.nr):
                # A soft cutoff is used instead of a hard truncation so high-R coefficients fade
                # out smoothly and do not inject ringing back into the interpolated bands.
                self._coeffs[:, :, ir] *= 0.5 * erfc((np.sqrt(r2[ir]) - rcut) / sigma)

        evals_pred = np.einsum("sbr,kr->skb", self._coeffs, sk).real
        self.mae = np.abs(eigenvalues - evals_pred).mean()

    def _find_stars(
        self,
        n_need: int,
        rmax: npt.NDArray[np.int_],
    ) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.float64], bool]:
        """Construct symmetry-unique reciprocal vectors used as star generators.

        Args:
            n_need: Requested number of generators.
            rmax: Search radius in reciprocal-lattice index units.

        Returns:
            Tuple with generator vectors, squared lengths, and sufficiency flag.

        """
        ranges = [range(-int(rmax[i]), int(rmax[i]) + 1) for i in range(3)]
        pts = np.array(list(itertools.product(*ranges)), dtype=int)
        r2 = np.einsum("ij,jk,ik->i", pts, self._recip_metric, pts)

        order = np.argsort(r2)
        pts = pts[order]
        r2 = r2[order]

        bounds = [0]
        for idx in range(1, len(r2)):
            left = r2[bounds[-1]]
            right = r2[idx]
            # Shells are grouped by nearly equal metric length so symmetry representatives are
            # chosen from physically equivalent reciprocal distances rather than raw indices.
            if abs(right - left) > max(left, right) * 1e-8 + 1e-30:
                bounds.append(idx)
        bounds.append(len(r2))

        generators = deque()
        for shell_idx in range(len(bounds) - 1):
            lo = bounds[shell_idx]
            hi = bounds[shell_idx + 1]
            if hi - lo == 1:
                generators.append(tuple(pts[lo]))
                continue

            seen = {tuple(pts[lo])}
            for idx in range(lo + 1, hi):
                vector = pts[idx]
                if all(tuple(rot @ vector) not in seen for rot in self._pg):
                    seen.add(tuple(vector))
            generators.extend(seen)

        generators = np.array(list(generators), dtype=int)
        gr2 = np.einsum("ij,jk,ik->i", generators, self._recip_metric, generators)
        order = np.argsort(gr2)
        generators = generators[order]
        gr2 = gr2[order]

        n_take = min(len(generators), n_need)
        return generators[:n_take], gr2[:n_take], len(generators) >= n_need

    def evaluate(self, kpoints_frac: npt.ArrayLike, kchunk: int = 4096) -> npt.NDArray[np.float64]:
        """Evaluate interpolated band energies on fractional k-points.

        Args:
            kpoints_frac: Fractional k-points to evaluate.
            kchunk: Chunk size for batched evaluation.

        Returns:
            Interpolated energies with shape ``(nspin, nk, nbands)``.

        """
        kpoints_frac = np.ascontiguousarray(np.atleast_2d(kpoints_frac), dtype=np.float64)
        nk = len(kpoints_frac)
        out = np.empty((self.nspin, nk, self.nbands), dtype=np.float64)

        for start in range(0, nk, kchunk):
            end = min(start + kchunk, nk)
            sk = nb_star_batch(kpoints_frac[start:end], self._rot_rpts_flat, self._npg, self.nr)
            for spin in range(self.nspin):
                out[spin, start:end] = nb_eval_energy_from_star(self._coeffs[spin], sk).T

        return out

    def eval_energy_velocity(
        self,
        kpoints_frac: npt.ArrayLike,
        spin: int,
        chunk_size: int = 4096,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Evaluate interpolated energies and group velocities for one spin channel.

        Args:
            kpoints_frac: Fractional k-points to evaluate.
            spin: Spin-channel index.
            chunk_size: Chunk size for batched evaluation.

        Returns:
            Tuple ``(e_all, vel_all)`` with shapes ``(nbands, nk)`` and ``(nbands, nk, 3)``.

        """
        kpoints_frac = np.ascontiguousarray(np.atleast_2d(kpoints_frac), dtype=np.float64)
        nk = len(kpoints_frac)
        e_all = np.empty((self.nbands, nk), dtype=np.float64)
        vel_all = np.empty((self.nbands, nk, 3), dtype=np.float64)
        coeffs_spin = self._coeffs[spin]

        for start in range(0, nk, chunk_size):
            end = min(start + chunk_size, nk)
            sk, dsf = nb_star_and_grad_batch(kpoints_frac[start:end], self._rot_rpts_flat, self._npg, self.nr)
            e_chunk, vel_chunk = nb_eval_energy_velocity_from_star(coeffs_spin, sk, dsf, self._frac_to_cart_t)
            e_all[:, start:end] = e_chunk
            vel_all[:, start:end] = vel_chunk

        return e_all, vel_all
