"""Numba-accelerated kernels for tetrahedron transport and SKW evaluation."""

import math

import numba as _nb
import numpy as np
import numpy.typing as npt

from quantum_macaroni.core.constants import HBAR, TWO_PI

# Minimal energy span used to treat a tetrahedron as effectively degenerate.
MIN_TETRAEDRON_ENERGY = 1e-14
# Minimal occupation-weight derivative contribution considered above numerical noise.
MIN_TETRAEDRON_WEIGHT = 1e-30
# kBT threshold below which the zero-temperature integration fallback is used.
MIN_KBT_ZERO_T = 1e-14


@_nb.njit(cache=True, fastmath=True)
def nb_occ_weights(ef: float, e1: float, e2: float, e3: float, e4: float) -> npt.NDArray[np.float64]:
    """Return tetrahedron corner occupation weights at a target energy.

    Args:
        ef: Target energy.
        e1: Sorted corner energy 1.
        e2: Sorted corner energy 2.
        e3: Sorted corner energy 3.
        e4: Sorted corner energy 4.

    Returns:
        Occupation weights for four tetrahedron corners.

    """
    w = np.zeros(4, dtype=np.float64)

    if ef <= e1:
        return w
    if ef >= e4:
        w[0] = 0.25
        w[1] = 0.25
        w[2] = 0.25
        w[3] = 0.25
        return w

    eps = 1e-14
    e21 = max(e2 - e1, eps)
    e31 = max(e3 - e1, eps)
    e41 = max(e4 - e1, eps)
    e32 = max(e3 - e2, eps)
    e42 = max(e4 - e2, eps)
    e43 = max(e4 - e3, eps)

    if ef < e2:
        x = ef - e1
        c = x * x * x / (4.0 * e21 * e31 * e41)
        w[0] = c * (4.0 - x * (1.0 / e21 + 1.0 / e31 + 1.0 / e41))
        w[1] = c * x / e21
        w[2] = c * x / e31
        w[3] = c * x / e41
    elif ef < e3:
        c1 = (ef - e1) * (ef - e1) / (4.0 * e41 * e31)
        c2 = (ef - e1) * (ef - e2) * (e3 - ef) / (4.0 * e41 * e32 * e31)
        c3 = (ef - e2) * (ef - e2) * (e4 - ef) / (4.0 * e42 * e32 * e41)
        c12 = c1 + c2
        c123 = c12 + c3
        w[0] = c1 + c12 * (e3 - ef) / e31 + c123 * (e4 - ef) / e41
        w[1] = c123 + (c2 + c3) * (e3 - ef) / e32 + c3 * (e4 - ef) / e42
        w[2] = c12 * (ef - e1) / e31 + (c2 + c3) * (ef - e2) / e32
        w[3] = c123 * (ef - e1) / e41 + c3 * (ef - e2) / e42
    else:
        x = e4 - ef
        c = x * x * x / (4.0 * e41 * e42 * e43)
        w[0] = 0.25 - c * x / e41
        w[1] = 0.25 - c * x / e42
        w[2] = 0.25 - c * x / e43
        w[3] = 0.25 - c * (4.0 - x * (1.0 / e41 + 1.0 / e42 + 1.0 / e43))

    return w


@_nb.njit(parallel=True, cache=True, fastmath=True)
def nb_transport_dos_flat(  # noqa: C901, PLR0912
    e_all: npt.NDArray[np.float64],
    vel_all: npt.NDArray[np.float64],
    tetrahedra: npt.NDArray[np.int32],
    tau: float,
    e_grid: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Accumulate flattened transport DOS tensor on an energy grid.

    Args:
        e_all: Band energies with shape ``(nbands, nk)``.
        vel_all: Velocities with shape ``(nbands, nk, 3)``.
        tetrahedra: Tetrahedron vertex indices with shape ``(ntet, 4)``.
        tau: Relaxation time in seconds.
        e_grid: Energy grid in eV.

    Returns:
        Flattened transport DOS with shape ``(ne, 9)``.

    """
    nbands = e_all.shape[0]
    ne = e_grid.shape[0]
    ntet = tetrahedra.shape[0]
    de = e_grid[1] - e_grid[0]
    half_de = 0.5 * de
    inv_de = 1.0 / de
    e_min = e_grid[0]

    tdos_bands = np.zeros((nbands, ne, 9), dtype=np.float64)

    for ib in _nb.prange(nbands):  # ty:ignore[not-iterable]
        acc = tdos_bands[ib]
        for it in range(ntet):
            i0 = tetrahedra[it, 0]
            i1 = tetrahedra[it, 1]
            i2 = tetrahedra[it, 2]
            i3 = tetrahedra[it, 3]

            ee0 = e_all[ib, i0]
            ee1 = e_all[ib, i1]
            ee2 = e_all[ib, i2]
            ee3 = e_all[ib, i3]

            if ee0 > ee1:
                ee0, ee1 = ee1, ee0
                i0, i1 = i1, i0
            if ee2 > ee3:
                ee2, ee3 = ee3, ee2
                i2, i3 = i3, i2
            if ee0 > ee2:
                ee0, ee2 = ee2, ee0
                i0, i2 = i2, i0
            if ee1 > ee3:
                ee1, ee3 = ee3, ee1
                i1, i3 = i3, i1
            if ee1 > ee2:
                ee1, ee2 = ee2, ee1
                i1, i2 = i2, i1

            # The tetrahedron weights assume ordered corner energies; sorting locally is cheaper
            # than carrying a pre-sorted connectivity for every band and k-mesh.
            if ee3 - ee0 < MIN_TETRAEDRON_ENERGY:
                continue

            ie_lo = int((ee0 - half_de - e_min) * inv_de)
            ie_lo = max(ie_lo, 0)

            ie_hi = int((ee3 + half_de - e_min) * inv_de) + 1
            if ie_hi >= ne:
                ie_hi = ne - 1

            for ie in range(ie_lo, ie_hi + 1):
                eps = e_grid[ie]
                eps_lo = eps - half_de
                eps_hi = eps + half_de
                if eps_hi <= ee0 or eps_lo >= ee3:
                    continue

                w_lo = nb_occ_weights(eps_lo, ee0, ee1, ee2, ee3)
                w_hi = nb_occ_weights(eps_hi, ee0, ee1, ee2, ee3)

                for ic in range(4):
                    # `ic` is the local corner index (0..3) inside the current tetrahedron.
                    # Map local tetrahedron corner order to the corresponding global k-point index.
                    if ic == 0:
                        ki = i0
                    elif ic == 1:
                        ki = i1
                    elif ic == 2:  # noqa: PLR2004
                        ki = i2
                    else:
                        ki = i3

                    dw = (w_hi[ic] - w_lo[ic]) * inv_de
                    # Tiny contributions are skipped to avoid spending time on numerical noise
                    # that is below the resolution implied by the chosen energy grid.
                    if dw < MIN_TETRAEDRON_WEIGHT:
                        continue

                    dw_tau = dw * tau
                    va0 = vel_all[ib, ki, 0]
                    va1 = vel_all[ib, ki, 1]
                    va2 = vel_all[ib, ki, 2]
                    acc[ie, 0] += dw_tau * va0 * va0
                    acc[ie, 1] += dw_tau * va0 * va1
                    acc[ie, 2] += dw_tau * va0 * va2
                    acc[ie, 3] += dw_tau * va1 * va0
                    acc[ie, 4] += dw_tau * va1 * va1
                    acc[ie, 5] += dw_tau * va1 * va2
                    acc[ie, 6] += dw_tau * va2 * va0
                    acc[ie, 7] += dw_tau * va2 * va1
                    acc[ie, 8] += dw_tau * va2 * va2

    tdos = np.zeros((ne, 9), dtype=np.float64)
    for ie in _nb.prange(ne):  # ty:ignore[not-iterable]
        for ib in range(nbands):
            for ab in range(9):
                tdos[ie, ab] += tdos_bands[ib, ie, ab]

    return tdos


@_nb.njit(cache=True, fastmath=True)
def nb_onsager_from_tdos_flat(
    tdos: npt.NDArray[np.float64],
    e_grid: npt.NDArray[np.float64],
    fermi: float,
    kbt: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Integrate transport DOS into flattened Onsager tensors.

    Args:
        tdos: Flattened transport DOS with shape ``(ne, 9)``.
        e_grid: Energy grid in eV.
        fermi: Fermi level in eV.
        kbt: Thermal energy in eV.

    Returns:
        Tuple ``(l0, l1, l2)`` where each element is a flattened 3x3 tensor.

    """
    ne = e_grid.shape[0]
    de = e_grid[1] - e_grid[0]

    l0 = np.zeros(9, dtype=np.float64)
    l1 = np.zeros(9, dtype=np.float64)
    l2 = np.zeros(9, dtype=np.float64)

    if kbt < MIN_KBT_ZERO_T:
        # In the zero-temperature limit the Fermi derivative collapses to a delta-like peak,
        # so the nearest grid point is a robust discrete approximation.
        idx = 0
        dmin = abs(e_grid[0] - fermi)
        for ie in range(1, ne):
            dist = abs(e_grid[ie] - fermi)
            if dist < dmin:
                dmin = dist
                idx = ie
        for ab in range(9):
            l0[ab] = tdos[idx, ab]
        return l0, l1, l2

    for ie in range(ne):
        energy = e_grid[ie]
        x = (energy - fermi) / kbt

        ch = math.cosh(0.5 * x)
        dfde = 1.0 / (4.0 * kbt * ch * ch)
        e_mu = energy - fermi
        w0 = dfde * de
        w1 = w0 * e_mu
        w2 = w1 * e_mu
        for ab in range(9):
            value = tdos[ie, ab]
            l0[ab] += value * w0
            l1[ab] += value * w1
            l2[ab] += value * w2

    return l0, l1, l2


@_nb.njit(parallel=True, cache=True, fastmath=True)
def nb_star_batch(
    kpoints: npt.NDArray[np.float64],
    rot_rpts_flat: npt.NDArray[np.float64],
    npg: int,
    nr: int,
) -> npt.NDArray[np.complex128]:
    """Evaluate star-function basis values on a k-point batch.

    Args:
        kpoints: Fractional k-points with shape ``(nk, 3)``.
        rot_rpts_flat: Flattened rotated star vectors.
        npg: Number of point-group operations.
        nr: Number of star vectors.

    Returns:
        Complex star basis values with shape ``(nk, nr)``.

    """
    nk = kpoints.shape[0]
    sk = np.zeros((nk, nr), dtype=np.complex128)
    inv_npg = 1.0 / npg

    for k in _nb.prange(nk):  # ty:ignore[not-iterable]
        kx = kpoints[k, 0]
        ky = kpoints[k, 1]
        kz = kpoints[k, 2]
        for ipg in range(npg):
            base = ipg * 3 * nr
            r0_base = base
            r1_base = base + nr
            r2_base = base + 2 * nr
            for r in range(nr):
                phase = TWO_PI * (
                    kx * rot_rpts_flat[r0_base + r] + ky * rot_rpts_flat[r1_base + r] + kz * rot_rpts_flat[r2_base + r]
                )
                sk[k, r] += complex(math.cos(phase), math.sin(phase))

        for r in range(nr):
            sk[k, r] *= inv_npg

    return sk


@_nb.njit(parallel=True, cache=True, fastmath=True)
def nb_star_and_grad_batch(
    kpoints: npt.NDArray[np.float64],
    rot_rpts_flat: npt.NDArray[np.float64],
    npg: int,
    nr: int,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """Evaluate star-function basis and gradients on a k-point batch.

    Args:
        kpoints: Fractional k-points with shape ``(nk, 3)``.
        rot_rpts_flat: Flattened rotated star vectors.
        npg: Number of point-group operations.
        nr: Number of star vectors.

    Returns:
        Tuple ``(sk, dsf)`` for basis values and gradients.

    """
    nk = kpoints.shape[0]
    sk = np.zeros((nk, nr), dtype=np.complex128)
    dsf = np.zeros((nk, 3, nr), dtype=np.complex128)
    inv_npg = 1.0 / npg
    grad_factor = complex(0.0, TWO_PI * inv_npg)

    for k in _nb.prange(nk):  # ty:ignore[not-iterable]
        kx = kpoints[k, 0]
        ky = kpoints[k, 1]
        kz = kpoints[k, 2]
        for ipg in range(npg):
            base = ipg * 3 * nr
            r0_base = base
            r1_base = base + nr
            r2_base = base + 2 * nr
            for r in range(nr):
                r0 = rot_rpts_flat[r0_base + r]
                r1 = rot_rpts_flat[r1_base + r]
                r2 = rot_rpts_flat[r2_base + r]
                phase = TWO_PI * (kx * r0 + ky * r1 + kz * r2)
                ph = complex(math.cos(phase), math.sin(phase))
                sk[k, r] += ph
                dsf[k, 0, r] += ph * r0
                dsf[k, 1, r] += ph * r1
                dsf[k, 2, r] += ph * r2

        for r in range(nr):
            sk[k, r] *= inv_npg
            dsf[k, 0, r] *= grad_factor
            dsf[k, 1, r] *= grad_factor
            dsf[k, 2, r] *= grad_factor

    return sk, dsf


@_nb.njit(parallel=True, cache=True, fastmath=True)
def nb_eval_energy_from_star(
    coeffs_spin: npt.NDArray[np.complex128],
    sk: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]:
    """Evaluate interpolated band energies from star coefficients.

    Args:
        coeffs_spin: Coefficients for one spin channel with shape ``(nbands, nr)``.
        sk: Star basis values with shape ``(nk, nr)``.

    Returns:
        Interpolated energies with shape ``(nbands, nk)``.

    """
    nbands = coeffs_spin.shape[0]
    nr = coeffs_spin.shape[1]
    nk = sk.shape[0]
    out = np.empty((nbands, nk), dtype=np.float64)

    for b in _nb.prange(nbands):  # ty:ignore[not-iterable]
        for k in range(nk):
            total = 0.0 + 0.0j
            for r in range(nr):
                total += coeffs_spin[b, r] * sk[k, r]
            out[b, k] = total.real

    return out


@_nb.njit(parallel=True, cache=True, fastmath=True)
def nb_eval_energy_velocity_from_star(
    coeffs_spin: npt.NDArray[np.complex128],
    sk: npt.NDArray[np.complex128],
    dsf: npt.NDArray[np.complex128],
    frac_to_cart_t: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Evaluate interpolated energies and velocities.

    Args:
        coeffs_spin: Coefficients for one spin channel with shape ``(nbands, nr)``.
        sk: Star basis values with shape ``(nk, nr)``.
        dsf: Star basis gradients with shape ``(nk, 3, nr)``.
        frac_to_cart_t: Transform from fractional to Cartesian reciprocal basis.

    Returns:
        Tuple ``(e_all, vel_all)``.

    """
    nbands = coeffs_spin.shape[0]
    nr = coeffs_spin.shape[1]
    nk = sk.shape[0]

    e_all = np.empty((nbands, nk), dtype=np.float64)
    vel_all = np.empty((nbands, nk, 3), dtype=np.float64)

    f00 = frac_to_cart_t[0, 0]
    f01 = frac_to_cart_t[0, 1]
    f02 = frac_to_cart_t[0, 2]
    f10 = frac_to_cart_t[1, 0]
    f11 = frac_to_cart_t[1, 1]
    f12 = frac_to_cart_t[1, 2]
    f20 = frac_to_cart_t[2, 0]
    f21 = frac_to_cart_t[2, 1]
    f22 = frac_to_cart_t[2, 2]

    for b in _nb.prange(nbands):  # ty:ignore[not-iterable]
        for k in range(nk):
            e_sum = 0.0 + 0.0j
            g0 = 0.0 + 0.0j
            g1 = 0.0 + 0.0j
            g2 = 0.0 + 0.0j

            for r in range(nr):
                coeff = coeffs_spin[b, r]
                e_sum += coeff * sk[k, r]
                g0 += coeff * dsf[k, 0, r]
                g1 += coeff * dsf[k, 1, r]
                g2 += coeff * dsf[k, 2, r]

            gf0 = g0.real
            gf1 = g1.real
            gf2 = g2.real
            e_all[b, k] = e_sum.real

            vel_all[b, k, 0] = (gf0 * f00 + gf1 * f01 + gf2 * f02) / HBAR * 1e-10
            vel_all[b, k, 1] = (gf0 * f10 + gf1 * f11 + gf2 * f12) / HBAR * 1e-10
            vel_all[b, k, 2] = (gf0 * f20 + gf1 * f21 + gf2 * f22) / HBAR * 1e-10

    return e_all, vel_all
