"""Boltzmann transport calculators and high-level orchestration entry point."""

import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from boltzpy.calculators.base import get_calculator, register_calculator, tensor_average
from boltzpy.core.constants import E_CHARGE, KB_EV, TWO_PI
from boltzpy.core.numerics import nb_onsager_from_tdos_flat, nb_transport_dos_flat
from boltzpy.interpolation import SKWInterpolator
from boltzpy.mesh import TetrahedronMesh
from boltzpy.parsers import DEFAULT_PARSER, ElectronicStructureParser, get_parser


def _validate_bound_parameter(name: str, value: float, *, allow_zero: bool = False) -> None:
    """Validate numeric parameter bounds for transport configuration."""
    if allow_zero:
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative")
        return
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")


class BoltzmannTransportCalculator:
    """Evaluate electrical and thermal transport tensors from an interpolator."""

    name = "boltzmann"

    def __init__(
        self,
        interpolator: SKWInterpolator,
        tau: float = 1e-14,
        chunk_size: int = 4096,
        energy_window_kbt_factor: float = 10.0,
        min_energy_window: float = 0.5,
        energy_step_kbt_divisor: float = 10.0,
        min_energy_step: float = 1e-4,
        low_temp_kbt_threshold: float = 1e-10,
        low_temp_energy_window: float = 0.5,
        low_temp_energy_step: float = 1e-3,
    ) -> None:
        """Initialize transport calculator.

        Args:
            interpolator: Interpolator used to evaluate energies and velocities.
            tau: Constant relaxation time in seconds.
            chunk_size: k-point chunk used in batched evaluations.
            energy_window_kbt_factor: Thermal-window multiplier for the integration range.
            min_energy_window: Lower bound for integration half-window in eV.
            energy_step_kbt_divisor: Thermal scaling divisor for energy-grid spacing.
            min_energy_step: Lower bound for energy-grid spacing in eV.
            low_temp_kbt_threshold: ``kBT`` threshold where low-temperature fallback is used.
            low_temp_energy_window: Integration half-window used in low-temperature fallback.
            low_temp_energy_step: Energy-grid spacing used in low-temperature fallback.

        """
        self.interp = interpolator
        self.tau = tau
        self.chunk_size = chunk_size
        self.tetra_mesh: TetrahedronMesh | None = None
        _validate_bound_parameter("energy_window_kbt_factor", energy_window_kbt_factor)
        _validate_bound_parameter("min_energy_window", min_energy_window)
        _validate_bound_parameter("energy_step_kbt_divisor", energy_step_kbt_divisor)
        _validate_bound_parameter("min_energy_step", min_energy_step)
        _validate_bound_parameter("low_temp_kbt_threshold", low_temp_kbt_threshold, allow_zero=True)
        _validate_bound_parameter("low_temp_energy_window", low_temp_energy_window)
        _validate_bound_parameter("low_temp_energy_step", low_temp_energy_step)
        self.energy_window_kbt_factor = energy_window_kbt_factor
        self.min_energy_window = min_energy_window
        self.energy_step_kbt_divisor = energy_step_kbt_divisor
        self.min_energy_step = min_energy_step
        self.low_temp_kbt_threshold = low_temp_kbt_threshold
        self.low_temp_energy_window = low_temp_energy_window
        self.low_temp_energy_step = low_temp_energy_step

    def _ensure_tetra_mesh(self, kpoint_mesh: tuple[int, int, int]) -> None:
        """Create or update tetrahedron mesh for current integration grid.

        Args:
            kpoint_mesh: Integration mesh dimensions.

        """
        mesh_array = np.array(kpoint_mesh, dtype=np.int32)
        if self.tetra_mesh is None or not np.array_equal(self.tetra_mesh.mesh, mesh_array):
            self.tetra_mesh = TetrahedronMesh(self.interp._lat, kpoint_mesh)

    def calculate_onsager_coefficients(  # noqa:
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Compute Onsager tensors $L_0$, $L_1$, and $L_2$.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            tuple of 3x3 Onsager tensors ``(l0, l1, l2)``.

        """
        if kchunk is None:
            kchunk = self.chunk_size

        self._ensure_tetra_mesh(kpoint_mesh)
        mesh = self.tetra_mesh
        if mesh is None:
            raise RuntimeError("tetrahedron mesh was not initialized")
        kbt = KB_EV * temperature

        if kbt > self.low_temp_kbt_threshold:
            # The grid is centered on the Fermi level because only that window contributes
            # materially to transport; tying spacing to kBT keeps the integration sharp enough
            # as temperature changes without exploding the grid at room temperature.
            e_window = max(self.energy_window_kbt_factor * kbt, self.min_energy_window)
            de = max(kbt / self.energy_step_kbt_divisor, self.min_energy_step)
        else:
            e_window = self.low_temp_energy_window
            de = self.low_temp_energy_step

        ne = int(2.0 * e_window / de) + 1
        e_grid = np.linspace(fermi_level - e_window, fermi_level + e_window, ne, dtype=np.float64)

        l0_total = np.zeros(9, dtype=np.float64)
        l1_total = np.zeros(9, dtype=np.float64)
        l2_total = np.zeros(9, dtype=np.float64)

        for spin in range(self.interp.nspin):
            e_all, vel_all = self.interp.eval_energy_velocity(mesh.full_kpoints, spin, kchunk)
            tdos = nb_transport_dos_flat(e_all, vel_all, mesh.tetrahedra, self.tau, e_grid)
            l0, l1, l2 = nb_onsager_from_tdos_flat(tdos, e_grid, fermi_level, kbt)
            l0_total += l0
            l1_total += l1
            l2_total += l2

        # A non-spin-polarized calculation stores one spin channel, but the physical transport
        # still carries spin degeneracy two, so we restore that factor here instead of in parsers.
        if self.interp.nspin == 1:
            l0_total *= 2.0
            l1_total *= 2.0
            l2_total *= 2.0

        norm = mesh.tetra_vol / (TWO_PI**3) * 1e30
        l0_total *= norm
        l1_total *= norm
        l2_total *= norm

        return l0_total.reshape(3, 3), l1_total.reshape(3, 3), l2_total.reshape(3, 3)

    def calculate_transport(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Return transport tensors for the given state point.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            tuple ``(sigma, seebeck, kappa, l0, l1, l2)``.

        """
        l0, l1, l2 = self.calculate_onsager_coefficients(fermi_level, temperature, kpoint_mesh, kchunk=kchunk)

        l0_inv = np.linalg.inv(l0)
        sigma = E_CHARGE * l0
        seebeck = -(l1 @ l0_inv) / temperature
        kappa = (E_CHARGE / temperature) * (l2 - l1 @ l0_inv @ l1)

        return sigma, seebeck, kappa, l0, l1, l2

    def calculate_conductivity(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """Return electrical conductivity tensor.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            Electrical conductivity tensor.

        """
        sigma, _, _, _, _, _ = self.calculate_transport(fermi_level, temperature, kpoint_mesh, kchunk=kchunk)
        return sigma

    def calculate_seebeck(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """Return Seebeck tensor.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            Seebeck tensor.

        """
        _, seebeck, _, _, _, _ = self.calculate_transport(fermi_level, temperature, kpoint_mesh, kchunk=kchunk)
        return seebeck

    def calculate_thermal_conductivity(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """Return electronic thermal conductivity tensor.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            Electronic thermal conductivity tensor.

        """
        _, _, kappa, _, _, _ = self.calculate_transport(fermi_level, temperature, kpoint_mesh, kchunk=kchunk)
        return kappa

    def calculate_zt(
        self,
        fermi_level: float,
        temperature: float,
        kpoint_mesh: tuple[int, int, int],
        lattice_thermal_conductivity: float = 0.0,
        kchunk: int | None = None,
    ) -> float:
        """Return isotropic thermoelectric figure of merit $ZT$.

        Args:
            fermi_level: Fermi level in eV.
            temperature: Temperature in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            lattice_thermal_conductivity: Lattice thermal conductivity in W/(m*K).
            kchunk: Optional chunk size override.

        Returns:
            Scalar isotropic figure of merit.

        """
        sigma, seebeck, kappa_el, _, _, _ = self.calculate_transport(
            fermi_level,
            temperature,
            kpoint_mesh,
            kchunk=kchunk,
        )

        sigma_avg = tensor_average(sigma)
        seebeck_avg = tensor_average(seebeck)
        kappa_el_avg = tensor_average(kappa_el)
        power_factor = seebeck_avg * seebeck_avg * sigma_avg
        return power_factor * temperature / (kappa_el_avg + lattice_thermal_conductivity)

    def _build_unified_energy_grid(
        self,
        mu_values: npt.NDArray[np.float64],
        temperatures: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Build a single energy grid covering all (mu, T) combinations.

        The grid is fine enough for the tightest thermal resolution (lowest T)
        and wide enough for the broadest window (highest T) shifted to every
        chemical potential value.
        """
        kbts = KB_EV * temperatures
        normal_mask = kbts > self.low_temp_kbt_threshold

        if np.any(normal_mask):
            normal_kbts = kbts[normal_mask]
            max_window = max(
                float(np.max(self.energy_window_kbt_factor * normal_kbts)),
                self.min_energy_window,
            )
            min_de = max(
                float(np.min(normal_kbts / self.energy_step_kbt_divisor)),
                self.min_energy_step,
            )
        else:
            max_window = self.low_temp_energy_window
            min_de = self.low_temp_energy_step

        if np.any(~normal_mask):
            max_window = max(max_window, self.low_temp_energy_window)
            min_de = min(min_de, self.low_temp_energy_step)

        e_lo = float(np.min(mu_values)) - max_window
        e_hi = float(np.max(mu_values)) + max_window

        ne = int((e_hi - e_lo) / min_de) + 1
        return np.linspace(e_lo, e_hi, ne, dtype=np.float64)

    def calculate_transport_scan(
        self,
        mu_values: npt.NDArray[np.float64],
        temperatures: npt.NDArray[np.float64],
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> dict[str, Any]:
        """Compute transport tensors for all (mu, T) pairs efficiently.

        Interpolation and transport DOS are computed **once**; only the cheap
        Onsager integration is repeated per (mu, T) pair.

        Args:
            mu_values: Absolute chemical potential values in eV (not shifts).
            temperatures: Temperature array in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            Dictionary with keys ``"sigma"``, ``"seebeck"``, ``"kappa"``
            (shape ``(n_mu, n_T, 3, 3)``) and ``"*_avg"`` (shape ``(n_mu, n_T)``).

        """
        if kchunk is None:
            kchunk = self.chunk_size

        self._ensure_tetra_mesh(kpoint_mesh)
        mesh = self.tetra_mesh
        if mesh is None:
            raise RuntimeError("tetrahedron mesh was not initialized")

        e_grid = self._build_unified_energy_grid(mu_values, temperatures)

        # compute TDOS once for all spins
        tdos_total = np.zeros((e_grid.shape[0], 9), dtype=np.float64)
        for spin in range(self.interp.nspin):
            e_all, vel_all = self.interp.eval_energy_velocity(
                mesh.full_kpoints,
                spin,
                kchunk,
            )
            tdos_total += nb_transport_dos_flat(
                e_all,
                vel_all,
                mesh.tetrahedra,
                self.tau,
                e_grid,
            )

        if self.interp.nspin == 1:
            tdos_total *= 2.0

        norm = mesh.tetra_vol / (TWO_PI**3) * 1e30

        n_mu = mu_values.shape[0]
        n_t = temperatures.shape[0]
        sigma_all = np.empty((n_mu, n_t, 3, 3), dtype=np.float64)
        seebeck_all = np.empty((n_mu, n_t, 3, 3), dtype=np.float64)
        kappa_all = np.empty((n_mu, n_t, 3, 3), dtype=np.float64)
        sigma_avg = np.empty((n_mu, n_t), dtype=np.float64)
        seebeck_avg = np.empty((n_mu, n_t), dtype=np.float64)
        kappa_avg = np.empty((n_mu, n_t), dtype=np.float64)

        for imu in range(n_mu):
            mu = float(mu_values[imu])
            for it in range(n_t):
                temp = float(temperatures[it])
                kbt = KB_EV * temp

                l0, l1, l2 = nb_onsager_from_tdos_flat(tdos_total, e_grid, mu, kbt)
                l0 = l0 * norm
                l1 = l1 * norm
                l2 = l2 * norm

                l0_mat = l0.reshape(3, 3)
                l1_mat = l1.reshape(3, 3)
                l2_mat = l2.reshape(3, 3)

                l0_inv = np.linalg.inv(l0_mat)
                sig = E_CHARGE * l0_mat
                see = -(l1_mat @ l0_inv) / temp
                kap = (E_CHARGE / temp) * (l2_mat - l1_mat @ l0_inv @ l1_mat)

                sigma_all[imu, it] = sig
                seebeck_all[imu, it] = see
                kappa_all[imu, it] = kap
                sigma_avg[imu, it] = tensor_average(sig)
                seebeck_avg[imu, it] = tensor_average(see)
                kappa_avg[imu, it] = tensor_average(kap)

        return {
            "sigma": sigma_all,
            "sigma_avg": sigma_avg,
            "seebeck": seebeck_all,
            "seebeck_avg": seebeck_avg,
            "kappa": kappa_all,
            "kappa_avg": kappa_avg,
        }


register_calculator(BoltzmannTransportCalculator.name, BoltzmannTransportCalculator)


def _prepare_temperature_array(temperature: float | npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Validate and normalize temperature input to a 1-D array."""
    temp_input = np.asarray(temperature, dtype=np.float64)
    temperatures = np.array([float(temp_input)], dtype=np.float64) if temp_input.ndim == 0 else np.ravel(temp_input)
    if temperatures.size == 0:
        raise ValueError("temperature array must not be empty")
    if np.any(temperatures <= 0.0):
        raise ValueError("temperature must be positive")
    return temperatures


def _prepare_mu_shifts(chemical_potential: float | npt.ArrayLike | None) -> npt.NDArray[np.float64]:
    """Validate and normalize chemical potential input to a 1-D array."""
    mu_input = np.asarray(chemical_potential, dtype=np.float64)
    mu_shifts = np.array([float(mu_input)], dtype=np.float64) if mu_input.ndim == 0 else np.ravel(mu_input)
    if mu_shifts.size == 0:
        raise ValueError("chemical_potential array must not be empty")
    return mu_shifts


def _compute_mu_scan(
    calc: BoltzmannTransportCalculator,
    fermi: float,
    mu_shifts: npt.NDArray[np.float64],
    temperatures: npt.NDArray[np.float64],
    kpoint_mesh: tuple[int, int, int],
    chunk_size: int,
    metadata: dict[str, Any],
) -> dict[float | str, Any]:
    """Sweep chemical potentials and temperatures, returning nested dict."""
    mu_abs = fermi + mu_shifts
    scan = calc.calculate_transport_scan(mu_abs, temperatures, kpoint_mesh, kchunk=chunk_size)

    results: dict[float | str, Any] = {}
    for imu, dmu in enumerate(mu_shifts):
        mu_key = float(dmu)
        t_dict: dict[float, dict[str, Any]] = {}
        for it, temp in enumerate(temperatures):
            t_dict[float(temp)] = {
                "sigma": scan["sigma"][imu, it],
                "sigma_avg": float(scan["sigma_avg"][imu, it]),
                "seebeck": scan["seebeck"][imu, it],
                "seebeck_avg": float(scan["seebeck_avg"][imu, it]),
                "kappa": scan["kappa"][imu, it],
                "kappa_avg": float(scan["kappa_avg"][imu, it]),
            }
        results[mu_key] = t_dict
    results["meta"] = metadata
    return results


def calculate_spin_polarized_transport(
    filepath: str | Path,
    temperature: float | npt.ArrayLike = 300.0,
    chemical_potential: float | npt.ArrayLike | None = None,
    tau: float = 1e-14,
    kpoint_mesh: tuple[int, int, int] = (20, 20, 20),
    lr_ratio: int = 5,
    band_window: tuple[float, float] | None = None,
    chunk_size: int = 4096,
    energy_window_kbt_factor: float = 10.0,
    min_energy_window: float = 0.5,
    energy_step_kbt_divisor: float = 10.0,
    min_energy_step: float = 1e-4,
    low_temp_kbt_threshold: float = 1e-10,
    low_temp_energy_window: float = 0.5,
    low_temp_energy_step: float = 1e-3,
    parser: str | ElectronicStructureParser = DEFAULT_PARSER,
    calculator: str = "boltzmann",
) -> dict[str, Any] | dict[float | str, Any]:
    """Run full parser -> interpolation -> transport workflow.

    Args:
        filepath: Path to electronic-structure file.
        temperature: Temperature in kelvin, scalar or array-like.
        chemical_potential: Chemical potential shift(s) relative to the Fermi
            level in eV.  Accepts a scalar or array-like.  When provided the
            result is a nested dictionary keyed by chemical potential and
            temperature values:
            ``{mu: {T: {"sigma": ..., "seebeck": ..., "kappa": ..., ...}}}``.
            When *None* (the default) the original flat-dict format is returned.
        tau: Constant relaxation time in seconds.
        kpoint_mesh: Integration mesh dimensions.
        lr_ratio: Interpolator star-vector ratio.
        band_window: Optional relative energy window around Fermi level.
        chunk_size: k-point chunk size for batched evaluation.
        energy_window_kbt_factor: Thermal-window multiplier for integration range.
        min_energy_window: Lower bound for integration half-window in eV.
        energy_step_kbt_divisor: Thermal scaling divisor for energy-grid spacing.
        min_energy_step: Lower bound for energy-grid spacing in eV.
        low_temp_kbt_threshold: ``kBT`` threshold where low-temperature fallback is used.
        low_temp_energy_window: Integration half-window used in low-temperature fallback.
        low_temp_energy_step: Energy-grid spacing used in low-temperature fallback.
        parser: Parser name or parser instance.
        calculator: Registered calculator name.

    Returns:
        When *chemical_potential* is *None*: dictionary with tensors and metadata
        in a uniform shape contract (``temperature`` has shape ``(nT,)``,
        ``sigma``/``seebeck``/``kappa`` have shape ``(nT, 3, 3)``, and ``*_avg``
        have shape ``(nT,)``).

        When *chemical_potential* is given: nested dictionary
        ``{mu: {T: {"sigma": 3x3, "seebeck": 3x3, "kappa": 3x3,
        "sigma_avg": float, "seebeck_avg": float, "kappa_avg": float}}}``
        with an extra ``"meta"`` key containing Fermi energy, jspins, parser and
        calculator info.

    Raises:
        ValueError: If requested band window does not include any bands.

    """
    parser_obj = get_parser(parser) if isinstance(parser, str) else parser

    calculator_cls = get_calculator(calculator)

    filepath = str(filepath)
    parsed = parser_obj.parse(filepath)
    fermi = parsed.fermi_energy
    eigenvalues = parsed.eigenvalues

    # for debugging
    # print(f"Parsing {Path(filepath).name} with parser '{parser_obj.name}'...")
    # print(f"  jspins={parsed.jspins}, nk={parsed.nk}, nbands={parsed.nbands}")
    # print(f"  E_Fermi = {fermi:.4f} eV")
    # print(f"  det(A) = {np.linalg.det(parsed.lattice):.2f} A^3")
    # print(f"  symops = {len(parsed.symops)}")

    if band_window is not None:
        emin, emax = band_window
        emin_abs = fermi + emin
        emax_abs = fermi + emax
        all_bands = eigenvalues.reshape(-1, parsed.nbands)
        band_min = all_bands.min(axis=0)
        band_max = all_bands.max(axis=0)
        # Filtering bands before interpolation keeps the linear solve smaller and avoids spending
        # most of the runtime on states that can never contribute near the chosen chemical window.
        mask = (band_max > emin_abs) & (band_min < emax_abs)
        indices = np.where(mask)[0]
        if len(indices) == 0:
            raise ValueError(f"No bands found in window [{emin}, {emax}] eV")
        b_lo = int(indices[0])
        b_hi = int(indices[-1]) + 1
        eigenvalues = eigenvalues[:, :, b_lo:b_hi]
        print(f"  band window [{emin}, {emax}] eV -> bands {b_lo}..{b_hi - 1} ({b_hi - b_lo})")

    print(f"SKW interpolation (lr_ratio={lr_ratio})...")
    t_skw = time.time()
    interp = SKWInterpolator(
        kpoints=parsed.kpoints,
        eigenvalues=eigenvalues,
        cell=parsed.lattice,
        symops=parsed.symops,
        time_reversal=True,
        lr_ratio=lr_ratio,
    )
    print(f"  MAE = {interp.mae:.6f} eV")
    print(f"  NR = {interp.nr}, NPG = {interp._npg}, dt = {time.time() - t_skw:.2f} s")

    temp_input = np.asarray(temperature, dtype=np.float64)
    is_scalar_temp = temp_input.ndim == 0
    temperatures = np.array([float(temp_input)], dtype=np.float64) if is_scalar_temp else np.ravel(temp_input)
    if temperatures.size == 0:
        raise ValueError("temperature array must not be empty")
    if np.any(temperatures <= 0.0):
        raise ValueError("temperature must be positive")

    # --- chemical potential grid -------------------------------------------
    use_mu_scan = chemical_potential is not None
    mu_shifts = _prepare_mu_shifts(chemical_potential) if use_mu_scan else np.array([0.0], dtype=np.float64)

    print(
        f"Transport (mesh={kpoint_mesh}, T={temperature} K, tau={tau:.1e} s, "
        f"chunk={chunk_size}, calculator='{calculator}')..."
    )
    if use_mu_scan:
        print(f"  Chemical potential shifts: {mu_shifts} eV relative to E_Fermi={fermi:.4f} eV")

    calc = calculator_cls(
        interp,
        tau=tau,
        chunk_size=chunk_size,
        energy_window_kbt_factor=energy_window_kbt_factor,
        min_energy_window=min_energy_window,
        energy_step_kbt_divisor=energy_step_kbt_divisor,
        min_energy_step=min_energy_step,
        low_temp_kbt_threshold=low_temp_kbt_threshold,
        low_temp_energy_window=low_temp_energy_window,
        low_temp_energy_step=low_temp_energy_step,
    )

    metadata = {
        "fermi_energy": fermi,
        "jspins": parsed.jspins,
        "parser": parser_obj.name,
        "calculator": calculator,
    }

    if use_mu_scan:
        return _compute_mu_scan(calc, fermi, mu_shifts, temperatures, kpoint_mesh, chunk_size, metadata)

    fermi_arr = np.array([fermi], dtype=np.float64)
    scan = calc.calculate_transport_scan(fermi_arr, temperatures, kpoint_mesh, kchunk=chunk_size)

    results = {
        "temperature": temperatures,
        "sigma": scan["sigma"][0],
        "sigma_avg": scan["sigma_avg"][0],
        "seebeck": scan["seebeck"][0],
        "seebeck_avg": scan["seebeck_avg"][0],
        "kappa": scan["kappa"][0],
        "kappa_avg": scan["kappa_avg"][0],
        **metadata,
    }

    return results
