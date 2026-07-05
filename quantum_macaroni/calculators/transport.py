"""Boltzmann transport calculators and high-level orchestration entry point."""

import hashlib
import time
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from quantum_macaroni.calculators.base import TransportTuple, get_calculator, register_calculator, tensor_average
from quantum_macaroni.checkpointing import (
    BaseSystemState,
    Checkpoint,
    CheckpointFormatError,
    CheckpointManager,
    ExecutionProgress,
    RuntimeParameters,
    stable_fingerprint,
    to_json_value,
)
from quantum_macaroni.checkpointing.state import ArrayMap, Metadata
from quantum_macaroni.core.constants import ANG3_TO_M3, E_CHARGE, KB_EV, TWO_PI
from quantum_macaroni.core.numerics import nb_onsager_from_tdos_flat, nb_transport_dos_flat
from quantum_macaroni.interpolation import SKWInterpolator
from quantum_macaroni.mesh import TetrahedronMesh
from quantum_macaroni.parsers import DEFAULT_PARSER, ElectronicStructureParser, ParserResult, get_parser

TRANSPORT_CHECKPOINT_WORKFLOW = "spin-polarized-transport"
FILE_DIGEST_CHUNK_SIZE = 1024 * 1024


def _validate_bound_parameter(name: str, value: float, *, allow_zero: bool = False) -> None:
    """Validate numeric parameter bounds for transport configuration."""
    if allow_zero:
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative")
        return
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")


@dataclass
class EnergyGridDefaults:
    """Shared default parameters for energy-grid construction."""

    energy_window_kbt_factor: float = 10.0
    min_energy_window: float = 0.5
    energy_step_kbt_divisor: float = 10.0
    min_energy_step: float = 1e-4
    low_temp_kbt_threshold: float = 1e-10
    low_temp_energy_window: float = 0.5
    low_temp_energy_step: float = 1e-3


DEFAULTS = EnergyGridDefaults()


class TransportWorkflowStage(IntEnum):
    """Ordered checkpoint stages for the spin-polarized transport workflow."""

    PARSED = 1
    INTERPOLATED = 2
    TRANSPORT_DOS = 3
    COMPLETED = 4


@dataclass(slots=True)
class TransportDOS:
    """Transport density-of-states payload reused across scan integrations.

    Attributes:
        energy_grid: Energy grid in eV.
        values: Flattened transport DOS with shape ``(ne, 9)``.
        norm: Brillouin-zone normalization factor applied during integration.
        kpoint_mesh: Integration mesh dimensions used to build the DOS.

    """

    energy_grid: npt.NDArray[np.float64]
    values: npt.NDArray[np.float64]
    norm: float
    kpoint_mesh: tuple[int, int, int]


def _onsager_to_transport(
    l0: npt.NDArray[np.float64],
    l1: npt.NDArray[np.float64],
    l2: npt.NDArray[np.float64],
    temperature: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Convert Onsager coefficients to physical transport tensors.

    Args:
        l0: 3x3 Onsager L0 tensor.
        l1: 3x3 Onsager L1 tensor.
        l2: 3x3 Onsager L2 tensor.
        temperature: Temperature in kelvin.

    Returns:
        Tuple ``(sigma, seebeck, kappa)`` of 3x3 transport tensors.

    """
    l0_inv = np.linalg.inv(l0)
    # The moment integrals use eV as their energy unit. Only one factor
    # of the elementary charge is needed for sigma and kappa, and no charge
    # factor is needed for Seebeck because eV per elementary charge is a volt.
    sigma = E_CHARGE * l0
    seebeck = -(l0_inv @ l1) / temperature
    kappa = (E_CHARGE / temperature) * (l2 - l1 @ l0_inv @ l1)
    return sigma, seebeck, kappa


class BoltzmannTransportCalculator:
    """Evaluate electrical and thermal transport tensors from an interpolator."""

    name = "boltzmann"

    def __init__(
        self,
        interpolator: SKWInterpolator,
        tau: float = 1e-14,
        chunk_size: int = 4096,
        energy_window_kbt_factor: float = DEFAULTS.energy_window_kbt_factor,
        min_energy_window: float = DEFAULTS.min_energy_window,
        energy_step_kbt_divisor: float = DEFAULTS.energy_step_kbt_divisor,
        min_energy_step: float = DEFAULTS.min_energy_step,
        low_temp_kbt_threshold: float = DEFAULTS.low_temp_kbt_threshold,
        low_temp_energy_window: float = DEFAULTS.low_temp_energy_window,
        low_temp_energy_step: float = DEFAULTS.low_temp_energy_step,
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
        _validate_bound_parameter("tau", tau)
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
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

    def _spin_degeneracy(self) -> float:
        """Return spin degeneracy factor (2 for non-spin-polarized, 1 otherwise)."""
        return 2.0 if self.interp.nspin == 1 else 1.0

    @staticmethod
    def _bz_norm(mesh: TetrahedronMesh) -> float:
        """Return Boltzmann normalization factor converting to SI units.

        Args:
            mesh: Tetrahedron mesh providing BZ volume.

        Returns:
            Normalization factor with units converting Å⁻³ to m⁻³.

        """
        return mesh.tetra_vol / (TWO_PI**3) * ANG3_TO_M3

    def _ensure_tetra_mesh(self, kpoint_mesh: tuple[int, int, int]) -> None:
        """Create or update tetrahedron mesh for current integration grid.

        Args:
            kpoint_mesh: Integration mesh dimensions.

        """
        mesh_array = np.array(kpoint_mesh, dtype=np.int32)
        if self.tetra_mesh is None or not np.array_equal(self.tetra_mesh.mesh, mesh_array):
            self.tetra_mesh = TetrahedronMesh(self.interp._lat, kpoint_mesh)

    def calculate_onsager_coefficients(  # noqa: PLR0913
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

        Raises:
            RuntimeError: If tetrahedron mesh was not initialized.

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
        spin_factor = self._spin_degeneracy()
        l0_total *= spin_factor
        l1_total *= spin_factor
        l2_total *= spin_factor

        norm = self._bz_norm(mesh)
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
    ) -> TransportTuple:
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

        sigma, seebeck, kappa = _onsager_to_transport(l0, l1, l2, temperature)

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

        Args:
            mu_values: Absolute chemical potential values in eV.
            temperatures: Temperature array in kelvin.

        Returns:
            Unified energy grid in eV.

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

    def build_transport_dos(
        self,
        mu_values: npt.NDArray[np.float64],
        temperatures: npt.NDArray[np.float64],
        kpoint_mesh: tuple[int, int, int],
        kchunk: int | None = None,
    ) -> TransportDOS:
        """Build transport DOS once for a chemical-potential/temperature scan.

        This is the expensive restart boundary: interpolation is evaluated over
        the full k-point mesh and converted into a transport density of states.

        Args:
            mu_values: Absolute chemical potential values in eV (not shifts).
            temperatures: Temperature array in kelvin.
            kpoint_mesh: Integration mesh dimensions.
            kchunk: Optional chunk size override.

        Returns:
            Transport DOS payload ready for repeated scan integrations.

        Raises:
            RuntimeError: If tetrahedron mesh was not initialized.

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

        spin_factor = self._spin_degeneracy()
        tdos_total *= spin_factor

        return TransportDOS(
            energy_grid=e_grid,
            values=tdos_total,
            norm=self._bz_norm(mesh),
            kpoint_mesh=kpoint_mesh,
        )

    def integrate_transport_dos_scan(
        self,
        transport_dos: TransportDOS,
        mu_values: npt.NDArray[np.float64],
        temperatures: npt.NDArray[np.float64],
    ) -> dict[str, Any]:
        """Integrate a transport DOS over all ``(mu, T)`` scan points.

        Args:
            transport_dos: Precomputed transport DOS payload.
            mu_values: Absolute chemical potential values in eV.
            temperatures: Temperature array in kelvin.

        Returns:
            Dictionary with keys ``"sigma"``, ``"seebeck"``, ``"kappa"``
            (shape ``(n_mu, n_T, 3, 3)``) and ``"*_avg"`` (shape ``(n_mu, n_T)``).

        """
        e_grid = transport_dos.energy_grid
        tdos_total = transport_dos.values
        norm = transport_dos.norm

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
                l0 *= norm
                l1 *= norm
                l2 *= norm

                l0_mat = l0.reshape(3, 3)
                l1_mat = l1.reshape(3, 3)
                l2_mat = l2.reshape(3, 3)

                sig, see, kap = _onsager_to_transport(l0_mat, l1_mat, l2_mat, temp)

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
        transport_dos = self.build_transport_dos(mu_values, temperatures, kpoint_mesh, kchunk=kchunk)
        return self.integrate_transport_dos_scan(transport_dos, mu_values, temperatures)


register_calculator(BoltzmannTransportCalculator.name, BoltzmannTransportCalculator)


def _transport_result_units() -> dict[str, str]:
    """Return units for transport output fields.

    Returns:
        Dictionary mapping field names to their physical units.

    """
    return {
        "temperature": "K",
        "chemical_potential_shift": "eV",
        "fermi_energy": "eV",
        "sigma": "S/m",
        "sigma_avg": "S/m",
        "seebeck": "V/K",
        "seebeck_avg": "V/K",
        "kappa": "W/(m*K)",
        "kappa_avg": "W/(m*K)",
    }


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


def _as_metadata(values: dict[str, Any]) -> Metadata:
    """Normalize a dictionary into checkpoint JSON metadata."""
    metadata = to_json_value(values)
    if not isinstance(metadata, dict):
        raise TypeError("checkpoint metadata must normalize to a dictionary")
    return metadata


def _file_digest(filepath: str | Path) -> str:
    """Return the SHA-256 digest of an input file."""
    digest = hashlib.sha256()
    with Path(filepath).open("rb") as file_obj:
        while True:
            chunk = file_obj.read(FILE_DIGEST_CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _energy_grid_settings(
    energy_window_kbt_factor: float,
    min_energy_window: float,
    energy_step_kbt_divisor: float,
    min_energy_step: float,
    low_temp_kbt_threshold: float,
    low_temp_energy_window: float,
    low_temp_energy_step: float,
) -> dict[str, float]:
    """Return checkpoint-compatible energy-grid settings."""
    return {
        "energy_window_kbt_factor": energy_window_kbt_factor,
        "min_energy_window": min_energy_window,
        "energy_step_kbt_divisor": energy_step_kbt_divisor,
        "min_energy_step": min_energy_step,
        "low_temp_kbt_threshold": low_temp_kbt_threshold,
        "low_temp_energy_window": low_temp_energy_window,
        "low_temp_energy_step": low_temp_energy_step,
    }


def _transport_workflow_fingerprints(
    filepath: str | Path,
    parser_name: str,
    calculator: str,
    temperatures: npt.NDArray[np.float64],
    mu_shifts: npt.NDArray[np.float64],
    tau: float,
    kpoint_mesh: tuple[int, int, int],
    lr_ratio: int,
    band_window: tuple[float, float] | None,
    energy_settings: dict[str, float],
    use_mu_scan: bool,
) -> dict[str, str]:
    """Return staged fingerprints for checkpoint compatibility."""
    resolved_path = Path(filepath).resolve()
    input_fingerprint = stable_fingerprint(
        {
            "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
            "stage": "input",
            "filepath": str(resolved_path),
            "file_sha256": _file_digest(resolved_path),
            "parser": parser_name,
        }
    )
    interpolation_fingerprint = stable_fingerprint(
        {
            "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
            "stage": "interpolation",
            "input": input_fingerprint,
            "lr_ratio": lr_ratio,
            "band_window": band_window,
        }
    )
    tdos_fingerprint = stable_fingerprint(
        {
            "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
            "stage": "transport_dos",
            "interpolation": interpolation_fingerprint,
            "calculator": calculator,
            "tau": tau,
            "kpoint_mesh": kpoint_mesh,
            "temperatures": temperatures,
            "mu_shifts": mu_shifts,
            "energy_settings": energy_settings,
        }
    )
    result_fingerprint = stable_fingerprint(
        {
            "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
            "stage": "completed",
            "transport_dos": tdos_fingerprint,
            "use_mu_scan": use_mu_scan,
        }
    )
    return {
        "input": input_fingerprint,
        "interpolation": interpolation_fingerprint,
        "transport_dos": tdos_fingerprint,
        "result": result_fingerprint,
    }


def _checkpoint_reached(checkpoint: Checkpoint, stage: TransportWorkflowStage) -> bool:
    """Return whether a checkpoint has reached a transport workflow stage."""
    return checkpoint.progress.step >= int(stage)


def _checkpoint_fingerprint_matches(checkpoint: Checkpoint, name: str, expected: str) -> bool:
    """Return whether a loaded checkpoint carries an expected staged fingerprint."""
    fingerprints = checkpoint.runtime_parameters.values.get("fingerprints")
    return isinstance(fingerprints, dict) and fingerprints.get(name) == expected


def _parser_result_arrays(parsed: ParserResult) -> ArrayMap:
    """Return checkpoint arrays for a parser result."""
    return {
        "parsed.kpoints": parsed.kpoints,
        "parsed.eigenvalues": parsed.eigenvalues,
        "parsed.lattice": parsed.lattice,
        "parsed.symops": parsed.symops,
    }


def _parser_result_metadata(parsed: ParserResult) -> dict[str, Any]:
    """Return checkpoint metadata for a parser result."""
    return {
        "fermi_energy": parsed.fermi_energy,
        "jspins": parsed.jspins,
        "nbands": parsed.nbands,
        "nk": parsed.nk,
    }


def _parser_result_from_checkpoint(checkpoint: Checkpoint) -> ParserResult:
    """Restore a parser result from checkpoint base state."""
    parsed_metadata = checkpoint.base_system.metadata.get("parsed")
    if not isinstance(parsed_metadata, dict):
        raise CheckpointFormatError("checkpoint is missing parsed metadata")

    arrays = checkpoint.base_system.arrays
    return ParserResult(
        kpoints=np.ascontiguousarray(arrays["parsed.kpoints"], dtype=np.float64),
        eigenvalues=np.ascontiguousarray(arrays["parsed.eigenvalues"], dtype=np.float64),
        fermi_energy=float(parsed_metadata["fermi_energy"]),
        jspins=int(parsed_metadata["jspins"]),
        nbands=int(parsed_metadata["nbands"]),
        nk=int(parsed_metadata["nk"]),
        lattice=np.ascontiguousarray(arrays["parsed.lattice"], dtype=np.float64),
        symops=np.ascontiguousarray(arrays["parsed.symops"], dtype=int),
    )


def _prefix_arrays(prefix: str, arrays: dict[str, npt.NDArray[np.generic]]) -> ArrayMap:
    """Prefix checkpoint array names with a namespace."""
    return {f"{prefix}{name}": array for name, array in arrays.items()}


def _unprefix_arrays(prefix: str, arrays: ArrayMap) -> ArrayMap:
    """Return checkpoint arrays from one namespace without their prefix."""
    return {name.removeprefix(prefix): array for name, array in arrays.items() if name.startswith(prefix)}


def _interpolator_from_checkpoint(checkpoint: Checkpoint) -> SKWInterpolator:
    """Restore an interpolator from checkpoint base state."""
    metadata = checkpoint.base_system.metadata.get("interpolator")
    if not isinstance(metadata, dict):
        raise CheckpointFormatError("checkpoint is missing interpolator metadata")
    return SKWInterpolator.from_state(metadata, _unprefix_arrays("skw.", checkpoint.base_system.arrays))


def _transport_dos_arrays(transport_dos: TransportDOS) -> ArrayMap:
    """Return checkpoint arrays for a transport DOS payload."""
    return {
        "tdos.energy_grid": transport_dos.energy_grid,
        "tdos.values": transport_dos.values,
    }


def _transport_dos_metadata(transport_dos: TransportDOS) -> dict[str, Any]:
    """Return checkpoint metadata for a transport DOS payload."""
    return {
        "norm": transport_dos.norm,
        "kpoint_mesh": transport_dos.kpoint_mesh,
    }


def _transport_dos_from_checkpoint(checkpoint: Checkpoint) -> TransportDOS:
    """Restore transport DOS from checkpoint progress state."""
    metadata = checkpoint.progress.metadata.get("transport_dos")
    if not isinstance(metadata, dict):
        raise CheckpointFormatError("checkpoint is missing transport DOS metadata")
    mesh_raw = metadata.get("kpoint_mesh")
    if not isinstance(mesh_raw, list | tuple) or len(mesh_raw) != 3:  # noqa: PLR2004
        raise CheckpointFormatError("checkpoint transport DOS mesh must have three dimensions")
    arrays = checkpoint.progress.arrays
    return TransportDOS(
        energy_grid=np.ascontiguousarray(arrays["tdos.energy_grid"], dtype=np.float64),
        values=np.ascontiguousarray(arrays["tdos.values"], dtype=np.float64),
        norm=float(metadata["norm"]),
        kpoint_mesh=(int(mesh_raw[0]), int(mesh_raw[1]), int(mesh_raw[2])),
    )


def _scan_result_arrays(scan: dict[str, Any]) -> ArrayMap:
    """Return checkpoint arrays for transport scan results."""
    return {
        "result.sigma": scan["sigma"],
        "result.sigma_avg": scan["sigma_avg"],
        "result.seebeck": scan["seebeck"],
        "result.seebeck_avg": scan["seebeck_avg"],
        "result.kappa": scan["kappa"],
        "result.kappa_avg": scan["kappa_avg"],
    }


def _scan_result_from_checkpoint(checkpoint: Checkpoint) -> dict[str, Any]:
    """Restore scan result arrays from checkpoint progress state."""
    arrays = checkpoint.progress.arrays
    return {
        "sigma": np.ascontiguousarray(arrays["result.sigma"], dtype=np.float64),
        "sigma_avg": np.ascontiguousarray(arrays["result.sigma_avg"], dtype=np.float64),
        "seebeck": np.ascontiguousarray(arrays["result.seebeck"], dtype=np.float64),
        "seebeck_avg": np.ascontiguousarray(arrays["result.seebeck_avg"], dtype=np.float64),
        "kappa": np.ascontiguousarray(arrays["result.kappa"], dtype=np.float64),
        "kappa_avg": np.ascontiguousarray(arrays["result.kappa_avg"], dtype=np.float64),
    }


def _build_transport_checkpoint(
    stage: TransportWorkflowStage,
    fingerprints: dict[str, str],
    runtime_config: dict[str, Any],
    temperatures: npt.NDArray[np.float64],
    mu_shifts: npt.NDArray[np.float64],
    parsed: ParserResult | None = None,
    interpolator: SKWInterpolator | None = None,
    transport_dos: TransportDOS | None = None,
    scan: dict[str, Any] | None = None,
    result_metadata: dict[str, Any] | None = None,
) -> Checkpoint:
    """Build a generic checkpoint object for the transport workflow."""
    base_arrays: ArrayMap = {}
    base_metadata: dict[str, Any] = {
        "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
        "fingerprints": fingerprints,
    }
    if parsed is not None:
        base_arrays.update(_parser_result_arrays(parsed))
        base_metadata["parsed"] = _parser_result_metadata(parsed)
    if interpolator is not None:
        interpolator_metadata, interpolator_arrays = interpolator.to_state()
        base_arrays.update(_prefix_arrays("skw.", interpolator_arrays))
        base_metadata["interpolator"] = interpolator_metadata

    progress_arrays: ArrayMap = {}
    progress_metadata: dict[str, Any] = {
        "stage": stage.name.lower(),
        "fingerprints": fingerprints,
    }
    if transport_dos is not None:
        progress_arrays.update(_transport_dos_arrays(transport_dos))
        progress_metadata["transport_dos"] = _transport_dos_metadata(transport_dos)
    if scan is not None:
        progress_arrays.update(_scan_result_arrays(scan))
        progress_metadata["result"] = {
            "metadata": result_metadata or {},
        }

    runtime_parameters = RuntimeParameters(
        values=_as_metadata(
            {
                "workflow": TRANSPORT_CHECKPOINT_WORKFLOW,
                "fingerprints": fingerprints,
                "config": runtime_config,
            }
        ),
        arrays={
            "temperatures": temperatures,
            "mu_shifts": mu_shifts,
        },
    )
    return Checkpoint(
        base_system=BaseSystemState(
            version=fingerprints["interpolation"],
            arrays=base_arrays,
            metadata=_as_metadata(base_metadata),
        ),
        runtime_parameters=runtime_parameters,
        progress=ExecutionProgress(
            step=int(stage),
            arrays=progress_arrays,
            metadata=_as_metadata(progress_metadata),
        ),
    )


def _save_transport_checkpoint(
    manager: CheckpointManager | None,
    checkpoint_path: str | Path | None,
    stage: TransportWorkflowStage,
    fingerprints: dict[str, str],
    runtime_config: dict[str, Any],
    temperatures: npt.NDArray[np.float64],
    mu_shifts: npt.NDArray[np.float64],
    parsed: ParserResult | None = None,
    interpolator: SKWInterpolator | None = None,
    transport_dos: TransportDOS | None = None,
    scan: dict[str, Any] | None = None,
    result_metadata: dict[str, Any] | None = None,
) -> None:
    """Persist a transport checkpoint when checkpointing is enabled."""
    if manager is None or checkpoint_path is None:
        return
    checkpoint = _build_transport_checkpoint(
        stage,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
        parsed=parsed,
        interpolator=interpolator,
        transport_dos=transport_dos,
        scan=scan,
        result_metadata=result_metadata,
    )
    manager.save_checkpoint(checkpoint, checkpoint_path)


def _load_transport_checkpoint(
    manager: CheckpointManager | None,
    checkpoint_path: str | Path | None,
    resume_checkpoint: bool,
) -> Checkpoint | None:
    """Load a transport checkpoint when resume is enabled and the file exists."""
    if manager is None or checkpoint_path is None or not resume_checkpoint:
        return None
    path = Path(checkpoint_path)
    if not path.exists():
        return None
    checkpoint = manager.load_checkpoint(path)
    workflow = checkpoint.runtime_parameters.values.get("workflow")
    if workflow != TRANSPORT_CHECKPOINT_WORKFLOW:
        raise CheckpointFormatError(f"checkpoint workflow {workflow!r} is not supported here")
    return checkpoint


def _format_temperature_scan(
    temperatures: npt.NDArray[np.float64],
    scan: dict[str, Any],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    """Format one-Fermi-level scan arrays into the flat public result contract."""
    return {
        "temperature": temperatures,  # Temperature grid in K, shape: (nT,)
        "sigma": scan["sigma"][0],  # Electrical conductivity tensor, shape: (nT, 3, 3)
        "sigma_avg": scan["sigma_avg"][0],  # Isotropic conductivity average, shape: (nT,)
        "seebeck": scan["seebeck"][0],  # Seebeck tensor, shape: (nT, 3, 3)
        "seebeck_avg": scan["seebeck_avg"][0],  # Isotropic Seebeck average, shape: (nT,)
        "kappa": scan["kappa"][0],  # Electronic thermal conductivity tensor, shape: (nT, 3, 3)
        "kappa_avg": scan["kappa_avg"][0],  # Isotropic thermal conductivity average, shape: (nT,)
        "units": _transport_result_units(),  # Unit map for transport outputs
        **metadata,
    }


def _format_mu_scan(
    mu_shifts: npt.NDArray[np.float64],
    temperatures: npt.NDArray[np.float64],
    scan: dict[str, Any],
    metadata: dict[str, Any],
) -> dict[float | str, Any]:
    """Format scan arrays into the nested chemical-potential result contract."""
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
    results["meta"] = {
        **metadata,
        "units": _transport_result_units(),
    }
    return results


def _compute_mu_scan(
    calc: BoltzmannTransportCalculator,
    fermi: float,
    mu_shifts: npt.NDArray[np.float64],
    temperatures: npt.NDArray[np.float64],
    kpoint_mesh: tuple[int, int, int],
    chunk_size: int,
    metadata: dict[str, Any],
) -> dict[float | str, Any]:
    """Sweep chemical potentials and temperatures, returning nested dict.

    Args:
        calc: Transport calculator instance.
        fermi: Fermi level in eV.
        mu_shifts: Chemical potential shifts relative to Fermi level in eV.
        temperatures: Temperature array in kelvin.
        kpoint_mesh: Integration mesh dimensions.
        chunk_size: k-point chunk size for batched evaluation.
        metadata: Dictionary with Fermi energy, jspins, parser, and calculator info.

    Returns:
        Nested dictionary keyed by chemical potential shift and temperature.

    """
    mu_abs = fermi + mu_shifts
    scan = calc.calculate_transport_scan(mu_abs, temperatures, kpoint_mesh, kchunk=chunk_size)
    return _format_mu_scan(mu_shifts, temperatures, scan, metadata)


def _completed_result_from_checkpoint(
    loaded_checkpoint: Checkpoint | None,
    checkpoint_path: str | Path | None,
    fingerprints: dict[str, str],
    use_mu_scan: bool,
    mu_shifts: npt.NDArray[np.float64],
    temperatures: npt.NDArray[np.float64],
) -> dict[str, Any] | dict[float | str, Any] | None:
    """Return completed checkpoint results when the result fingerprint matches."""
    if (
        loaded_checkpoint is None
        or not _checkpoint_reached(loaded_checkpoint, TransportWorkflowStage.COMPLETED)
        or not _checkpoint_fingerprint_matches(loaded_checkpoint, "result", fingerprints["result"])
    ):
        return None

    result_entry = loaded_checkpoint.progress.metadata.get("result")
    if not isinstance(result_entry, dict) or not isinstance(result_entry.get("metadata"), dict):
        raise CheckpointFormatError("checkpoint is missing completed result metadata")
    scan = _scan_result_from_checkpoint(loaded_checkpoint)
    result_metadata = dict(result_entry["metadata"])
    print(f"Loaded completed transport result from checkpoint: {checkpoint_path}")
    if use_mu_scan:
        return _format_mu_scan(mu_shifts, temperatures, scan, result_metadata)
    return _format_temperature_scan(temperatures, scan, result_metadata)


def _load_or_parse_electronic_structure(
    parser_obj: ElectronicStructureParser,
    filepath: str | Path,
    loaded_checkpoint: Checkpoint | None,
    checkpoint_manager: CheckpointManager | None,
    checkpoint_path: str | Path | None,
    fingerprints: dict[str, str],
    runtime_config: dict[str, Any],
    temperatures: npt.NDArray[np.float64],
    mu_shifts: npt.NDArray[np.float64],
) -> ParserResult:
    """Load parser state from checkpoint or parse the input file."""
    if (
        loaded_checkpoint is not None
        and _checkpoint_reached(loaded_checkpoint, TransportWorkflowStage.PARSED)
        and _checkpoint_fingerprint_matches(loaded_checkpoint, "input", fingerprints["input"])
    ):
        print(f"Loaded parsed state from checkpoint: {checkpoint_path}")
        return _parser_result_from_checkpoint(loaded_checkpoint)

    parsed = parser_obj.parse(filepath)
    _save_transport_checkpoint(
        checkpoint_manager,
        checkpoint_path,
        TransportWorkflowStage.PARSED,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
        parsed=parsed,
    )
    return parsed


def _select_interpolation_eigenvalues(
    parsed: ParserResult,
    fermi: float,
    band_window: tuple[float, float] | None,
) -> npt.NDArray[np.float64]:
    """Return eigenvalues after applying the optional relative band window."""
    eigenvalues = parsed.eigenvalues
    if band_window is None:
        return eigenvalues

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
    print(f"  band window [{emin}, {emax}] eV -> bands {b_lo}..{b_hi - 1} ({b_hi - b_lo})")
    return eigenvalues[:, :, b_lo:b_hi]


def _load_or_fit_interpolator(
    parsed: ParserResult,
    fermi: float,
    lr_ratio: int,
    band_window: tuple[float, float] | None,
    loaded_checkpoint: Checkpoint | None,
    checkpoint_manager: CheckpointManager | None,
    checkpoint_path: str | Path | None,
    fingerprints: dict[str, str],
    runtime_config: dict[str, Any],
    temperatures: npt.NDArray[np.float64],
    mu_shifts: npt.NDArray[np.float64],
) -> SKWInterpolator:
    """Load a fitted interpolator from checkpoint or fit a new one."""
    if (
        loaded_checkpoint is not None
        and _checkpoint_reached(loaded_checkpoint, TransportWorkflowStage.INTERPOLATED)
        and _checkpoint_fingerprint_matches(loaded_checkpoint, "interpolation", fingerprints["interpolation"])
    ):
        print(f"Loaded SKW interpolation from checkpoint: {checkpoint_path}")
        return _interpolator_from_checkpoint(loaded_checkpoint)

    eigenvalues = _select_interpolation_eigenvalues(parsed, fermi, band_window)
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
    _save_transport_checkpoint(
        checkpoint_manager,
        checkpoint_path,
        TransportWorkflowStage.INTERPOLATED,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
        parsed=parsed,
        interpolator=interp,
    )
    return interp


def _calculate_scan_with_checkpoint(
    calc: Any,
    mu_abs: npt.NDArray[np.float64],
    temperatures: npt.NDArray[np.float64],
    kpoint_mesh: tuple[int, int, int],
    chunk_size: int,
    parsed: ParserResult,
    interp: SKWInterpolator,
    loaded_checkpoint: Checkpoint | None,
    checkpoint_manager: CheckpointManager | None,
    checkpoint_path: str | Path | None,
    fingerprints: dict[str, str],
    runtime_config: dict[str, Any],
    mu_shifts: npt.NDArray[np.float64],
) -> tuple[dict[str, Any], TransportDOS | None]:
    """Return scan arrays, using transport-DOS checkpoints when available."""
    transport_dos: TransportDOS | None = None
    if (
        loaded_checkpoint is not None
        and _checkpoint_reached(loaded_checkpoint, TransportWorkflowStage.TRANSPORT_DOS)
        and _checkpoint_fingerprint_matches(loaded_checkpoint, "transport_dos", fingerprints["transport_dos"])
    ):
        transport_dos = _transport_dos_from_checkpoint(loaded_checkpoint)
        print(f"Loaded transport DOS from checkpoint: {checkpoint_path}")

    if transport_dos is None and checkpoint_manager is None:
        return calc.calculate_transport_scan(mu_abs, temperatures, kpoint_mesh, kchunk=chunk_size), None

    if transport_dos is None:
        transport_dos = calc.build_transport_dos(mu_abs, temperatures, kpoint_mesh, kchunk=chunk_size)
        _save_transport_checkpoint(
            checkpoint_manager,
            checkpoint_path,
            TransportWorkflowStage.TRANSPORT_DOS,
            fingerprints,
            runtime_config,
            temperatures,
            mu_shifts,
            parsed=parsed,
            interpolator=interp,
            transport_dos=transport_dos,
        )

    return calc.integrate_transport_dos_scan(transport_dos, mu_abs, temperatures), transport_dos


def calculate_spin_polarized_transport(
    filepath: str | Path,
    temperature: float | npt.ArrayLike = 300.0,
    chemical_potential: float | npt.ArrayLike | None = None,
    tau: float = 1e-14,
    kpoint_mesh: tuple[int, int, int] = (20, 20, 20),
    lr_ratio: int = 5,
    band_window: tuple[float, float] | None = None,
    chunk_size: int = 4096,
    energy_window_kbt_factor: float = DEFAULTS.energy_window_kbt_factor,
    min_energy_window: float = DEFAULTS.min_energy_window,
    energy_step_kbt_divisor: float = DEFAULTS.energy_step_kbt_divisor,
    min_energy_step: float = DEFAULTS.min_energy_step,
    low_temp_kbt_threshold: float = DEFAULTS.low_temp_kbt_threshold,
    low_temp_energy_window: float = DEFAULTS.low_temp_energy_window,
    low_temp_energy_step: float = DEFAULTS.low_temp_energy_step,
    parser: str | ElectronicStructureParser = DEFAULT_PARSER,
    calculator: str = "boltzmann",
    checkpoint_path: str | Path | None = None,
    resume_checkpoint: bool = True,
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
        checkpoint_path: Optional ``.npz`` checkpoint path used to persist and resume workflow state.
        resume_checkpoint: Whether to resume from ``checkpoint_path`` when it already exists.

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
    temperatures = _prepare_temperature_array(temperature)

    # --- chemical potential grid -------------------------------------------
    use_mu_scan = chemical_potential is not None
    mu_shifts = _prepare_mu_shifts(chemical_potential) if use_mu_scan else np.array([0.0], dtype=np.float64)

    energy_settings = _energy_grid_settings(
        energy_window_kbt_factor,
        min_energy_window,
        energy_step_kbt_divisor,
        min_energy_step,
        low_temp_kbt_threshold,
        low_temp_energy_window,
        low_temp_energy_step,
    )
    runtime_config: dict[str, Any] = {
        "filepath": str(Path(filepath).resolve()),
        "parser": parser_obj.name,
        "calculator": calculator,
        "tau": tau,
        "kpoint_mesh": kpoint_mesh,
        "lr_ratio": lr_ratio,
        "band_window": band_window,
        "chunk_size": chunk_size,
        "use_mu_scan": use_mu_scan,
        "energy_settings": energy_settings,
    }
    checkpoint_manager = CheckpointManager() if checkpoint_path is not None else None
    fingerprints = (
        _transport_workflow_fingerprints(
            filepath,
            parser_obj.name,
            calculator,
            temperatures,
            mu_shifts,
            tau,
            kpoint_mesh,
            lr_ratio,
            band_window,
            energy_settings,
            use_mu_scan,
        )
        if checkpoint_manager is not None
        else {}
    )
    loaded_checkpoint = _load_transport_checkpoint(checkpoint_manager, checkpoint_path, resume_checkpoint)

    completed_result = _completed_result_from_checkpoint(
        loaded_checkpoint,
        checkpoint_path,
        fingerprints,
        use_mu_scan,
        mu_shifts,
        temperatures,
    )
    if completed_result is not None:
        return completed_result

    parsed = _load_or_parse_electronic_structure(
        parser_obj,
        filepath,
        loaded_checkpoint,
        checkpoint_manager,
        checkpoint_path,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
    )

    fermi = parsed.fermi_energy
    metadata = {
        "fermi_energy": fermi,
        "jspins": parsed.jspins,
        "parser": parser_obj.name,
        "calculator": calculator,
    }

    interp = _load_or_fit_interpolator(
        parsed,
        fermi,
        lr_ratio,
        band_window,
        loaded_checkpoint,
        checkpoint_manager,
        checkpoint_path,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
    )

    print(
        f"Transport (mesh={kpoint_mesh}, T={temperature} K, tau={tau:.1e} s, "
        f"chunk={chunk_size}, calculator='{calculator}')..."
    )
    if use_mu_scan:
        print(f"  Chemical potential shifts: {mu_shifts} eV relative to E_Fermi={fermi:.4f} eV")

    calc: Any = calculator_cls(
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

    scan, transport_dos = _calculate_scan_with_checkpoint(
        calc,
        fermi + mu_shifts,
        temperatures,
        kpoint_mesh,
        chunk_size,
        parsed,
        interp,
        loaded_checkpoint,
        checkpoint_manager,
        checkpoint_path,
        fingerprints,
        runtime_config,
        mu_shifts,
    )

    _save_transport_checkpoint(
        checkpoint_manager,
        checkpoint_path,
        TransportWorkflowStage.COMPLETED,
        fingerprints,
        runtime_config,
        temperatures,
        mu_shifts,
        parsed=parsed,
        interpolator=interp,
        transport_dos=transport_dos,
        scan=scan,
        result_metadata=metadata,
    )

    if use_mu_scan:
        return _format_mu_scan(mu_shifts, temperatures, scan, metadata)
    return _format_temperature_scan(temperatures, scan, metadata)
