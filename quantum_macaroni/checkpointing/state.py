"""State containers for restartable computations."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TypeAlias

import numpy as np
import numpy.typing as npt

JsonScalar: TypeAlias = str | int | float | bool | None
JsonValue: TypeAlias = JsonScalar | list["JsonValue"] | dict[str, "JsonValue"]
Metadata: TypeAlias = dict[str, JsonValue]
ArrayMap: TypeAlias = dict[str, npt.NDArray[np.generic]]


@dataclass(slots=True)
class BaseSystemState:
    """Heavy foundational state reused across restarted or branched runs.

    Attributes:
        version: User-defined version or fingerprint for the foundational system.
        arrays: Large numerical arrays that define the base system.
        metadata: Small JSON-compatible descriptors for the base system.

    """

    version: str = "1"
    arrays: ArrayMap = field(default_factory=dict)
    metadata: Metadata = field(default_factory=dict)


@dataclass(slots=True)
class RuntimeParameters:
    """Run-specific parameters that may be replaced when branching a checkpoint.

    Attributes:
        values: Small JSON-compatible runtime settings.
        arrays: Numerical runtime arrays such as boundary conditions.

    """

    values: Metadata = field(default_factory=dict)
    arrays: ArrayMap = field(default_factory=dict)


@dataclass(slots=True)
class ExecutionProgress:
    """Intermediate computation progress needed to resume a run.

    Attributes:
        step: Current iteration, time step, epoch, or work-unit index.
        arrays: Numerical arrays that evolve during execution.
        metadata: Small JSON-compatible progress descriptors.

    """

    step: int = 0
    arrays: ArrayMap = field(default_factory=dict)
    metadata: Metadata = field(default_factory=dict)


@dataclass(slots=True)
class Checkpoint:
    """Complete persisted snapshot split into base, runtime, and progress state.

    Attributes:
        base_system: Heavy foundational state loaded exactly as saved.
        runtime_parameters: Run-specific settings for this branch.
        progress: Intermediate computation state at the checkpointed step.

    """

    base_system: BaseSystemState
    runtime_parameters: RuntimeParameters
    progress: ExecutionProgress

    def with_runtime_parameters(self, runtime_parameters: RuntimeParameters) -> Checkpoint:
        """Return a branch that keeps base/progress state and swaps runtime parameters.

        Args:
            runtime_parameters: Runtime settings for the new branch.

        Returns:
            New checkpoint view with replaced runtime parameters.

        """
        return replace(self, runtime_parameters=runtime_parameters)
