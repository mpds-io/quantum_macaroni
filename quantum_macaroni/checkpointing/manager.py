"""Atomic NumPy checkpoint serialization."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Final, cast

import numpy as np
import numpy.typing as npt

from quantum_macaroni.checkpointing.state import (
    ArrayMap,
    BaseSystemState,
    Checkpoint,
    ExecutionProgress,
    JsonValue,
    Metadata,
    RuntimeParameters,
)

CHECKPOINT_FORMAT: Final = "quantum_macaroni.checkpoint.v1"
MANIFEST_KEY: Final = "__manifest__"
BASE_ARRAY_PREFIX: Final = "base/arrays/"
RUNTIME_ARRAY_PREFIX: Final = "runtime/arrays/"
PROGRESS_ARRAY_PREFIX: Final = "progress/arrays/"
OBJECT_DTYPE_KIND: Final = "O"


class CheckpointFormatError(ValueError):
    """Raised when a checkpoint archive is malformed or unsupported."""


class CheckpointManager:
    """Save and load decoupled computation state using atomic NumPy archives."""

    def __init__(self, root: str | Path = ".") -> None:
        """Initialize a manager rooted at ``root``.

        Args:
            root: Base directory used for relative checkpoint paths.

        """
        self.root = Path(root)

    def save_checkpoint(self, checkpoint: Checkpoint, filepath: str | Path) -> Path:
        """Save a checkpoint atomically and return the resolved path.

        Args:
            checkpoint: State snapshot to persist.
            filepath: Destination ``.npz`` path. Relative paths are resolved under ``root``.

        Returns:
            Resolved checkpoint path.

        Raises:
            TypeError: If an array has object dtype.
            ValueError: If metadata is not JSON-serializable.

        """
        target = self._resolve_path(filepath)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = _checkpoint_to_payload(checkpoint)
        _write_npz_atomically(target, payload)
        return target

    def load_checkpoint(
        self,
        filepath: str | Path,
        runtime_parameters: RuntimeParameters | None = None,
    ) -> Checkpoint:
        """Load a checkpoint, optionally replacing its runtime parameters.

        Args:
            filepath: Source checkpoint path. Relative paths are resolved under ``root``.
            runtime_parameters: Optional runtime parameters for a new branch.

        Returns:
            Loaded checkpoint. If ``runtime_parameters`` is provided, the returned checkpoint
            keeps the persisted base/progress state and uses the replacement parameters.

        Raises:
            CheckpointFormatError: If the archive is malformed or unsupported.

        """
        target = self._resolve_path(filepath)
        with np.load(target, allow_pickle=False) as archive:
            checkpoint = _checkpoint_from_archive(archive)

        if runtime_parameters is None:
            return checkpoint
        return checkpoint.with_runtime_parameters(runtime_parameters)

    def load_for_branch(self, filepath: str | Path, runtime_parameters: RuntimeParameters) -> Checkpoint:
        """Load a checkpoint and inject new runtime parameters for a branch.

        Args:
            filepath: Source checkpoint path. Relative paths are resolved under ``root``.
            runtime_parameters: Runtime settings for the resumed branch.

        Returns:
            Branched checkpoint with persisted base/progress and new runtime settings.

        """
        return self.load_checkpoint(filepath, runtime_parameters=runtime_parameters)

    def _resolve_path(self, filepath: str | Path) -> Path:
        """Return an absolute checkpoint path under the manager root when needed."""
        path = Path(filepath)
        if path.is_absolute():
            return path
        return self.root / path


def _checkpoint_to_payload(checkpoint: Checkpoint) -> dict[str, npt.NDArray[np.generic]]:
    """Convert a checkpoint into arrays accepted by ``numpy.savez_compressed``."""
    base_arrays = _normalize_arrays(checkpoint.base_system.arrays)
    runtime_arrays = _normalize_arrays(checkpoint.runtime_parameters.arrays)
    progress_arrays = _normalize_arrays(checkpoint.progress.arrays)

    manifest: dict[str, JsonValue] = {
        "format": CHECKPOINT_FORMAT,
        "created_utc": datetime.now(UTC).isoformat(),
        "base": {
            "version": checkpoint.base_system.version,
            "array_names": sorted(base_arrays),
            "metadata": checkpoint.base_system.metadata,
        },
        "runtime": {
            "array_names": sorted(runtime_arrays),
            "values": checkpoint.runtime_parameters.values,
        },
        "progress": {
            "step": checkpoint.progress.step,
            "array_names": sorted(progress_arrays),
            "metadata": checkpoint.progress.metadata,
        },
    }

    payload: dict[str, npt.NDArray[np.generic]] = {MANIFEST_KEY: _manifest_to_array(manifest)}
    payload.update(_prefix_arrays(BASE_ARRAY_PREFIX, base_arrays))
    payload.update(_prefix_arrays(RUNTIME_ARRAY_PREFIX, runtime_arrays))
    payload.update(_prefix_arrays(PROGRESS_ARRAY_PREFIX, progress_arrays))
    return payload


def _checkpoint_from_archive(archive: Any) -> Checkpoint:
    """Build a checkpoint from a loaded ``numpy.load`` archive."""
    if MANIFEST_KEY not in archive.files:
        raise CheckpointFormatError("checkpoint is missing manifest")

    manifest = _manifest_from_array(archive[MANIFEST_KEY])
    if manifest.get("format") != CHECKPOINT_FORMAT:
        raise CheckpointFormatError("unsupported checkpoint format")

    base_manifest = _require_manifest_mapping(manifest, "base")
    runtime_manifest = _require_manifest_mapping(manifest, "runtime")
    progress_manifest = _require_manifest_mapping(manifest, "progress")

    base_system = BaseSystemState(
        version=_require_string(base_manifest, "version"),
        arrays=_load_arrays(archive, BASE_ARRAY_PREFIX, _require_string_sequence(base_manifest, "array_names")),
        metadata=_require_metadata(base_manifest, "metadata"),
    )
    runtime_parameters = RuntimeParameters(
        values=_require_metadata(runtime_manifest, "values"),
        arrays=_load_arrays(archive, RUNTIME_ARRAY_PREFIX, _require_string_sequence(runtime_manifest, "array_names")),
    )
    progress = ExecutionProgress(
        step=_require_int(progress_manifest, "step"),
        arrays=_load_arrays(archive, PROGRESS_ARRAY_PREFIX, _require_string_sequence(progress_manifest, "array_names")),
        metadata=_require_metadata(progress_manifest, "metadata"),
    )
    return Checkpoint(base_system=base_system, runtime_parameters=runtime_parameters, progress=progress)


def _manifest_to_array(manifest: Mapping[str, JsonValue]) -> npt.NDArray[np.uint8]:
    """Encode manifest JSON into a byte array stored inside the NumPy archive."""
    try:
        encoded = json.dumps(manifest, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError("checkpoint metadata must be JSON-serializable without NaN or Infinity") from exc
    return np.frombuffer(encoded, dtype=np.uint8)


def _manifest_from_array(array: npt.NDArray[np.generic]) -> dict[str, Any]:
    """Decode manifest JSON from a byte array."""
    try:
        decoded = np.asarray(array, dtype=np.uint8).tobytes().decode("utf-8")
        manifest = json.loads(decoded)
    except (TypeError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CheckpointFormatError("checkpoint manifest is not valid JSON") from exc

    if not isinstance(manifest, dict):
        raise CheckpointFormatError("checkpoint manifest must be a JSON object")
    return manifest


def _normalize_arrays(arrays: Mapping[str, npt.NDArray[np.generic]]) -> ArrayMap:
    """Validate array names and return contiguous non-object arrays."""
    normalized: ArrayMap = {}
    for name, array in arrays.items():
        _validate_array_name(name)
        array_value = np.ascontiguousarray(array)
        if array_value.dtype.kind == OBJECT_DTYPE_KIND:
            raise TypeError(f"checkpoint array '{name}' has object dtype, which would require pickle")
        normalized[name] = array_value
    return normalized


def _validate_array_name(name: str) -> None:
    """Validate an array key before placing it in the archive namespace."""
    if not name:
        raise ValueError("checkpoint array names must be non-empty")
    if name.startswith("__") or "/" in name or "\\" in name:
        raise ValueError(f"unsupported checkpoint array name: {name!r}")


def _prefix_arrays(prefix: str, arrays: Mapping[str, npt.NDArray[np.generic]]) -> dict[str, npt.NDArray[np.generic]]:
    """Return archive array mapping with prefixed keys."""
    return {f"{prefix}{name}": array for name, array in arrays.items()}


def _load_arrays(archive: Any, prefix: str, names: Sequence[str]) -> ArrayMap:
    """Load named arrays from an archive prefix."""
    arrays: ArrayMap = {}
    for name in names:
        archive_key = f"{prefix}{name}"
        if archive_key not in archive.files:
            raise CheckpointFormatError(f"checkpoint is missing array payload {archive_key!r}")
        arrays[name] = np.ascontiguousarray(archive[archive_key])
    return arrays


def _require_manifest_mapping(manifest: Mapping[str, Any], key: str) -> dict[str, Any]:
    """Return a required manifest child object."""
    value = manifest.get(key)
    if not isinstance(value, dict):
        raise CheckpointFormatError(f"checkpoint manifest entry {key!r} must be an object")
    return value


def _require_metadata(manifest: Mapping[str, Any], key: str) -> Metadata:
    """Return required JSON metadata object."""
    value = manifest.get(key)
    if not isinstance(value, dict):
        raise CheckpointFormatError(f"checkpoint manifest entry {key!r} must be an object")
    return value


def _require_string(manifest: Mapping[str, Any], key: str) -> str:
    """Return a required string manifest value."""
    value = manifest.get(key)
    if not isinstance(value, str):
        raise CheckpointFormatError(f"checkpoint manifest entry {key!r} must be a string")
    return value


def _require_int(manifest: Mapping[str, Any], key: str) -> int:
    """Return a required integer manifest value."""
    value = manifest.get(key)
    if not isinstance(value, int):
        raise CheckpointFormatError(f"checkpoint manifest entry {key!r} must be an integer")
    return value


def _require_string_sequence(manifest: Mapping[str, Any], key: str) -> tuple[str, ...]:
    """Return a required sequence of string manifest values."""
    value = manifest.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise CheckpointFormatError(f"checkpoint manifest entry {key!r} must be a list of strings")
    return tuple(value)


def _write_npz_atomically(target: Path, payload: Mapping[str, npt.NDArray[np.generic]]) -> None:
    """Write an ``.npz`` payload to a temporary file before replacing ``target``."""
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=target.parent,
            prefix=f".{target.name}.",
            suffix=".tmp",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)
            np.savez_compressed(temp_file, **cast(dict[str, Any], payload))
            temp_file.flush()
            os.fsync(temp_file.fileno())

        os.replace(temp_path, target)
        _fsync_directory(target.parent)
    except BaseException:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise


def _fsync_directory(path: Path) -> None:
    """Best-effort fsync for directory metadata after atomic replacement."""
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY

    try:
        directory_fd = os.open(path, flags)
    except OSError:
        return

    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
