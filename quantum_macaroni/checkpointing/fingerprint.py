"""Stable JSON fingerprints for checkpoint compatibility checks."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from quantum_macaroni.checkpointing.state import JsonValue


def to_json_value(value: Any) -> JsonValue:
    """Return a JSON-compatible representation of common scientific values.

    Args:
        value: Python, NumPy, or path-like value to normalize.

    Returns:
        JSON-compatible value.

    Raises:
        ValueError: If a floating-point value is NaN or infinite.
        TypeError: If the value cannot be represented as JSON.

    """
    normalized: JsonValue
    if value is None or isinstance(value, str | bool | int):
        normalized = value
    elif isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("JSON metadata cannot contain NaN or Infinity")
        normalized = value
    elif isinstance(value, np.generic):
        normalized = to_json_value(value.item())
    elif isinstance(value, np.ndarray):
        normalized = to_json_value(value.tolist())
    elif isinstance(value, Path):
        normalized = str(value)
    elif isinstance(value, Mapping):
        normalized = {
            str(key): to_json_value(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    elif isinstance(value, Sequence) and not isinstance(value, bytes | bytearray):
        normalized = [to_json_value(item) for item in value]
    else:
        raise TypeError(f"unsupported JSON metadata value of type {type(value).__name__}")
    return normalized


def stable_fingerprint(payload: Mapping[str, Any]) -> str:
    """Return a deterministic SHA-256 fingerprint for a JSON-like payload.

    Args:
        payload: Mapping containing values normalized by :func:`to_json_value`.

    Returns:
        Hex SHA-256 digest.

    """
    normalized = to_json_value(payload)
    encoded = json.dumps(normalized, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
