## Why

The codebase has accumulated technical debt: duplicated logic across `transport.py`, a misspelled constant (`TETRAEDRON`), magic numbers without names, incomplete docstrings, legacy `Tuple` imports, hardcoded element tables that duplicate a dependency, commented-out debug code, and missing input validation. Cleaning these up now improves maintainability, reduces bug surface, and aligns the code with the conventions documented in `AGENTS.md`.

## What Changes

- Extract duplicated Onsager-to-transport conversion into a shared helper `_onsager_to_transport()`
- Remove duplicated temperature/chemical-potential normalization — call existing `_prepare_temperature_array()` and `_prepare_mu_shifts()` instead of inlining
- Extract shared default parameters for `BoltzmannTransportCalculator` and `calculate_spin_polarized_transport` into a module-level defaults dataclass or namespace
- Fix typo: rename `MIN_TETRAEDRON_ENERGY` → `MIN_TETRAHEDRON_ENERGY`, `MIN_TETRAEDRON_WEIGHT` → `MIN_TETRAHEDRON_WEIGHT`
- Replace legacy `from typing import Tuple` with `tuple[...]` in `tetrahedron.py`
- Replace hardcoded `_sym2z()` with `ase.data.chemical_symbols` (ase is already a dependency)
- Extract magic numbers into named constants: `ANG3_TO_M3 = 1e30` (converts per Å³ → per m³; the lattice is in Ångströms, not Bohr), `ANG_TO_M = 1e-10`, SKW regularization defaults, shell-grouping tolerances
- Add missing `Args`/`Returns`/`Raises` sections to docstrings in `calculators/base.py`, `calculators/transport.py`, `parsers/fleur_outxml.py`
- Add input validation for `tau`, `kpoint_mesh`, `chunk_size`, `lr_ratio`, and `spin` index
- Add graceful error messages for missing XML elements in `FleurOutxmlParser.parse`
- Remove commented-out debug prints from `transport.py`
- Use `TransportTuple` return type alias where the full 6-tuple is spelled out

## Capabilities

### New Capabilities
- `code-quality`: Refactoring, naming fixes, docstring completion, input validation, magic-number extraction, and dead code removal across the `quantum_macaroni` package

### Modified Capabilities

## Impact

- `quantum_macaroni/core/numerics.py` — constant renames, local `eps` → reference module constant
- `quantum_macaroni/core/constants.py` — new named constants for unit conversions
- `quantum_macaroni/core/symmetry.py` — minor type-hint update (`np.int_` → `np.int64`)
- `quantum_macaroni/mesh/tetrahedron.py` — remove `from typing import Tuple`
- `quantum_macaroni/interpolation/skw.py` — extract SKW defaults to named constants, improve private attribute naming
- `quantum_macaroni/parsers/fleur_outxml.py` — replace `_sym2z` with `ase.data.chemical_symbols`, add error messages for missing XML elements
- `quantum_macaroni/parsers/base.py` — add docstring sections to registry functions
- `quantum_macaroni/calculators/base.py` — add docstring sections to registry functions and `tensor_average`
- `quantum_macaroni/calculators/transport.py` — extract helpers, remove duplication, add validation, remove dead code, complete docstrings
- **No public API changes** — all changes are internal refactoring and documentation improvements