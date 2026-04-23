## Context

`quantum_macaroni` has grown organically with duplicated logic (Onsager conversion, temperature normalization, default parameters), naming inconsistencies (`TETRAEDRON` typo, legacy `Tuple` import), magic numbers without names, incomplete docstrings, and missing input validation. The `AGENTS.md` file now establishes conventions but the codebase itself does not yet follow them consistently.

The change is cross-cutting: it touches `core/`, `parsers/`, `mesh/`, `interpolation/`, and `calculators/`. No new features are added—this is purely internal quality improvement.

## Goals / Non-Goals

**Goals:**
- Eliminate code duplication in `transport.py` by extracting shared helpers
- Align all code with `AGENTS.md` conventions (type hints, docstrings, naming)
- Replace magic numbers with named constants
- Add input validation for public-facing parameters
- Remove dead code
- Make error messages informative rather than `AttributeError` crashes

**Non-Goals:**
- Splitting `transport.py` into separate modules (too disruptive for a refactor change)
- Introducing a generic `Registry[T]` class (only 2 registries, not worth the abstraction)
- Changing the public API surface (all changes are internal)
- Modifying Numba kernel internals beyond constant reuse and typo fix
- Adding new dependencies beyond what already exists

## Decisions

1. **Onsager-to-transport helper as a module-level function** — Extract `_onsager_to_transport(l0, l1, l2, temperature)` returning `(sigma, seebeck, kappa)`. Alternative: method on `BoltzmannTransportCalculator`. Rationale: the conversion is a pure mathematical function with no state; a free function is simpler and callable from both `calculate_transport` and `calculate_transport_scan`.

2. **Shared defaults via a `dataclass`** — Create `EnergyGridDefaults` dataclass holding the 7 energy-grid parameters. Both `BoltzmannTransportCalculator.__init__` and `calculate_spin_polarized_transport` accept `EnergyGridDefaults` or individual overrides. Alternative: module-level constants. Rationale: dataclass groups related parameters, makes defaults explicit, and avoids duplicating 7 default values.

3. **Replace `_sym2z()` with `ase.data.chemical_symbols`** — `ase` is already a dependency. Alternative: keep inline table. Rationale: removes 100 lines, uses a maintained upstream source, and `ase.data.chemical_symbols` is a stable API.

4. **Name magic constants after their physics meaning** — `ANG3_TO_M3 = 1e30` (converts per Å³ → per m³; the lattice is in Ångströms, not Bohr — 1/Bohr³ ≈ 6.748e30, which is not 1e30), `ANG_TO_M = 1e-10`, `SKW_DAMPING_C1 = 0.25`, `SKW_DAMPING_C2 = 0.25`, `SHELL_GROUPING_RTOL = 1e-8`, `SHELL_GROUPING_FLOOR = 1e-30`. Rationale: names convey intent; values are easier to audit; wrong unit labels cause silent bugs.

5. **Fix typo as a rename with backward-compatible alias** — Rename `MIN_TETRAEDRON_ENERGY` → `MIN_TETRAHEDRON_ENERGY` and `MIN_TETRAEDRON_WEIGHT` → `MIN_TETRAHEDRON_WEIGHT`. Add deprecated aliases for one release cycle. Alternative: silent rename. Rationale: these are module-level constants that may be imported by external code.

6. **Validation via existing `_validate_bound_parameter`** — Extend `_validate_bound_parameter` usage to `tau`, `chunk_size`, `kpoint_mesh`, and `lr_ratio`. Alternative: `pydantic` or `__post_init__`. Rationale: the validation helper already exists and is used; adding calls is minimal and consistent.

7. **Graceful XML error handling** — Add explicit `None` checks after `root.find(...)` with descriptive `ValueError` messages. Alternative: `try/except AttributeError`. Rationale: explicit checks give better error messages and are more readable.

## Risks / Trade-offs

- [Risk: Typo rename breaks external imports of `MIN_TETRAEDRON_*`] → Mitigation: add deprecated aliases that log a warning for one release
- [Risk: Replacing `_sym2z` with `ase.data.chemical_symbols` changes element ordering] → Mitigation: `ase.data.chemical_symbols` returns list indexed by atomic number (1-indexed with empty string at 0), which is functionally identical to the dict we build
- [Risk: Refactoring `transport.py` may introduce subtle numerical differences] → Mitigation: extract helpers as exact copies first, then optimize; run `pytest` after each change
- [Risk: Dataclass defaults change adds import dependency] → Mitigation: `dataclasses` is stdlib; no new dependency