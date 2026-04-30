## 1. Constants and Naming Fixes

- [x] 1.1 Rename `MIN_TETRAEDRON_ENERGY` → `MIN_TETRAHEDRON_ENERGY` and `MIN_TETRAEDRON_WEIGHT` → `MIN_TETRAHEDRON_WEIGHT` in `core/numerics.py`; add deprecated aliases for backward compatibility; update all references
- [x] 1.2 Add named constants `ANG3_TO_M3 = 1e30` (converts per Å³ → per m³; lattice is in Ångströms, not Bohr) and `ANG_TO_M = 1e-10` to `core/constants.py`; replace magic numbers in `calculators/transport.py` and `core/numerics.py`
- [x] 1.3 Add named constants `SKW_DAMPING_C1 = 0.25`, `SKW_DAMPING_C2 = 0.25`, `SHELL_GROUPING_RTOL = 1e-8`, `SHELL_GROUPING_FLOOR = 1e-30` in `interpolation/skw.py`; replace inline magic numbers
- [x] 1.4 Replace `from typing import Tuple` with built-in `tuple[int, int, int]` in `mesh/tetrahedron.py`; remove the import

## 2. Deduplication and Extraction

- [x] 2.1 Extract `_onsager_to_transport(l0, l1, l2, temperature)` helper from `calculators/transport.py`; refactor `calculate_transport` and `calculate_transport_scan` to call it
- [x] 2.2 Replace inline temperature normalization in `calculate_spin_polarized_transport` with a call to `_prepare_temperature_array()`
- [x] 2.3 Create `EnergyGridDefaults` dataclass with the 7 energy-grid default parameters; update `BoltzmannTransportCalculator.__init__` and `calculate_spin_polarized_transport` to reference shared defaults
- [x] 2.4 Extract spin-degeneracy factor and normalization into a shared approach (method or constant) to avoid duplicating `if nspin == 1: *= 2.0` and `norm = tetra_vol / (TWO_PI**3) * ANG3_TO_M3`

## 3. Replace Hardcoded Element Table

- [x] 3.1 Replace `_sym2z()` in `parsers/fleur_outxml.py` with `ase.data.chemical_symbols`; update `structure_from_outxml` to use the ASE data source; remove the `_sym2z` function

## 4. Input Validation

- [x] 4.1 Add `_validate_bound_parameter("tau", tau)` call in `BoltzmannTransportCalculator.__init__`
- [x] 4.2 Add validation for `chunk_size > 0`, `lr_ratio > 0`, and `spin` index bounds in `SKWInterpolator.__init__` and `eval_energy_velocity`
- [x] 4.3 Add validation for positive `kpoint_mesh` dimensions in `TetrahedronMesh.__init__`
- [x] 4.4 Add `None` checks with descriptive `ValueError` messages for missing XML elements in `FleurOutxmlParser.parse` and `structure_from_outxml` (`<magnetism>`, `<FermiEnergy>`, `<bravaisMatrix>`)

## 5. Documentation and Cleanup

- [x] 5.1 Complete missing docstring sections (`Args`, `Returns`, `Raises`) in `calculators/base.py` for `register_calculator`, `get_calculator`, `available_calculators`, `tensor_average`
- [x] 5.2 Complete missing docstring sections in `calculators/transport.py` for `_build_unified_energy_grid`, `_transport_result_units`, `_compute_mu_scan`
- [x] 5.3 Complete missing docstring sections in `parsers/fleur_outxml.py` for `_sym2z` (or its replacement) and add `Raises` to `FleurOutxmlParser.parse`
- [x] 5.4 Add `Raises` sections to `BoltzmannTransportCalculator.calculate_onsager_coefficients` and `calculate_transport_scan`
- [x] 5.5 Remove commented-out debug `print` statements from `calculate_spin_polarized_transport`
- [x] 5.6 Use `TransportTuple` return type alias in `BoltzmannTransportCalculator.calculate_transport` instead of inlining the 6-tuple

## 6. Verify

- [x] 6.1 Run `ruff check .` and fix any violations
- [x] 6.2 Run `ruff format .` to ensure consistent formatting
- [x] 6.3 Run `pytest` to confirm no regressions