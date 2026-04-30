## ADDED Requirements

### Requirement: Duplicated Onsager-to-transport conversion extracted to helper
The Onsager coefficient to physical transport tensor conversion SHALL be implemented in a single shared function `_onsager_to_transport(l0, l1, l2, temperature)` in `calculators/transport.py`. Both `calculate_transport` and `calculate_transport_scan` SHALL call this helper instead of inlining the conversion.

#### Scenario: Onsager conversion is not duplicated
- **WHEN** the codebase is searched for `l0_inv = np.linalg.inv(l0)` patterns
- **THEN** only one implementation exists, inside `_onsager_to_transport`

### Requirement: Duplicated temperature normalization removed
`calculate_spin_polarized_transport` SHALL call `_prepare_temperature_array()` instead of reimplementing the same logic inline.

#### Scenario: Temperature validation is centralized
- **WHEN** `calculate_spin_polarized_transport` receives a temperature argument
- **THEN** it delegates to `_prepare_temperature_array()` rather than containing its own `np.asarray` + validation logic

### Requirement: Shared default parameters for energy grid
The 7 energy-grid parameters (`energy_window_kbt_factor`, `min_energy_window`, `energy_step_kbt_divisor`, `min_energy_step`, `low_temp_kbt_threshold`, `low_temp_energy_window`, `low_temp_energy_step`) SHALL be defined in a single place. `BoltzmannTransportCalculator.__init__` and `calculate_spin_polarized_transport` SHALL reference these shared defaults.

#### Scenario: Default parameter values are not duplicated
- **WHEN** a default energy-grid parameter value needs to change
- **THEN** it only needs to be updated in one location

### Requirement: Misspelled tetrahedron constants renamed
`MIN_TETRAEDRON_ENERGY` SHALL be renamed to `MIN_TETRAHEDRON_ENERGY` and `MIN_TETRAEDRON_WEIGHT` SHALL be renamed to `MIN_TETRAHEDRON_WEIGHT` in `core/numerics.py`. Deprecated aliases with the old names SHALL be provided temporarily.

#### Scenario: Correctly spelled constants are used
- **WHEN** a developer references `MIN_TETRAHEDRON_ENERGY`
- **THEN** it resolves to the same value as the old `MIN_TETRAEDRON_ENERGY`

### Requirement: Legacy Tuple import replaced with built-in tuple
`from typing import Tuple` SHALL be removed from `mesh/tetrahedron.py` and `Tuple[int, int, int]` SHALL be replaced with `tuple[int, int, int]`.

#### Scenario: No typing.Tuple imports in mesh module
- **WHEN** `ruff check` is run on `mesh/tetrahedron.py`
- **THEN** no `from typing import Tuple` import exists and `tuple[...]` generic syntax is used

### Requirement: Hardcoded element table replaced with ASE data
The `_sym2z()` function in `parsers/fleur_outxml.py` SHALL be replaced with `ase.data.chemical_symbols`, which is already available through the `ase` dependency.

#### Scenario: Element table is not hardcoded
- **WHEN** `FleurOutxmlParser.parse` resolves atomic numbers from element symbols
- **THEN** it uses `ase.data.chemical_symbols` rather than a locally defined mapping function

### Requirement: Magic numbers extracted to named constants
The following magic numbers SHALL be replaced with named constants in appropriate modules:

- `1e30` unit conversion factor → `ANG3_TO_M3` in `core/constants.py` (converts per Å³ → per m³; the lattice is in Ångströms, not Bohr)
- `1e-10` angstrom-to-meter conversion → `ANG_TO_M` in `core/constants.py`
- `0.25` SKW damping parameters → `SKW_DAMPING_C1`, `SKW_DAMPING_C2` in `interpolation/skw.py`
- `1e-8` / `1e-30` shell-grouping tolerances → named constants in `interpolation/skw.py`
- `6.0` tetrahedra-per-cube → use `TetrahedronMesh._TETRAHEDRA_CONFIGS.shape[1]` or a named constant

#### Scenario: Named constants convey physical meaning
- **WHEN** a developer reads `mesh.tetra_vol / (TWO_PI**3) * ANG3_TO_M3`
- **THEN** the intent (unit conversion from per Å³ to per m³) is clear from the constant name

### Requirement: Missing docstring sections completed
All public functions and methods across `calculators/base.py`, `calculators/transport.py`, and `parsers/fleur_outxml.py` SHALL have Google-style docstrings with `Args`, `Returns`, and `Raises` sections where applicable.

#### Scenario: Docstrings are complete
- **WHEN** `ruff check --select D` is run
- **THEN** no missing-docstring-section violations are reported for these modules

### Requirement: Input validation added for key parameters
The following parameters SHALL be validated:

- `tau` in `BoltzmannTransportCalculator.__init__` — must be positive
- `chunk_size` — must be positive
- `kpoint_mesh` dimensions in `TetrahedronMesh` — must be positive integers
- `lr_ratio` in `SKWInterpolator.__init__` — must be positive
- `spin` in `SKWInterpolator.eval_energy_velocity` — must be in range

#### Scenario: Invalid tau raises ValueError
- **WHEN** `BoltzmannTransportCalculator(interpolator, tau=-1e-14)` is called
- **THEN** a `ValueError` is raised with a descriptive message

#### Scenario: Invalid kpoint_mesh raises ValueError
- **WHEN** `TetrahedronMesh(lattice, (0, 20, 20))` is called
- **THEN** a `ValueError` is raised indicating mesh dimensions must be positive

### Requirement: Graceful XML error messages in FleurOutxmlParser
`FleurOutxmlParser.parse` and its helper functions SHALL check for missing XML elements (`<magnetism>`, `<FermiEnergy>`, `<bravaisMatrix>`) and raise `ValueError` with descriptive messages rather than crashing with `AttributeError` on `None`.

#### Scenario: Missing magnetism element
- **WHEN** a FLEUR XML file lacks the `<magnetism>` element
- **THEN** a `ValueError` is raised with message "Missing required element <magnetism> in out.xml"

### Requirement: Dead commented-out code removed
The commented-out debug `print` statements in `calculate_spin_polarized_transport` (lines 573-578 of transport.py) SHALL be removed.

#### Scenario: No commented-out debug prints in transport module
- **WHEN** the `calculators/transport.py` file is searched for `# print(` patterns
- **THEN** no commented-out debug prints exist