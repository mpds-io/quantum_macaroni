# quantum_macaroni

quantum_macaroni is a modular Boltzmann-transport workflow for electronic-structure data.

Current pipeline:
- Parser plugins for electronic-structure outputs (default parser: Fleur out.xml).
- SKW interpolation of band energies.
- Tetrahedron k-space integration mesh.
- Transport-property calculators (default: Boltzmann transport calculator).

## Features

- Plugin architecture for parsers and calculators.
- Transport tensors and isotropic averages.
- Temperature sweep with scalar or array input.
- Chemical-potential sweep relative to Fermi level.
- CLI interface with JSON output.

## Requirements

- Python 3.11+
- Core dependencies:
	- ase
	- lxml
	- numba
	- numpy
	- scipy

Dependencies are defined in [pyproject.toml](pyproject.toml).

## Installation

Use your preferred environment manager.

Example with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Quick Start (CLI)

Main entry point is [main.py](main.py).

Minimal run:

```bash
python main.py examples/PbTe-nospin/out-nospin.xml
```

Run with temperature and chemical-potential sweeps:

```bash
python main.py examples/PbTe-nospin/out-nospin.xml \
	--temperature 300 900 7 \
	--chemical-potential -0.5 0.5 11 \
	--kmesh 80 80 80 \
	--lr-ratio 20 \
	--band-window -3 3 \
	--output transport_results.json
```

### CLI Argument Rules for Temperature and Chemical Potential

Both arguments accept either:
- one number
- or three numbers: start, stop, number_of_points

Examples:
- `--temperature 300`
- `--temperature 300 900 7`
- `--chemical-potential 0.0`
- `--chemical-potential -0.3 0.3 13`

For three-number form, the third value must be a positive integer (number of points).

## CLI Reference

```text
python main.py FILEPATH [options]

Options:
	--temperature T [T ...]                one value or (start stop npoints)
	--chemical-potential MU [MU ...]       one value or (start stop npoints), eV shift from E_F
	--tau FLOAT                            relaxation time in seconds (default: 1e-14)
	--kmesh NX NY NZ                       k-point mesh (default: 80 80 80)
	--lr-ratio INT                         SKW interpolator star-vector ratio (default: 20)
	--band-window EMIN EMAX                band window relative to E_F in eV (default: -3 3)
	--chunk-size INT                       chunk size for batched evaluation (default: 4096)
	--parser {available_parsers}           parser plugin (default: fleur-outxml)
	--calculator {available_calculators}   calculator plugin (default: boltzmann)
	--output PATH                          output JSON file path (default: transport_results.json)
```

Available parser/calculator names come from the runtime registries in
[quantum_macaroni/parsers/__init__.py](quantum_macaroni/parsers/__init__.py) and
[quantum_macaroni/calculators/__init__.py](quantum_macaroni/calculators/__init__.py).

## Output JSON Format

The CLI stores calculation output to a JSON file (default: transport_results.json).

When chemical potential is provided, structure is:

```json
{
	"-0.5": {
		"300.0": {
			"sigma": [[...], [...], [...]],
			"sigma_avg": 0.0,
			"seebeck": [[...], [...], [...]],
			"seebeck_avg": 0.0,
			"kappa": [[...], [...], [...]],
			"kappa_avg": 0.0
		}
	},
	"0.0": {
		"300.0": {
			"sigma": [[...], [...], [...]],
			"sigma_avg": 0.0,
			"seebeck": [[...], [...], [...]],
			"seebeck_avg": 0.0,
			"kappa": [[...], [...], [...]],
			"kappa_avg": 0.0
		}
	},
	"meta": {
		"fermi_energy": 0.0,
		"jspins": 1,
		"parser": "fleur-outxml",
		"calculator": "boltzmann"
	}
}
```

Note: JSON keys are strings, so numeric keys for chemical potential and temperature are serialized as strings.

## Python API

Public API is exported from [quantum_macaroni/__init__.py](quantum_macaroni/__init__.py).

Main high-level function:
- `calculate_spin_polarized_transport`

Example:

```python
import numpy as np
from quantum_macaroni import calculate_spin_polarized_transport

result = calculate_spin_polarized_transport(
		"examples/PbTe-nospin/out-nospin.xml",
		temperature=np.linspace(300.0, 900.0, 7),
		chemical_potential=np.linspace(-0.5, 0.5, 11),
		tau=1e-14,
		kpoint_mesh=(80, 80, 80),
		lr_ratio=20,
		band_window=(-3.0, 3.0),
		chunk_size=4096,
)
```

For backward-compatible example script, see [examples/PbTe-nospin/boltz.py](examples/PbTe-nospin/boltz.py).

## Project Layout

```text
quantum_macaroni/
	calculators/      transport calculators and registry
	core/             constants and numerics
	interpolation/    SKW interpolator
	mesh/             tetrahedron mesh
	parsers/          parser interfaces and implementations
examples/
	PbTe-nospin/      sample input and usage script
main.py             CLI entry point
```
