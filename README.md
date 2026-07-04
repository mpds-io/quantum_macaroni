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
- Restartable checkpoints for parsed input, fitted interpolation, transport DOS, and completed scans.
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

By default, the CLI writes restartable checkpoint state to `transport_state.npz`
and resumes from it on later compatible runs.

Run with temperature and chemical-potential sweeps:

```bash
python main.py examples/PbTe-nospin/out-nospin.xml \
	--temperature 300 900 7 \
	--chemical-potential -0.5 0.5 11 \
	--kmesh 80 80 80 \
	--lr-ratio 20 \
	--band-window -3 3 \
	--checkpoint transport_state.npz \
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
	--checkpoint PATH                      restartable .npz checkpoint (default: transport_state.npz)
	--no-checkpoint                        disable checkpoint writing and resume
	--no-resume                            ignore an existing checkpoint while writing a fresh one
```

Available parser/calculator names come from the runtime registries in
[quantum_macaroni/parsers/__init__.py](quantum_macaroni/parsers/__init__.py) and
[quantum_macaroni/calculators/__init__.py](quantum_macaroni/calculators/__init__.py).

## Checkpointing

The CLI saves restartable workflow state to `transport_state.npz` by default as the run advances. A later run with the same input and compatible configuration resumes automatically from the latest valid stage:

- parsed electronic structure
- fitted SKW interpolation
- transport density of states
- completed transport scan

Checkpoints are compressed NumPy archives with a JSON manifest and non-object NumPy arrays, written atomically to avoid partial files. Compatibility is checked with staged fingerprints over the input file content and calculation settings. Use `--checkpoint PATH` to choose a different file, `--no-resume` to ignore an existing checkpoint and replace it during the new run, or `--no-checkpoint` to disable checkpointing.

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
		checkpoint_path="transport_state.npz",
)
```

For backward-compatible example script, see [examples/PbTe-nospin/boltz.py](examples/PbTe-nospin/boltz.py).

## Checkpointing API

The checkpointing layer separates restartable computations into three explicit containers:
- `BaseSystemState`: heavy foundational arrays and metadata that should be reused exactly.
- `RuntimeParameters`: run-specific values and optional runtime arrays that can be replaced when branching.
- `ExecutionProgress`: current step and evolving arrays needed to resume work.

Checkpoints are stored as compressed NumPy `.npz` archives. Array payloads stay in NumPy binary form, while small manifest metadata is stored as JSON bytes inside the archive. Writes are atomic: the manager writes a same-directory temporary file, fsyncs it, and then replaces the target path.

```python
import numpy as np
from quantum_macaroni import BaseSystemState, Checkpoint, CheckpointManager, ExecutionProgress, RuntimeParameters

manager = CheckpointManager("runs")
checkpoint = Checkpoint(
		base_system=BaseSystemState(arrays={"matrix": np.eye(3)}),
		runtime_parameters=RuntimeParameters(values={"seed": 7}),
		progress=ExecutionProgress(step=10, arrays={"state": np.ones(3)}),
)

manager.save_checkpoint(checkpoint, "step-10.npz")
branched = manager.load_for_branch(
		"step-10.npz",
		runtime_parameters=RuntimeParameters(values={"seed": 99}),
)
```

For a complete run-save-restart-branch workflow, see [examples/checkpoint_branching.py](examples/checkpoint_branching.py):

```bash
python -m examples.checkpoint_branching --checkpoint /tmp/quantum_macaroni_branching_example.npz
```

## Project Layout

```text
quantum_macaroni/
	calculators/      transport calculators and registry
	checkpointing/    decoupled checkpoint state and persistence
	core/             constants and numerics
	interpolation/    SKW interpolator
	mesh/             tetrahedron mesh
	parsers/          parser interfaces and implementations
examples/
	PbTe-nospin/      sample input and usage script
main.py             CLI entry point
```
