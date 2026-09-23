# Quantum Macaroni ⚛️🍝

This toolbox is a modular Boltzmann-transport workflow for post-processing the electronic-structure data.

Current pipeline:
- Parser plugins for FLEUR out.xml and CRYSTAL outp electronic-structure outputs (default: FLEUR).
- SKW interpolation of band energies.
- Tetrahedron k-space integration mesh.
- Transport-property calculators (default: Boltzmann transport calculator).


### Features

- Plugin architecture for parsers and calculators.
- Transport tensors and isotropic averages.
- Temperature sweep with scalar or array input.
- Chemical-potential sweep relative to Fermi level.
- Restartable checkpoints for parsed input, fitted interpolation, transport DOS, and completed scans.
- CLI interface with JSON output.


### Requirements

- Python 3.11+:
- ase
- lxml
- numba
- numpy
- scipy


## Installation

Use your preferred environment manager. An example with pip reads:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```


## Quick Start (CLI)

Main entry point is [main.py](main.py). Minimal run:

```bash
python main.py examples/PbTe-nospin/out-nospin.xml
```

CRYSTAL properties output is detected automatically:

```bash
python main.py examples/outp/mno2afm.outp
```

Run transport with an accompanying CRYSTAL DOS file:

```bash
python main.py system.outp --dos-file system.DOSS --output results.json
# Equivalent positional form:
python main.py system.outp system.DOSS --output results.json
```

The outp file supplies bands, lattice, and symmetry. By default, the DOSS file
supplies the Fermi energy for band selection and chemical-potential shifts.
Use `--fermi-source outp` to use the outp reference instead. The output includes
the parsed DOS under `dos`, alongside transport results. Both files must describe
the same system with the same energy convention; their spin-channel counts are
checked. Changing the selected Fermi energy invalidates affected checkpoint stages.

The properties file must include direct lattice vectors (`COORPRT`), Cartesian
symmetry matrices (`SYMMOPS`), a k-point coordinate table, eigenvalues, and a numeric
Fermi energy. Alpha and beta eigenvalues are kept as two spin channels, including
AFM systems with zero total spin. Non-spin-polarized output is read as one channel;
the transport calculator supplies its factor of two for spin degeneracy. The parser
converts atomic-unit energies to eV and uses the printed Fermi energy as the reference
for chemical-potential shifts. It accepts one properties dataset per file.

E and T symmetry labels are handled as two- and three-state multiplets. When
the file prints one energy per multiplet, the parser expands it into two or three
band entries before interpolation. When each component is already printed,
it retains those entries without multiplying them again. `NUMBER OF AO` resolves
the encoding; otherwise it must be unambiguous from band counts across k-points.
Ambiguous or inconsistent counts raise an error instead of guessing. This orbital
multiplicity is separate from spin degeneracy.

The same parser is available in Python as `CrystalOutpParser()` or through
`calculate_spin_polarized_transport(path, parser="crystal-outp")`.

The supplied AFM example prints a Fermi energy of -1 Hartree (about -27.21 eV)
and reports `SPIN LOCKING: NO ENERGY GAP COMPUTED`. This value is preserved;
choose chemical-potential shifts using your intended transport reference. Empty
transport windows produce zero tensors. For singular conductivity tensors,
Seebeck components in nonconducting directions are reported as zero by convention.

CRYSTAL electronic DOS files can also be exported directly from the CLI:

```bash
python main.py system.DOSS --output dos.json
```

`DOSS.DAT` and `fort.25` inputs are also detected automatically. DOS-only runs
default to `dos_results.json` and do not create transport checkpoints. The JSON
contains energies relative to the Fermi energy in eV, absolute energies in eV,
and DOS with shape `(spin, energy, projection)` in states/eV/cell.
Use `--parser crystal-doss` or `--parser crystal-outp` to select a format explicitly.

The DOS parser is also available in Python:

```python
from quantum_macaroni import CrystalDOSSParser

dos = CrystalDOSSParser().parse("DOSS.DAT")  # Also accepts .DOSS and fort.25 files
energy = dos.energies                       # E - E_F, in eV; shape (nenergy,)
alpha = dos.dos[0]                          # states/eV/cell; shape (nenergy, nprojections)
if dos.jspins == 2:
    beta = dos.dos[1]
absolute_energy = dos.absolute_energies
```

`dos.fermi_energy` stores the absolute Fermi energy in eV. Restricted DOSS files
return one channel; polarized files return alpha then beta. Projection columns
retain their original order, including any total-DOS column. The parser reverses
CRYSTAL's beta plotting sign and converts states/Hartree to states/eV without an
extra spin factor. It supports combined `fort.25` files containing both BAND and
DOSS blocks and validates matching energy grids across spins and projections.
The reader parses CRYSTAL's text records directly using Python and NumPy.
Electronic DOS returns a `DOSResult`, separate from the band parser registry: DOSS files do not provide
the band velocities required to calculate transport. Phonon DOS is unsupported.

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

```
python main.py FILEPATH [DOSS_FILEPATH] [options]

Options:
	--dos-file PATH, --doss PATH            accompanying CRYSTAL DOSS file for outp input
	--fermi-source {doss,outp}              combined-run Fermi-energy reference (default: doss)
	--temperature T [T ...]                one value or (start stop npoints)
	--chemical-potential MU [MU ...]       one value or (start stop npoints), eV shift from E_F
	--tau FLOAT                            relaxation time in seconds (default: 1e-14)
	--kmesh NX NY NZ                       k-point mesh (default: 80 80 80)
	--lr-ratio INT                         SKW interpolator star-vector ratio (default: 20)
	--band-window EMIN EMAX                band window relative to E_F in eV (default: -3 3)
	--chunk-size INT                       chunk size for batched evaluation (default: 4096)
	--parser {available_parsers,crystal-doss} input parser (default: automatic detection)
	--calculator {available_calculators}   calculator plugin (default: boltzmann)
	--output PATH                          JSON path (default: transport_results.json or dos_results.json)
	--checkpoint PATH                      restartable .npz checkpoint (default: transport_state.npz)
	--no-checkpoint                        disable checkpoint writing and resume
	--no-resume                            ignore an existing checkpoint while writing a fresh one
```

Available parser/calculator names come from the runtime registries in
`parsers/__init__.py` and `calculators/__init__.py`. The CLI additionally provides
the `crystal-doss` mode for DOS export.

## Checkpointing

The CLI saves restartable workflow state to `transport_state.npz` by default as the run advances. A later run with the same input and compatible configuration resumes automatically from the latest valid stage:

- parsed electronic structure
- fitted SKW interpolation
- transport density of states
- completed transport scan

Checkpoints are compressed NumPy archives with a JSON manifest and non-object NumPy arrays, written atomically to avoid partial files. Compatibility is checked with staged fingerprints over the input file content and calculation settings. Use `--checkpoint PATH` to choose a different file, `--no-resume` to ignore an existing checkpoint and replace it during the new run, or `--no-checkpoint` to disable checkpointing.

## Output JSON Format

The CLI stores calculation output to a JSON file (default: `transport_results.json`).

If a chemical potential is provided, its structure is:

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

NB: since JSON keys must be strings, the numeric keys for chemical potential and temperature are serialized.

## Python API

Public API is exported from `quantum_macaroni/__init__.py`.

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

### File Layout

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
main.py               CLI entry point
```
