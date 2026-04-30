# AGENTS.md

Guidelines for AI coding agents and contributors working on `quantum_macaroni`.

## 1. Project Structure & Module Organization

```
quantum_macaroni/
├── __init__.py              # Public API re-exports
├── core/
│   ├── constants.py         # Physical constants & conversion factors
│   ├── numerics.py          # Numba-accelerated kernels
│   └── symmetry.py         # Point-group construction utilities
├── parsers/
│   ├── base.py              # Parser protocol, ParserResult, registry
│   └── fleur_outxml.py      # FLEUR out.xml parser plugin
├── mesh/
│   └── tetrahedron.py       # k-space mesh & tetrahedral decomposition
├── interpolation/
│   └── skw.py               # SKW band-structure interpolation
└── calculators/
    ├── base.py              # Calculator protocol, registry, tensor_average
    └── transport.py         # Boltzmann transport calculator
```

### Subpackage roles

| Subpackage | Purpose | Key exports |
|---|---|---|
| `core` | Low-level constants, numerics, symmetry | `BOHR_TO_ANG`, `HBAR`, `E_CHARGE`, `KB_EV`, `HTR_TO_EV`, `TWO_PI` |
| `parsers` | Parse DFT output into normalized `ParserResult` | `ElectronicStructureParser`, `ParserResult`, `FleurOutxmlParser`, `DEFAULT_PARSER` |
| `mesh` | k-space meshing for Brillouin-zone integration | `TetrahedronMesh` |
| `interpolation` | Band-structure interpolation algorithms | `SKWInterpolator` |
| `calculators` | Transport tensor computation | `BoltzmannTransportCalculator`, `calculate_spin_polarized_transport` |

### Plugin registries

Parsers and calculators use a registration pattern. To add a new parser:

```python
from quantum_macaroni.parsers.base import ElectronicStructureParser, register_parser

class MyParser:
    name = "my-parser"

    def parse(self, filepath, iteration="last"):
        ...

register_parser(MyParser())
```

To add a new calculator:

```python
from quantum_macaroni.calculators.base import TransportCalculator, register_calculator

class MyCalculator:
    name = "my-calc"

    def calculate_transport(self, fermi_level, temperature, kpoint_mesh, kchunk=None):
        ...

register_calculator(MyCalculator.name, MyCalculator)
```

Entry points: `quantum_macaroni.__init__` re-exports the public API. New plugins should also be wired into their subpackage `__init__.py`.

## 2. Security & Configuration Tips

- Never commit secrets, API keys, tokens, or credentials to the repository
- Environment files (`.env`, `.envrc`) and virtual environments (`.venv`) are excluded via `.gitignore`
- `.pypirc` (PyPI credentials) is also excluded
- Do not hard-code file paths or system-specific configuration; accept them as parameters
- Dependencies are managed through `pyproject.toml` and locked via `uv.lock` — never install packages ad-hoc
- When adding dependencies, use `uv add <package>` and verify the lockfile updates

## 3. Automation & Agent Workflow

Always run these commands after making code changes:

```bash
ruff check .            # Lint all Python files
ruff format .           # Auto-format all Python files
ty check                # Type-check the project
pytest                  # Run the test suite
```

Run `ruff check --fix .` to auto-fix lint issues where possible. Run `ruff format --check .` to verify formatting without writing.

Recommended workflow:
1. Make code changes
2. Run `ruff check .` and fix any violations
3. Run `ruff format .` to ensure consistent formatting
4. Run `ty check` to verify type correctness
5. Run `pytest` to confirm no regressions
6. Commit only after all checks pass

## 4. Commit & Pull Request Guidelines

### Commit messages

Use conventional commit format:

```
<type>(<scope>): <description>

[optional body]
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`

Examples:
```
feat(parsers): add VASP parser plugin
fix(numerics): fix tetrahedron weight edge case at equal energies
docs(core): document physical constant conventions
chore: update dependencies
```

### Branch naming

Use kebab-case with a type prefix: `<type>/<short-description>`
- `feat/add-vasp-parser`
- `fix/tetra-weight-bug`
- `docs/update-agents-md`

### Pull requests

- Include a summary of changes and motivation
- Reference related issues when applicable
- Ensure all lint, type-check, and test commands pass before requesting review
- Keep PRs focused on a single concern

## 5. Coding Style & Naming Conventions

### ruff configuration (source: `ruff.toml`)

| Setting | Value |
|---|---|
| Line length | 120 |
| Target version | Python 3.11 (`py311`) |
| Indent width | 4 spaces |
| Quote style | Double quotes |
| Indent style | Spaces |

**Active lint rules:** `E` (pycodestyle errors), `F` (pyflakes), `B` (flake8-bugbear), `SIM` (flake8-simplify), `I` (isort), `C901` (complexity), `NPY` (NumPy), `N` (pep8-naming), `PL` (Pylint), `D` (pydocstyle)

**Ignored rules:** `E501` (line length enforced by formatter), `PLR0913` (too many arguments), `PLR0914` (too many locals), `PLR0915` (too many statements), `PLR0917` (too many positional arguments)

**Dummy variable pattern:** `^(_+|(_+[a-zA-Z0-9_]*[a-zA-Z0-9]+?)$` — underscore-prefixed variables are allowed as unused.

### Docstrings

Use Google-style docstrings with `Args`, `Returns`, and `Raises` sections:

```python
def calculate_transport(self, fermi_level: float, temperature: float) -> npt.NDArray[np.float64]:
    """Return transport tensors for the given state point.

    Args:
        fermi_level: Fermi level in eV.
        temperature: Temperature in kelvin.

    Returns:
        Electrical conductivity tensor with shape (3, 3).

    Raises:
        ValueError: If temperature is non-positive.

    """
```

Module-level docstrings are single-line summaries.

### Type hints

- Use `numpy.typing` aliases: `npt.NDArray[np.float64]`, `npt.NDArray[np.int_]`, `npt.ArrayLike`
- Use `|` union syntax (Python 3.10+): `str | Path`, `int | None`
- Use `tuple[...]` generic syntax (not `Tuple` from typing), except when importing `Tuple` is needed for legacy compatibility
- Prefer concrete NumPy dtypes in annotations: `np.float64`, `np.int32`, `np.complex128`

### NumPy / SciPy / Numba conventions

- Numba kernels are prefixed with `nb_` and decorated with `@_nb.njit(cache=True, fastmath=True)` or `@_nb.njit(parallel=True, cache=True, fastmath=True)`
- Use `_nb.prange` for parallel loops in Numba; add `# ty:ignore[not-iterable]` on the same line
- Suppress complexity violations on Numba kernels with `# noqa: C901` on the `def` line
- Import `numba` as `_nb` to indicate it is an internal acceleration detail
- Use `np.ascontiguousarray` before passing arrays to Numba functions
- Physical constants are centralized in `core/constants.py` — never redefine them locally