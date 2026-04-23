## 1. Create AGENTS.md with Project Structure & Security sections

- [x] 1.1 Create `AGENTS.md` at repository root with the "Project Structure & Module Organization" section documenting `quantum_macaroni/` layout, subpackage roles (`core`, `parsers`, `mesh`, `interpolation`, `calculators`), public API exports, and plugin registry pattern
- [x] 1.2 Add "Security & Configuration Tips" section covering secrets handling (never commit keys/credentials), `.env`/`.venv`/`.pypirc` exclusion via `.gitignore`, and dependency safety

## 2. Add Automation, Commit/PR, and Coding Style sections

- [x] 2.1 Add "Automation & Agent Workflow" section documenting `ruff check .`, `ruff format .`, `ty check`, and `pytest` commands
- [x] 2.2 Add "Commit & Pull Request Guidelines" section with conventional commit format, branch naming (kebab-case), and PR description expectations
- [x] 2.3 Add "Coding Style & Naming Conventions" section summarizing `ruff.toml` settings (line-length 120, py311 target, active rules E/F/B/SIM/I/C901/NPY/N/PL/D, ignored rules E501/PLR0913/PLR0914/PLR0915/PLR0917, double quotes, space indents), docstring conventions (Google-style with Args/Returns/Raises), type hints (`npt.NDArray`, `npt.ArrayLike`, `|` union syntax), and NumPy/SciPy/numba patterns

## 3. Verify

- [x] 3.1 Confirm all five required sections exist in `AGENTS.md`
- [x] 3.2 Run `ruff check AGENTS.md` or verify file is excluded from lint and contains no broken markdown