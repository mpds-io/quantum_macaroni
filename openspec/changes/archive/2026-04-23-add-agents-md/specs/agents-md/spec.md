## ADDED Requirements

### Requirement: AGENTS.md exists at repository root
The repository SHALL contain an `AGENTS.md` file at the top level that documents conventions for AI agents and human contributors.

#### Scenario: File presence
- **WHEN** the repository is checked out
- **THEN** an `AGENTS.md` file exists at the repository root

### Requirement: Project Structure & Module Organization section
`AGENTS.md` SHALL include a "Project Structure & Module Organization" section that documents the package layout, the role of each subpackage (`core`, `parsers`, `mesh`, `interpolation`, `calculators`), and the public API surface exported from `quantum_macaroni.__init__`.

#### Scenario: Subpackage documentation
- **WHEN** a reader consults the Project Structure section
- **THEN** each of the five subpackages is described with its purpose and key exports

### Requirement: Security & Configuration Tips section
`AGENTS.md` SHALL include a "Security & Configuration Tips" section that covers secrets handling, `.env` exclusion, and `.gitignore` conventions.

#### Scenario: Secrets guidance
- **WHEN** a reader consults the Security section
- **THEN** the document instructs never to commit secrets, API keys, or credentials and references `.gitignore` entries for `.env`, `.venv`, and `.pypirc`

### Requirement: Automation & Agent Workflow section
`AGENTS.md` SHALL include an "Automation & Agent Workflow" section that documents the exact commands for linting, formatting, type checking, and testing.

#### Scenario: Lint and format commands
- **WHEN** a reader consults the Automation section
- **THEN** the exact `ruff check`, `ruff format`, `ty`, and `pytest` commands are listed

### Requirement: Commit & Pull Request Guidelines section
`AGENTS.md` SHALL include a "Commit & Pull Request Guidelines" section that specifies commit message format, branch naming conventions, and PR process.

#### Scenario: Commit message format
- **WHEN** a reader consults the Commit & PR section
- **THEN** a conventional commit message format is documented with examples

### Requirement: Coding Style & Naming Conventions section
`AGENTS.md` SHALL include a "Coding Style & Naming Conventions" section that summarizes `ruff.toml` settings (line length, target version, rule selection, ignored rules, quote style, indent style) and documents docstring conventions, type hint patterns, and NumPy/SciPy conventions used in the codebase.

#### Scenario: Ruff settings summary
- **WHEN** a reader consults the Coding Style section
- **THEN** all active ruff lint rules, ignored rules, line-length (120), target version (py311), quote style (double), and indent style (space) are listed