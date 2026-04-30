## Why

The project lacks an `AGENTS.md` file that documents repository conventions for AI coding agents and human contributors. Without it, agents have no reference for code style, project structure, security practices, commit conventions, or automation rules—leading to inconsistent code generation, wrong formatting choices, and missed lint/typecheck steps.

## What Changes

- Add a new `AGENTS.md` file at the repository root with five sections:
  1. **Project Structure & Module Organization** — package layout, subpackage roles, public API surface
  2. **Security & Configuration Tips** — secrets handling, environment files, dependency safety
  3. **Automation & Agent Workflow** — lint/format/typecheck commands, test runner, CI expectations
  4. **Commit & Pull Request Guidelines** — message format, branch naming, PR process
  5. **Coding Style & Naming Conventions** — ruff configuration summary, docstring style, type hints, NumPy conventions

## Capabilities

### New Capabilities
- `agents-md`: Repository guidelines document for AI agents and contributors covering structure, security, automation, commit/PR, and coding style

### Modified Capabilities

## Impact

- Adds one new documentation file (`AGENTS.md`) at the repository root
- No code changes, no API changes, no dependency changes
- Improves onboarding for automated agents and human contributors