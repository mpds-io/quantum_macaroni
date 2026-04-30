## Context

`quantum_macaroni` is a Python package for Boltzmann transport calculations from DFT output. It has five subpackages (`core`, `parsers`, `mesh`, `interpolation`, `calculators`) with plugin registries for parsers and calculators. The project uses `ruff` for linting/formatting (configured in `ruff.toml`), `uv` for dependency management, and `pytest` for testing. There is currently no onboarding document describing conventions for agents or contributors.

## Goals / Non-Goals

**Goals:**
- Provide a single, authoritative reference at `AGENTS.md` that AI agents and contributors can read to produce code consistent with the project's conventions
- Cover all five required sections: structure, security, automation, commit/PR, and coding style
- Derive all guidelines from actual codebase patterns and `ruff.toml` configuration so the document is grounded in reality

**Non-Goals:**
- Replacing inline code documentation or docstrings
- Establishing new conventions beyond what the codebase already follows
- Adding CI/CD pipeline configuration (only document existing commands)

## Decisions

1. **Single file at repository root** — `AGENTS.md` lives at the top level so agents discover it immediately. Alternatives: per-directory `AGENTS.md` (too scattered) or inside `docs/` (not discovered). Rationale: root-level placement is the de facto standard for agent instruction files.

2. **Sections mirror the user's requested structure** — The five sections (Project Structure, Security, Automation, Commit & PR, Coding Style) map directly to the request. No additional sections are added to keep the document focused.

3. **Guidelines derived from existing `ruff.toml`** — All lint/format settings (line length 120, Python 3.11 target, rule selection, ignored rules, quote style, indent style) are reflected verbatim so agents don't need to read the config file separately.

4. **Registry pattern documented** — The parser and calculator plugin registries (`register_parser`/`get_parser`, `register_calculator`/`get_calculator`) are a key architectural pattern that agents must follow when adding new parsers or calculators.

## Risks / Trade-offs

- [Risk: AGENTS.md becomes stale] → Mitigation: reference `ruff.toml` and `pyproject.toml` as source-of-truth; note that agents should re-read those files if conventions change
- [Risk: Document is too long to fit in agent context window] → Mitigation: keep each section concise and actionable, avoid duplicating information available elsewhere