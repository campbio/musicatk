# 1. Record architecture decisions

- **Status:** Accepted
- **Date:** 2026-09-18
- **Deciders:** Campbell Lab maintainers

## Context

musicatk has accumulated structural decisions whose rationale exists only in
individual memory or in pull-request comments. The clearest example is the
v2.0.0 change that folded the `musica_result` class into `musica`: the reasoning
is not recoverable from the code or from `git log`, so it gets re-litigated.

The problem sharpens with AI agents in the loop. An agent reads the current state
of the code and cannot distinguish a deliberate constraint from an accident, so
it may "fix" something that was decided on purpose — and do it confidently.

## Decision

We record architecturally significant decisions as Architecture Decision Records
in `dev/adr/`, using the MADR-style template in `dev/adr/template.md`. The
criteria for "architecturally significant", the proposal process, and the
append-only rule are documented in `dev/adr/README.md`.

Records live in `dev/` — never in `docs/`, which is pkgdown build output and is
regenerated wholesale.

## Consequences

- Structural changes now require an ADR before implementation. This slows them
  down deliberately.
- `AGENTS.md` points agents at `dev/adr/`, so the decision history is
  discoverable by both humans and tools.
- The index table in `dev/adr/README.md` must be updated with each new ADR.
  Nothing enforces this automatically.
- ADR 0001 is self-justifying: the decision to keep records is itself the first
  record.

## Alternatives considered

- **Keep decisions in GitHub issues and PR discussions.** This is the status quo
  and is what failed — discussion threads are not indexed, not versioned with the
  code, and are effectively invisible to agents reading a checkout.
- **A single DECISIONS.md changelog.** Simpler, but merge conflicts on a
  shared file and no natural place for the context/alternatives of any one
  decision. One file per decision keeps the diff clean.
- **Store ADRs under `docs/adr/`.** Rejected: `docs/` is pkgdown output for this
  package, and `pkgdown::clean_site()` would delete them.
