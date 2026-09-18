# Architecture Decision Records

Decisions about musicatk that are hard to reverse, or that a future contributor
would otherwise question and re-litigate, are recorded here as numbered ADRs.

## When an ADR is required

Write one before implementing any change that:

- alters the object model (S4 class design, slots, accessor conventions);
- adds, removes, or changes a DESCRIPTION dependency;
- splits, merges, or reorganizes source files;
- changes the public API or breaks backward compatibility;
- changes the build, test, release, or deployment machinery.

**Not** for routine choices — a new plotting argument, a bug fix, a refactor
inside one function, or anything a reviewer could evaluate from the diff alone.
The test is whether someone six months from now would ask "why is it like this?"

## Process

1. Propose the decision in a GitHub issue.
2. Draft the ADR from `template.md` with status **proposed**. The `adr-author`
   skill helps here.
3. The maintainer approves; status becomes **accepted** and the ADR is numbered.
4. Implement.

ADRs are append-only. An accepted ADR is never edited or deleted — a later ADR
supersedes it, and the old one's status changes to **superseded by NNNN** with
both kept in place. The record of what we used to think is the point.

## Index

| # | Title | Status | Date |
|---|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted | 2026-09-18 |
