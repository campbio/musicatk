# Contributing to musicatk

Thanks for contributing. This guide covers both human contributors and AI coding
agents; the authoritative instructions for agents are in [AGENTS.md](AGENTS.md).

## Setup

```r
BiocManager::install("musicatk", dependencies = TRUE)
```

musicatk requires R >= 4.4.0 and depends on large BSgenome/TxDb annotation
packages. The first install is slow. There is no compilation step.

## Branch and PR model

Campbell Lab members with write access branch directly on `campbio/musicatk`.
External contributors fork and branch on their fork.

1. Branch from `devel`. Never commit directly to `devel` or `master`.
2. Make your change. Run `make test` after every edit.
3. Run `make check` before opening the PR.
4. Add a `NEWS.md` entry for any user-facing change.
5. Open the PR against `devel` and fill in the template.

## Canonical commands

Use the Makefile. It is the single answer to "how do I run this", for humans,
agents, and CI alike.

| Command | Purpose |
|---|---|
| `make test` | fast testthat loop — after every change |
| `make check` | full `R CMD check` — before every PR |
| `make bioccheck` | BiocCheck on the built tarball |
| `make docs` | regenerate `man/` + `NAMESPACE` from roxygen |
| `make lint` | lintr across `R/` and `inst/shiny` |
| `make site-check` | verify the pkgdown reference index |
| `make app` | launch the Shiny app |

Run `make help` for the full list.

## Things that are generated — never hand-edit

- `man/` and `NAMESPACE` — owned by roxygen2. Edit the roxygen comments, then
  `make docs`.
- `docs/` — pkgdown output. Edit roxygen comments, vignettes, or `_pkgdown.yml`.
  Nothing hand-written is ever stored there — maintainer docs live in `dev/`.

## Architectural changes

Structural refactors — splitting files, changing DESCRIPTION dependencies,
redesigning S4 classes — need an approved ADR before implementation. Open a
GitHub issue proposing it, then record the outcome in `dev/adr/` using
`dev/adr/template.md`.

## AI agent tooling

Lab members should install the shared skill sets once, at user level:

1. [Bioconductor official skills](https://github.com/Bioconductor/ai-agent-skills)
   — `bioc-pkg-dev`, `build-check-bioccheck`, `analyze-r-package`,
   `improve-code-coverage`, `security-audit-r-package`, `update-r-news`,
   `adr-author`, `bioc-howto`. Reference from `~/.claude/CLAUDE.md` per that
   repo's `instructions/claude.md`.
2. [Posit r-lib plugin](https://github.com/posit-dev/skills) —
   `/plugin marketplace add posit-dev/skills` then
   `/plugin install r-lib@posit-dev-skills`.
3. [Sean Davis's Bioconductor skills](https://github.com/seandavi/ai-agent-skills)
   — complementary workflows; keep what adds value over #1.

Workflow habits that need no install: plan mode for non-trivial changes,
`/code-review` before PRs, `/security-review` before releases, `/simplify`
after features land.
