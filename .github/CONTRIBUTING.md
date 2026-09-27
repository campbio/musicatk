# Contributing to musicatk

Thanks for contributing. This guide covers both human contributors and AI coding
agents; the authoritative instructions for agents are in [AGENTS.md](../AGENTS.md).

## Setup

```r
BiocManager::install("musicatk", dependencies = TRUE)
```

musicatk requires R >= 4.4.0 and depends on large BSgenome/TxDb annotation
packages. The first install is slow. There is no compilation step.

## Repository layout

Alongside the standard R package directories, the repository carries a set of
files supporting AI-assisted development. Each location has one job:

```
musicatk/
├── AGENTS.md          canonical agent instructions — the file to edit
├── CLAUDE.md          one line: @AGENTS.md
├── GEMINI.md          one line: @AGENTS.md
├── Makefile           canonical commands (make help)
├── .lintr             lint rules; covers inst/shiny
├── .claude/
│   └── settings.json  what agents may and may not run (shared, committed)
├── .github/
│   ├── CONTRIBUTING.md  SECURITY.md  PULL_REQUEST_TEMPLATE.md
│   └── workflows/       CI
├── dev/               maintainer-facing, not shipped in the package
│   ├── adr/           architecture decision records (README.md = index)
│   ├── RELEASE.md     Bioconductor release checklist
│   ├── ROADMAP.md     where the package is going
│   ├── AUDIT.md       periodic dependency/deprecation audit
│   ├── agent-log.md   findings awaiting triage into issues
│   └── hooks/         load-standards.sh (session start), lint-changed.sh
└── R/ man/ tests/ vignettes/ inst/ docs/     standard package structure
```

Three rules explain the placement:

- **Root** holds only what agent harnesses auto-load, plus the `Makefile`.
  `AGENTS.md` is the single source of truth and is self-contained — `CLAUDE.md`
  and `GEMINI.md` just import it, so there is one file to keep current.
- **`.github/`** holds the community health files GitHub surfaces on its own.
- **`dev/`** holds maintainer documentation. One `.Rbuildignore` entry (`^dev$`)
  keeps all of it out of the package tarball.

`docs/` is pkgdown build output. Nothing is ever stored there by hand — that is
why decision records live in `dev/adr/` rather than `docs/adr/`.

Instructions tell agents what to know; `.claude/settings.json` controls what they
are permitted to do. Anything that must never happen belongs in the settings file,
not in prose.

## Branch and PR model

Campbell Lab members with write access branch directly on `campbio/musicatk`.
External contributors fork and branch on their fork.

`master` is the GitHub default branch, but it is an automatic copy of the
current Bioconductor release, kept in sync by a GitHub Action. Never commit to
it, branch from it, or open a PR against it.

1. Branch from `devel`. Never commit directly to `devel` or `master`.
2. Make your change. Run `make test` after every edit.
3. Run `make check-full` and `make bioccheck` before opening the PR.
4. Add a `NEWS.md` entry for any user-facing change.
5. Open the PR against `devel` (not the suggested `master`) and fill in the
   template.

## Canonical commands

Use the Makefile. It is the single answer to "how do I run this", for humans,
agents, and CI alike.

| Command | Purpose |
|---|---|
| `make test` | fast testthat loop — after every change |
| `make check` | quick `R CMD check` (no vignettes) — any time |
| `make check-full` | full `R CMD check` — before every PR |
| `make bioccheck` | BiocCheck on the built tarball — before every PR |
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
