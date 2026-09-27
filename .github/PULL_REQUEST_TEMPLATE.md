<!-- Base branch: devel. Use RELEASE_X_Y only for an approved release fix.
     Never target main/master, which is updated automatically. -->

## What changed and why

<!-- Link the issue if there is one. -->

## How it was tested

## Checklist

- [ ] Tests added or updated, `make test` passes, and `make coverage`
      didn't drop
- [ ] `make check-full` and `make bioccheck` pass with no new errors or warnings
- [ ] `make lint` clean for the files I touched
- [ ] `make docs` run, if roxygen comments changed
- [ ] New exports added to `_pkgdown.yml`, and `make site-check` passes
- [ ] NEWS.md updated for user-facing changes
- [ ] Version bumped (z) if this will be pushed to Bioconductor
- [ ] Plan review and `/code-review` run; findings fixed or answered
- [ ] Related issue linked
- [ ] Shiny only: verified with a screenshot of the running app (`make app`)

## ADR

<!-- Link to dev/adr/NNNN-*.md, or "N/A — not architectural".
     See dev/adr/README.md for when one is required. -->

## Scientific correctness

<!-- REQUIRES HUMAN JUDGMENT — do not let an agent tick this for you.
     If this PR changes numerical results, signature inference, exposure
     prediction, or plotting semantics, state who verified the output is
     scientifically correct and how. Passing tests is not sufficient.
     Otherwise write "No change to results". -->

## Generated content

<!-- If an AI agent wrote part of this PR, say which parts, so reviewers can
     weight their attention accordingly. -->
