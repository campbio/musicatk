## Summary

<!-- What changes and why. Link the issue if there is one. -->

## Checklist

- [ ] `make test` passes
- [ ] `make check` passes (no new WARNINGs or NOTEs)
- [ ] `make lint` clean for the files I touched
- [ ] `NEWS.md` updated, if this is a user-facing change
- [ ] `make docs` run, if I changed roxygen comments
- [ ] New exported functions added to the `_pkgdown.yml` reference index
      (`make site-check` passes)
- [ ] `/code-review` run on this branch
- [ ] ADR linked below if this is an architectural or dependency change
- [ ] Shiny only: verified with a screenshot of the running app (`make app`)

## ADR

<!-- Link to dev/adr/NNNN-*.md, or "N/A — not architectural".
     See dev/adr/README.md for when one is required. -->

## Scientific correctness

<!-- REQUIRES HUMAN JUDGMENT — do not let an agent tick this for you.
     If this PR changes numerical results, signature inference, exposure
     prediction, or plotting semantics, state who verified the output is
     scientifically correct and how. Passing tests is not sufficient. -->

## Generated content

<!-- If an AI agent wrote part of this PR, say which parts, so reviewers can
     weight their attention accordingly. -->
