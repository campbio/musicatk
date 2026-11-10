# 0002. Import withr to seed exposure prediction

- **Status:** Accepted
- **Date:** 2026-11-09
- **Deciders:** Joshua D. Campbell (maintainer), 2026-11-09

## Context

`predict_exposure(algorithm = "lda")` starts each sample by assigning its
mutations to signatures at random with `sample()`, so two runs on the same
data give different exposures (#152). `discover_signatures()` already takes
`seed = 1` and passes it to topicmodels or NMF, which seed internally.
Prediction has no such argument, and users are told in the vignette to wrap
calls in `withr::with_seed()` themselves.

The maintainer chose to give `predict_exposure()` and `auto_predict_grid()`
a `seed = 1` argument (2026-11-09, PR F of #208). Package code must then set
the RNG seed and restore the caller's RNG state afterwards. BiocCheck reports
any `set.seed()` call in `R/` as a WARNING ("Remove set.seed usage"), and
WARNINGs block a release, so `set.seed()` with a hand-written restore is not
an option.

withr is in Suggests today (tests, vignette). Moving it to Imports is a
DESCRIPTION dependency change, which needs an ADR. #204 (PR L) separately
reviews unused and single-use Imports.

## Decision

We will move withr from Suggests to Imports and use
`withr::with_seed(seed, ...)` wherever package code needs a seeded random
step, starting with LDA exposure prediction. A `seed = NULL` argument means
"use the caller's RNG state", as in `discover_signatures()`.

## Consequences

- Prediction is reproducible by default, and the caller's RNG state is left
  as it was, with no BiocCheck warning.
- No new install footprint: withr is already a hard dependency of ggplot2
  and testthat, both of which musicatk installs.
- One more Import to keep; #204 should list withr as used.
- Other stochastic functions (`create_umap()`, `compare_k_vals()`
  bootstraps) can adopt the same pattern later; that is follow-up work, not
  part of this decision.

## Alternatives considered

- **`set.seed()` plus a manual restore of `.Random.seed`.** Same behaviour,
  but BiocCheck reports a WARNING for every `set.seed()` call.
- **A deterministic starting assignment in `lda_posterior()`** (each
  sample's mutations split evenly across signatures). Removes the
  randomness without any dependency, but gives no `seed` argument and
  changes the algorithm's start; not the option the maintainer chose.
- **No seeding; document `withr::with_seed()` for users.** The current
  state, which leaves results irreproducible by default.
