# musicatk Release Checklist

Bioconductor releases twice a year, roughly **April** and **October**. The
package freeze lands a few weeks before each; check the official schedule at
https://bioconductor.org/developers/release-schedule/ each cycle and note the
current freeze date here.

## Version numbering

Bioconductor uses an even/odd `x.y.z` scheme:

- **devel** branch carries an **odd** `y` (e.g. 2.3.z)
- **release** branch carries an **even** `y` (e.g. 2.2.z)
- `z` increments with each change pushed to a branch; it resets on a `y` bump.

The Bioconductor build system performs the `y` bumps at release time. You bump
`z` yourself with every change that goes to the Bioconductor git server.

## Before the freeze

- [ ] `make check-full` clean — no ERRORs, no WARNINGs, NOTEs triaged.
- [ ] `make bioccheck` clean — this runs both `BiocCheck()` on the built tarball
      and `BiocCheckGitClone()` on the checkout.
- [ ] Triage every remaining WARNING/NOTE into a GitHub issue with a fix plan.
      Fix what you can, re-run, and PR the fixes to `devel`.
- [ ] `make test` passes; coverage has not regressed.
- [ ] `NEWS.md` has a section for this version listing user-facing changes.
- [ ] `make docs` run; `man/` and `NAMESPACE` are current.
- [ ] `make site-check` passes — every exported function appears in the
      `_pkgdown.yml` reference index.
- [ ] `/security-review` run on the release branch.
- [ ] Tarball size checked: `make bioccheck` prints it. BiocCheck requires
      **under 10 MB**, with **no single file over 5 MB**; musicatk has
      historically carried avoidable payload in `vignettes/` (see
      dev/ROADMAP.md, technical debt).

## Syncing to Bioconductor

The `bioc` remote must point at the Bioconductor git server:

```bash
git remote set-url bioc git@git.bioconductor.org:packages/musicatk.git
git remote -v   # verify before pushing
```

Then:

```bash
git fetch bioc
git checkout devel && git merge bioc/devel     # reconcile
make check-full && make bioccheck              # re-verify after the merge
git push bioc devel
```

## After the release

- [ ] Confirm the package appears in the new release on bioconductor.org.
- [ ] Check the nightly build report —
      https://bioconductor.org/checkResults/ — for both release and devel.
      This is the authoritative status; CI BiocCheck is only an early warning.
- [ ] `make site-deploy` to publish the updated pkgdown site.
- [ ] Open the next cycle's tracking issue.
