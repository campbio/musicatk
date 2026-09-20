# Security Policy

## Reporting a vulnerability

Report suspected vulnerabilities to the package maintainer, Joshua Campbell
(camp@bu.edu), rather than opening a public issue. Include the affected version,
reproduction steps, and impact. You can expect an acknowledgement within a week.

For issues in Bioconductor infrastructure rather than this package, contact the
Bioconductor core team via https://bioconductor.org/about/.

## Scope

musicatk is an analysis library. It reads genomic variant files (VCF, MAF),
downloads annotation from Bioconductor AnnotationHub/BSgenome, and runs a local
Shiny app. The security-relevant surface is therefore:

- Parsing untrusted VCF/MAF input.
- File paths supplied to `extract_variants_*` and the Shiny upload handlers.
- The Shiny app when run on a shared or networked host rather than locally.

## Rules for automated tools and AI agents

These apply to any agent operating in this repository:

- Never read, log, echo, or commit credentials, tokens, API keys, `.Renviron`,
  `.Rprofile` secrets, or `~/.ssh` material.
- Never commit absolute local paths or machine-specific configuration.
- Never transmit patient, clinical, or unpublished genomic data to any external
  service. Test fixtures in `inst/extdata/` are public TCGA samples and are the
  only data that may appear in the repository.
- Never add a network call to package code without an ADR.
- `/security-review` is run before each Bioconductor release (see dev/RELEASE.md).
