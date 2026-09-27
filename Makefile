# Standard targets expected by r-bioc-dev-standards, plus musicatk extras.
# Run `make help` to list them. Recipe lines must start with a tab.
# Which extra targets agents may run is listed in AGENTS.md.

.PHONY: help docs test test-one check check-full bioccheck lint coverage \
	site-check build style site site-deploy app test-app clean

help:  ## List targets
	@grep -E '^[a-zA-Z_-]+:.*## ' $(MAKEFILE_LIST) | \
	  awk 'BEGIN {FS = ":.*## "}; {printf "  %-12s %s\n", $$1, $$2}'

docs:  ## Regenerate man/*.Rd and NAMESPACE from roxygen comments
	Rscript -e 'devtools::document()'

test:  ## Run the full test suite
	Rscript -e 'devtools::test(stop_on_failure = TRUE)'

test-one:  ## Run matching test files: make test-one FILTER=<pattern>
	@test -n "$(FILTER)" || { echo "Usage: make test-one FILTER=<pattern>"; exit 1; }
	Rscript -e 'devtools::test(filter = "$(FILTER)", stop_on_failure = TRUE)'

check:  ## Quick R CMD check: skips vignettes and the PDF manual
	Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--ignore-vignettes"), build_args = "--no-build-vignettes", error_on = "warning", check_dir = tempdir())'

check-full:  ## Full check: rebuilds vignettes, runs \donttest examples; needs all Suggests
	Rscript -e 'devtools::check(document = FALSE, vignettes = TRUE, run_dont_test = TRUE, force_suggests = TRUE, error_on = "warning", check_dir = tempdir())'

bioccheck:  ## BiocCheckGitClone on the repo, then BiocCheck on a built tarball
	Rscript -e 'BiocCheck::BiocCheckGitClone(".", "quit-with-status" = TRUE)'
	@tmp="$$(mktemp -d)"; \
	  cd "$$tmp" && R CMD build "$(CURDIR)" && \
	  Rscript -e 'tb <- list.files(pattern = "[.]tar[.]gz$$")[1]; message(sprintf("Tarball %s: %.2f MB (limit 10 MB)", tb, file.size(tb) / 1e6)); BiocCheck::BiocCheck(tb, "quit-with-status" = TRUE)'; \
	  status=$$?; echo "BiocCheck output kept in $$tmp"; exit $$status

lint:  ## Run lintr on the package (reports only; changes nothing)
	Rscript -e 'print(lintr::lint_package())'

coverage:  ## Print test coverage, overall and per file
	Rscript -e 'print(covr::package_coverage())'

site-check:  ## Check the pkgdown reference index lists every export (no site build)
	Rscript -e 'pkgdown::check_pkgdown()'

# ---- musicatk extras --------------------------------------------------------

build:  ## Build the source tarball in the repo root
	R CMD build .

style:  ## Apply styler to the whole package (people only; not in CI)
	Rscript -e 'styler::style_pkg()'

site:  ## Full local pkgdown build (people only; slow)
	Rscript -e 'pkgdown::build_site()'

site-deploy:  ## Maintainer action: build locally, push to gh-pages
	Rscript -e 'pkgdown::deploy_to_branch()'

app:  ## Launch the Shiny app for visual verification
	Rscript -e 'shiny::runApp(system.file("shiny", package = "musicatk"))'

test-app:  ## shinytest2 smoke suite (placeholder until tests/app/ exists)
	@if [ -d tests/app ]; then \
		Rscript -e 'testthat::test_dir("tests/app")'; \
		else \
		echo "tests/app/ not present yet; see AGENTS.md (Tests)."; \
		fi

clean:  ## Remove build artifacts (people only)
	rm -f musicatk_*.tar.gz
	rm -rf musicatk.Rcheck
