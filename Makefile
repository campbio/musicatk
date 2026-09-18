# Canonical entry points for musicatk development.
# Humans, agents, and CI all use these targets. Non-obvious flags live here.

.PHONY: help test check bioccheck docs lint style site site-deploy app test-app clean

help:  ## list targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-14s %s\n", $$1, $$2}'

test:  ## fast loop — run after every change
	Rscript -e 'devtools::test()'

check:  ## full check — run before opening a PR
	_R_CHECK_FORCE_SUGGESTS_=false Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual"), error_on = "warning")'

build:  ## build the source tarball
	R CMD build .

bioccheck:  ## build tarball first — matches how the Bioc build system runs it
	R CMD build . && Rscript -e 'BiocCheck::BiocCheck(Sys.glob("*.tar.gz")); BiocCheck::BiocCheckGitClone(".")'

docs:  ## regenerate man/ and NAMESPACE from roxygen comments
	Rscript -e 'devtools::document()'

lint:  ## lintr over R/ and inst/shiny
	Rscript -e 'lintr::lint_package()'

style:  ## apply styler (does not run in CI; run before committing)
	Rscript -e 'styler::style_pkg()'

site:  ## local full pkgdown build (slow — never use as routine verification)
	Rscript -e 'pkgdown::build_site()'

site-deploy:  ## maintainer action: build locally, push to gh-pages
	Rscript -e 'pkgdown::deploy_to_branch()'

site-check:  ## cheap structural check of the pkgdown reference index
	Rscript -e 'pkgdown::check_pkgdown()'

app:  ## launch the Shiny app for visual verification
	Rscript -e 'shiny::runApp(system.file("shiny", package = "musicatk"))'

test-app:  ## shinytest2 smoke suite (PLACEHOLDER — tests/app/ does not exist yet)
	@if [ -d tests/app ]; then \
		Rscript -e 'testthat::test_dir("tests/app")'; \
		else \
		echo "tests/app/ not present yet — see AGENTS.md (package-specific notes)."; \
		fi

clean:  ## remove build artifacts
	rm -f musicatk_*.tar.gz
	rm -rf musicatk.Rcheck
