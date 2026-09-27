# Makefile for an r-bioc-dev-standards package. Run `make help` to list
# targets. Recipe lines must start with a tab.
#
# The standard targets (test, check, bioccheck, ...) live in the shared
# standards.mk in r-bioc-dev-standards, not here. It's downloaded to a cache
# on first use and refreshed at every Claude session start.

# Package settings. Uncomment to change a default (see standards.mk).
# FORCE_SUGGESTS = FALSE

# Only FILTER may be set on the make command line. Settings allow
# `make test-one` with any arguments, and a command-line variable could
# change which makefile, shell, or source is used. Environment variables
# still work, e.g. R_BIOC_STANDARDS_REF=<branch> make help.
ifneq ($(findstring $$,$(value FILTER)),)
  $(error FILTER may contain only letters, digits, '.', '_' and '-')
endif
ifneq ($(filter-out FILTER=%,$(MAKEOVERRIDES)),)
  $(error Only FILTER=<pattern> may be set on the make command line; set other variables in the environment or above)
endif

R_BIOC_STANDARDS_REF ?= v1
R_BIOC_STANDARDS_BASE ?= https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)
STANDARDS_MK := $(or $(XDG_CACHE_HOME),$(HOME)/.cache)/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)/standards.mk

include $(STANDARDS_MK)

$(STANDARDS_MK):
	@mkdir -p "$(@D)"
	curl -fsSL --max-time 30 "$(R_BIOC_STANDARDS_BASE)/shared/standards.mk" -o "$@.tmp"
	@mv -f "$@.tmp" "$@"

# Extra targets for this package go below, each listed in AGENTS.md. To
# replace a standard target, define it here; make warns that it overrides
# the shared recipe.

.PHONY: build style site site-deploy app test-app clean

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
