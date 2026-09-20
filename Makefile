R_PKG := 02_code/r

.PHONY: help install test check cli gui deps clean-builds

help:
	@echo "nanoamp project targets:"
	@echo "  make install      Install the R package (R CMD INSTALL $(R_PKG))"
	@echo "  make test         Run testthat tests"
	@echo "  make check        Build and R CMD check into 05_builds/r"
	@echo "  make cli          Run 'nanoamp doctor' from the repository CLI"
	@echo "  make gui          Launch the Shiny GUI"
	@echo "  make deps         Fetch bundled external tools where possible"
	@echo "  make clean-builds Remove 05_builds/r contents"

install:
	R CMD INSTALL $(R_PKG)

test:
	Rscript -e 'devtools::test("$(R_PKG)", reporter = "summary", stop_on_failure = TRUE)'

check:
	mkdir -p 05_builds/r
	R CMD build $(R_PKG) --no-build-vignettes
	mv nanoamp_*.tar.gz 05_builds/r/
	cd 05_builds/r && R CMD check --no-manual --no-build-vignettes nanoamp_*.tar.gz

cli:
	sh 02_code/cli/nanoamp doctor

gui:
	Rscript 02_code/gui/run_gui.R

deps:
	bash 03_dependence/fetch_dependencies.sh

clean-builds:
	rm -rf 05_builds/r/*
