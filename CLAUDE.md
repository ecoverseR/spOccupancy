# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

spOccupancy is an R package that fits Bayesian occupancy models (Polya-Gamma augmentation, NNGP/GP spatial models). Each model is an R wrapper in `R/<Model>.R` plus a Gibbs sampler in `src/<Model>.cpp`, written against the plain R C API (no Rcpp).

## Commands

- Dev iteration: `devtools::load_all()` / `devtools::install()`
- Before a CRAN release: `R CMD build . && R CMD check --as-cran spOccupancy_<version>.tar.gz`
- Tests: run only the file for the model you touched, e.g. `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-PGOcc.R")'`. The full suite is slow (26 files of MCMC fits). Don't run it unless asked.
- `.Rbuildignore` excludes `tests/` and `vignettes/`, so `R CMD check` on a tarball runs no tests. Tests only run from source via devtools/testthat.

## Gotchas

- **No roxygen.** `NAMESPACE` and `man/*.Rd` are hand-written. Never run `devtools::document()`/`roxygenise()`; it will clobber NAMESPACE. Edit them by hand.
- `.Call` uses string names and positional args (often 35–65). `R_useDynamicSymbols` is FALSE, so every routine must be declared in `src/spOccupancy.h` and registered in `src/init.cpp` `CallEntries[]`. The arg count in init.cpp must equal both the C signature and the number of args passed from R.
- `ppcOcc()`, `waicOcc()`, `updateMCMC()` and `getSVCSamples()` are not S3 generics. They branch on hard-coded `class(object) %in% c(...)` lists, including a validity check near the top. A new model class must be added to those lists.
- A new exported function must also go in `_pkgdown.yml` reference sections.
- Before `.Call`, R code sets `storage.mode(x) <- "double"/"integer"` on every argument, and model dimensions are packed into an integer `consts` vector.
- `tests/`, `vignettes/` (pkgdown-only), `local_testing/` and `background_docs/` are not shipped in the package.

## Conventions

- R: dot.separated argument names (`occ.formula`, `n.omp.threads`, `n.burn`), and priors/inits named like `beta.normal`, `sigma.sq.psi.ig`. Use 2-space indent. Class checks use `class(x) %in% c(...)` to match existing code.
- C++: follow the `src/PGOcc.cpp` header boilerplate (`USE_FC_LEN_T`, `R_NO_REMAP` so use the `Rf_` prefix, `FCONE` on BLAS/LAPACK char args). SEXP params use the `_r` suffix and C locals are camelCase. Use `R_alloc` for scratch memory, and `PROTECT` with an `nProtect` counter followed by `UNPROTECT(nProtect)`. Use `GetRNGstate`/`PutRNGstate`, `R_CheckUserInterrupt()` in sampler loops and `Rprintf` for output. Guard OpenMP with `#ifdef _OPENMP`.
- `NEWS.md`: `# spOccupancy x.y.z` heading, then `+ ` bullets, linking GitHub issues/PRs.
- For multi-file changes (e.g. wiring up a new model), propose a plan before editing. Use the `/add-model` skill for the full checklist.
