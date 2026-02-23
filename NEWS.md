<!-- README.md is generated from README.Rmd. Please edit that file -->

# New Features in zoomgrid 1.1.0

## CRAN compliance improvements

- Removed external resources in vignettes to comply with CRAN policies.
  All vignette figures and assets are now bundled within the package
  source.

## Parallel backend refactoring

- Replaced `parallel::mcmapply()` with the `future` and `future.apply`
  framework.
- Parallel execution is now cross-platform and works correctly on
  Windows.
- The `future` ecosystem is listed under Suggests, and parallel
  execution is enabled only when the required packages are available.
- Improved core detection and worker selection logic for safer
  execution.

## Bug fixes

- Fixed a potential out-of-bounds indexing issue in `build_subgrids()`
  related to endpoint index handling.
- Improved internal index validation to prevent runtime boundary errors.

## User-facing output

- Improved console output in `grid_search()` and `grid_search_check()`
  for clearer progress reporting.
- Added compact CLI-style headers and runtime summaries when
  `silent = FALSE`.

# New Features in zoomgrid 1.0.0

- All the functions
- Documentation including README, LICENSE and etc.
