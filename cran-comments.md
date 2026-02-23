## Resubmission

This is a resubmission of zoomgrid (previously archived).

The package has been revised to address the issues raised in the
previous CRAN checks and to improve overall stability and compliance.

## Main changes

1. Vignette compliance

   External resources have been removed from the vignette.
   All figures and assets are now bundled within the package source.

2. Parallel backend refactoring

   The previous use of parallel::mcmapply() has been replaced with the
   future / future.apply framework.

   - Parallel execution now works cross-platform, including Windows.
   - The packages 'future' and 'future.apply' are listed in Suggests
     and are only required when parallel = TRUE.
   - Worker selection uses future::availableCores().

3. Bug fixes

   - Fixed a potential boundary indexing issue in build_subgrids().
   - Improved input validation for parallel and cores arguments.
   - Removed unintended modification of RNG state in
     grid_search_check().

4. Documentation updates

   - Updated package-level documentation to reflect the new
     parallel backend.
   - Corrected minor typographical issues.

## R CMD check results

0 errors | 0 warnings | 0 notes

The package was tested on:

- macOS (local)
- Windows (win-builder)
- R-devel and R-release

Thank you for reviewing the package.