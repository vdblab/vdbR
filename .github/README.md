# GitHub Actions
The default `r-lib/actions/setup-r-dependencies@v2` uses `pqk`, which seems to be struggling to find the bioconductor dependencies despite specifying `BiocViews:` in the Description file.

The renv.lock file is used as an alternative.
