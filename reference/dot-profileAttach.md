# Attach profiling state without deleting existing raw records

This is mainly a safeguard for worker processes that enter a profiled
helper after the main process has already created the profile directory.

## Usage

``` r
.profileAttach(pathProject)
```
