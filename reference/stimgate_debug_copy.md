# Copy the latest stimgate debug file to the working directory

Copies the most recent debug file created by stimgate to the current
working directory. The copied file uses the same filename as the source.

## Usage

``` r
stimgate_debug_copy(pathDir = getwd())
```

## Arguments

- pathDir:

  character. Directory to copy the debug file to. Default is
  [`getwd()`](https://rdrr.io/r/base/getwd.html) (i.e. the working
  directory).

## Value

character Path to the copied file (invisibly), or NULL if no debug file
exists.
