# Print debug message conditionally

Writes debug output directly to pathProject/debug/debug.txt when
STIMGATE_DEBUG is enabled.

## Usage

``` r
.debug(msg, val = NULL)
```

## Arguments

- msg:

  character Message to print.

- val:

  object Optional value to append to message. Default: NULL.

## Value

logical invisibly TRUE if message was written, FALSE otherwise.
