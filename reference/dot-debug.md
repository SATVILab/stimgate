# Print debug message conditionally

Writes debug output to the most recent stimgate debug file when
STIMGATE_DEBUG is enabled.

## Usage

``` r
.debug(msg, val = NULL)
```

## Arguments

- msg:

  character Message to print

- val:

  object Optional value to append to message. Default is NULL.

## Value

logical invisibly TRUE if message was written, FALSE otherwise
