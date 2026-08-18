# Get markers and channels

From a cytometry object (e.g. flowFrame or flowSet), either get a
character vector of markers or channels (getChnl and getMarker), or get
a named vector that converts between channel names and marker names
(e.g. chnlToMarker).

## Usage

``` r
chnlLab(data)
```

## Arguments

- data:

  object of class flowFrame, flowSet. Channel and corresponding marker
  names are drawn from here.

## Value

A named character vector.

## Details

Note that chnlLab is equivalent to chnlToMarker, and markerLab is
equivalent to markerToChnl.

## Examples

``` r
# \donttest{
if (requireNamespace("flowCore", quietly = TRUE)) {
  exprs <- matrix(
    seq_len(8),
    ncol = 2,
    dimnames = list(NULL, c("FSC-A", "FL1-H"))
  )
  ff <- flowCore::flowFrame(exprs)

  # Get channel to marker mapping
  chnlLab(ff)
}
#>   FSC-A   FL1-H 
#> "FSC-A" "FL1-H" 
# }
```
