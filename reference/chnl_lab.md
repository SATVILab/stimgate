# Get markers and channels

From a cytometry object (e.g. flowFrame or flowSet), either get a
character vector of markers or channels (get_chnl and get_marker), or
get a named vector that converts between channel names and marker names
(e.g. chnl_to_marker).

## Usage

``` r
chnl_lab(data)
```

## Arguments

- data:

  object of class flowFrame, flowSet. Channel and corresponding marker
  names are drawn from here.

## Value

A named character vector.

## Details

Note that chnl_lab is equivalent to chnl_to_marker, and marker_lab is
equivalent to marker_to_chnl.

## Examples

``` r
# \donttest{
# Create example flowFrame-like data structure
data(GvHD, package = "flowCore") 
fs <- GvHD[1:2]

# Get channel to marker mapping
chnl_lab(fs[[1]])
#>               FSC-H               SSC-H               FL1-H               FL2-H 
#>        "FSC-Height"        "SSC-Height"         "CD15 FITC"           "CD45 PE" 
#>               FL3-H               FL2-A               FL4-H                Time 
#>        "CD14 PerCP"             "FL2-A"          "CD33 APC" "Time (51.20 sec.)" 
# }
```
