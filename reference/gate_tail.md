# Constructs a cutpoint for a flowFrame by using a derivative of the kernel density estimate

We determine a gating cutpoint using either the first or second
derivative of the kernel density estimate (KDE) of the `x`.

## Usage

``` r
.cytokine_cutpoint(
  x,
  num_peaks = 1,
  ref_peak = 1,
  method = c("first_deriv", "second_deriv"),
  tol = 0.01,
  adjust = 1,
  side = "right",
  strict = TRUE,
  plot = FALSE,
  auto_tol = FALSE,
  ...
)
```

## Arguments

- x:

  a `numeric` vector used as input data

- num_peaks:

  the number of peaks expected to see. This effectively removes any
  peaks that are artifacts of smoothing

- ref_peak:

  After `num_peaks` are found, this argument provides the index of the
  reference population from which a gate will be obtained. By default,
  the peak farthest to the left is used.

- method:

  the method used to select the cutpoint. See details.

- tol:

  the tolerance value

- adjust:

  the scaling adjustment applied to the bandwidth used in the first
  derivative of the kernel density estimate

- side:

  character specifying the side of the gate, either `'right'` (default)
  or `'left'`

- strict:

  `logical` when the actual number of peaks detected is less than
  `ref_peak`. an error is reported by default. But if `strict` is set to
  FALSE, then the reference peak will be reset to the peak of the far
  right.

- plot:

  logical specifying whether to plot the peaks found

- auto_tol:

  when TRUE, it tries to set the tolerance automatically.

- ...:

  additional arguments passed to `.deriv_density`

## Value

the cutpoint along the x-axis

## Details

By default, we compute the first derivative of the kernel density
estimate. Next, we determine the lowest valley from the derivative,
which corresponds to the density's mode for cytokines. We then contruct
a gating cutpoint as the value less than the tolerance value `tol` in
magnitude and is also greater than the lowest valley.

Alternatively, if the `method` is selected as `second_deriv`, we select
a cutpoint from the second derivative of the KDE. Specifically, we
choose the cutpoint as the largest peak of the second derivative of the
KDE density which is greater than the reference peak.
