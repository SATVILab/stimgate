# Manage axis limits

Manage axis limits. Fix axis limits to be equal between x- and y-axes,
and/or expand axis coordinates. The primary use of `axisLimits` is
forcing the x- and y-axes to have the same limits "automatically" (i.e.
by inspecting the `ggplot` object, thus not requiring the user to
manually calculate limits to pass to
[`ggplot2::expand_limits`](https://ggplot2.tidyverse.org/reference/expand_limits.html)).

## Usage

``` r
axisLimits(p, limitsExpand = NULL, limitsEqual = FALSE)

axis_limits(p, limits_expand = NULL, limits_equal = FALSE)
```

## Arguments

- p:

  object of class 'ggplot' Limits are adjusted for this plot.

- limitsExpand:

  list or NULL If not NULL, then it is (effectively) passed onto
  ggplot2::expand_limits to ensure that certain values are included in
  the plot (such as, for example, 0 if that is the minimum value
  possible but it may not be plotted). If not named, then must consist
  of one numeric vector that will then force all values in the numeric
  value to be included in the plot. If named, then must have names x
  and/or y, with the elements again being numeric vectors that must be
  included in plot. Default is NULL.

- limitsEqual:

  logical If TRUE, then the ranges on the x- and y-axes must be equal.
  Effectively applied after expand_grid is applied. Default is FALSE.

- limits_expand:

  list or NULL Alias for `limitsExpand`. Default is NULL.

- limits_equal:

  logical Alias for `limitsEqual`. Default is FALSE.

## Value

A ggplot object with adjusted axis limits.

## Examples

``` r
data("cars", package = "datasets")
library(ggplot2)
p <- ggplot(cars, aes(speed, dist)) +
  geom_point()

axisLimits(
  p,
  limitsEqual = TRUE
)


# both axes
axisLimits(
  p,
  limitsExpand = list(200)
)

# x only
axisLimits(
  p,
  limitsExpand = list(x = 75)
)

# y only
axisLimits(
  p,
  limitsExpand = list(y = 200)
)

# lower and upper expansion
axisLimits(
  p,
  limitsExpand = list(
    y = c(-50, 200),
    x = c(-10, 75)
  )
)


# note that when fixing range and expanding, range is fixed
# after expansions are applied, so effectively the larger
# expansions apply to both.
# compare the following output to the previous output:
axisLimits(
  p,
  limitsExpand = list(
    y = c(-50, 200),
    x = c(-10, 75)
  ),
  limitsEqual = TRUE
)
```
