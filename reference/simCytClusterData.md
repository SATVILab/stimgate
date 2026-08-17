# Generate simulated data from mixture component

Sample from multivariate normal or t distribution.

## Usage

``` r
simCytClusterData(mixtureType, clusterNumber, nCell, muVec, sigmaMat)
```

## Arguments

- mixtureType:

  Character. "gaussianOnly", "tOnly", or "tPlusGauss".

- clusterNumber:

  Integer. Used for alternating distributions in "tPlusGauss".

- nCell:

  Integer. Number of samples.

- muVec:

  Numeric vector. Mean vector.

- sigmaMat:

  Numeric matrix. Covariance matrix.

## Value

Numeric matrix of sampled data (nCell x length(muVec)).
