# Simulate cytometric data for a single stimulation condition

Simulate flow cytometric data for a single stimulation condition.

## Usage

``` r
simCytCondition(
  nMarker,
  nCell,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVec,
  probExact = FALSE,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
)
```

## Arguments

- nMarker:

  Integer. Number of markers/dimensions.

- nCell:

  Integer. Total number of cells to simulate.

- transformationFunc:

  Function. Transformation to apply to simulated data.

- mixtureType:

  Character. Type of mixture distribution.

- meanExprMat:

  Numeric matrix. Cluster mean vectors.

- probVec:

  Character vector. Cluster labels.

- probExact:

  Logical. If TRUE, use exact probabilities for cluster assignment; if
  FALSE, sample from a multinomial distribution. Default is `FALSE`.

- clusterPerturbationSd:

  Numeric. Cluster-level perturbation SD.

## Value

A list with `conditionMatrix` and `conditionLabels`.
