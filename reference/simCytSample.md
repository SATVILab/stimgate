# Simulate all cytokine combinations for a set of stimulation conditions for a single biological sample

Simulate a set of stimulation conditions (e.g., stimulated and
unstimulated) for a single biological sample, where the first condition
is always the unstimulated condition.

## Usage

``` r
simCytSample(
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns,
  probExact,
  probResponseVecByStimCondition = NULL,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
)
```

## Arguments

- nMarker:

  Integer. Number of markers/dimensions.

- nCondition:

  Integer. Number of stimulation conditions (must be \>= 2).

- nCluster:

  Integer. Number of clusters (must be a power of 2 between 2 and 1024).

- nCellByCondition:

  Integer or numeric vector. Number of cells per condition. If a single
  value, it is recycled for all conditions.

- transformationFunc:

  Function. Transformation to apply to simulated data (e.g., identity or
  gamma transformation).

- mixtureType:

  Character. Type of mixture distribution: "gaussianOnly", "tOnly", or
  "tPlusGauss".

- meanExprMat:

  Numeric matrix. Cluster mean vectors (nCluster x nMarker).

- clusterLabelVec:

  Character vector. Labels for each cluster (length nCluster).

- probVecUns:

  Numeric vector. Probability distribution for unstimulated condition
  (length nCluster, sums to 1).

- probExact:

  Logical. If TRUE, use exact probabilities for cluster assignment; if
  FALSE, sample from a multinomial distribution. Default is `FALSE`.

- probResponseVecByStimCondition:

  List. Probability response vectors for each stimulated condition
  (length nCondition - 1, each of length nCluster).

- conditionPerturbationSd:

  Numeric. Standard deviation of condition-level perturbations to
  cluster means.

- clusterPerturbationSd:

  Numeric. Standard deviation of cluster-level perturbations within each
  condition.

## Value

A list of length nCondition with named elements "unstim", "stim1", etc.
Each element contains:

- `conditionMatrix`: Numeric matrix of simulated data (nCell x nMarker).

- `conditionLabels`: Character vector of cluster labels for each cell.
