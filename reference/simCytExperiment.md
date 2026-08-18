# Simulate a set of stimulation conditions for multiple biological samples

Simulate stimulation conditions (e.g., stimulated and unstimulated) for
multiple biological samples, where each sample has one unstimulated
condition and one or more stimulated conditions. The outputs are
returned as lists of `flowFrame` objects representing each
sample-condition combination, and matching lists of cellular cluster
labels.

## Usage

``` r
simCytExperiment(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns,
  probExact = FALSE,
  probResponseVecByStimCondition = NULL,
  samplePerturbationSd = 0,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2
)
```

## Arguments

- nSample:

  Integer. Number of biological samples to simulate.

- nMarker:

  Integer. Number of markers/dimensions.

- nCondition:

  Integer. Number of conditions per sample. The first condition is
  unstimulated and the rest are stimulated.

- nCluster:

  Integer. Number of clusters. Must equal `2^nMarker`.

- nCellByCondition:

  Numeric or integer vector. Number of cells per condition. If length is
  1, the value is recycled across all conditions.

- transformationFunc:

  Function. Transformation applied marker-wise to simulated expression
  values.

- mixtureType:

  Character. Mixture distribution used for simulation. Supported values
  are "gaussianOnly", "tOnly", and "tPlusGauss".

- meanExprMat:

  Numeric matrix. Baseline cluster means with dimensions
  `nCluster x nMarker`.

- clusterLabelVec:

  Character vector. Cluster labels of length `nCluster`.

- probVecUns:

  Numeric vector. Baseline cluster probabilities for the unstimulated
  condition. Must have length `nCluster` and sum to 1.

- probExact:

  Logical. If TRUE, use exact probabilities for cluster assignment; if
  FALSE, sample from a multinomial distribution. Default is `FALSE`.

- probResponseVecByStimCondition:

  NULL or list. If provided, must be a list of length `nCondition - 1`,
  where each element is a numeric vector of length `nCluster`. Each
  vector is added to `probVecUns` to construct the stimulated-condition
  cluster probabilities.

- samplePerturbationSd:

  Numeric. Standard deviation of sample-level perturbations added to
  cluster means. Default is 0.

- conditionPerturbationSd:

  Numeric. Standard deviation of condition-level perturbations added to
  cluster means within each sample. Default is 0.

- clusterPerturbationSd:

  Numeric. Standard deviation of cluster-level perturbations applied
  during cell-level simulation. Default is 0.

- covEvMin:

  Numeric. Minimum eigenvalue for cluster covariance matrices. Default
  is 1.

- covEvMax:

  Numeric. Maximum eigenvalue for cluster covariance matrices. Default
  is 2.

## Value

A list with two elements:

- `flowFrameList`: A named list of
  [`flowCore::flowFrame`](https://rdrr.io/pkg/flowCore/man/flowFrame-class.html)
  objects.

- `labelsList`: A named list of character vectors of per-cell cluster
  labels. Names of list elements are formatted as `sampleXXxUnstim`,
  `sampleXXX_stim1`, etc.
