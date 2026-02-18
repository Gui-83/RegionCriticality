# Brain Criticality Documentation

The folder `brain-criticality` contains functions to analyze criticality in avalanches, using the lifetime (`T`) and size (`S`) of a set of avalanches.

--- 

## Get Data from `regions` Class  

To compute the avalanche data from a `regions` class:  

```matlab
session = "/mnt/cortex-data-311/Rat386-20180921/Rat386-20180921.xml";
R = regions(session);
R = R.loadSpikes();
R = R.computeAvalanches();
```

You can also compute avalanches with pca reduction step :

```matlab
session = "/mnt/cortex-data-311/Rat386-20180921/Rat386-20180921.xml";
R = regions(session);
R = R.loadSpikes();
R = R.computeAvalanches(dopc=true, var=0.9);
```

Here the `var` (percent) argument specified how much of the variance you want to keep. If you want to compare it to the remaining pc, meaning keep only low variance pc that explains `1-var` percent of the variance : 

```matlab
R = R.computeAvalanches(dopc=true, var=0.9, first=false);
```

The relevant variables are stored in the `regions_array` member of the `R` object.

## Compute Power-Law Fit

Use the `swipeMLEDiscretePowerLawBounds` function to compute the power-law fit. This function takes six arguments:

- **`x`**: The data to fit the power law (e.g., `T` or `S`).
- **`significanceLevel`**: The minimum `p`-value for the goodness-of-fit test. Data is not considered power-law distributed if the `p`-value is below this threshold.
- **`dicoStep`**: The number of dichotomy steps during optimization of the power-law parameter $α$. The parameter $α$ is accurate to $C × 2^{-dicoStep}$, where $C` defaults to `5` (assuming $α$ lies between `1` and `6`).
- **`min_decade`**: The minimum number of decades the fit range must span.
- **`base`**: A parameter for grid search to optimize the fit range. The search tests `xmin` and `xmax` as $x_{\text{min}} = \text{base}^i$ and $x_{\text{max}} = x_{\text{min}} + \text{base}^j$ for integers $i` and `j`.
- **`verbose`**: A boolean for verbose output.

**Outputs:**

- **`alpha`**: The power-law parameter.
- **`xmin`**: The minimum value of the fitted range.
- **`xmax`**: The maximum value of the fitted range.
- **`p`**: The goodness-of-fit `p`-value.

**Example:**

```matlab
[alpha, xmin, xmax, p] = swipeMLEDiscretePowerLawBounds(T, 0.05, 10, 1, 3, true);
```

---

## Plotting Power-Law Fit

To plot the result of `swipeMLEDiscretePowerLawBounds`, use the `plotDiscretePowerLawDensityFit` function to overlay the empirical PDF and fitted PDF on a log-log scale:

```matlab
plotDiscretePowerLawDensityFit(T, xmin, xmax, alpha);
```

---

## Compute Area vs. Lifetime Relationship

To compute the `area`:

```matlab
A = getArea(T, S);
```

This calculates the `area` for each lifetime value between $1$ and $\max(T)$. If no avalanche with lifetime `t` exists, then `A(t) = 0`.

Fit a linear model:

```matlab
lm = fitPowerFunction(A);
```

Here, `lm` is the output of MATLAB’s built-in `fitlm` function.

---

## Plot Area vs. Lifetime Relationship

To plot the linear model:

```matlab
plotPowerLawFunction(lm);
```

---

## Compute Shape Collapse

To compute shape collapse, use `aval_timeDependentSize`, a one-dimensional list where avalanches are described by their time-dependent size, separated by `0`. 

**Example:** `[0 1 1 5 0 1 0 3 1 3 0]` represents two avalanches.

1. Separate avalanches and compute their lengths:

    ```matlab
    [ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);
    ```

2. Compute the average time-dependent size for each lifetime:

    ```matlab
    ST_AVG = AvalAverageSizeTimeDependent(ST, lengthST);
    ```

3. Interpolate to homogenize shapes:

    ```matlab
    [x, shape, T] = transformCollapseShape(ST_AVG);
    ```

    - **`x`**: A linearly spaced list (size `1×100`) from `0` (start) to `1` (end).
    - **`shape`**: Interpolated average avalanche shapes (size `100×number_of_unique_lifetimes`).
    - **`T`**: List of lifetimes for avalanches in `shape`.

4. Filter by lifetime threshold (e.g., `4`):

    ```matlab
    [T, shape] = lifetimeThresholdCollapseShapeTransformed(T, shape, 4);
    ```

5. Fit $\gamma$ value:

    ```matlab
    gam = fitCollapseShape(T, shape, 10);
    ```

---

## Plot Shape Collapse

To visualize the shape collapse:

```matlab
plotScalledShapeCollapseTransform(x, T, shape, gam);
```

---

## Compute Branching Ratio

To compute the branching ratio, use the `profile` list (available in the region class as an attribute):

```matlab
[m, r, lm] = branchingRatio(profile, 10);
```

- **`m`**: The branching ratio.
- **`r`**: The lagged branching ratio list.
- **`lm`**: The log-linear model fitted to `r`.

---

## Plot Branching Ratio

To plot the result:

```matlab
plotBranchingRatioFit(m, r, lm);
```

---

## Compute Slow Decay of Autocorrelation

1. Separate avalanches and compute their lengths:

    ```matlab
    [ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);
    ```

2. Choose a reference time `T_ref` in $(0, 1)$ and compute:

    ```matlab
    [~, decay, mean_autocorr, lm] = autocorrelationDecay(ST, lengthST, 1/3);
    ```

---

## Plot Slow Decay of Autocorrelation

To plot the decay:

```matlab
plotAutocorrelationDecay(decay, mean_autocorr, lm, 1/3);
```

## Compute overall analysis

It is possible to do all the analysis at once provided `S`, `T`, `A`, `ST`, `lenghtST` and `profile` mentioned above :

```matlab
[alpha, Txmin, Txmax, Tp, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = analysis(S, T, A, ST, lengthST, profile, T_ref=1/3)
```

## Plot overall analysis

```matlab
plotAnalysis(T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm, 1/3)
```

## Compute overall analysis from an instance of regions

provided a instance of the region class you can perfom the analysis :

```matlab
session = "/mnt/cortex-data-311/Rat386-20180921/Rat386-20180921.xml";
R = regions(session);
R = R.loadSpikes();
R = R.computeAvalanches();
[T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = rawRegionAnalysis(R, regions=11);
```
It is very importent to provide a region number number (it can be a list of number) to `rawRegionAnalysis` because the default case raise an error. That is because all the regions data was stored in a region with the number `0` but that got deleted with the latest update of `regions`, see `regionsBridge/rawRegion/getter/getRawRegionProfile.m`.

## Plot overall analysis from an instance of regions

```matlab
session = "/mnt/cortex-data-311/Rat386-20180921/Rat386-20180921.xml";
R = regions(session);
R = R.loadSpikes();
R = R.computeAvalanches();
plotRawRegionAnalysis(R, regions=11);
```