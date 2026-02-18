# Code Base Structure  

## Folder `regionsBridge`  

The functions in this folder take a `regions` instance as input and extract data from it. One deprecated feature is the separation of data with respect to the state (due to too many changes from the previous version). The code still works but will always include all available data regardless of the state. However, it can still separate regions. To fix this, the function in `regionsBridge/rawRegion/getter` must be changed.  

Note that the functions in this folder are only wrappers for functions in the `analysis` folder. These functions are primarily for convenience in extracting data from a `regions` instance.  

### Folder `rawRegion`  

- `rawRegionSwipePhases`: Computes the gamma estimations over all combinations of `phase`, `states`, and `region`. It is useful for identifying conditions under which criticality arises. (Mostly deprecated because states cannot be differentiated.)  
- `rawRegionAnalysis`: Computes all relevant indicators of criticality, as explained in the `README.md` file. Since the default ID value of `0` (which designated "all regions") no longer exists, you must provide a list of region IDs using the `regions` parameter. Additionally, since states cannot be differentiated, the `state` parameter is unnecessary—it will automatically include all available states in the `region` instance.  
- `rawRegionSwipe`: Computes gamma, branching ratio, and slow decay of autocorrelation estimations over all combinations of `states` and `region`. It is useful for identifying conditions under which criticality arises. (Mostly deprecated because states cannot be differentiated.)  
- `rawRegionBootstrapAnalysis`: Performs a bootstrap analysis on the three estimators of $\gamma$.  

### Folder `ICARegion` and `ICA`  

These functions were previously used for detecting avalanches within the ICA framework. However, due to a lack of results and interest, their development was discontinued. Do not use those functions.

## Folder `plot`  

This folder contains all functions used for generating plots.  

- `plotCriticalityScorePhases`: Plots the criticality score from `rawRegionSwipePhases`. It takes the results of `rawRegionSwipePhases` and visualizes them to identify conditions under which criticality arises. (Mostly generates random matrices, though.) This is an auxiliary function of `plotRawRegionSwipePhases`, which is simpler to use.  
- `plotCriticalityScore`: An auxiliary function of `plotCriticalityScorePhases`. `plotCriticalityScorePhases` plots the same graph for each `phase`.  

### Folder `ShapeCollapse`  

- `plotShapeCollapseTransformed`: Plots all avalanche shapes provided by `AvalAverageSizeTimeDependent`. These are interpolated avalanche shapes so that every avalanche has the same number of sample points regardless of its lifetime.  
- `plotScaledShapeCollapseTransform`: Similar to `plotShapeCollapseTransformed`, but scales the shapes according to the provided \(\gamma\) parameter. It also plots the resulting mean shape. If the \(\gamma\) value is correct, the resulting shape should be approximately quadratic, with all avalanches collapsing onto the same function.  

### Folder `Raw`  

Contains functions for plotting statistics computed on raw data (i.e., without ICA or PCA reduction) within a `Region` class instance. These functions are wrappers for those in `regionsBridge/rawRegion`, meaning they take the results from those functions and visualize them. They require a `regions` instance as input, which is convenient for data handling but has the same deprecation issue as `regionsBridge/rawRegion`.  

- `plotRawRegionSwipePhases`: Plots the criticality score, which represents the precision of the three \(\gamma\) estimates. The precision is defined as:  

  $
  p = \frac{1}{\operatorname{var}(\gamma_1, \gamma_2, \gamma_3)}
  $

  where each \(\gamma_i\) is an estimation of the criticality parameter \(\gamma\) using a different method:  
  - One from the scaling relation,  
  - One from the Area vs. Lifetime relation,  
  - One from the shape collapse.  

  This function uses the results of `rawRegionSwipePhases`.  

- `plotRawRegionAnalysis`: Plots the results of `rawRegionAnalysis`. Like `rawRegionAnalysis`, you must provide region numbers, and the `state` parameter is deprecated. This function plots power-law fits, the Area vs. Lifetime relation, shape collapse, branching ratio, and slow decay of autocorrelation.  
- `plotRawRegionSwipe`: Plots power-law fits, the Area vs. Lifetime relation, shape collapse, branching ratio, and slow decay of autocorrelation parameters for each combination of state and region. (Remember, state acquisition is deprecated, meaning it doesn't work—fix the getter function in `regionsBridge/rawRegion`.)  
- `plotRawRegionBranchingRatio`: Plots the linear regression of the delayed branching ratio computed by `branchingRatio`.  
- `plotRawRegionBootstrapAnalysis`: Plots the PDF of the three estimators of \(\gamma\) and their p-values against a normal distribution.  

### Folder `powerLaw`  

These functions generate power-law fit plots from raw data.  

- `plotDiscretePowerLawDensityFit`: Plots the empirical PDF and the fitted PDF, given the raw data and the parameters of the fitted PDF.  
- `plotPowerLawFunction`: Plots the raw data and the fitted power-law function, given the `fitlm` MATLAB object. This function is used to visualize the Area vs. Lifetime relationship.  

## Folder `compute`  

Contains processing functions. This folder is mostly independent of the `region` class and operates directly on raw data.  

- `swipeMLEDiscretePowerLawBounds`: Fits the power-law distribution to the data.  
- `preprocess`: Removes outliers (i.e., occurrences that appear only once) from raw data.  
- `MLEDiscretePowerLawDensityAprox`: Uses an analytical approximation to estimate the \(\alpha\) parameter (does not take `xmax` into account).  

### Folder `timeDependentSize`  

Contains functions for analyzing avalanche size over time.  

- `separateAvalSizeTimeDependent`: Extracts avalanches from a profile list where each avalanche is separated by a `0`.  
- `AvalAverageSizeTimeDependent`: Computes the mean time-dependent avalanche size for each avalanche lifetime based on the results of `separateAvalSizeTimeDependent`.  

### Folder `shapeCollapse`  

- `transformCollapseShape`: Takes the result of `AvalAverageSizeTimeDependent` and interpolates the avalanche shapes to a uniform length.  
- `scaleCollapseShape`: Scales the shapes by $T^{1-\gamma}$ to collapse all avalanches onto the same shape (if the \(\gamma\) parameter is correct).  
- `lifetimeThresholdCollapseShapeTransformed`: Removes all shapes with a lifetime below the given threshold.  
- `fitCollapseShape`: Finds the optimal \(\gamma\) parameter by minimizing the variability of the resulting shape.  

### Folder `goodnessFit`  

Contains functions for statistical goodness-of-fit tests.  

- `KolmogorovTest`: Computes the Kolmogorov goodness-of-fit test.  
- `ChiSquaredTest`: Computes the chi-squared goodness-of-fit test (not used).  

### Folder `DiscreteBoundedPowerLawDensityMLE`  

- `DiscreteBoundedPowerLawMLE`: Estimates the optimal \(\alpha\) parameter, given the raw data and `xmin` and `xmax`.  
- `DiscreteBoundedPowerLawLikelihoodDerivative`: Computes the derivative of the log-likelihood. The algorithm attempts to find the zero-crossing of this function, assuming the likelihood is concave.  

### Folder `analysis`  

- `regionSwipeAnalysis`: Performs a sweep over all possible state and region combinations to estimate the three \(\gamma\) values. This assigns a criticality score to each state-region combination.  
- `analysis`: Performs the overall criticality analysis as described in the `README.md` file.  