# Power Law Fit Method  

This is the method used in `swipeMLEDiscretePowerLawBounds` to compute the power law fit. See Klaus et al., 2011, Statistical Analyses Support Power Law Distributions Found in Neuronal Avalanches.

## Maximum Likelihood Estimator  

Given the raw data, `xmin`, and `xmax`, the $\alpha$ parameter is computed by maximizing the likelihood function. In this case, it is much easier to find the point at which the derivative of the likelihood function is $0$ because this derivative is concave, meaning the dichotomy algorithm can be used. The dichotomy algorithm allows us to obtain exponentially more precise estimates of $\gamma$ with each iteration.  

The MLE is computed by the `DiscreteBoundedPowerLawMLE` function.  

## Grid Search for Range Parameters  

The only way to find the optimal `xmin` and `xmax` is, unfortunately, to perform a grid search over the parameter space. This parameter space is too large for an exhaustive grid search. To resolve this issue, we take advantage of the logarithmic nature of the power law and search through logarithmically spaced values of `xmin` and `xmax$, meaning we check:  

$
x_{\text{min}} = \alpha + b^i, \quad x_{\text{max}} = x_{\text{min}} + b^j
$

where $b$ is the `base` parameter, and $i$ and $j$ iterate over all possible integers such that $x_{\text{min}}$ and $x_{\text{max}}$ always remain within the range of the data. Additionally, we enforce the condition:  

$
\log(x_{\text{max}}/x_{\text{min}}) > d
$

where $d$ is the `min_decade` parameter. Here, $\alpha$ is the smallest value present in the raw data.  

### Finding the Maximum Value of $i$  

The function `upper_i` provides the maximum value $i_{\max}$ such that:  

$
\alpha + b^{i_{\max}} = \beta - 10^d
$

where $\beta$ is the maximum value in the raw data, and $d$ is the minimum number of decades required for the range.  

### Determining the Range of $j$  

The function `bound_j` determines the valid range of $j$ such that:  

$
x_{\text{min}} + 10^d \leq x_{\text{min}} + b^j \leq \beta
$

This ensures that there are at least $d$ decades between $x_{\text{min}}$ and $x_{\text{max}}$ while also preventing $x_{\text{max}}$ from exceeding the bounds of the data.  

## Goodness of Fit Test  

For each fit, the Kolmogorov-Smirnov goodness-of-fit test is performed. If the p-value of this test is above the specified threshold `significanceLevel`, then the corresponding range parameters become candidates for the optimal fit. The candidate with the largest range is selected. 