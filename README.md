# grLasso
Fast Lasso Regularization Paths for Data with a Large number of Health Care Providers}

## Overview
The increasing availability of large-scale datasets collected by national registries has provided extraordinary opportunities for biomedical research. However, as the number of health care providers included escalates along with the sample size, many commonly used statistical tools for high-dimensional variable selection that are designed for general purposes are no longer computationally efficient, especially when the magnitude of the number of providers is comparable with or far exceeds the dimension of risk factors.

Here, we propose an efficient statistical tool for conducting variable selection on data containing a large number of center effects. Our proposed solution is motivated by Block Ascent Newton introduced by Dr.He(2013) and can be seen as an extension of the (block) coordinate descent algorithm introduced by Dr.Breheny(2015) and Dr.Friedman(2010). For GLM problems, center effects and penalized parameters of risk factors are updated by an iterative procedure. In the outer layer of each iteration, the unpenalized center effects are updated by Newton's method, then the coefficients of risk factors are updated using the subdifferential-based (block) coordinate descent method based on the updated center effect in the inner layer. For discrete survival models, an additional middle layer was used to update the required estimate of baseline hazard over different time points without data expansion. Because of the special structure of data, the efficiency of using the Newton method to update the center effect and baseline hazard is much higher than that of the coordinate descent method, consequently, our proposed algorithm can markedly reduce the processing time and memory required. 


## Usage:

The code has been compressed into a temporary R package. You can download the file named "ppSurv.zip" (for discrete survival models) and "ppLasso.zip" (for GLMs) to your local computer, unzip it and run the following code to install it:

```{r }
install.packages(".../ppLasso", repos = NULL, type = "source")
install.packages(".../ppSurv", repos = NULL, type = "source")
```

## Simulation:

### 1. Running time:

Below we compare the running speed of our statistical tool to existing efficient statistical tools. Under each simulation setting, we simulated 20 data and show the mean running time and standard deviation. The running time is measured based on the 10-fold cross-validation across the whole regularization path. Experiments were conducted on Intel Xeon Gold 6230 processor with base frequency 2.10 GHz and Memory 192GB. Both ppSurv and pplasso were implemented using Rcpp and RcppArmadillo

(1) Time to convergence under the different number of centers for logistic regression model without (left) or with (right) grouped variable. We simulate the number of transplant centers ranging from 50 to 400, and include 50 risk factors with 10 of them having non-zero effects.

![Runtime_ppLasso](https://user-images.githubusercontent.com/93620754/207072629-38c87d5a-963f-488f-89c4-5f45dddc0e53.png)


(2) Time to convergence under the different number of centers for discrete survival logistic models with 30 (left) and 50 (right) discrete time points. The baseline hazard follows the Weibull distribution. Due to the limitation of computing power, we only show the comparison between algorithms under a moderate size of data (100 centers) as an illustrating example.

![Runtime_ppSurv](https://user-images.githubusercontent.com/93620754/207072650-54c2c38a-9ff3-4f07-b17d-8b81e45bd7fc.png)


### 2. Regularization Coefficient Paths:

#### (1) Regularization path of model without unpenalized groups

![Regularization_Path_unpenalize](https://user-images.githubusercontent.com/93620754/207072941-47c1f5ab-c7d7-4a28-9df8-9af9ec0b05aa.png)


#### (2) Regularization path of model with unpenalized groups

![Regularization_Path_penalize](https://user-images.githubusercontent.com/93620754/207072975-d2106f33-f030-4bb7-ac1a-88470ada2787.png)


### 3.  Cross entropy Loss by cross validation:

Here, we use the "Birthwt" data provided with the grpreg package to show how cross-validation can be used to select the best regularization parameter. Currently, we use cross entropy loss as our selection criterion.

[You can find a detailed description of this data here](https://www.rdocumentation.org/packages/grpreg/versions/3.4.0/topics/Birthwt)

![Cross_Entropy_Loss](https://user-images.githubusercontent.com/93620754/207073019-37e6da80-f07b-4f4a-a766-442e3ed41dc9.png)


