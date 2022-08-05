# grLasso
Fast Lasso Regularization Paths for Logistic Regression Models with Grouped Covariates

## Overview
Provier profiling has been recognized as a useful tool in monitoring health care quality an improving medical cost-effectiveness. However, as the number of providers increase rapidly, existing packages are no longer able to handle such a large amount of data. For example, both glmnet and grpreg need an unpredictable amount time to convert provider information into thousands of dummy variables. And it can take hours or even days to solve such high-dimensional data by using decent algorithms.

Therefore, we combine the advantages of block ascent Newton method and coordinate descent algorithm and propose a new algorithm to deal with such computational challenges. Simulaiton studies show that our algorithm performs as well as grpreg package in low-dimension provider settings, and as the provider dimensions increase, our algorithm will be much faster than existing packages.


## Usage:

All the code here has been compressed into a R package. You can download the file named "TmpGrlasso.zip" to your local computer, unzip it and run the following code to install it:

```{r }
install.packages(".../TmpGrlasso", repos = NULL, type="source")
```

## Simulation:

### 1. Fixed Effects Model vs. Pooled Model:

<img src="https://drive.google.com/uc?export=view&id=1sIYLQXqYhRbNY493rsnNMF3_b33D0Cgp"  width=75% height=75%>


### 2. Regularization Coefficient Paths:

#### (1) Regularization path of model without unpenalized groups

Ten risk factors were randomly assigned to three groups. One group has true zero effect on the outcome.

<img src="https://drive.google.com/uc?export=view&id=1qL_72lY9ajn-qdN7y3iCI9DRKEdQczZW"  width=75% height=75%>

#### (2) Regularization path of model with unpenalized groups

Ten risk factors were randomly assigned to four groups. Among the four groups of variables, the true effect of one group on outcome was 0. Among the three groups with a non-zero real effect on outcome, we randomly chose one group to not be penalized.

<img src="https://drive.google.com/uc?export=view&id=1UlxAT7nXw5qB3O4d0DOYIpympKiGxD-H"  width=75% height=75%>


### 3.  Cross entropy Loss by cross validation:

Here, we use the "Birthwt" data provided with the grpreg package to show how cross-validation can be used to select the best regularization parameter. Currently, we use cross entropy loss as our selection criterion.


[You can find a detailed description of this data here](https://www.rdocumentation.org/packages/grpreg/versions/3.4.0/topics/Birthwt)

<img src="https://drive.google.com/uc?export=view&id=1XgdVbr3fjXp97nqeQA4EkKqWLAGrZhhK"  width=75% height=75%>


## Comparison with grpreg:


### With only one provider (intercept)

In this section, we will compare the running speed and Estimation difference of our algorithm with grpreg.

The number of risk factors ranges from 4 to 50. For each number of risk factors, 10 different datasets were generated. Sample size of each dataset follows *Poisson*(4000). 

#### (1) Running Speed:

<img src="https://drive.google.com/uc?export=view&id=1MV9T0r8lOzQlz1GsWG59zSVTKyUIPRMX"  width=80% height=80%>

#### (2) Estimation difference:

<img src="https://drive.google.com/uc?export=view&id=13kKehyF41Xcw6ASH4ZWso2i2UHlB0zZe"  width=80% height=80%>

### With multiple providers

In this section, we will compare the running speed and model accuracy of our algorithm with grpreg.

The number of providers ranges from 50 to 1100. For each number of providers, we simulated 50 risk factors and randomly assigned them to 10 different groups (each group had 5 risk factors). For the sake of sparsity, only the first 2 groups had non-zero real effects, and the remaining 8 groups had zero real effects on the outcome.

Since the data is randomly simulated, the convergence speed of the algorithms may vary under different data, so the time of a single run may not increase strictly with the increase of provider size. Therefore, we randomly generated 10 datasets for each number of providers, which contained the same number of risk factors but with different effects. We will finally evaluate the performance of the model under the corresponding number of providers by mean and standard deviation.

#### (1) Running Speed:

As we mentioned above, the running time is reported by the mean and its standard deviation.

We only provide a naive example here, with provider counts from 50 to 1100. Even with relatively “small” number of providers, our model is already significantly faster than grpreg package. In real world data, the number of providers may reach thousands or tens of thousands, and then the running time of grpreg will be unacceptable.

<img src="https://drive.google.com/uc?export=view&id=1ugDFN2w3HPGB-1UztibVq0zkMAOszFoM"  width=80% height=80%>

<img src="https://drive.google.com/uc?export=view&id=1Gf_UVuOhchyexQCvXYJHuB3ufmgxbmXu"  width=80% height=80%>


*Note that the running time of the algorithm excludes the time required by grpreg to convert high-dimensional provider data into dummy variables. But we need to keep in mind that as the number of providers grows, data conversion can be incredibly time and memory consuming.*

#### (2) Model Accuracy:

The accuracy of model is evaluated by root mean square error (RMSE), root model error (RME), cross entropy loss and prediction error rate.

##### Root Mean Square Error

<img src="https://drive.google.com/uc?export=view&id=13n3E81ig_7iUSnB63tAtgjKq5ig_gVlu"  width=80% height=80%>

##### Root Model Error

<img src="https://drive.google.com/uc?export=view&id=110jNt6p9NklXbIaWO8AmzsJDRxH_BU5_"  width=80% height=80%>

##### Cross Entropy Loss

<img src="https://drive.google.com/uc?export=view&id=1YeoMEhWGbR6mpAPLHQyHLevfsU0RaCLi"  width=80% height=80%>

##### Prediction Error Rate

<img src="https://drive.google.com/uc?export=view&id=1t3yD_7XnZn5_8l-Ha_ARN3Xx6jogYLZV"  width=80% height=80%>


In terms of model accuracy, the figure above shows that our algorithm performs as well as grpreg for different numbers of providers.

