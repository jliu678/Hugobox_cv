# Negative Binomial Deviance

Negative Binomial deviance is a measure of the goodness-of-fit for a Negative Binomial Generalized Linear Model (GLM). It quantifies how well a statistical model with a Negative Binomial error distribution fits the observed data, especially when dealing with overdispersed count data.

## What is Deviance in GLMs?

In the context of Generalized Linear Models (GLMs), the deviance serves as a generalization of the residual sum of squares from ordinary linear regression. It's a key statistic used to assess how well a model explains the variability in the response variable.

Formally, the deviance is defined as twice the difference between the log-likelihood of a saturated model and the log-likelihood of the current (fitted) model.

**Saturated Model**: A theoretical model that perfectly fits the observed data. It has as many parameters as data points, meaning each observation's predicted value is exactly equal to its observed value. This model represents the maximum possible log-likelihood.

**Current Model**: The statistical model you've built, with its estimated parameters.

The formula for deviance {{< math >}} $D$ {{< /math >}} is generally:

{{< math >}} 
$$D = 2[\log L(\text{saturated model}) - \log L(\text{current model})] \quad (1)$$
{{< /math >}}

A smaller deviance indicates a better fit of the current model to the data, as it means the log-likelihood of the current model is closer to the maximum possible log-likelihood of the saturated model.

## Why Negative Binomial Deviance?

The Negative Binomial distribution is often used for count data (non-negative integers like 0,1,2,…), especially when the data exhibits overdispersion.

**Overdispersion**: This occurs when the observed variance in count data is greater than its mean. The Poisson distribution assumes that its mean and variance are equal. If this assumption is violated (i.e., data is overdispersed), using a Poisson model can lead to underestimated standard errors and incorrect inferences. The Negative Binomial distribution is more flexible because it includes an additional dispersion parameter (often denoted as {{< math >}} $\theta$ {{< /math >}} or {{< math >}} $\alpha$ {{< /math >}}) that allows the variance to be greater than the mean.

Just as Poisson deviance is derived from the Poisson log-likelihood, Negative Binomial deviance is derived from the Negative Binomial log-likelihood. When fitting a Negative Binomial GLM, the model parameters (including the dispersion parameter) are estimated by maximizing the Negative Binomial log-likelihood, which is equivalent to minimizing the Negative Binomial deviance.

## Formula for Negative Binomial Deviance

The unit deviance for a single observation {{< math >}} $y$ {{< /math >}} with predicted mean {{< math >}} $\mu$ {{< /math >}} and dispersion parameter {{< math >}} $\theta$ {{< /math >}} (where variance is {{< math >}} $\mu + \mu^2/\theta$ {{< /math >}}) is:

{{< math >}} 
$$d(y,\mu;\theta) = 2\left[y\log\left(\frac{y}{\mu}\right) - (y+\theta^{-1})\log\left(\frac{1+\theta y}{1+\theta\mu}\right)\right] \quad (2)$$
{{< /math >}}

The total Negative Binomial deviance for a dataset is the sum of these unit deviances over all observations:

{{< math >}} 
$$D_{NB} = \sum_i d(y_i,\mu_i;\theta) = \sum_i 2\left[y_i\log\left(\frac{y_i}{\mu_i}\right) - (y_i+\theta^{-1})\log\left(\frac{1+\theta y_i}{1+\theta\mu_i}\right)\right] \quad (3)$$
{{< /math >}}

### Key Components:

- {{< math >}} $y_i$ {{< /math >}}: The observed count for observation i
- {{< math >}} $\mu_i$ {{< /math >}}: The mean count predicted by the model for observation i
- {{< math >}} $\theta$ {{< /math >}} (or {{< math >}} $\alpha$ {{< /math >}}): The dispersion parameter of the Negative Binomial distribution. This parameter accounts for the extra variance beyond what is explained by the mean

As {{< math >}} $\theta \to \infty$ {{< /math >}} (or {{< math >}} $\alpha \to 0$ {{< /math >}}), the Negative Binomial distribution approaches the Poisson distribution, and the Negative Binomial deviance approaches the Poisson deviance.

Conversely, smaller values of {{< math >}} $\theta$ {{< /math >}} (or larger {{< math >}} $\alpha$ {{< /math >}}) indicate greater overdispersion.

## Use in Practice

Negative Binomial deviance is crucial for:

**Model Fitting**: It's the objective function that optimization algorithms (e.g., iterative reweighted least squares) minimize to find the best-fitting parameters for a Negative Binomial GLM.

**Model Comparison**: Differences in deviance between nested Negative Binomial models can be used to perform likelihood ratio tests to assess the significance of adding or removing terms from a model.

**Goodness-of-Fit Assessment**: While not a formal test, the ratio of residual deviance to degrees of freedom (often called dispersion statistic) should ideally be close to 1 for a good fit. Values significantly greater than 1 suggest remaining overdispersion or model misspecification.

In fields like single-cell RNA sequencing, where UMI counts are often overdispersed, Negative Binomial models are frequently preferred over Poisson models, and minimizing Negative Binomial deviance is the standard approach for fitting such models.

## Derivation of Negative Binomial Deviance

To derive the Negative Binomial deviance formula, we follow the general definition of deviance in Generalized Linear Models (GLMs). The deviance is twice the difference between the log-likelihood of the saturated model (a model that perfectly fits the data) and the log-likelihood of the current (fitted) model.

We'll use a common parameterization of the Negative Binomial distribution where its probability mass function (PMF) is given by:

{{< math >}} 
$$P(Y=y|\mu,\theta) = \frac{\Gamma(y+\theta^{-1})}{y!\Gamma(\theta^{-1})} \left(\frac{\theta^{-1}}{\theta^{-1}+\mu}\right)^{\theta^{-1}} \left(\frac{\mu}{\theta^{-1}+\mu}\right)^y \quad (4)$$
{{< /math >}}

where:
- {{< math >}} $y$ {{< /math >}} is the observed count ({{< math >}} $y=0,1,2,\ldots$ {{< /math >}})
- {{< math >}} $\mu$ {{< /math >}} is the mean of the distribution, {{< math >}} $E[Y]=\mu$ {{< /math >}}
- {{< math >}} $\theta$ {{< /math >}} is the dispersion parameter ({{< math >}} $\theta > 0$ {{< /math >}}). The variance is {{< math >}} $\text{Var}[Y]=\mu+\theta\mu^2$ {{< /math >}}. (Note: some parameterizations use {{< math >}} $\alpha=1/\theta$ {{< /math >}}, so variance is {{< math >}} $\mu+\mu^2/\alpha$ {{< /math >}}). For this derivation, we use {{< math >}} $\theta$ {{< /math >}}.

### 1. Log-Likelihood for a Single Observation

First, let's write the log-likelihood for a single observation {{< math >}} $y_i$ {{< /math >}} with its corresponding mean {{< math >}} $\mu_i$ {{< /math >}}:

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = \log\left(\frac{\Gamma(y_i+\theta^{-1})}{y_i!\Gamma(\theta^{-1})} \left(\frac{\theta^{-1}}{\theta^{-1}+\mu_i}\right)^{\theta^{-1}} \left(\frac{\mu_i}{\theta^{-1}+\mu_i}\right)^{y_i}\right) \quad (5)$$
{{< /math >}}

Using logarithm properties ({{< math >}} $\log(ab)=\log a+\log b$ {{< /math >}}), we get:

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = \log(\Gamma(y_i+\theta^{-1})) - \log(y_i!) - \log(\Gamma(\theta^{-1})) + \theta^{-1}\log(\theta^{-1}) - \theta^{-1}\log(\theta^{-1}+\mu_i) + y_i\log(\mu_i) - y_i\log(\theta^{-1}+\mu_i) \quad (6)$$
{{< /math >}}

Rearranging terms to group by {{< math >}} $\mu_i$ {{< /math >}} and {{< math >}} $y_i$ {{< /math >}}:

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = y_i\log(\mu_i) - y_i\log(\theta^{-1}+\mu_i) - \theta^{-1}\log(\theta^{-1}+\mu_i) + \log(\Gamma(y_i+\theta^{-1})) - \log(y_i!) - \log(\Gamma(\theta^{-1})) + \theta^{-1}\log(\theta^{-1}) \quad (7)$$
{{< /math >}}

This can be rewritten as:

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = y_i\log\left(\frac{\mu_i}{\theta^{-1}+\mu_i}\right) - \theta^{-1}\log(\theta^{-1}+\mu_i) + C_1(y_i,\theta) \quad (8)$$
{{< /math >}}

where {{< math >}} $C_1(y_i,\theta) = \log(\Gamma(y_i+\theta^{-1})) - \log(y_i!) - \log(\Gamma(\theta^{-1})) + \theta^{-1}\log(\theta^{-1})$ {{< /math >}} are terms that depend on {{< math >}} $y_i$ {{< /math >}} and {{< math >}} $\theta$ {{< /math >}} but not {{< math >}} $\mu_i$ {{< /math >}}.

We can also write {{< math >}} $\log(\theta^{-1}+\mu_i) = \log(1/\theta+\mu_i) = \log\left(\frac{1+\theta\mu_i}{\theta}\right) = \log(1+\theta\mu_i) - \log(\theta)$ {{< /math >}}.

Substituting this into the log-likelihood:

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = y_i\log\left(\frac{\mu_i\theta}{1+\theta\mu_i}\right) - \theta^{-1}(\log(1+\theta\mu_i) - \log(\theta)) + C_1(y_i,\theta) \quad (9)$$
{{< /math >}}

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = y_i\log(\mu_i) + y_i\log(\theta) - y_i\log(1+\theta\mu_i) - \theta^{-1}\log(1+\theta\mu_i) + \theta^{-1}\log(\theta) + C_1(y_i,\theta) \quad (10)$$
{{< /math >}}

{{< math >}} 
$$\log L(y_i|\mu_i,\theta) = y_i\log(\mu_i) - (y_i+\theta^{-1})\log(1+\theta\mu_i) + C_2(y_i,\theta) \quad (11)$$
{{< /math >}}

where {{< math >}} $C_2(y_i,\theta) = C_1(y_i,\theta) + y_i\log(\theta) + \theta^{-1}\log(\theta)$ {{< /math >}} are terms that do not depend on {{< math >}} $\mu_i$ {{< /math >}}.

### 2. Log-Likelihood of the Saturated Model

For the saturated model, the predicted mean {{< math >}} $\mu_i^S$ {{< /math >}} is exactly equal to the observed value {{< math >}} $y_i$ {{< /math >}}.

So, substitute {{< math >}} $\mu_i = y_i$ {{< /math >}} into the log-likelihood function:

{{< math >}} 
$$\log L(y_i|y_i,\theta) = y_i\log(y_i) - (y_i+\theta^{-1})\log(1+\theta y_i) + C_2(y_i,\theta) \quad (12)$$
{{< /math >}}

### 3. Deriving the Deviance

The deviance for a single observation {{< math >}} $d(y_i,\mu_i;\theta)$ {{< /math >}} is defined as {{< math >}} $2 \times (\log L(y_i|y_i,\theta) - \log L(y_i|\mu_i,\theta))$ {{< /math >}}.

{{< math >}} 
$$d(y_i,\mu_i;\theta) = 2[(y_i\log(y_i) - (y_i+\theta^{-1})\log(1+\theta y_i) + C_2(y_i,\theta)) - (y_i\log(\mu_i) - (y_i+\theta^{-1})\log(1+\theta\mu_i) + C_2(y_i,\theta))] \quad (13)$$
{{< /math >}}

The {{< math >}} $C_2(y_i,\theta)$ {{< /math >}} terms cancel out:

{{< math >}} 
$$d(y_i,\mu_i;\theta) = 2[y_i\log(y_i) - y_i\log(\mu_i) - (y_i+\theta^{-1})\log(1+\theta y_i) + (y_i+\theta^{-1})\log(1+\theta\mu_i)] \quad (14)$$
{{< /math >}}

Rearranging the terms to match the target formula:

{{< math >}} 
$$d(y_i,\mu_i;\theta) = 2[y_i(\log(y_i) - \log(\mu_i)) - (y_i+\theta^{-1})(\log(1+\theta y_i) - \log(1+\theta\mu_i))] \quad (15)$$
{{< /math >}}

{{< math >}} 
$$d(y_i,\mu_i;\theta) = 2\left[y_i\log\left(\frac{y_i}{\mu_i}\right) - (y_i+\theta^{-1})\log\left(\frac{1+\theta y_i}{1+\theta\mu_i}\right)\right] \quad (16)$$
{{< /math >}}

Finally, the total Negative Binomial deviance {{< math >}} $D_{NB}$ {{< /math >}} is the sum of these unit deviances over all observations {{< math >}} $i$ {{< /math >}}:

{{< math >}} 
$$D_{NB} = \sum_i d(y_i,\mu_i;\theta) = \sum_i 2\left[y_i\log\left(\frac{y_i}{\mu_i}\right) - (y_i+\theta^{-1})\log\left(\frac{1+\theta y_i}{1+\theta\mu_i}\right)\right] \quad (17)$$
{{< /math >}}

This derivation shows how the Negative Binomial deviance arises directly from the log-likelihoods of the saturated and fitted models under the Negative Binomial distribution. Minimizing this deviance is equivalent to maximizing the likelihood of the Negative Binomial GLM.