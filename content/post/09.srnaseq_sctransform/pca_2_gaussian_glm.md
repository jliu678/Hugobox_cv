# PCA and Maximum Likelihood Equivalence

When PCA minimizes the squared reconstruction error, it's simultaneously finding the parameters ({{< math >}} $u_i'$ {{< /math >}} and {{< math >}} $v_j$ {{< /math >}}) that are most likely to have generated the observed data {{< math >}} $z_{ij}$ {{< /math >}}, assuming the data follows a Gaussian distribution with a constant variance.

Let's break down why this equivalence holds:

## Gaussian Model Likelihood

The model states that each data point {{< math >}} $z_{ij}$ {{< /math >}} is drawn from a Normal distribution:

{{< math >}} 
$$z_{ij} \sim N(\mu_{ij}, \sigma^2) \quad (1)$$
{{< /math >}}

where {{< math >}} $\mu_{ij} = u_i' v_j$ {{< /math >}} and {{< math >}} $\sigma^2$ {{< /math >}} is a constant variance across all data points.

The **likelihood function** for a single data point {{< math >}} $z_{ij}$ {{< /math >}} given this model is:

{{< math >}} 
$$L(z_{ij} | u_i', v_j, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(z_{ij} - u_i' v_j)^2}{2\sigma^2}} \quad (2)$$
{{< /math >}}

Assuming all data points are independent, the total likelihood for the entire dataset (matrix Z) is the product of individual likelihoods:

{{< math >}} 
$$L(Z | U, V, \sigma^2) = \prod_{i,j} \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(z_{ij} - u_i' v_j)^2}{2\sigma^2}} \quad (3)$$
{{< /math >}}

## Maximizing the Likelihood

To maximize this likelihood function, it's often easier to maximize its logarithm (the **log-likelihood**), which doesn't change the location of the maximum:

{{< math >}} 
$$\log L = \sum_{i,j}\left[-\frac{1}{2}\log(2\pi\sigma^2) - \frac{(z_{ij} - u_i' v_j)^2}{2\sigma^2}\right] \quad (4)$$
{{< /math >}}

{{< math >}} 
$$\log L = C - \frac{1}{2\sigma^2}\sum_{i,j}(z_{ij} - u_i' v_j)^2 \quad (5)$$
{{< /math >}}

where {{< math >}} $C = -\frac{IJ}{2}\log(2\pi\sigma^2)$ {{< /math >}} is a constant with respect to {{< math >}} $u_i'$ {{< /math >}} and {{< math >}} $v_j$ {{< /math >}}.

To **maximize** this log-likelihood, we need to **minimize** the term:

{{< math >}} 
$$\sum_{i,j}(z_{ij} - u_i' v_j)^2 \quad (6)$$
{{< /math >}}

## The Equivalence

This is precisely the **PCA objective function** we discussed earlier:

{{< math >}} 
$$\min_{u,v} \sum_{i,j}(z_{ij} - u_i' v_j)^2 \quad (7)$$
{{< /math >}}

Therefore, minimizing the sum of squared reconstruction errors (the PCA objective) is mathematically equivalent to maximizing the likelihood of observing the data Z under the assumption that each data point {{< math >}} $z_{ij}$ {{< /math >}} is drawn from a Gaussian distribution with a mean equal to its low-rank approximation ({{< math >}} $u_i' v_j$ {{< /math >}}) and a constant variance ({{< math >}} $\sigma^2$ {{< /math >}}).

This connection highlights PCA's probabilistic interpretation, showing that it's not just a geometric method but also a form of **statistical modeling**. It implicitly assumes that the noise or error in the data is Gaussian and independent for each element. This provides a stronger theoretical foundation for using PCA when data are believed to have such underlying characteristics.