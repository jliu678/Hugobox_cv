# Moment Generating Function (MGF) of a Poisson Random Variable

## 🔷 Claim

If {{< math >}}$Y_j \sim \text{Poisson}(\lambda_j)${{< /math >}}, then the MGF of {{< math >}}$Y_j${{< /math >}}, evaluated at {{< math >}}$t_j${{< /math >}}, is:

{{< math >}}
$$ M_{Z_j}(t_j) = E[e^{t_j Y_j}] = \exp(\lambda_j(e^{t_j} - 1)) $$
{{< /math >}} {#eq:poisson-mgf-claim}

## 🔶 Step-by-Step Proof

### Step 1: Definition of MGF

The moment generating function (MGF) of a discrete random variable {{< math >}}$Y_j${{< /math >}} is:

{{< math >}}
$$ M_{Z_j}(t_j) = E[e^{t_j Y_j}] = \sum_{y=0}^{\infty} e^{t_j y} \cdot P(Y_j = y) $$
{{< /math >}} {#eq:mgf-definition}

For {{< math >}}$Y_j \sim \text{Poisson}(\lambda_j)${{< /math >}}, the PMF is:

{{< math >}}
$$ P(Y_j = y) = \frac{e^{-\lambda_j} \lambda_j^y}{y!} $$
{{< /math >}} {#eq:poisson-pmf}

Plug this into the sum:

{{< math >}}
$$ M_{Z_j}(t_j) = \sum_{y=0}^{\infty} e^{t_j y} \cdot \frac{e^{-\lambda_j} \lambda_j^y}{y!} $$
{{< /math >}} {#eq:mgf-substitution}

### Step 2: Simplify the Expression

Combine the exponentials:

{{< math >}}
$$ M_{Z_j}(t_j) = e^{-\lambda_j} \sum_{y=0}^{\infty} \frac{(\lambda_j e^{t_j})^y}{y!} $$
{{< /math >}} {#eq:mgf-combine-exp}

This is the Taylor series for the exponential function:

{{< math >}}
$$ \sum_{y=0}^{\infty} \frac{x^y}{y!} = e^x $$
{{< /math >}} {#eq:taylor-series}

with {{< math >}}$x = \lambda_j e^{t_j}${{< /math >}}. So:

{{< math >}}
$$ M_{Z_j}(t_j) = e^{-\lambda_j} \cdot e^{\lambda_j e^{t_j}} = \exp(-\lambda_j + \lambda_j e^{t_j}) $$
{{< /math >}} {#eq:mgf-taylor-apply}

### Final Answer

{{< math >}}
$$ M_{Z_j}(t_j) = \exp(\lambda_j(e^{t_j} - 1)) $$
{{< /math >}} {#eq:poisson-mgf-final}

✅ **Proved.**

## Key Insights

This result shows that the MGF of a Poisson random variable has a particularly elegant exponential form. The parameter {{< math >}}$\lambda_j${{< /math >}} appears both as a multiplicative factor and within the exponential transformation {{< math >}}$(e^{t_j} - 1)${{< /math >}}, which is fundamental to many applications in probability theory and statistics.