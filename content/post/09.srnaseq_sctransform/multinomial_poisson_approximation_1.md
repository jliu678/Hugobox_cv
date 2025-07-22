# 🧠 Why Multinomial is Inconvenient for Some Models like GLMs

The multinomial distribution models counts across multiple categories under the constraint that the total sum is fixed (i.e., all outcomes must add up to {{< math >}}$n${{< /math >}}). While this is natural for things like sequencing counts or proportions, it creates challenges for generalized linear models (GLMs) and related statistical tools.

## 🚫 Challenges of the Multinomial in GLMs

### Dependent Outcomes
In a multinomial distribution, the counts across categories are not independent—if one goes up, another must go down to maintain the fixed total.

This inherent correlation complicates modeling, especially when trying to assign independent predictors (like gene expression or experimental covariates) to each category.

### No Canonical Link Function
GLMs rely on exponential family distributions and convenient link functions (like log for Poisson, logit for Binomial). The multinomial doesn't have a natural link function for modeling each category separately.

### Parameter Constraints
The multinomial probabilities {{< math >}}$\pi_j${{< /math >}} must sum to 1, which introduces additional constraints in the optimization. That makes estimation harder compared to unconstrained GLM fitting.

## ✅ Why Poisson is Used Instead (Poisson Trick)

To make things easier, statisticians often approximate the multinomial model by treating the counts as independent Poisson random variables:

{{< math >}}
$$ y_{ij} \sim \text{Poisson}(\mu_{ij}) $$
{{< /math >}} {#eq:poisson-approx}

Why is this helpful?

### No Fixed Sum Constraint
Poisson counts are independent and not constrained to sum to anything. This simplifies modeling.

### Easier GLM Fitting
Poisson is part of the exponential family with a convenient log link. This fits naturally into the GLM framework.

### Good Approximation for Large Counts
When the total count {{< math >}}$n_i${{< /math >}} is large, and probabilities {{< math >}}$\pi_j${{< /math >}} are small, the Poisson approximation to the multinomial is mathematically valid:

{{< math >}}
$$ y_{ij} \sim \text{Multinomial}(n_i, \pi_{ij}) \approx \text{Poisson}(n_i \pi_{ij}) $$
{{< /math >}} {#eq:multinomial-poisson-approx}

This works particularly well in high-throughput sequencing data, where total counts per sample are large.

## 📌 Summary

- **Multinomial**: Realistic but hard to work with in GLMs (due to dependence and constraints).
- **Poisson**: Simpler, flexible, and closely approximates multinomial under common conditions.

Therefore, many models "pretend" the counts are Poisson, even when they come from a multinomial.

# 🔍 Prove multinomial can be well approximated by poisson

We want to show:

If

{{< math >}}
$$ (y_{i1}, y_{i2}, \ldots, y_{iJ}) \sim \text{Multinomial}(n_i, \pi_{i1}, \pi_{i2}, \ldots, \pi_{iJ}) $$
{{< /math >}} {#eq:multinomial-goal}

then for large {{< math >}}$n_i${{< /math >}} and small {{< math >}}$\pi_{ij}${{< /math >}}, each {{< math >}}$y_{ij}${{< /math >}} can be well-approximated by

{{< math >}}
$$ y_{ij} \sim \text{Poisson}(\lambda_{ij} = n_i \pi_{ij}) $$
{{< /math >}} {#eq:poisson-goal}

and furthermore, the joint distribution becomes approximately a product of independent Poisson distributions.

We'll break this down into:

## 1: Multinomial Definition

The multinomial probability mass function (PMF) is:

{{< math >}}
$$ P(y_{i1}, \ldots, y_{iJ}) = \frac{n_i!}{y_{i1}! \cdots y_{iJ}!} \prod_{j=1}^J \pi_{ij}^{y_{ij}} $$
{{< /math >}} {#eq:multinomial-pmf}

subject to {{< math >}}$\sum_j y_{ij} = n_i${{< /math >}}

So the counts are dependent — they must sum to {{< math >}}$n_i${{< /math >}}.

## 2: Poisson Joint Distribution

Let's consider independent Poisson random variables:

{{< math >}}
$$ z_{i1}, \ldots, z_{iJ} \sim \text{Poisson}(\lambda_{ij} = n_i \pi_{ij}) $$
{{< /math >}} {#eq:independent-poisson}

The joint distribution of independent Poissons is:

{{< math >}}
$$ P(z_{i1}, \ldots, z_{iJ}) = \prod_{j=1}^J \frac{e^{-\lambda_{ij}} \lambda_{ij}^{z_{ij}}}{z_{ij}!} $$
{{< /math >}} {#eq:poisson-joint}

**Note**: These do not require that {{< math >}}$\sum_j z_{ij} = n_i${{< /math >}}.

## 3: Conditional Distribution of Poissons Given Sum is Multinomial

Let's now condition the independent Poisson variables on the event that their sum equals {{< math >}}$n_i${{< /math >}}. That is:

{{< math >}}
$$ (y_{i1}, \ldots, y_{iJ}) \mid \left(\sum_{j=1}^J z_{ij} = n_i\right) $$
{{< /math >}} {#eq:conditional-setup}

This conditional distribution is exactly the multinomial:

{{< math >}}
$$ (y_{i1}, \ldots, y_{iJ}) \mid \sum_j y_{ij} = n_i \sim \text{Multinomial}(n_i, \pi_{ij}) $$
{{< /math >}} {#eq:conditional-multinomial}

with {{< math >}}$\pi_{ij} = \frac{\lambda_{ij}}{\sum_k \lambda_{ik}} = \frac{n_i \pi_{ij}}{n_i} = \pi_{ij}${{< /math >}}

Let's formally prove that the multinomial distribution arises as the conditional distribution of independent Poisson random variables conditioned on their sum.

### Overview of Formal Proof

Let {{< math >}}$Z_1, Z_2, \ldots, Z_J \sim \text{i.i.d. Poisson}(\lambda_j)${{< /math >}}

Then:

{{< math >}}
$$ (Y_1, \ldots, Y_J) := (Z_1, \ldots, Z_J) \mid \sum_{j=1}^J Z_j = n \text{ follows } \text{Multinomial}(n; \pi_1, \ldots, \pi_J) $$
{{< /math >}} {#eq:goal}

with

{{< math >}}
$$ \pi_j = \frac{\lambda_j}{\sum_{k=1}^J \lambda_k} $$
{{< /math >}} {#eq:pi-definition}

We'll prove this by:

1. Writing the joint distribution of independent Poisson variables.
2. Conditioning this joint distribution on the total count {{< math >}}$\sum Z_j = n${{< /math >}}.
3. Showing the resulting conditional PMF matches the multinomial PMF.

### Step 1: Joint PMF of Independent Poissons

Let {{< math >}}$Z_1, \ldots, Z_J${{< /math >}} be independent with {{< math >}}$Z_j \sim \text{Poisson}(\lambda_j)${{< /math >}}, so:

{{< math >}}
$$ P(Z_1 = y_1, \ldots, Z_J = y_J) = \prod_{j=1}^J \frac{e^{-\lambda_j} \lambda_j^{y_j}}{y_j!} $$
{{< /math >}} {#eq:joint-poisson}

This can be rewritten as:

{{< math >}}
$$ P(Z_1 = y_1, \ldots, Z_J = y_J) = \left(\prod_{j=1}^J \frac{\lambda_j^{y_j}}{y_j!}\right) \cdot \exp\left(-\sum_{j=1}^J \lambda_j\right) $$
{{< /math >}} {#eq:joint-poisson-factored}

Let {{< math >}}$\mathbf{y} = (y_1, \ldots, y_J)${{< /math >}}, and suppose {{< math >}}$\sum_j y_j = n${{< /math >}}. Then the probability of this particular configuration conditional on the sum being {{< math >}}$n${{< /math >}} is:

{{< math >}}
$$ P(\mathbf{Y} = \mathbf{y} \mid \sum_j Z_j = n) = \frac{P(Z_1 = y_1, \ldots, Z_J = y_J)}{P(\sum_j Z_j = n)} $$
{{< /math >}} {#eq:conditional-definition}

### Step 2: Express the Conditional PMF

So:

{{< math >}}
$$ P(\mathbf{Y} = \mathbf{y} \mid \sum Z_j = n) = \frac{\left(\prod_{j=1}^J \frac{\lambda_j^{y_j}}{y_j!}\right) e^{-\sum \lambda_j}}{P(\sum Z_j = n)} $$
{{< /math >}} {#eq:conditional-numerator}

We now calculate the denominator {{< math >}}$P(\sum Z_j = n)${{< /math >}}.

Let {{< math >}}$S = \sum_j Z_j${{< /math >}}. Since the sum of independent Poisson random variables is also Poisson with parameter {{< math >}}$\Lambda = \sum_j \lambda_j${{< /math >}}, we have:

{{< math >}}
$$ S \sim \text{Poisson}\left(\sum_{j=1}^J \lambda_j\right) \Rightarrow P(S = n) = \frac{e^{-\sum \lambda_j} \left(\sum \lambda_j\right)^n}{n!} $$
{{< /math >}} {#eq:sum-poisson}

### Step 3: Substitute and Simplify

Substitute back into equation {{< math >}}$\eqref{eq:conditional-numerator}${{< /math >}}:

{{< math >}}
$$
\begin{align}
P(\mathbf{Y} = \mathbf{y} \mid \sum Z_j = n) &= \frac{\left(\prod_{j=1}^J \frac{\lambda_j^{y_j}}{y_j!}\right) e^{-\sum \lambda_j}}{\frac{e^{-\sum \lambda_j} \left(\sum \lambda_j\right)^n}{n!}} \\[0.5em]
&= \frac{n!}{\prod_j y_j!} \cdot \prod_j \left(\frac{\lambda_j}{\sum_k \lambda_k}\right)^{y_j}
\end{align}
$$
{{< /math >}} {#eq:conditional-simplified}

Let {{< math >}}$\pi_j = \frac{\lambda_j}{\sum_k \lambda_k}${{< /math >}}, then:

{{< math >}}
$$ P(\mathbf{Y} = \mathbf{y} \mid \sum Z_j = n) = \frac{n!}{\prod_j y_j!} \cdot \prod_j \pi_j^{y_j} $$
{{< /math >}} {#eq:final-result}

### Step 4: Recognition

Equation {{< math >}}$\eqref{eq:final-result}${{< /math >}} is exactly the PMF of the multinomial distribution:

{{< math >}}
$$ \text{Multinomial}(n; \pi_1, \ldots, \pi_J) $$
{{< /math >}} {#eq:multinomial-result}

where the multinomial coefficient is:

{{< math >}}
$$ \binom{n}{y_1, \ldots, y_J} = \frac{n!}{y_1! \cdots y_J!} $$
{{< /math >}} {#eq:multinomial-coefficient}

### 🎉 Conclusion

We have proven that conditioning independent Poisson random variables {{< math >}}$Z_j \sim \text{Poisson}(\lambda_j)${{< /math >}} on their sum {{< math >}}$\sum_j Z_j = n${{< /math >}} yields a multinomial distribution with:

- **Number of trials**: {{< math >}}$n${{< /math >}}
- **Success probabilities**: {{< math >}}$\pi_j = \frac{\lambda_j}{\sum_k \lambda_k}${{< /math >}}

This result provides a fundamental connection between Poisson and multinomial distributions and explains why multinomial models are natural for count data with fixed totals.

### Biological Interpretation

In the context of single-cell RNA sequencing:

- **{{< math >}}$Z_j${{< /math >}}**: Represents the "true" number of transcripts for gene {{< math >}}$j${{< /math >}} that would be captured under ideal conditions
- **{{< math >}}$\lambda_j${{< /math >}}**: The expected capture rate for gene {{< math >}}$j${{< /math >}} (reflects true expression level)
- **Conditioning on {{< math >}}$\sum Z_j = n${{< /math >}}**: Accounts for the fact that we observe a fixed library size {{< math >}}$n${{< /math >}} per cell
- **{{< math >}}$\pi_j${{< /math >}}**: The relative expression proportion of gene {{< math >}}$j${{< /math >}} within the cell

This mathematical framework justifies using multinomial models for scRNA-seq count data with fixed library sizes.