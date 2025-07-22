# Multinomial Modeling Framework for Single-Cell RNA Sequencing Data

## 1. Observed Data

**{{< math >}}$y_{ij}${{< /math >}}**: This is the observed UMI count for **cell (or droplet)** {{< math >}}$i${{< /math >}} and **gene** {{< math >}}$j${{< /math >}}.

For example, {{< math >}}$y_{i3} = 10${{< /math >}} means 10 UMIs (counts) for gene 3 in cell {{< math >}}$i${{< /math >}}.

## 2. Total UMIs

**{{< math >}}$n_i = \sum_j y_{ij}${{< /math >}}**: This is the **total number of UMI counts in cell** {{< math >}}$i${{< /math >}} **across all genes** {{< math >}}$j${{< /math >}}.

This is a fixed number for each cell — we observed this when doing the sequencing.

## 3. True Underlying Proportions

**{{< math >}}$\pi_{ij}${{< /math >}}**: This is the **true (but unknown) relative abundance** (proportion) of gene {{< math >}}$j${{< /math >}} in cell {{< math >}}$i${{< /math >}}.

So, for each cell {{< math >}}$i${{< /math >}}, the vector {{< math >}}$\pi_i = (\pi_{i1}, \dots, \pi_{iJ})${{< /math >}} represents a probability distribution over the genes.

It satisfies:

{{< math >}}
$$ \sum_j \pi_{ij} = 1 \quad \text{and} \quad \pi_{ij} \geq 0 $$
{{< /math >}}

## 4. Multinomial Modeling

The vector of observed counts for cell {{< math >}}$i${{< /math >}}:

{{< math >}}
$$ \mathbf{y}_i = (y_{i1}, \dots, y_{iJ})^\top $$
{{< /math >}}

is assumed to follow a **multinomial distribution** with:
- **Number of trials**: {{< math >}}$n_i${{< /math >}} (total UMIs for cell {{< math >}}$i${{< /math >}})
- **Probabilities**: {{< math >}}$(\pi_{i1}, \dots, \pi_{iJ})${{< /math >}}

## 5. Multinomial Probability Mass Function

The **density (or PMF)** of a multinomial distribution is:

{{< math >}}
$$ f(\mathbf{y}_i) = \frac{n_i!}{y_{i1}!\cdots y_{iJ}!} \prod_j \pi_{ij}^{y_{ij}} $$
{{< /math >}} {#eq:multinomial-pmf}

This is what the original formula means:

{{< math >}}
$$ f(\mathbf{y}_i) = \binom{n_i}{y_{i1}, \dots, y_{iJ}} \prod_j \pi_{ij}^{y_{ij}} $$
{{< /math >}} {#eq:multinomial-coefficient}

Where:
- {{< math >}}$\binom{n_i}{y_{i1}, \dots, y_{iJ}} = \frac{n_i!}{y_{i1}!\cdots y_{iJ}!}${{< /math >}} is the **multinomial coefficient**, and
- The product over {{< math >}}$j${{< /math >}} is the probability part.

## Interpretation

This model says:

- For each cell {{< math >}}$i${{< /math >}}, we fix the total number of molecules {{< math >}}$n_i${{< /math >}}, and assume they are **randomly allocated** to genes {{< math >}}$j${{< /math >}} according to the **true relative abundances** {{< math >}}$\pi_{ij}${{< /math >}}.
- The variability in observed UMI counts {{< math >}}$y_{ij}${{< /math >}} is purely from this **sampling process**, not from changes in {{< math >}}$n_i${{< /math >}} (which is fixed).

## Key Statistical Insights

The multinomial framework treats each cell as having a fixed "library size" (total UMIs {{< math >}}$n_i${{< /math >}}) that gets randomly allocated across genes according to their true underlying proportions {{< math >}}$\pi_{ij}${{< /math >}}. This elegantly separates two sources of variation:

- **Technical sampling noise**: The randomness in which molecules get captured/sequenced
- **Biological signal**: The true gene expression proportions {{< math >}}$\pi_{ij}${{< /math >}} that we want to infer


# Marginal Binomial Distribution from Multinomial Modeling

## 🔁 Recap of Setup

From before:
- For a fixed cell {{< math >}}$i${{< /math >}}:

{{< math >}}
$$ \mathbf{y}_i = (y_{i1}, \dots, y_{iJ}) \sim \text{Multinomial}(n_i, \pi_{i1}, \dots, \pi_{iJ}) $$
{{< /math >}} {#eq:multinomial-setup}

- That is, out of the total {{< math >}}$n_i${{< /math >}} UMI counts in cell {{< math >}}$i${{< /math >}}, each count is assigned to gene {{< math >}}$j${{< /math >}} with probability {{< math >}}$\pi_{ij}${{< /math >}}.

## 🧠 Focusing on One Gene at a Time

The **marginal distribution** (i.e. just looking at {{< math >}}$y_{ij}${{< /math >}}, ignoring the other genes) is:

{{< math >}}
$$ y_{ij} \sim \text{Binomial}(n_i, \pi_{ij}) $$
{{< /math >}} {#eq:marginal-binomial}

### Why Binomial?

Because in a multinomial distribution, each **individual component** (e.g., {{< math >}}$y_{ij}${{< /math >}}) has a binomial marginal distribution:
- You're drawing {{< math >}}$n_i${{< /math >}} total UMIs.
- Each has a {{< math >}}$\pi_{ij}${{< /math >}} chance of being gene {{< math >}}$j${{< /math >}}.

## 📈 Marginal Mean

{{< math >}}
$$ \mathbb{E}[y_{ij}] = n_i \pi_{ij} = \mu_{ij} $$
{{< /math >}} {#eq:marginal-mean}

This is intuitive: on average, gene {{< math >}}$j${{< /math >}} accounts for {{< math >}}$\pi_{ij}${{< /math >}} proportion of {{< math >}}$n_i${{< /math >}} total counts.

## 📊 Marginal Variance

{{< math >}}
$$ \text{Var}(y_{ij}) = n_i \pi_{ij}(1 - \pi_{ij}) = \mu_{ij} - \frac{1}{n_i} \mu_{ij}^2 $$
{{< /math >}} {#eq:marginal-variance}

Notice:
- This variance is **less than** a Poisson variance {{< math >}}$\mu_{ij}${{< /math >}}, because of the {{< math >}}$- \mu_{ij}^2/n_i${{< /math >}} term.
- That's a key trait of multinomial sampling: the **counts are negatively correlated** across genes, and there's a **fixed total** {{< math >}}$n_i${{< /math >}}, which "compresses" the variance.

## 0️⃣ Probability of a Zero Count

The probability that **gene** {{< math >}}$j${{< /math >}} gets **zero counts** in cell {{< math >}}$i${{< /math >}}:

{{< math >}}
$$ \mathbb{P}(y_{ij} = 0) = (1 - \pi_{ij})^{n_i} = \left(1 - \frac{\mu_{ij}}{n_i} \right)^{n_i} $$
{{< /math >}} {#eq:zero-probability}

This approximation shows:
- Even if the true expression (i.e., {{< math >}}$\mu_{ij}${{< /math >}}) is small but nonzero, you still might see a zero count — **due to random sampling**.
- This is one cause of **dropouts** in scRNA-seq data — not necessarily true zero expression, but stochastic sampling from a multinomial process with small probabilities.

## 📌 Summary

| Quantity | Expression | Interpretation |
|----------|------------|----------------|
| **Marginal dist.** | {{< math >}}$y_{ij} \sim \text{Binomial}(n_i, \pi_{ij})${{< /math >}} | Counts for one gene follow binomial |
| **Mean** | {{< math >}}$\mu_{ij} = n_i \pi_{ij}${{< /math >}} | Expected count |
| **Variance** | {{< math >}}$\mu_{ij} - \frac{1}{n_i} \mu_{ij}^2${{< /math >}} | Lower than Poisson |
| **Zero prob.** | {{< math >}}$\left(1 - \frac{\mu_{ij}}{n_i} \right)^{n_i}${{< /math >}} | Chance of dropout |

## Key Insights

The marginal binomial framework reveals several important properties of single-cell RNA-seq data:

1. **Underdispersion**: Variance is lower than Poisson due to the fixed library size constraint
2. **Negative correlation**: Genes compete for the fixed pool of {{< math >}}$n_i${{< /math >}} UMIs
3. **Dropout mechanism**: Zero counts can arise from sampling even when true expression is non-zero


# 📐 Derivation of Covariance in Multinomial Distribution

The known formula for the multinomial is:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \begin{cases} 
n_i \pi_{ij}(1-\pi_{ij}) & \text{if } j=k \\
-n_i \pi_{ij} \pi_{ik} & \text{if } j \neq k 
\end{cases} $$
{{< /math >}} {#eq:multinomial-covariance}

Here's the intuition for the off-diagonal case {{< math >}}$j \neq k${{< /math >}}:

## Think of each UMI draw

You draw {{< math >}}$n_i${{< /math >}} UMIs, one at a time. Each draw is:
- gene {{< math >}}$j${{< /math >}} with probability {{< math >}}$\pi_{ij}${{< /math >}}
- gene {{< math >}}$k${{< /math >}} with probability {{< math >}}$\pi_{ik}${{< /math >}}
- and so on

The total counts {{< math >}}$y_{ij}${{< /math >}} and {{< math >}}$y_{ik}${{< /math >}} are just sums over these draws.

**Now the key point:**
If one UMI goes to gene {{< math >}}$j${{< /math >}}, it cannot go to gene {{< math >}}$k${{< /math >}} (because it's assigned to only one gene). This makes gene counts negatively dependent.

Mathematically, this competition causes:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = -n_i \pi_{ij} \pi_{ik} $$
{{< /math >}} {#eq:negative-covariance}

## 🎯 Why Exactly That Formula?

Because each of the {{< math >}}$n_i${{< /math >}} trials is:
- independent,
- assigns probability {{< math >}}$\pi_{ij}${{< /math >}} to gene {{< math >}}$j${{< /math >}},

and across trials, the number of co-occurrences is:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \sum_{\text{trials}} \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) $$
{{< /math >}} {#eq:sum-covariances}

where {{< math >}}$Z_{ij}^{(t)} = 1${{< /math >}} if the {{< math >}}$t${{< /math >}}-th UMI goes to gene {{< math >}}$j${{< /math >}}, else 0.

Then:

{{< math >}}
$$ \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) = -\pi_{ij}\pi_{ik} $$
{{< /math >}} {#eq:single-trial-covariance}

Sum over {{< math >}}$n_i${{< /math >}} trials: gives {{< math >}}$-n_i \pi_{ij} \pi_{ik}${{< /math >}}

## Why the covariance between {{< math >}}$y_{ij}${{< /math >}} and {{< math >}}$y_{ik}${{< /math >}} is the sum of covariances between individual trials:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \sum_{t=1}^{n_i} \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) $$
{{< /math >}} {#eq:covariance-decomposition}

Let's break this down.

## 🧠 Setup: Breaking the Multinomial into Bernoulli Trials

The multinomial distribution is just the sum of {{< math >}}$n_i${{< /math >}} independent categorical trials. Each trial assigns one UMI to one of {{< math >}}$J${{< /math >}} genes.

Define the indicator variables:

{{< math >}}
$$ Z_{ij}^{(t)} = 1 \text{ if the } t\text{-th UMI is assigned to gene } j, \text{ and 0 otherwise.} $$
{{< /math >}} {#eq:indicator-variables}

So each {{< math >}}$y_{ij} = \sum_{t=1}^{n_i} Z_{ij}^{(t)}${{< /math >}}

These {{< math >}}$Z_{ij}^{(t)}${{< /math >}} are not random draws over {{< math >}}$t${{< /math >}} — we're thinking of the {{< math >}}$t${{< /math >}}-th draw as a Bernoulli trial that outputs a gene according to probabilities {{< math >}}$\pi_{i1}, \ldots, \pi_{iJ}${{< /math >}}.

## 📘 Covariance Between Sums

Now suppose:

{{< math >}}
$$ y_{ij} = \sum_{t=1}^{n_i} Z_{ij}^{(t)}, \quad y_{ik} = \sum_{s=1}^{n_i} Z_{ik}^{(s)} $$
{{< /math >}} {#eq:sums-definition}

Then:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \text{Cov}\left(\sum_{t=1}^{n_i} Z_{ij}^{(t)}, \sum_{s=1}^{n_i} Z_{ik}^{(s)}\right) = \sum_{t=1}^{n_i} \sum_{s=1}^{n_i} \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(s)}) $$
{{< /math >}} {#eq:bilinearity-application}

Now separate the cases:

1. **When {{< math >}}$t \neq s${{< /math >}}**: the draws are independent, so the covariance is 0.

2. **When {{< math >}}$t = s${{< /math >}}**: you have:

{{< math >}}
$$ \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) = \mathbb{E}[Z_{ij}^{(t)} Z_{ik}^{(t)}] - \mathbb{E}[Z_{ij}^{(t)}] \mathbb{E}[Z_{ik}^{(t)}] $$
{{< /math >}} {#eq:covariance-definition}

Since {{< math >}}$Z_{ij}^{(t)} Z_{ik}^{(t)} = 0${{< /math >}} (because a UMI can go to either gene {{< math >}}$j${{< /math >}} or gene {{< math >}}$k${{< /math >}}, but not both), the product is always 0. Therefore:

{{< math >}}
$$ \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) = 0 - \pi_{ij}\pi_{ik} = -\pi_{ij}\pi_{ik} $$
{{< /math >}} {#eq:single-covariance-result}

So the only nonzero terms in the double sum come when {{< math >}}$t = s${{< /math >}}, and there are {{< math >}}$n_i${{< /math >}} such terms. So:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \sum_{t=1}^{n_i} (-\pi_{ij}\pi_{ik}) = -n_i \pi_{ij} \pi_{ik} $$
{{< /math >}} {#eq:final-covariance}

## Why is below decomposition valid?

{{< math >}}
$$ \text{Cov}\left(\sum_{t=1}^{n_i} Z_{ij}^{(t)}, \sum_{s=1}^{n_i} Z_{ik}^{(s)}\right) = \sum_{t=1}^{n_i} \sum_{s=1}^{n_i} \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(s)}) $$
{{< /math >}} {#eq:decomposition-question}


## ✅ Yes — This Identity is Always Valid

This follows from the **bilinearity of covariance**. Specifically, for any two finite sums of random variables:

{{< math >}}
$$ \text{Cov}\left(\sum_a X_a, \sum_b Y_b\right) = \sum_a \sum_b \text{Cov}(X_a, Y_b) $$
{{< /math >}} {#eq:bilinearity-general}

This is a general and well-known property of covariance.

## 🔍 Proof Sketch (using bilinearity)

Let's recall the definition of covariance:

{{< math >}}
$$ \text{Cov}(X,Y) = \mathbb{E}[XY] - \mathbb{E}[X]\mathbb{E}[Y] $$
{{< /math >}} {#eq:covariance-basic-definition}

So for two sums:

{{< math >}}
$$ \text{Cov}\left(\sum_t X_t, \sum_s Y_s\right) = \mathbb{E}\left[\left(\sum_t X_t\right)\left(\sum_s Y_s\right)\right] - \mathbb{E}\left[\sum_t X_t\right]\mathbb{E}\left[\sum_s Y_s\right] $$
{{< /math >}} {#eq:covariance-sums-expanded}

Use distributivity of expectation and multiplication:

{{< math >}}
$$ = \sum_t \sum_s \mathbb{E}[X_t Y_s] - \sum_t \mathbb{E}[X_t] \sum_s \mathbb{E}[Y_s] = \sum_t \sum_s (\mathbb{E}[X_t Y_s] - \mathbb{E}[X_t]\mathbb{E}[Y_s]) = \sum_t \sum_s \text{Cov}(X_t, Y_s) $$
{{< /math >}} {#eq:bilinearity-proof}

## 🔁 Apply It to Your Case

So for your multinomial indicator variables:
- {{< math >}}$y_{ij} = \sum_{t=1}^{n_i} Z_{ij}^{(t)}${{< /math >}}
- {{< math >}}$y_{ik} = \sum_{s=1}^{n_i} Z_{ik}^{(s)}${{< /math >}}

Then:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = \text{Cov}\left(\sum_t Z_{ij}^{(t)}, \sum_s Z_{ik}^{(s)}\right) = \sum_t \sum_s \text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(s)}) $$
{{< /math >}} {#eq:final-application}

And you then use:
- {{< math >}}$\text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(s)}) = 0${{< /math >}} if {{< math >}}$t \neq s${{< /math >}} (independent trials)
- {{< math >}}$\text{Cov}(Z_{ij}^{(t)}, Z_{ik}^{(t)}) = -\pi_{ij}\pi_{ik}${{< /math >}} if {{< math >}}$t = s${{< /math >}}

So:

{{< math >}}
$$ \text{Cov}(y_{ij}, y_{ik}) = n_i \cdot (-\pi_{ij}\pi_{ik}) = -n_i \pi_{ij} \pi_{ik} $$
{{< /math >}} {#eq:final-result}


# $ \text{Var}(y_{ij}) = ? $

## 🔁 Background Setup

Let:

- {{< math >}}$y_{ij}${{< /math >}} = number of UMIs in cell {{< math >}}$i${{< /math >}} that map to gene {{< math >}}$j${{< /math >}}
- {{< math >}}$n_i${{< /math >}} = total number of UMI counts (draws)
- {{< math >}}$\pi_{ij}${{< /math >}} = probability that a UMI maps to gene {{< math >}}$j${{< /math >}} in cell {{< math >}}$i${{< /math >}}

Then:

{{< math >}}
$$ \mathbf{y}_i = (y_{i1}, \ldots, y_{iJ}) \sim \text{Multinomial}(n_i, \pi_{i1}, \ldots, \pi_{iJ}) $$
{{< /math >}} {#eq:multinomial-distribution}

By definition:

- {{< math >}}$y_{ij}${{< /math >}} is the number of times gene {{< math >}}$j${{< /math >}} is chosen among {{< math >}}$n_i${{< /math >}} draws.
- Each draw is independent and has probability {{< math >}}$\pi_{ij}${{< /math >}} of selecting gene {{< math >}}$j${{< /math >}}.

## 🧱 Step 1: Express {{< math >}}$y_{ij}${{< /math >}} as a sum of indicator variables

Let's define:

{{< math >}}
$$ Z_{ij}^{(t)} = \begin{cases} 
1 & \text{if the } t\text{-th UMI maps to gene } j \\
0 & \text{otherwise}
\end{cases} $$
{{< /math >}} {#eq:indicator-variable}

Then:

{{< math >}}
$$ y_{ij} = \sum_{t=1}^{n_i} Z_{ij}^{(t)} $$
{{< /math >}} {#eq:sum-of-indicators}

So now, we've written {{< math >}}$y_{ij}${{< /math >}} as a sum of {{< math >}}$n_i${{< /math >}} Bernoulli trials.

Each {{< math >}}$Z_{ij}^{(t)} \sim \text{Bernoulli}(\pi_{ij})${{< /math >}}, and all {{< math >}}$Z_{ij}^{(t)}${{< /math >}} are independent.

## 🧮 Step 2: Use variance rules

When you have a sum of independent random variables:

{{< math >}}
$$ \text{Var}\left(\sum_t X_t\right) = \sum_t \text{Var}(X_t) $$
{{< /math >}} {#eq:variance-sum-rule}

Since the {{< math >}}$Z_{ij}^{(t)}${{< /math >}} are independent:

{{< math >}}
$$ \text{Var}(y_{ij}) = \sum_{t=1}^{n_i} \text{Var}(Z_{ij}^{(t)}) $$
{{< /math >}} {#eq:variance-decomposition}

Each {{< math >}}$Z_{ij}^{(t)}${{< /math >}} has:

{{< math >}}
$$ \text{Var}(Z_{ij}^{(t)}) = \pi_{ij}(1 - \pi_{ij}) $$
{{< /math >}} {#eq:bernoulli-variance}

So:

{{< math >}}
$$ \text{Var}(y_{ij}) = n_i \cdot \pi_{ij}(1 - \pi_{ij}) $$
{{< /math >}} {#eq:final-variance}

## 📊 Final Result

The variance of the UMI count for gene {{< math >}}$j${{< /math >}} in cell {{< math >}}$i${{< /math >}} under the multinomial model is:

{{< math >}}
$$ \boxed{\text{Var}(y_{ij}) = n_i \pi_{ij}(1 - \pi_{ij})} $$
{{< /math >}} {#eq:boxed-result}

This shows that the variance depends on:
- **Library size** ({{< math >}}$n_i${{< /math >}}): Larger libraries have higher variance
- **Gene expression level** ({{< math >}}$\pi_{ij}${{< /math >}}): Maximum variance occurs at {{< math >}}$\pi_{ij} = 0.5${{< /math >}}
- **Sampling noise**: The {{< math >}}$(1 - \pi_{ij})${{< /math >}} term captures the binomial-like uncertainty