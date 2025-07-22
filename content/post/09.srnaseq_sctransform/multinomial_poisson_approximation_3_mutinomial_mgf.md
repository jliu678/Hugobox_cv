# Moment Generating Function of the Multinomial Distribution

## Multinomial Distribution Setup

Let {{< math >}}$\mathbf{X} = (X_1, X_2, \ldots, X_k) \sim \text{Multinomial}(n; p_1, p_2, \ldots, p_k)${{< /math >}}

where:
- {{< math >}}$X_i${{< /math >}} is the count of outcomes in category {{< math >}}$i${{< /math >}}
- {{< math >}}$n = \sum_{i=1}^k X_i${{< /math >}}
- {{< math >}}$\sum_{i=1}^k p_i = 1${{< /math >}}, {{< math >}}$p_i \geq 0${{< /math >}}
- {{< math >}}$X_i \in \mathbb{N}${{< /math >}}

## Goal

Derive the MGF of {{< math >}}$\mathbf{X}${{< /math >}}, defined as:

{{< math >}}
$$ M_{\mathbf{X}}(\mathbf{t}) = \mathbb{E}[e^{\mathbf{t}^\top \mathbf{X}}] = \mathbb{E}[e^{t_1 X_1 + \cdots + t_k X_k}] $$
{{< /math >}} {#eq:mgf-definition}

## Step 1: Probability Mass Function

The PMF of {{< math >}}$\mathbf{X} \sim \text{Mult}(n, \mathbf{p})${{< /math >}} is:

{{< math >}}
$$ P(X_1 = x_1, \ldots, X_k = x_k) = \frac{n!}{x_1! \cdots x_k!} p_1^{x_1} \cdots p_k^{x_k} $$
{{< /math >}} {#eq:multinomial-pmf}

for all {{< math >}}$x_i \geq 0${{< /math >}}, {{< math >}}$\sum x_i = n${{< /math >}}.

## Step 2: Plug into MGF Definition

We plug into the definition of the MGF:

{{< math >}}
$$ M_{\mathbf{X}}(\mathbf{t}) = \sum_{\mathbf{x}: \sum x_i = n} e^{\sum_{i=1}^k t_i x_i} \cdot \frac{n!}{x_1! \cdots x_k!} p_1^{x_1} \cdots p_k^{x_k} $$
{{< /math >}} {#eq:mgf-expansion}

Note that:

{{< math >}}
$$ e^{\sum_{i=1}^k t_i x_i} = \prod_{i=1}^k e^{t_i x_i} = \prod_{i=1}^k (e^{t_i})^{x_i} $$
{{< /math >}} {#eq:exponential-product}

So the MGF becomes:

{{< math >}}
$$ M_{\mathbf{X}}(\mathbf{t}) = \sum_{\mathbf{x}: \sum x_i = n} \frac{n!}{x_1! \cdots x_k!} (p_1 e^{t_1})^{x_1} \cdots (p_k e^{t_k})^{x_k} $$
{{< /math >}} {#eq:mgf-rearranged}

## Step 3: Recognize the Multinomial Generating Formula

This is exactly the multinomial expansion of:

{{< math >}}
$$ \left(\sum_{i=1}^k p_i e^{t_i}\right)^n $$
{{< /math >}} {#eq:multinomial-expansion}

Hence:

{{< math >}}
$$ M_{\mathbf{X}}(\mathbf{t}) = \left(\sum_{i=1}^k p_i e^{t_i}\right)^n $$
{{< /math >}} {#eq:mgf-result}

---

# The Multinomial Expansion

## 🔷 Theoretical Statement

Let {{< math >}}$y_1 + y_2 + \cdots + y_J = n${{< /math >}}, and {{< math >}}$a_1, \ldots, a_J \in \mathbb{R}_+${{< /math >}}. The multinomial identity says:

{{< math >}}
$$ \sum_{\sum_j y_j = n} \frac{n!}{y_1! \cdots y_J!} \prod_{j=1}^J a_j^{y_j} = \left(\sum_{j=1}^J a_j\right)^n $$
{{< /math >}} {#eq:multinomial-identity}

This is the multinomial expansion of {{< math >}}$(a_1 + \cdots + a_J)^n${{< /math >}}.

It generalizes the binomial identity:

{{< math >}}
$$ (a + b)^n = \sum_{k=0}^n \binom{n}{k} a^k b^{n-k} $$
{{< /math >}} {#eq:binomial-identity}

## 💡 Intuition

Imagine expanding {{< math >}}$(a_1 + a_2 + \cdots + a_J)^n${{< /math >}}. Each term in the expansion results from choosing one of {{< math >}}$a_1, \ldots, a_J${{< /math >}} {{< math >}}$n${{< /math >}} times, and multiplying the chosen ones together.

For example, choosing {{< math >}}$a_1${{< /math >}} twice and {{< math >}}$a_2${{< /math >}} once gives {{< math >}}$a_1^2 a_2${{< /math >}}, and this can happen in several different orders.

The number of such orderings is {{< math >}}$\frac{n!}{y_1! y_2! \cdots y_J!}${{< /math >}}

So the total coefficient on {{< math >}}$a_1^{y_1} a_2^{y_2} \cdots a_J^{y_J}${{< /math >}} is exactly this multinomial coefficient.

## 🔢 Numerical Example: J = 2, n = 3

Let's take:
- {{< math >}}$a_1 = 2${{< /math >}}
- {{< math >}}$a_2 = 1${{< /math >}}  
- {{< math >}}$n = 3${{< /math >}}

We'll verify:

{{< math >}}
$$ (2 + 1)^3 = \sum_{y_1 + y_2 = 3} \frac{3!}{y_1! y_2!} \cdot 2^{y_1} \cdot 1^{y_2} $$
{{< /math >}} {#eq:numerical-example}

**Left side:**
{{< math >}}$(2 + 1)^3 = 3^3 = 27${{< /math >}}

**Right side:** all combinations of {{< math >}}$y_1 + y_2 = 3${{< /math >}}

| {{< math >}}$y_1${{< /math >}} | {{< math >}}$y_2${{< /math >}} | Multinomial Coef {{< math >}}$\frac{3!}{y_1! y_2!}${{< /math >}} | {{< math >}}$2^{y_1} \cdot 1^{y_2}${{< /math >}} | Term Value |
|---|---|---|---|---|
| 0 | 3 | {{< math >}}$\frac{6}{1 \cdot 6} = 1${{< /math >}} | {{< math >}}$1 \cdot 1 = 1${{< /math >}} | 1 |
| 1 | 2 | {{< math >}}$\frac{6}{1 \cdot 2} = 3${{< /math >}} | {{< math >}}$2 \cdot 1 = 2${{< /math >}} | 6 |
| 2 | 1 | {{< math >}}$\frac{6}{2 \cdot 1} = 3${{< /math >}} | {{< math >}}$4 \cdot 1 = 4${{< /math >}} | 12 |
| 3 | 0 | {{< math >}}$\frac{6}{6 \cdot 1} = 1${{< /math >}} | {{< math >}}$8 \cdot 1 = 8${{< /math >}} | 8 |

**Total:** {{< math >}}$1 + 6 + 12 + 8 = 27${{< /math >}} ✅

Matches the left-hand side exactly.

## 🧠 Application in Probability: MGF of Multinomial

This identity is also the reason why the moment generating function (MGF) of a multinomial is:

{{< math >}}
$$ M_{\mathbf{Y}}(\mathbf{t}) = \mathbb{E}\left[e^{\sum_j t_j Y_j}\right] = \left(\sum_j \pi_j e^{t_j}\right)^n $$
{{< /math >}} {#eq:mgf-application}

Because:

{{< math >}}
$$
\begin{align}
M_{\mathbf{Y}}(\mathbf{t}) &= \sum_{\sum y_j = n} P(\mathbf{Y} = \mathbf{y}) \cdot e^{\sum_j t_j y_j} \\
&= \sum_{\sum y_j = n} \frac{n!}{\prod y_j!} \prod \pi_j^{y_j} \cdot e^{t_j y_j} \\
&= \sum_{\sum y_j = n} \frac{n!}{\prod y_j!} \prod (\pi_j e^{t_j})^{y_j} \\
&= \left(\sum_j \pi_j e^{t_j}\right)^n
\end{align}
$$
{{< /math >}} {#eq:mgf-derivation}

## Intuition Behind the Multinomial Count

After {{< math >}}$a_i${{< /math >}} and {{< math >}}$j_i${{< /math >}} have been assigned, we get {{< math >}}$n${{< /math >}} elements and each is from {{< math >}}$a.${{< /math >}} to fill the {{< math >}}$n${{< /math >}} positions.

Think of distributing labels (or outcomes, like the terms {{< math >}}$a_i${{< /math >}}) across {{< math >}}$n${{< /math >}} positions. If all entries were distinct, the number of permutations would simply be:

{{< math >}}$n!${{< /math >}}

But in the multinomial setting, some labels are repeated—there are {{< math >}}$y_1${{< /math >}} of {{< math >}}$a_1${{< /math >}}, {{< math >}}$y_2${{< /math >}} of {{< math >}}$a_2${{< /math >}}, etc. So, many permutations are indistinguishable due to identical labels.

To correct for this overcounting, we divide by the factorial of each group:

{{< math >}}
$$ \text{Number of distinct assignments} = \frac{n!}{y_1! y_2! \cdots y_J!} $$
{{< /math >}} {#eq:multinomial-coefficient}