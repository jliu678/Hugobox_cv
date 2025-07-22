# Rigorous Justification: Poisson Limit of the Multinomial Distribution

## Theorem Statement

We rigorously justify that:

{{< math >}}
$$ \text{Multinomial}\left(n, \pi_1 = \frac{\lambda_1}{n}, \ldots, \pi_J = \frac{\lambda_J}{n}\right) \xrightarrow{d} \prod_{j=1}^J \text{Poisson}(\lambda_j) \quad \text{as } n \to \infty $$
{{< /math >}} {#eq:poisson-limit-theorem}

This is known as the **Poisson limit of the multinomial**.

We'll prove this using moment generating functions (MGFs) and asymptotic analysis.

## 🧮 Step-by-Step Proof via MGF (Moment Generating Function)

We'll use the fact that convergence in distribution can be shown by proving convergence of the joint moment generating function (MGF), provided the limit MGF characterizes a unique distribution.

### 🔧 Setup

Let 

{{< math >}}
$$ \mathbf{Y}^{(n)} = (Y_1^{(n)}, \ldots, Y_J^{(n)}) \sim \text{Multinomial}\left(n, \frac{\lambda_1}{n}, \ldots, \frac{\lambda_J}{n}\right) $$
{{< /math >}} {#eq:multinomial-setup}

Let {{< math >}}$\mathbf{t} = (t_1, \ldots, t_J) \in \mathbb{R}^J${{< /math >}}. The joint MGF of {{< math >}}$\mathbf{Y}^{(n)}${{< /math >}} is:

{{< math >}}
$$ M_{\mathbf{Y}^{(n)}}(\mathbf{t}) = \mathbb{E}\left[\exp\left(\sum_{j=1}^J t_j Y_j^{(n)}\right)\right] $$
{{< /math >}} {#eq:mgf-definition}

### 📐 Step 1: MGF of the Multinomial

The MGF of the multinomial distribution with total {{< math >}}$n${{< /math >}} and probabilities {{< math >}}$\pi_j${{< /math >}} is:

{{< math >}}
$$ M_{\mathbf{Y}^{(n)}}(\mathbf{t}) = \left(\sum_{j=1}^J \pi_j e^{t_j}\right)^n $$
{{< /math >}} {#eq:multinomial-mgf}

In our case, {{< math >}}$\pi_j = \frac{\lambda_j}{n}${{< /math >}}, so:

{{< math >}}
$$ M_{\mathbf{Y}^{(n)}}(\mathbf{t}) = \left(\sum_{j=1}^J \frac{\lambda_j}{n} e^{t_j}\right)^n = \left(\frac{1}{n}\sum_{j=1}^J \lambda_j e^{t_j}\right)^n $$
{{< /math >}} {#eq:multinomial-mgf-substituted}

Let's factor out the {{< math >}}$n${{< /math >}}-dependence:

{{< math >}}
$$ M_{\mathbf{Y}^{(n)}}(\mathbf{t}) = \left(1 + \frac{1}{n}\left(\sum_j \lambda_j e^{t_j} - n\right)\right)^n $$
{{< /math >}} {#eq:mgf-factored}

Note that {{< math >}}$\sum_j \lambda_j = n${{< /math >}}. 



Using the fundamental limit {{< math >}}$(1 + \frac{x}{n})^n \to e^x${{< /math >}} as {{< math >}}$n \to \infty${{< /math >}}:

{{< math >}}
$$ M_{\mathbf{Y}^{(n)}}(\mathbf{t}) \to \exp\left(\sum_{j=1}^J \lambda_j (e^{t_j} - 1)\right) \quad \text{as } n \to \infty $$
{{< /math >}} {#eq:limiting-mgf}

### 📌 Step 2: Limiting MGF Matches Independent Poissons

Now, recall the MGF of a Poisson random variable {{< math >}}$Z_j \sim \text{Poisson}(\lambda_j)${{< /math >}}:

{{< math >}}
$$ M_{Z_j}(t_j) = \exp(\lambda_j (e^{t_j} - 1)) $$
{{< /math >}} {#eq:poisson-mgf}

So the joint MGF of independent Poisson variables {{< math >}}$Z_j \sim \text{Poisson}(\lambda_j)${{< /math >}} is:

{{< math >}}
$$ M_{\mathbf{Z}}(\mathbf{t}) = \prod_{j=1}^J \exp(\lambda_j (e^{t_j} - 1)) = \exp\left(\sum_{j=1}^J \lambda_j (e^{t_j} - 1)\right) $$
{{< /math >}} {#eq:independent-poisson-mgf}

This matches the limit of the multinomial MGF from equation {{< math >}}$(11)${{< /math >}}.

### ✅ Step 3: Conclude Distributional Convergence

Since:

1. The MGF of {{< math >}}$\mathbf{Y}^{(n)}${{< /math >}} converges pointwise to that of {{< math >}}$\mathbf{Z} = (Z_1, \ldots, Z_J)${{< /math >}}, a vector of independent Poisson variables;

2. The limit MGF corresponds to a proper distribution;

3. MGFs uniquely determine distribution (in a neighborhood of zero);

We conclude:

{{< math >}}
$$ \mathbf{Y}^{(n)} \xrightarrow{d} \mathbf{Z} \sim \prod_{j=1}^J \text{Poisson}(\lambda_j) \quad \text{as } n \to \infty $$
{{< /math >}} {#eq:convergence-conclusion}

That is, the multinomial converges in distribution to a product of independent Poisson distributions, as {{< math >}}$n \to \infty${{< /math >}} and {{< math >}}$\pi_j = \frac{\lambda_j}{n}${{< /math >}}.

## 🔁 Why This Works: Intuition

- The multinomial enforces the sum {{< math >}}$\sum y_j = n${{< /math >}}, which creates negative correlations.

- But when {{< math >}}$n \to \infty${{< /math >}} and each {{< math >}}$\pi_j \to 0${{< /math >}}, individual events are rare and approximately independent — this is the **law of rare events**.

- Hence the joint distribution behaves like independent Poissons.

## 📚 Summary

| Quantity | Multinomial | Poisson Approximation |
|----------|-------------|----------------------|
| **Mean** | {{< math >}}$n\pi_j = \lambda_j${{< /math >}} | {{< math >}}$\lambda_j${{< /math >}} |
| **Variance** | {{< math >}}$n\pi_j(1-\pi_j) \to \lambda_j${{< /math >}} | {{< math >}}$\lambda_j${{< /math >}} |
| **Covariance** ({{< math >}}$j \neq k${{< /math >}}) | {{< math >}}$-n\pi_j\pi_k \to 0${{< /math >}} | {{< math >}}$0${{< /math >}} |
| **MGF** | {{< math >}}$\left(\sum_j \pi_j e^{t_j}\right)^n \to e^{\sum_j \lambda_j (e^{t_j} - 1)}${{< /math >}} | {{< math >}}$e^{\sum_j \lambda_j (e^{t_j} - 1)}${{< /math >}} |

The key insight is that as the number of trials increases while the individual success probabilities decrease proportionally, the multinomial distribution converges to a product of independent Poisson distributions. This is fundamental to understanding why Poisson models work well for count data with large sample sizes and small event probabilities.