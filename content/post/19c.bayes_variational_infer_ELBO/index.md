---
title: 🎓 Math Derivation for ELBO/KL in Bayesian Inference and VAEs
summary: Marginal likelihood is often intractable in Bayesian Inference. Thus ELBO is used as a lower bound to the log marginal likelihood. I derives the ELBO and KL divergence in Bayesian inference, and shows how they are used in Variational Autoencoders (VAEs).
date: 2025-04-24
authors:
  - admin
tags:
  - Neural Network Variational Autoencoders
  - ELBO
  - KL Divergence
  - Bayesian
  - Marginal Distribution
image:
  caption: 'Image credit: [**Shawn Day on Unsplash**](https://unsplash.com)'
---
## 🔍 Marginal likelihood is often intractable in Bayesian Inference
Bayesian inference and Variational Autoencoders (VAEs), a marriage of deep learning and Bayesian inference, are powerful and foundational for probabilistic modeling used in computational biology. But rely on the concept of marginal likelihood, which is often intractable to compute directly. This is because both frameworks typically involve integrating over latent variables or parameters, which can lead to high-dimensional integrals that are computationally intractable.

Specifically, in Bayesian inference, we want the posterior distribution:

{{< math >}}
$$p(\theta \mid D) = \frac{p(D \mid \theta) \, p(\theta)}{p(D)} \tag{1}$$
{{< /math >}}

But the denominator—the evidence or marginal likelihood—is:

{{< math >}}
$$p(D) = \int p(D \mid \theta) \, p(\theta) \, d\theta \tag{2}$$
{{< /math >}}

This integral is often intractable.

## 💡 Variational Inference: The Workaround

We introduce a simpler variational distribution {{< math >}}$q(\theta)${{< /math >}} to approximate {{< math >}}$p(\theta \mid D)${{< /math >}}, and we try to make {{< math >}}$q${{< /math >}} close to the true posterior.

We measure closeness using KL divergence:

{{< math >}}
$$\text{KL}(q(\theta) \| p(\theta \mid D)) = \int q(\theta) \log \frac{q(\theta)}{p(\theta \mid D)} d\theta \tag{3}$$
{{< /math >}}

This is hard to compute directly because it involves {{< math >}}$p(D)${{< /math >}}, so we rearrange terms.

We can rewrite {{< math >}}$\log p(D)${{< /math >}} as:

{{< math >}}
$$\log p(D) = \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] + \text{KL}(q(\theta) \| p(\theta \mid D)) \tag{4}$$
{{< /math >}}

Thus:

{{< math >}}
$$\log p(D) = \text{ELBO}(q) + \text{KL}(q(\theta) \| p(\theta \mid D)) \tag{5}$$
{{< /math >}}

Since the KL divergence is always ≥ 0:

{{< math >}}
$$\text{ELBO}(q) \leq \log p(D) \tag{6}$$
{{< /math >}}

That's why it's called a lower bound.

## 🧮 Derive ELBO and KL

Let {{< math >}}$q(\theta)${{< /math >}} be any distribution over {{< math >}}$\theta${{< /math >}} such that its support covers that of {{< math >}}$p(\theta \mid D)${{< /math >}}. We'll exploit a classic trick: insert {{< math >}}$q(\theta)${{< /math >}} into the log marginal likelihood using expectation and apply properties of KL divergence.

### Step 1: Start with log evidence

We take the logarithm of {{< math >}}$p(D)${{< /math >}}, and "multiply and divide" inside by {{< math >}}$q(\theta)${{< /math >}}:

{{< math >}}
$$\log p(D) = \log \int \frac{q(\theta)}{q(\theta)} p(D \mid \theta) p(\theta) d\theta \tag{7}$$
{{< /math >}}

{{< math >}}
$$= \log \int q(\theta) \cdot \frac{p(D \mid \theta) p(\theta)}{q(\theta)} d\theta \tag{8}$$
{{< /math >}}

{{< math >}}
$$= \log \mathbb{E}_{q(\theta)} \left[ \frac{p(D \mid \theta) p(\theta)}{q(\theta)} \right] \tag{9}$$
{{< /math >}}

This is Jensen's inequality time.

### Step 2: Apply Jensen's Inequality

{{< math >}}
$$\log \mathbb{E}_{q(\theta)} \left[ \frac{p(D \mid \theta) p(\theta)}{q(\theta)} \right] \geq \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D \mid \theta) p(\theta)}{q(\theta)} \right] \tag{10}$$
{{< /math >}}

That gives us the ELBO:

{{< math >}}
$$\text{ELBO}(q) = \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] \tag{11}$$
{{< /math >}}

So:

{{< math >}}
$$\log p(D) \geq \text{ELBO}(q) \tag{12}$$
{{< /math >}}

But we can go further — let's rewrite {{< math >}}$\log p(D)${{< /math >}} exactly in terms of ELBO + KL divergence.

### Step 3: Add and Subtract the Same Quantity

We now write:

{{< math >}}
$$\log p(D) = \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] + \left( \log p(D) - \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] \right) \tag{13}$$
{{< /math >}}

Now we observe that the term in parentheses is exactly the KL divergence between {{< math >}}$q(\theta)${{< /math >}} and the true posterior:

{{< math >}}
$$\text{KL}(q(\theta) \| p(\theta \mid D)) = \mathbb{E}_{q(\theta)} \left[ \log \frac{q(\theta)}{p(\theta \mid D)} \right] \tag{14}$$
{{< /math >}}

But recall:

{{< math >}}
$$p(\theta \mid D) = \frac{p(D, \theta)}{p(D)} \Rightarrow \log p(\theta \mid D) = \log p(D, \theta) - \log p(D) \tag{15}$$
{{< /math >}}

Then:

{{< math >}}
$$\log \frac{q(\theta)}{p(\theta \mid D)} = \log \frac{q(\theta)}{p(D, \theta)} + \log p(D) \tag{16}$$
{{< /math >}}

Take expectation over {{< math >}}$q(\theta)${{< /math >}}:

{{< math >}}
$$\text{KL}(q(\theta) \| p(\theta \mid D)) = -\mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] + \log p(D) \tag{17}$$
{{< /math >}}

Rearranged:

{{< math >}}
$$\log p(D) = \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] + \text{KL}(q(\theta) \| p(\theta \mid D)) \tag{18}$$
{{< /math >}}

### Definition of Expectation used above
Note the above derivations used multiple times the definition that **expectation of a function {{< math >}}$f(\theta)${{< /math >}} under a probability distribution {{< math >}}$q(\theta)${{< /math >}}** is:

{{< math >}}
$$\mathbb{E}_{q(\theta)}[f(\theta)] = \int q(\theta) \, f(\theta) \, d\theta \tag{19}$$
{{< /math >}}

## 📐 ELBO Expression used in VAEs

{{< math >}}
$$\text{ELBO}(q) = \mathbb{E}_{q(\theta)}[\log p(D \mid \theta)] - \text{KL}(q(\theta) \| p(\theta)) \tag{20}$$
{{< /math >}}


This is the most widely used form in variational inference and VAEs. It comes from expanding the joint $p(D, \theta)$, and interpreting the ELBO as a trade-off between reconstruction and regularization.

**Interpretation:**

- The first term encourages {{< math >}}$q(\theta)${{< /math >}} to explain the data well.
- The second term encourages {{< math >}}$q(\theta)${{< /math >}} to stay close to the prior.

Its derivation start from:

### 1. ELBO–KL decomposition:

{{< math >}}
$$\log p(D) = \text{ELBO}(q) + \text{KL}(q(\theta) \| p(\theta \mid D)) \tag{5}$$
{{< /math >}}

This is always true by the definition of the Kullback-Leibler divergence and Jensen's inequality. Rearranging:

{{< math >}}
$$\text{ELBO}(q) = \log p(D) - \text{KL}(q(\theta) \| p(\theta \mid D)) \tag{22}$$
{{< /math >}}

### 2. Definition of ELBO via expected joint:

Alternatively, ELBO is often defined as:

{{< math >}}
$$\text{ELBO}(q) = \mathbb{E}_{q(\theta)} \left[ \log \frac{p(D, \theta)}{q(\theta)} \right] = \mathbb{E}_{q(\theta)}[\log p(D, \theta)] - \mathbb{E}_{q(\theta)}[\log q(\theta)] \tag{23}$$
{{< /math >}}

Now recall:

{{< math >}}
$$\log p(D, \theta) = \log p(D \mid \theta) + \log p(\theta) \tag{24}$$
{{< /math >}}

So:

{{< math >}}
$$\text{ELBO}(q) = \mathbb{E}_{q(\theta)}[\log p(D \mid \theta)] + \mathbb{E}_{q(\theta)}[\log p(\theta)] - \mathbb{E}_{q(\theta)}[\log q(\theta)] \tag{25}$$
{{< /math >}}

Group terms:

{{< math >}}
$$\text{ELBO}(q) = \mathbb{E}_{q(\theta)}[\log p(D \mid \theta)] - \text{KL}(q(\theta) \| p(\theta)) \tag{20}$$
{{< /math >}}

