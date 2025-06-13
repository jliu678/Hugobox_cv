---
title: 📊 Math Derivation in Bayesian Inference-- Bypassing Marginals Using Conjugate Priors 
summary: Conjugate priors are a powerful concept in Bayesian statistics that allow us to simplify the process of updating our beliefs about a parameter given new data.  
date: 2025-04-12
authors:
  - admin
tags:
  - Bayesian
  - Conjugate Priors
  - Marginal Distribution
  - Math derivation
image:
  caption: 'Image credit: [**Karan Suthar on Unsplash**](https://unsplash.com)'
---
## Introduction
Conjugate priors are a powerful concept in Bayesian statistics that allow us to simplify the process of updating our beliefs about a parameter given new data. Using conjugate priors allows you to bypass computing marginal distributions because the posterior has the same functional form as the prior, making the normalization constant analytically tractable.

#### Standard Bayes' Theorem

{{< math >}} $$P(\theta|D) = \frac{P(D|\theta)P(\theta)}{P(D)}$$ {{< /math >}}

where the marginal likelihood is:

{{< math >}} $$P(D) = \int P(D|\theta)P(\theta)d\theta$$ {{< /math >}}

#### The Problem

Computing {{< math >}} $P(D)$ {{< /math >}} often requires intractable integrals, especially in high dimensions.

#### Conjugate Prior Solution

When the prior {{< math >}} $P(\theta)$ {{< /math >}} is conjugate to the likelihood {{< math >}} $P(D|\theta)$ {{< /math >}}, we have:

1. **Prior:** {{< math >}} $P(\theta) \propto f(\theta; \alpha_0)$ {{< /math >}} (some parametric form)
2. **Likelihood:** {{< math >}} $P(D|\theta) \propto g(\theta; D)$ {{< /math >}}
3. **Posterior:** {{< math >}} $P(\theta|D) \propto f(\theta; \alpha_n)$ {{< /math >}} (same form as prior)

where {{< math >}} $\alpha_n = h(\alpha_0, D)$ {{< /math >}} is a simple update function. Specifically, since we know the posterior belongs to the same family as the prior, we can write:

{{< math >}} $$P(\theta|D) = \frac{f(\theta; \alpha_n)}{Z(\alpha_n)}$$ {{< /math >}}

where {{< math >}} $Z(\alpha_n)$ {{< /math >}} is the normalization constant for the known distribution family, which has a closed form.

## Example: Beta-Binomial Conjugacy

- **Prior:** {{< math >}} $P(\theta) = \text{Beta}(\alpha, \beta)$ {{< /math >}}
- **Likelihood:** {{< math >}} $P(D|\theta) = \text{Binomial}(n, \theta)$ {{< /math >}} with {{< math >}} $s$ {{< /math >}} successes
- **Posterior:** {{< math >}} $P(\theta|D) = \text{Beta}(\alpha + s, \beta + n - s)$ {{< /math >}}

We never need to compute:

{{< math >}} $$P(D) = \int_0^1 \binom{n}{s}\theta^s(1-\theta)^{n-s} \cdot \frac{\theta^{\alpha-1}(1-\theta)^{\beta-1}}{B(\alpha,\beta)} d\theta$$ {{< /math >}}

Instead, we directly get the posterior parameters and use the known Beta normalization. This explains why conjugate priors are so powerful-- they let us skip the hardest part of Bayesian computation. Let me break down what's happening:

#### The Hard Way (Without Conjugate Priors)

To get the posterior, we'd normally need to compute this integral:

{{< math >}} $$P(D) = \int_0^1 \binom{n}{s}\theta^s(1-\theta)^{n-s} \cdot \frac{\theta^{\alpha-1}(1-\theta)^{\beta-1}}{B(\alpha,\beta)} d\theta$$ {{< /math >}}

This integral combines:
- {{< math >}} $\binom{n}{s}\theta^s(1-\theta)^{n-s}$ {{< /math >}} (binomial likelihood)
- {{< math >}} $\frac{\theta^{\alpha-1}(1-\theta)^{\beta-1}}{B(\alpha,\beta)}$ {{< /math >}} (Beta prior)

#### Why This Integral is Nasty

Even though this particular integral has a closed form, in general such integrals:
- May not have analytical solutions
- Require numerical integration (expensive, approximate)
- Get exponentially harder in higher dimensions

#### The Conjugate Prior Shortcut

Instead of computing that integral, we use the conjugate relationship:

1. **Recognize the pattern:**
   {{< math >}} $$\theta^s(1-\theta)^{n-s} \times \theta^{\alpha-1}(1-\theta)^{\beta-1} = \theta^{(s+\alpha)-1}(1-\theta)^{(n-s+\beta)-1} \quad (1)$${{< /math >}}

   This looks like a Beta distribution with parameters {{< math >}} $(\alpha + s, \beta + n - s)$ {{< /math >}}
   
   We derive below in [next section](#math-derivation) how this leads to the marginal likelihood:

   {{< math >}} $$P(D) = \frac{B(\alpha+s, \beta+n-s)}{B(\alpha, \beta)} \binom{n}{s} \quad (2)$$ {{< /math >}}

2. **Write the posterior directly:**

   We substitute equation (1) and (2) in Standard Bayes' Theorem to get the posterior distribution:

   {{< math >}} $$P(\theta|D) = \text{Beta}(\alpha + s, \beta + n - s) = \frac{\theta^{\alpha+s-1}(1-\theta)^{\beta+n-s-1}}{B(\alpha+s, \beta+n-s)}$$ {{< /math >}}



But for getting the posterior distribution, we don't even need {{< math >}} $P(D)$ {{< /math >}} - we just need the updated parameters!

#### Interim Summary

- *Without conjugacy:* Solve a potentially intractable integral
- *With conjugacy:* Simple parameter update: {{< math >}} $(\alpha, \beta) \rightarrow (\alpha + s, \beta + n - s)$ {{< /math >}}

The normalization "just works" because we're staying within the same distributional family where normalization constants are known.


## Math Derivation<a id="math-derivation"></a>

Below is a **step-by-step derivation of the marginal likelihood integral for the Beta-Binomial conjugate prior relationship**. This shows how we can transform a potentially difficult integral into a simple ratio of Beta functions, which is the essence of why conjugate priors are so powerful in Bayesian inference.

### Starting Point: The Marginal Likelihood Integral

{{< math >}} $$P(D) = \int_0^1 \binom{n}{s}\theta^s(1-\theta)^{n-s} \cdot \frac{\theta^{\alpha-1}(1-\theta)^{\beta-1}}{B(\alpha,\beta)} d\theta \quad (1)$$ {{< /math >}}

### Step 1: Factor Out Constants

The binomial coefficient {{< math >}} $\binom{n}{s}$ {{< /math >}} and {{< math >}} $\frac{1}{B(\alpha,\beta)}$ {{< /math >}} don't depend on {{< math >}} $\theta$ {{< /math >}}, so we can pull them outside the integral:

{{< math >}} $$P(D) = \binom{n}{s} \cdot \frac{1}{B(\alpha,\beta)} \int_0^1 \theta^s(1-\theta)^{n-s} \cdot \theta^{\alpha-1}(1-\theta)^{\beta-1} d\theta \quad (2)$$ {{< /math >}}

### Step 2: Combine the Powers

Using the exponent rule {{< math >}} $x^a \cdot x^b = x^{a+b}$ {{< /math >}}:

{{< math >}} $$P(D) = \binom{n}{s} \cdot \frac{1}{B(\alpha,\beta)} \int_0^1 \theta^{s+\alpha-1}(1-\theta)^{n-s+\beta-1} d\theta \quad (3)$$ {{< /math >}}

### Step 3: Recognize the Beta Function Integral

The integral {{< math >}} $\int_0^1 \theta^{a-1}(1-\theta)^{b-1} d\theta$ {{< /math >}} is exactly the definition of the Beta function {{< math >}} $B(a,b)$ {{< /math >}}.

In our case, we have:
- {{< math >}} $a = s + \alpha$ {{< /math >}}
- {{< math >}} $b = n - s + \beta$ {{< /math >}}

So:
{{< math >}} $$\int_0^1 \theta^{s+\alpha-1}(1-\theta)^{n-s+\beta-1} d\theta = B(s+\alpha, n-s+\beta) = B(\alpha+s, \beta+n-s) \quad (4)$$ {{< /math >}}

### Step 4: Substitute Back

{{< math >}} $$P(D) = \binom{n}{s} \cdot \frac{1}{B(\alpha,\beta)} \cdot B(\alpha+s, \beta+n-s) \quad (5)$$ {{< /math >}}

### Final Result

Rearranging:
{{< math >}} $$P(D) = \frac{B(\alpha+s, \beta+n-s)}{B(\alpha,\beta)} \binom{n}{s} \quad (6)$$ {{< /math >}}

### Conclusion

This is why conjugate priors are so computationally elegant: they transform intractable integrals into simple ratios of known functions. Specifically for the Beta-Binomial case:
1. **We avoided numerical integration** - instead of computing a potentially difficult integral, we used the known relationship between integrals and Beta functions
2. **Both Beta functions have closed forms** - {{< math >}} $B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$ {{< /math >}}
3. **This gives us the exact marginal likelihood** - no approximation needed

