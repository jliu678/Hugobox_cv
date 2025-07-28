# NOTE
the derivation did not show Dirichlet-Multinomial is Equivalent to Negative Binomial without condition!

# Dirichlet-Multinomial Equivalent to Negative Binomial Conditioned on Sum

A Dirichlet-Multinomial distribution can indeed be equivalently characterized as a collection of independent Negative Binomial distributions with the same success probability ({{< math >}}$p${{< /math >}}), when conditioned on their sum. This property is a fundamental result in compound probability distributions and is widely used in statistical modeling, particularly for overdispersed count data.

Here's the proof, demonstrating the equivalence of the joint probability mass functions (PMFs) under the specified conditions.

## Proof of Equivalence

We aim to show that the joint distribution of {{< math >}}$(x_1, \ldots, x_K)${{< /math >}} from independent Negative Binomials, when conditioned on their sum {{< math >}}$x_{\cdot} = \sum x_k${{< /math >}}, is equivalent to the Dirichlet-Multinomial distribution.

Let's define the two scenarios:

### Scenario 1: Independent Negative Binomials Conditioned on Their Sum

Let {{< math >}}$x_k \sim \text{NB}(r_k, p)${{< /math >}} independently for {{< math >}}$k = 1, \ldots, K${{< /math >}}.

The PMF for {{< math >}}$x_k${{< /math >}} (number of failures before {{< math >}}$r_k${{< /math >}} successes with success probability {{< math >}}$p${{< /math >}}) is:

{{< math >}}
$$P(x_k) = \binom{x_k + r_k - 1}{x_k} p^{r_k} (1-p)^{x_k} \quad (1)$$
{{< /math >}}

Since the {{< math >}}$x_k${{< /math >}} are independent, their joint PMF is:

{{< math >}}
$$P(x_1, \ldots, x_K) = \prod_{k=1}^K \left[\binom{x_k + r_k - 1}{x_k} p^{r_k} (1-p)^{x_k}\right] \quad (2)$$
{{< /math >}}

Let {{< math >}}$r_{\cdot} = \sum_{k=1}^K r_k${{< /math >}} and {{< math >}}$x_{\cdot} = \sum_{k=1}^K x_k${{< /math >}}.

{{< math >}}
$$P(x_1, \ldots, x_K) = \left(\prod_{k=1}^K \binom{x_k + r_k - 1}{x_k}\right) p^{r_{\cdot}} (1-p)^{x_{\cdot}} \quad (3)$$
{{< /math >}}

Now, let {{< math >}}$x_{\cdot} = \sum_{k=1}^K x_k${{< /math >}}. A key property of Negative Binomial distributions is that the sum of independent Negative Binomial random variables with the same success probability {{< math >}}$p${{< /math >}} is also a Negative Binomial random variable.

So, {{< math >}}$x_{\cdot} \sim \text{NB}(r_{\cdot}, p)${{< /math >}}.

The PMF for {{< math >}}$x_{\cdot}${{< /math >}} is:

{{< math >}}
$$P(x_{\cdot}) = \binom{x_{\cdot} + r_{\cdot} - 1}{x_{\cdot}} p^{r_{\cdot}} (1-p)^{x_{\cdot}} \quad (4)$$
{{< /math >}}

We are interested in the conditional distribution of {{< math >}}$(x_1, \ldots, x_K)${{< /math >}} given their sum {{< math >}}$x_{\cdot}${{< /math >}}:

{{< math >}}
$$P(x_1, \ldots, x_K \mid x_{\cdot}) = \frac{P(x_1, \ldots, x_K \text{ and } x_{\cdot})}{P(x_{\cdot})} \quad (5)$$
{{< /math >}}

Since {{< math >}}$x_{\cdot}${{< /math >}} is simply the sum of {{< math >}}$x_1, \ldots, x_K${{< /math >}}, the numerator {{< math >}}$P(x_1, \ldots, x_K \text{ and } x_{\cdot})${{< /math >}} is just {{< math >}}$P(x_1, \ldots, x_K)${{< /math >}} itself (provided {{< math >}}$\sum x_k = x_{\cdot}${{< /math >}}).

{{< math >}}
$$P(x_1, \ldots, x_K \mid x_{\cdot}) = \frac{\left(\prod_{k=1}^K \binom{x_k + r_k - 1}{x_k}\right) p^{r_{\cdot}} (1-p)^{x_{\cdot}}}{\binom{x_{\cdot} + r_{\cdot} - 1}{x_{\cdot}} p^{r_{\cdot}} (1-p)^{x_{\cdot}}} \quad (6)$$
{{< /math >}}

The terms {{< math >}}$p^{r_{\cdot}}${{< /math >}} and {{< math >}}$(1-p)^{x_{\cdot}}${{< /math >}} cancel out, leaving:

{{< math >}}
$$P(x_1, \ldots, x_K \mid x_{\cdot}) = \frac{\prod_{k=1}^K \binom{x_k + r_k - 1}{x_k}}{\binom{x_{\cdot} + r_{\cdot} - 1}{x_{\cdot}}} \quad (7)$$
{{< /math >}}

We can express the binomial coefficients using Gamma functions: {{< math >}}$\binom{n}{k} = \frac{\Gamma(n+1)}{\Gamma(k+1)\Gamma(n-k+1)}${{< /math >}}.

More generally, {{< math >}}$\binom{x+r-1}{x} = \frac{\Gamma(x+r)}{\Gamma(x+1)\Gamma(r)}${{< /math >}}.

So,

{{< math >}}
$$P(x_1, \ldots, x_K \mid x_{\cdot}) = \frac{\Gamma(x_{\cdot} + 1)\Gamma(r_{\cdot})}{\Gamma(x_{\cdot} + r_{\cdot})} \prod_{k=1}^K \frac{\Gamma(x_k + r_k)}{\Gamma(x_k + 1)\Gamma(r_k)} \quad (8)$$
{{< /math >}}

Rearranging the terms:

{{< math >}}
$$P(x_1, \ldots, x_K \mid x_{\cdot}) = \frac{x_{\cdot}!}{\prod_{k=1}^K x_k!} \cdot \frac{\Gamma(r_{\cdot})}{\Gamma(x_{\cdot} + r_{\cdot})} \cdot \prod_{k=1}^K \frac{\Gamma(x_k + r_k)}{\Gamma(r_k)} \quad (9)$$
{{< /math >}}

This is the PMF of a Dirichlet-Multinomial distribution with parameters {{< math >}}$n = x_{\cdot}${{< /math >}} (total count) and {{< math >}}$\alpha = (r_1, \ldots, r_K)${{< /math >}} (concentration parameters).

### Scenario 2: Dirichlet-Multinomial Distribution

A random vector {{< math >}}$(y_1, \ldots, y_K)${{< /math >}} follows a Dirichlet-Multinomial distribution with parameters {{< math >}}$n${{< /math >}} (total trials) and {{< math >}}$\alpha = (\alpha_1, \ldots, \alpha_K)${{< /math >}} if its PMF is:

{{< math >}}
$$P(y_1, \ldots, y_K) = \frac{n!}{y_1! \cdots y_K!} \cdot \frac{\Gamma(\sum \alpha_k)}{\Gamma(n + \sum \alpha_k)} \cdot \prod_{k=1}^K \frac{\Gamma(y_k + \alpha_k)}{\Gamma(\alpha_k)} \quad (10)$$
{{< /math >}}

In our context, if we set {{< math >}}$n = x_{\cdot}${{< /math >}} and {{< math >}}$\alpha_k = r_k${{< /math >}}, then the PMF becomes:

{{< math >}}
$$P(y_1, \ldots, y_K) = \frac{x_{\cdot}!}{y_1! \cdots y_K!} \cdot \frac{\Gamma(\sum r_k)}{\Gamma(x_{\cdot} + \sum r_k)} \cdot \prod_{k=1}^K \frac{\Gamma(y_k + r_k)}{\Gamma(r_k)} \quad (11)$$
{{< /math >}}

Since {{< math >}}$r_{\cdot} = \sum r_k${{< /math >}} and {{< math >}}$y_{\cdot} = \sum y_k = x_{\cdot}${{< /math >}}:

{{< math >}}
$$P(y_1, \ldots, y_K) = \frac{x_{\cdot}!}{y_1! \cdots y_K!} \cdot \frac{\Gamma(r_{\cdot})}{\Gamma(x_{\cdot} + r_{\cdot})} \cdot \prod_{k=1}^K \frac{\Gamma(y_k + r_k)}{\Gamma(r_k)} \quad (12)$$
{{< /math >}}

## Conclusion

Comparing the PMF derived from conditioning independent Negative Binomials on their sum (equation 9) with the PMF of the Dirichlet-Multinomial distribution (equation 12), we see they are identical when {{< math >}}$y_k = x_k${{< /math >}}.

This proves that a Dirichlet-Multinomial distribution is equivalent to a collection of independent Negative Binomial distributions with the same scale parameter ({{< math >}}$p${{< /math >}}) conditioned on their sum.

This equivalence is powerful because it provides two different generative processes for the same multivariate count data distribution, which can be useful for modeling and inference. It highlights that the overdispersion inherent in the Dirichlet-Multinomial can arise from underlying independent Negative Binomial processes, and the constraint of a fixed total sum then gives it the specific Dirichlet-Multinomial form.



# Probability Intersection and Implication: P(A and B) = P(A) When A Implies B

The numerator {{< math >}}$P(x_1, \ldots, x_K \text{ and } x_{\cdot})${{< /math >}} is equivalent to {{< math >}}$P(x_1, \ldots, x_K)${{< /math >}} because the event "{{< math >}}$(x_1, \ldots, x_K)${{< /math >}} occurs" implicitly **implies** that "their sum {{< math >}}$x_{\cdot} = \sum x_k${{< /math >}} occurs." In probability, if event A implies event B, then the probability of A and B occurring is simply the probability of A occurring, i.e., {{< math >}}$P(A \cap B) = P(A)${{< /math >}}.

## Detailed Explanation

Let's break this down:

### 1. Define Event A
Let {{< math >}}$A${{< /math >}} be the event that the specific values {{< math >}}$(x_1, \ldots, x_K)${{< /math >}} occur. This means {{< math >}}$X_1 = x_1, X_2 = x_2, \ldots, X_K = x_K${{< /math >}}. The probability of this event is:

{{< math >}}
$$P(A) = P(x_1, \ldots, x_K) \quad (1)$$
{{< /math >}}

### 2. Define Event B
Let {{< math >}}$B${{< /math >}} be the event that their sum {{< math >}}$x_{\cdot}${{< /math >}} occurs. This means {{< math >}}$\sum_{k=1}^K X_k = x_{\cdot}${{< /math >}}. The probability of this event is:

{{< math >}}
$$P(B) = P(x_{\cdot}) \quad (2)$$
{{< /math >}}

### 3. The Implication
If event {{< math >}}$A${{< /math >}} occurs (i.e., we have specific values for {{< math >}}$x_1, \ldots, x_K${{< /math >}}), then it automatically follows that their sum {{< math >}}$\sum_{k=1}^K x_k${{< /math >}} must equal some value, which we've defined as {{< math >}}$x_{\cdot}${{< /math >}}. Therefore, event {{< math >}}$A${{< /math >}} (the specific values) **implies** event {{< math >}}$B${{< /math >}} (the specific sum). In set notation, the set of outcomes corresponding to event {{< math >}}$A${{< /math >}} is a subset of the set of outcomes corresponding to event {{< math >}}$B${{< /math >}} (i.e., {{< math >}}$A \subseteq B${{< /math >}}).

This implication can be written as:

{{< math >}}
$$A \Rightarrow B \quad (3)$$
{{< /math >}}

### 4. Probability of Intersection
The probability {{< math >}}$P(A \text{ and } B)${{< /math >}}, or {{< math >}}$P(A \cap B)${{< /math >}}, represents the probability that both events {{< math >}}$A${{< /math >}} and {{< math >}}$B${{< /math >}} occur. Since {{< math >}}$A${{< /math >}} implies {{< math >}}$B${{< /math >}}, whenever {{< math >}}$A${{< /math >}} occurs, {{< math >}}$B${{< /math >}} must also occur. Thus, the outcomes where {{< math >}}$A${{< /math >}} and {{< math >}}$B${{< /math >}} both occur are exactly the outcomes where {{< math >}}$A${{< /math >}} occurs. Therefore:

{{< math >}}
$$P(A \cap B) = P(A) \quad (4)$$
{{< /math >}}

## Applying This to Our Specific Case

{{< math >}}$P(x_1, \ldots, x_K \text{ and } x_{\cdot})${{< /math >}} means "the vector {{< math >}}$(X_1, \ldots, X_K)${{< /math >}} takes on the values {{< math >}}$(x_1, \ldots, x_K)${{< /math >}} **AND** the sum of these values is {{< math >}}$x_{\cdot}${{< /math >}}."

Since it's necessarily true that if {{< math >}}$(X_1, \ldots, X_K)${{< /math >}} takes on values {{< math >}}$(x_1, \ldots, x_K)${{< /math >}}, then their sum is {{< math >}}$\sum x_k = x_{\cdot}${{< /math >}}, the second part of the "AND" condition provides no additional restriction or information. It's redundant.

Hence:

{{< math >}}
$$P(x_1, \ldots, x_K \text{ and } x_{\cdot}) = P(x_1, \ldots, x_K) \quad (5)$$
{{< /math >}}

## Mathematical Formulation

More formally, we can express this relationship as:

{{< math >}}
$$P\left(\bigcap_{k=1}^K \{X_k = x_k\} \cap \left\{\sum_{k=1}^K X_k = x_{\cdot}\right\}\right) = P\left(\bigcap_{k=1}^K \{X_k = x_k\}\right) \quad (6)$$
{{< /math >}}

This equality holds because:

{{< math >}}
$$\bigcap_{k=1}^K \{X_k = x_k\} \subseteq \left\{\sum_{k=1}^K X_k = x_{\cdot}\right\} \quad (7)$$
{{< /math >}}

whenever {{< math >}}$x_{\cdot} = \sum_{k=1}^K x_k${{< /math >}}.

## Conclusion

This fundamental principle of probability theory—that {{< math >}}$P(A \cap B) = P(A)${{< /math >}} when {{< math >}}$A \Rightarrow B${{< /math >}}—is what allows us to simplify the numerator in the conditional probability calculation, making the derivation of the Dirichlet-Multinomial equivalence proof both elegant and mathematically rigorous.