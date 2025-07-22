# Derivation of the Poisson Approximation to the Multinomial Distribution

## 🔢 Setup

Let:
{{< math >}}$$X = (X_1, X_2, \ldots, X_k) \sim \text{Multinomial}(n, p)$${{< /math >}}

where {{< math >}}$\sum_{i=1}^k X_i = n${{< /math >}}

- Each {{< math >}}$X_i${{< /math >}} counts the number of times category {{< math >}}$i${{< /math >}} occurs out of {{< math >}}$n${{< /math >}} trials
- {{< math >}}$\sum_{i=1}^k p_i = 1${{< /math >}}, and {{< math >}}$\lambda_i = np_i${{< /math >}}

We are interested in the regime where:
- {{< math >}}$n \to \infty${{< /math >}}, and
- {{< math >}}$p_i \to 0${{< /math >}}, such that {{< math >}}$\lambda_i = np_i${{< /math >}} stays constant.

## 🧮 Multinomial PMF

The probability mass function of the multinomial distribution is:

{{< math >}}$$P(X_1 = x_1, \ldots, X_k = x_k) = \frac{n!}{x_1! x_2! \cdots x_k!} p_1^{x_1} p_2^{x_2} \cdots p_k^{x_k}$${{< /math >}}

subject to {{< math >}}$\sum_{i=1}^k x_i = n${{< /math >}}

## 🔁 Approximate Using Stirling and Limits

Set {{< math >}}$\lambda_i = np_i${{< /math >}}, so {{< math >}}$p_i = \frac{\lambda_i}{n}${{< /math >}}

Substitute:

{{< math >}}$$P(X = x) = \frac{n!}{x_1! x_2! \cdots x_k!} \prod_{i=1}^k \left(\frac{\lambda_i}{n}\right)^{x_i}$${{< /math >}}

{{< math >}}$$= \frac{n!}{n^{\sum x_i} x_1! x_2! \cdots x_k!} \prod_{i=1}^k \lambda_i^{x_i}$${{< /math >}}

{{< math >}}$$= \frac{n!}{n^n \prod x_i!} \prod_{i=1}^k \lambda_i^{x_i}$${{< /math >}}

Now apply Stirling's approximation:
{{< math >}}$$n! \approx n^n e^{-n} \sqrt{2\pi n}$${{< /math >}}

So:
{{< math >}}$$P(X = x) \approx \frac{n^n e^{-n} \sqrt{2\pi n}}{n^n \prod x_i!} \prod \lambda_i^{x_i} = e^{-n} \sqrt{2\pi n} \cdot \frac{1}{\prod x_i!} \prod \lambda_i^{x_i}$${{< /math >}}

Also recall: {{< math >}}$\sum \lambda_i = n${{< /math >}}, so:
{{< math >}}$$e^{-n} = \prod e^{-\lambda_i}$${{< /math >}}

Thus:
{{< math >}}$$P(X = x) \approx \prod_{i=1}^k \left(\frac{e^{-\lambda_i} \lambda_i^{x_i}}{x_i!}\right)$${{< /math >}}

Which is just:
{{< math >}}$$\prod_{i=1}^k \text{Poisson}(x_i; \lambda_i)$${{< /math >}}

## ✅ Conclusion

{{< math >}}$$\text{Multinomial}(n, p_i = \lambda_i/n) \xrightarrow{n \to \infty} \text{Independent Poisson}(\lambda_i)$${{< /math >}}

The multinomial distribution converges to a product of independent Poisson distributions when the number of trials grows large while the individual probabilities shrink proportionally, keeping the expected counts {{< math >}}$\lambda_i = np_i${{< /math >}} fixed.