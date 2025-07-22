# Dirichlet-Multinomial ≈ Independent Negative Binomials

## Step 1: Dirichlet-Multinomial PMF

Let {{< math >}} $\mathbf{y} = (y_1, \dots, y_J) \sim \text{DirMult}(n, \boldsymbol\alpha)$ {{< /math >}}, where {{< math >}} $\sum_j y_j = n$ {{< /math >}}, and {{< math >}} $\boldsymbol\alpha = (\alpha_1, \dots, \alpha_J)$ {{< /math >}}, with {{< math >}} $\alpha_0 = \sum_j \alpha_j$ {{< /math >}}. Then:

{{< math >}}
$$P(\mathbf{y}) = \frac{n!}{\prod_j y_j!} \cdot \frac{\Gamma(\alpha_0)}{\Gamma(n + \alpha_0)} \cdot \prod_j \frac{\Gamma(y_j + \alpha_j)}{\Gamma(\alpha_j)} \tag{1}$$
{{< /math >}}

Taking logs:

{{< math >}}
$$\log P(\mathbf{y}) = \log n! - \sum_j \log y_j! + \log \Gamma(\alpha_0) - \log \Gamma(n + \alpha_0) + \sum_j \log \Gamma(y_j + \alpha_j) - \sum_j \log \Gamma(\alpha_j) \tag{2}$$
{{< /math >}}

## Step 2: Stirling's Approximation

Recall the full Stirling approximation:

{{< math >}}
$$\log \Gamma(z) \approx z \log z - z + \tfrac{1}{2} \log(2\pi) + \tfrac{1}{2} \log z \tag{3}$$
{{< /math >}}

In particular:

{{< math >}}
$$\log n! = \log \Gamma(n+1) \approx (n+1) \log(n+1) - (n+1) + \tfrac{1}{2} \log(2\pi) + \tfrac{1}{2} \log(n+1) \tag{4}$$
{{< /math >}}

But in large {{< math >}} $n$ {{< /math >}}, we approximate:

{{< math >}}
$$\boxed{ \log n! \approx n \log n - n + \tfrac{1}{2} \log(2\pi n) } \tag{5}$$
{{< /math >}}

This comes from:

{{< math >}}
$$\log \Gamma(n+1) = \log n! \approx \left(n + \tfrac{1}{2} \right) \log n - n + \tfrac{1}{2} \log(2\pi) = n \log n - n + \tfrac{1}{2} \log(2\pi n) \tag{6}$$
{{< /math >}}

## Step 3: Apply Stirling to All Terms

Apply the approximation to every term:

**3.1** {{< math >}} $\log n! \approx n \log n - n + \tfrac{1}{2} \log(2\pi n)$ {{< /math >}}

**3.2** {{< math >}} $\sum_j \log y_j! \approx \sum_j \left[ y_j \log y_j - y_j + \tfrac{1}{2} \log(2\pi y_j) \right]$ {{< /math >}}

**3.3** {{< math >}} $\log \Gamma(n + \alpha_0) \approx (n + \alpha_0) \log(n + \alpha_0) - (n + \alpha_0) + \tfrac{1}{2} \log(2\pi(n + \alpha_0))$ {{< /math >}}

**3.4** {{< math >}} $\sum_j \log \Gamma(y_j + \alpha_j) \approx \sum_j \left[ (y_j + \alpha_j) \log(y_j + \alpha_j) - (y_j + \alpha_j) + \tfrac{1}{2} \log(2\pi(y_j + \alpha_j)) \right]$ {{< /math >}}

**3.5** {{< math >}} $\sum_j \log \Gamma(\alpha_j) \approx \sum_j \left[ \alpha_j \log \alpha_j - \alpha_j + \tfrac{1}{2} \log(2\pi \alpha_j) \right]$ {{< /math >}}

## Step 4: Combine and Cancel Terms

Plug everything back into {{< math >}} $\log P(\mathbf{y})$ {{< /math >}} and group terms.

Let's expand and combine:

{{< math >}}
$$
\begin{align}
\log P(\mathbf{y}) &\approx [n \log n - n + \tfrac{1}{2} \log(2\pi n)] - \sum_j [y_j \log y_j - y_j + \tfrac{1}{2} \log(2\pi y_j)] \\
&\quad + \log \Gamma(\alpha_0) - [(n + \alpha_0) \log(n + \alpha_0) - (n + \alpha_0) + \tfrac{1}{2} \log(2\pi(n + \alpha_0))] \\
&\quad + \sum_j \left[(y_j + \alpha_j) \log(y_j + \alpha_j) - (y_j + \alpha_j) + \tfrac{1}{2} \log(2\pi(y_j + \alpha_j)) \right] \\
&\quad - \sum_j \left[ \alpha_j \log \alpha_j - \alpha_j + \tfrac{1}{2} \log(2\pi \alpha_j) \right] \tag{7}
\end{align}
$$
{{< /math >}}

### Cancel Constants Directly

1. **Cancel additive constants from factorial terms**: Since {{< math >}} $\sum_j y_j = n$ {{< /math >}}, we have:

{{< math >}}
$$-n + \sum_j y_j = 0 \quad \text{and} \quad -\sum_j \alpha_j + \alpha_0 = 0 \tag{8}$$
{{< /math >}}

So all standalone linear {{< math >}} $-n$ {{< /math >}}, {{< math >}} $-y_j$ {{< /math >}}, {{< math >}} $-\alpha_j$ {{< /math >}}, {{< math >}} $-\alpha_0$ {{< /math >}} terms cancel.

2. **Group and cancel all {{< math >}} $\log(2\pi \cdot)$ {{< /math >}} terms**: These include:

{{< math >}}
$$\tfrac{1}{2} \log(2\pi n) - \sum_j \tfrac{1}{2} \log(2\pi y_j) - \tfrac{1}{2} \log(2\pi(n + \alpha_0)) + \sum_j \tfrac{1}{2} \log(2\pi(y_j + \alpha_j)) - \sum_j \tfrac{1}{2} \log(2\pi \alpha_j) \tag{9}$$
{{< /math >}}

This full expression is of order:

{{< math >}}
$$\boxed{o(n)} \quad \text{(i.e., grows slower than linearly in } n \text{)} \tag{10}$$
{{< /math >}}

✅ **What does asymptotically {{< math >}} $o(n)$ {{< /math >}} mean?** It means these terms **vanish relative to {{< math >}} $n$ {{< /math >}}** as {{< math >}} $n \to \infty$ {{< /math >}}, i.e.

{{< math >}}
$$\frac{o(n)}{n} \to 0 \text{ as } n \to \infty \tag{11}$$
{{< /math >}}

So for large {{< math >}} $n$ {{< /math >}}, these have **negligible** effect on likelihood comparisons and can be dropped.

## Step 5: Main Terms After Cancellation

After canceling the constants, we are left with:

{{< math >}}
$$
\begin{align}
\log P(\mathbf{y}) &\approx n \log n - \sum_j y_j \log y_j - (n + \alpha_0) \log(n + \alpha_0) \\
&\quad + \sum_j (y_j + \alpha_j) \log(y_j + \alpha_j) - \sum_j \alpha_j \log \alpha_j + \log \Gamma(\alpha_0) \tag{12}
\end{align}
$$
{{< /math >}}

## Step 6: Expand Log Terms via Taylor Approximation

Focus on the key expansion:
{{< math >}}
$$(y_j + \alpha_j) \log(y_j + \alpha_j) \tag{13}$$
{{< /math >}}

Use Taylor expansion at {{< math >}} $y_j \gg \alpha_j$ {{< /math >}}:

{{< math >}}
$$\log(y_j + \alpha_j) = \log y_j + \frac{\alpha_j}{y_j} - \frac{1}{2} \left( \frac{\alpha_j}{y_j} \right)^2 + \cdots \tag{14}$$
{{< /math >}}

Then multiply:

{{< math >}}
$$
\begin{align}
(y_j + \alpha_j) \log(y_j + \alpha_j) &= y_j \log y_j + \alpha_j \log y_j + \alpha_j \\
&\quad + \underbrace{\alpha_j \cdot \frac{\alpha_j}{y_j} - \tfrac{1}{2} \cdot \frac{\alpha_j^2}{y_j}}_{\text{order } o(1)} + \cdots \tag{15}
\end{align}
$$
{{< /math >}}

So:

{{< math >}}
$$\boxed{ (y_j + \alpha_j) \log(y_j + \alpha_j) = y_j \log y_j + \alpha_j \log y_j + \alpha_j + o(1) } \tag{16}$$
{{< /math >}}

## Step 7: Approximate the Full Log Likelihood

Apply the expansion to the sum:

{{< math >}}
$$\sum_j (y_j + \alpha_j) \log(y_j + \alpha_j) \approx \sum_j \left[ y_j \log y_j + \alpha_j \log y_j + \alpha_j \right] + o(1) \tag{17}$$
{{< /math >}}

Now plug into the full log PMF:

{{< math >}}
$$
\begin{align}
\log P(\mathbf{y}) &\approx n \log n - \sum_j y_j \log y_j - (n + \alpha_0) \log(n + \alpha_0) \\
&\quad + \sum_j \left[ y_j \log y_j + \alpha_j \log y_j + \alpha_j \right] \\
&\quad - \sum_j \alpha_j \log \alpha_j + \log \Gamma(\alpha_0) \tag{18}
\end{align}
$$
{{< /math >}}

Cancel {{< math >}} $y_j \log y_j$ {{< /math >}} terms:

{{< math >}}
$$
\begin{align}
\log P(\mathbf{y}) &\approx - (n + \alpha_0) \log(n + \alpha_0) + n \log n + \sum_j \alpha_j \log y_j \\
&\quad + \sum_j \alpha_j - \sum_j \alpha_j \log \alpha_j + \log \Gamma(\alpha_0) + o(1) \tag{19}
\end{align}
$$
{{< /math >}}

## Final Log-PMF Approximation

{{< math >}}
$$\log P(\mathbf{y}) \approx \sum_j \left[ \alpha_j \log y_j \right] + \text{constant} \tag{20}$$
{{< /math >}}

But this form suggests a connection to the **log-pmf of NB**. Let's formalize that connection in the next steps.