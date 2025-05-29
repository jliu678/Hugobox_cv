---
title: Extending the Steady-State Model with Offset
math: true
date: '2025-05-28'
---

# Extending the Steady-State Model with Offset

## 🧩 Original Model (No Offset)

Previously, we assumed:

{{< math >}}
$$u_i \approx \gamma_0 \cdot s_i$$
{{< /math >}}

But this forces the regression line to go through the origin, which isn't always biologically realistic.

## ✅ Generalized Model (With Offset)

Now we model:

{{< math >}}
$$u_i \approx \gamma_0 \cdot s_i + o$$
{{< /math >}}

Where:

- {{< math >}}$\gamma_0${{< /math >}} is still the slope (steady-state ratio)
- {{< math >}}$o${{< /math >}} is a constant offset (intercept), modeling basal transcription

## 🔍 Least Squares Solution with Offset

To solve for both {{< math >}}$\gamma_0${{< /math >}} and {{< math >}}$o${{< /math >}}, we use ordinary least squares (OLS) for linear regression.

Given data {{< math >}}$u=(u_1,\ldots,u_n)${{< /math >}}, {{< math >}}$s=(s_1,\ldots,s_n)${{< /math >}}, we estimate:

{{< math >}}
$$\hat{u} = \gamma_0 s + o$$
{{< /math >}}

OLS gives:


### 🎯 Slope (Steady-State Ratio)

{{< math >}}
$$\gamma_0 = \frac{\text{Cov}(u,s)}{\text{Var}(s)}$$
{{< /math >}}

Where:

{{< math >}}
$$\text{Cov}(u,s) = \frac{1}{n}\sum_{i=1}^n (u_i-\bar{u})(s_i-\bar{s})$$
{{< /math >}}

{{< math >}}
$$\text{Var}(s) = \frac{1}{n}\sum_{i=1}^n (s_i-\bar{s})^2$$
{{< /math >}}

### 🧷 Offset

{{< math >}}
$$o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}

This centers the regression line at the mean of the data points.
