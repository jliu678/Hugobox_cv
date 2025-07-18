---
title: "SAVER Model Structure"
date: 2025-07-10
draft: True
---

##
more in another md file

## ✅ Goal:

{{< math >}} 
$$
 \hat{\lambda}_{gc} = \frac{Y_{gc} + \hat{\alpha}_{gc}}{s_c + \hat{\beta}_{gc}} = \frac{s_c}{s_c + \hat{\beta}_{gc}} \cdot \frac{Y_{gc}}{s_c} + \frac{\hat{\beta}_{gc}}{s_c + \hat{\beta}_{gc}} \cdot \hat{\mu}_{gc} 
$$
{{< /math >}}

**Gamma distribution gives:**

{{< math >}} 
$$ \hat{\alpha}_{gc} = \hat{\beta}_{gc} \hat{\mu}_{gc} \quad \text{(1)} $$
{{< /math >}}

## 🔁 Step-by-Step Proof:

We start with:

{{< math >}} 
$$ \hat{\lambda}_{gc} = \frac{Y_{gc} + \hat{\alpha}_{gc}}{s_c + \hat{\beta}_{gc}} \quad \text{(2)} $$
{{< /math >}}

Substitute {{< math >}} $\hat{\alpha}_{gc} = \hat{\beta}_{gc} \hat{\mu}_{gc}$ {{< /math >}} from equation (1):

{{< math >}} 
$$ = \frac{Y_{gc} + \hat{\beta}_{gc} \hat{\mu}_{gc}}{s_c + \hat{\beta}_{gc}} \quad \text{(3)} $$
{{< /math >}}

Split the numerator:

{{< math >}} 
$$ = \frac{Y_{gc}}{s_c + \hat{\beta}_{gc}} + \frac{\hat{\beta}_{gc} \hat{\mu}_{gc}}{s_c + \hat{\beta}_{gc}} \quad \text{(4)} $$
{{< /math >}}

Now multiply and divide the first term by {{< math >}} $s_c$ {{< /math >}}:

{{< math >}} 
$$ = \left(\frac{s_c}{s_c + \hat{\beta}_{gc}} \cdot \frac{Y_{gc}}{s_c}\right) + \left(\frac{\hat{\beta}_{gc}}{s_c + \hat{\beta}_{gc}} \cdot \hat{\mu}_{gc}\right) \quad \text{(5)} $$
{{< /math >}}

Which gives:

{{< math >}} 
$$ \hat{\lambda}_{gc} = \frac{s_c}{s_c + \hat{\beta}_{gc}} \cdot \frac{Y_{gc}}{s_c} + \frac{\hat{\beta}_{gc}}{s_c + \hat{\beta}_{gc}} \cdot \hat{\mu}_{gc} \quad \text{(6)} $$
{{< /math >}}

✅ **Proved.**

## 🔍 Summary:

This is a weighted average of:

- {{< math >}} $\frac{Y_{gc}}{s_c}$ {{< /math >}}: observed normalized expression
- {{< math >}} $\hat{\mu}_{gc}$ {{< /math >}}: predicted expression from prior

With weights proportional to data confidence ({{< math >}} $s_c$ {{< /math >}}) and prior confidence ({{< math >}} $\hat{\beta}_{gc}$ {{< /math >}}).