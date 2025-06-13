---
title: 🧮 Math Derivation for scRNAseq Imputation (SAVER)
summary: Mathematical models used in SAVER, a tool of scRNAseq imputation, could well depict the structure inherent to the scRNAseq datasets as suggested by its superior performance. Here we derive SAVER's Poisson–gamma mixture model (also known as negative binomial model) and its Bayesian framework that leverage conjugate priors to estimate the posterior distribution of gene expression levels.
date: 2025-04-24
authors:
  - admin
tags:
  - scRNAseq Poisson–Gamma Mixture Model
  - Imputation
  - Conjugate Priors
  - Bayesian
  - Marginal Distribution
image:
  caption: 'Image credit: [**Round Icons on Unsplash**](https://unsplash.com)'
---
## Introduction

Single-cell RNASeq could have a large amount of zero values, representing either missing data or no expression. Imputation approaches to deal with this issue have the risk of generating false positive or irreproducible results. Model-based imputation generate fewer false-positives compared with data smoothing based methods (MAGIC and knn-smooth), but this varied greatly depending on how well the model described the datasets.

[SAVER](https://www.nature.com/articles/s41592-018-0033-z#Sec2) was the least likely to generate false or irreproducible results in [a benchmark of common imputation methods](https://f1000research.com/articles/7-1740/v1). It suggests the mathematical model used in SAVER could well depict the structure inherent to the scRNAseq datasets. Here we will derive the Poisson–gamma mixture model (also known as negative binomial model) and its Bayesian framework used in SAVER that leverage conjugate priors to estimate the posterior distribution of gene expression levels.

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