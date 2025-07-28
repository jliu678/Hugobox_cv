# Pearson Residual Analysis for Rare Cell Populations in scRNA-seq

## Background Context

We're analyzing scRNA-seq UMI counts, which are often modeled using Poisson or Negative Binomial distributions.

**Given:**
- A dataset of {{< math >}} $n$ {{< /math >}} total cells
- A rare cell population {{< math >}} $P$ {{< /math >}} of size {{< math >}} $m$ {{< /math >}} (with {{< math >}} $m \ll n$ {{< /math >}})
- These {{< math >}} $m$ {{< /math >}} rare cells express a marker gene at {{< math >}} $\text{Poisson}(\lambda)$ {{< /math >}}
- The remaining {{< math >}} $n - m$ {{< /math >}} cells have zero expression for this gene
- Sequencing depth is equal across cells

## Step 1: Average Expression Across All Cells

Since only {{< math >}} $m$ {{< /math >}} of the {{< math >}} $n$ {{< /math >}} cells express the gene (with mean {{< math >}} $\lambda$ {{< /math >}}), and the rest are zero:

{{< math >}}
$$\text{Expected mean expression across all cells} = \frac{1}{n}(m \cdot \lambda + (n-m) \cdot 0) = \frac{m\lambda}{n} \quad (1)$$
{{< /math >}}

## Step 2: Pearson Residual for a Cell in Population P

The Pearson residual for a count {{< math >}} $x$ {{< /math >}} is defined as:

{{< math >}}
$$\text{residual} = \frac{\text{observed} - \text{expected}}{\sqrt{\text{variance}}} \quad (2)$$
{{< /math >}}

In this case:

- **Observed expression** in population P cells = {{< math >}} $\lambda$ {{< /math >}} (on average, from {{< math >}} $\text{Poisson}(\lambda)$ {{< /math >}})
- **Expected** = average across all cells = {{< math >}} $\frac{m\lambda}{n}$ {{< /math >}}
- **Variance** under Poisson = same as the mean = {{< math >}} $\frac{m\lambda}{n}$ {{< /math >}}

Therefore, the residual is:

{{< math >}}
$$\text{Residual} = \frac{\lambda - \frac{m\lambda}{n}}{\sqrt{\frac{m\lambda}{n}}} \quad (3)$$
{{< /math >}}

Simplifying the numerator:

{{< math >}}
$$\lambda - \frac{m\lambda}{n} = \lambda\left(1 - \frac{m}{n}\right) = \lambda \cdot \frac{n-m}{n} \quad (4)$$
{{< /math >}}

Substituting back:

{{< math >}}
$$\text{Residual} = \frac{\lambda \cdot \frac{n-m}{n}}{\sqrt{\frac{m\lambda}{n}}} = \lambda \cdot \frac{n-m}{n} \cdot \frac{1}{\sqrt{\frac{m\lambda}{n}}} \quad (5)$$
{{< /math >}}

Further simplification:

{{< math >}}
$$\text{Residual} = \frac{n-m}{n} \cdot \lambda \cdot \sqrt{\frac{n}{m\lambda}} = \frac{n-m}{n} \cdot \sqrt{\frac{\lambda n}{m}} \quad (6)$$
{{< /math >}}

If {{< math >}} $m \ll n$ {{< /math >}}, then {{< math >}} $\frac{n-m}{n} \approx 1$ {{< /math >}}, so:

{{< math >}}
$$\text{Residual} \approx \sqrt{\frac{\lambda n}{m}} \quad (7)$$
{{< /math >}}

## Conclusion

The Pearson residual for rare population cells is approximately:

{{< math >}}
$$\boxed{\sqrt{\frac{\lambda n}{m}}} \quad (8)$$
{{< /math >}}

### Key Insights:

- If {{< math >}} $\lambda$ {{< /math >}} is **high** (strong marker gene expression)
- If {{< math >}} $n$ {{< /math >}} is **large** (many total cells)  
- And {{< math >}} $m$ {{< /math >}} is **small** (rare population)

➡️ **The residual becomes large**, making the marker gene stand out clearly in residual-based dimension reduction methods (e.g., Pearson PCA or GLM-PCA).

This mathematical framework explains why Pearson residuals are particularly effective for identifying rare cell populations in single-cell RNA sequencing data, as the residual magnitude scales with both the expression strength of the marker gene and the rarity of the population expressing it.