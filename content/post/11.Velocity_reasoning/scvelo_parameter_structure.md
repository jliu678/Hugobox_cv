# scVelo parameter structure

## Gene-specific parameters (same for all cells)

* `.var['fit_alpha']` - Transcription rate for each gene ({{< math >}} $\alpha_g$ {{< /math >}})
* `.var['fit_beta']` - Splicing rate for each gene ({{< math >}} $\beta_g$ {{< /math >}})  
* `.var['fit_gamma']` - Degradation rate for each gene ({{< math >}} $\gamma_g$ {{< /math >}})
* `.var['fit_t_']` - Switching time for each gene ({{< math >}} $t_{switch,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{genes},)$ {{< /math >}} - one value per gene

## Cell-and-gene-specific data

* `.layers['fit_t']` - Latent time varies by both cell AND gene ({{< math >}} $t_{c,g}$ {{< /math >}})
* `.layers['velocity']` - Velocity varies by both cell AND gene ({{< math >}} $v_{c,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{cells} \times n_{genes})$ {{< /math >}} - different value for each cell-gene combination

## Mathematical representation

For a given gene {{< math >}} $g$ {{< /math >}} and cell {{< math >}} $c$ {{< /math >}}:

{{< math >}}
$$
\begin{align}
\alpha_g &= \text{transcription rate (gene-specific)} \\
\beta_g &= \text{splicing rate (gene-specific)} \\
\gamma_g &= \text{degradation rate (gene-specific)} \\
t_{c,g} &= \text{latent time (cell- and gene-specific)} \\
v_{c,g} &= \text{velocity (cell- and gene-specific)}
\end{align}
$$
{{< /math >}}

## Why this makes sense biologically

* Each gene has **intrinsic kinetic properties** - the rates {{< math >}} $\alpha_g$, $\beta_g$, $\gamma_g$ {{< /math >}} describe how fast the gene transcribes, splices, and degrades
* But different cells are at **different stages** of the process for each gene - the latent time {{< math >}} $t_{c,g}$ {{< /math >}} varies
* The **rates are gene properties**, the **timing/stage varies by cell**

This means:
- If you want to compare kinetic rates between genes → look at `.var`
- If you want to see where different cells are in the process for a specific gene → look at `.layers['fit_t']`

## Data access patterns

```python
# Same gene has same intrinsic kinetic rates across all cells
gene_alpha = adata.var['fit_alpha']['Ins1']  # α_Ins1 (same for all cells)
gene_beta = adata.var['fit_beta']['Ins1']    # β_Ins1 (same for all cells)  
gene_gamma = adata.var['fit_gamma']['Ins1']  # γ_Ins1 (same for all cells)

# But latent time varies per cell for the same gene
ins1_latent_times = adata.layers['fit_t'][:, gene_idx]  # t_c,Ins1 (different per cell)
# Cell 1 might be at t=0.2, Cell 2 at t=0.8 for the same gene
```