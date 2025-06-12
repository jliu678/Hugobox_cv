# 🔢 Model Structure and Fitting

## 1. Observation Model

{{< math >}}$$Y_{gc} \sim \text{Poisson}(s_c \lambda_{gc})$${{< /math >}}

- {{< math >}}$Y_{gc}${{< /math >}}: Observed count for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$c${{< /math >}}
- {{< math >}}$s_c${{< /math >}}: Size factor for cell {{< math >}}$c${{< /math >}}
- {{< math >}}$\lambda_{gc}${{< /math >}}: True (latent) expression level

## 2. Prior on Expression

{{< math >}}$$\lambda_{gc} \sim \text{Gamma}(\alpha_{gc}, \beta_{gc})$${{< /math >}}

This introduces a **Gamma prior** on the Poisson rate, forming a **Poisson-Gamma compound model**, which is equivalent to a **Negative Binomial** model for marginal counts.

# A. Prior Mean Estimation
To estimate the prior mean {{< math >}}$\mu_{gc}${{< /math >}}, the model uses a **Poisson generalized linear model (GLM)** with a **log link** and **LASSO penalty**:

{{< math >}}$$\log E\left(\frac{s_c Y_{gc}}{Y_{g'c}}\right) = \log \mu_{gc} = \gamma_{g0} + \sum_{g' \neq g} \gamma_{gg'} \log\left[\frac{s_c Y_{g'c} + 1}{s_c}\right]$${{< /math >}}

### a. 🔍 Explanation of Terms

- **{{< math >}}$Y_{gc}${{< /math >}}**: Observed expression count for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$c${{< /math >}}.
- **{{< math >}}$s_c${{< /math >}}**: Size factor for cell {{< math >}}$c${{< /math >}}, used for normalization.
- **{{< math >}}$\mu_{gc}${{< /math >}}**: Prior mean expression for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$c${{< /math >}}, estimated from the model.
- **{{< math >}}$\gamma_{g0}${{< /math >}}**: Intercept term for gene {{< math >}}$g${{< /math >}}.
- **{{< math >}}$\gamma_{gg'}${{< /math >}}**: Regression coefficient for the predictor gene {{< math >}}$g'${{< /math >}} when modeling gene {{< math >}}$g${{< /math >}}.
- **{{< math >}}$\log\left[\frac{s_c Y_{g'c} + 1}{s_c}\right]${{< /math >}}**: Log-normalized expression of gene {{< math >}}$g'${{< /math >}} in cell {{< math >}}$c${{< /math >}}, with a +1 pseudocount to avoid log(0).

### b. 📉 Loss Function

The **penalized loss function** minimized during fitting is [derived later](#poisson_lasso_loss_function) and is expressed as:

{{< math >}}$$L(\gamma) = -\sum_c \left[Y_{gc} \log(\mu_{gc}) - \mu_{gc}\right] + \lambda \sum_{g' \neq g} |\gamma_{gg'}|$${{< /math >}}

Where:
- The first term is the **Poisson negative log-likelihood**.
- The second term is the **L1 penalty** on the regression coefficients (excluding the intercept {{< math >}}$\gamma_{g0}${{< /math >}}).

### c. ⚙️ Fitting Procedure

- **Tool Used**: `glmnet` R package (version 2.0–5)
- **Penalty**: LASSO (L1 regularization)
- **Model Selection**: The penalty parameter {{< math >}}$\lambda${{< /math >}} is chosen via **fivefold cross-validation**, selecting the model with the **lowest cross-validation error** 
[(see later for details)](#lasso_cross_validation).

## 🧠 How It Relates to the Empirical Bayes-like Technique

- This model is used to **predict the prior mean** {{< math >}}$\mu_{gc}${{< /math >}} for each gene in each cell.
- The **Poisson GLM with LASSO penalty** selects a sparse set of predictor genes {{< math >}}$g'${{< /math >}} that best explain the expression of gene {{< math >}}$g${{< /math >}}, reflecting biological sparsity (i.e., gene interactions are limited).
- The **LASSO penalty** (via the `glmnet` package) ensures that only a subset of genes with strong predictive power have nonzero coefficients {{< math >}}$\gamma_{gg'}${{< /math >}}.
- The **predicted** {{< math >}}$\mu_{gc}${{< /math >}} from this model is then treated as the **prior mean** in a Bayesian framework, enabling **shrinkage** and **denoising** of gene expression estimates.

# A.1 Derive Loss Function of the Poisson-LASSO Model<a id="poisson_lasso_loss_function"></a>

## 🔢 Step 1: Poisson Probability Mass Function

For a Poisson-distributed random variable {{< math >}} $Y_c \sim \text{Poisson}(\lambda_c)$ {{< /math >}}, the probability of observing {{< math >}} $Y_c = y_c$ {{< /math >}} is:

{{< math >}}
$$P(Y_c = y_c) = \frac{e^{-\lambda_c} \lambda_c^{y_c}}{y_c!}$$
{{< /math >}}

## 🧮 Step 2: Log-Likelihood for All Observations

Assuming we have {{< math >}} $C$ {{< /math >}} independent observations {{< math >}} $\{Y_1, Y_2, \ldots, Y_C\}$ {{< /math >}}, the likelihood is the product of individual probabilities:

{{< math >}}
$$L(\lambda) = \prod_{c=1}^C \frac{e^{-\lambda_c} \lambda_c^{Y_c}}{Y_c!}$$
{{< /math >}}

Taking the log of the likelihood:

{{< math >}}
$$\log L(\lambda) = \sum_{c=1}^C \left[-\lambda_c + Y_c \log(\lambda_c) - \log(Y_c!)\right]$$
{{< /math >}}

## 🔍 Step 3: Negative Log-Likelihood (NLL)

The negative log-likelihood (which we minimize during model fitting) is:

{{< math >}} 
$$\mathcal{L}(\lambda) = -\log L(\lambda) = \sum_{c=1}^C [\lambda_c - Y_c \log(\lambda_c) + \log(Y_c!)]$$ 
{{< /math >}}

In Poisson regression, we model {{< math >}} $\lambda_c = \exp(\eta_c)$ {{< /math >}}, where {{< math >}} $\eta_c$ {{< /math >}} is a linear predictor (e.g., from LASSO regression). So:

{{< math >}} 
$$\lambda_c = \exp\left(\gamma_0 + \sum_j \gamma_j x_{jc}\right)$$ 
{{< /math >}}

Substituting into the NLL:

{{< math >}} 
$$\mathcal{L}(\gamma) = \sum_{c=1}^C [\exp(\eta_c) - Y_c \eta_c + \log(Y_c!)]$$ 
{{< /math >}}

The term **{{< math >}} $\log(Y_c!)$ {{< /math >}} is constant with respect to the parameters {{< math >}} $\gamma$ {{< /math >}}, so it's often omitted during optimization.**

## ➕ Adding LASSO Penalty

To enforce sparsity in the coefficients {{< math >}} $\gamma$ {{< /math >}}, we add an L1 penalty:

{{< math >}} 
$$\mathcal{L}_{\text{penalized}}(\gamma) = \sum_{c=1}^C [\exp(\eta_c) - Y_c \eta_c] + \lambda \sum_j |\gamma_j|$$ 
{{< /math >}}

This is the Poisson LASSO loss function used in your model.

## Key Components Summary

- **Observation Model**: {{< math >}} $Y_{gc} \sim \text{Poisson}(s_c \lambda_{gc})$ {{< /math >}}
- **Size Factor**: {{< math >}} $s_c$ {{< /math >}} normalizes for library size differences
- **Expression Prior**: {{< math >}} $\lambda_{gc} \sim \text{Gamma}(\alpha_{gc}, \beta_{gc})$ {{< /math >}}
- **GLM Framework**: Uses log-link function with LASSO regularization
- **Compound Distribution**: Poisson-Gamma yields Negative Binomial for overdispersed counts


# A.2 LASSO Cross-Validation<a id="lasso_cross_validation"></a>

## 🔧 What is the Penalty Parameter {{< math >}} $\lambda$ {{< /math >}}?

In LASSO regression, the penalty parameter {{< math >}} $\lambda$ {{< /math >}} controls the strength of the L1 regularization:

- **Large {{< math >}} $\lambda$ {{< /math >}}**: More shrinkage → more coefficients set to zero → simpler model.
- **Small {{< math >}} $\lambda$ {{< /math >}}**: Less shrinkage → more coefficients retained → more complex model.

## 🔁 What is Fivefold Cross-Validation?

Cross-validation is a technique to evaluate how well a model generalizes to unseen data. In fivefold cross-validation:

1. The dataset is split into 5 equal parts (folds).
2. The model is trained on 4 folds and tested on the remaining fold.
3. This process is repeated 5 times, each time using a different fold as the test set.
4. The average error across the 5 test sets is computed.

## 📉 How is {{< math >}} $\lambda$ {{< /math >}} Chosen?

1. A range of {{< math >}} $\lambda$ {{< /math >}} values is tested.
2. For each {{< math >}} $\lambda$ {{< /math >}}, fivefold cross-validation is performed.
3. The cross-validation error (e.g., deviance or mean squared error) is computed for each {{< math >}} $\lambda$ {{< /math >}}.
4. The {{< math >}} $\lambda$ {{< /math >}} with the lowest average cross-validation error is selected.

This ensures that the chosen {{< math >}} $\lambda$ {{< /math >}} balances model complexity and predictive accuracy, avoiding both underfitting and overfitting.

## ✅ Why It Matters

In your Poisson LASSO regression model, this process ensures that the selected genes (nonzero coefficients) are reliable predictors of the target gene's expression, based on how well they generalize across different subsets of the data.

## Mathematical Framework

The LASSO objective function can be written as:

{{< math >}}
$$\min_{\beta} \left\{ \frac{1}{2n} \|y - X\beta\|_2^2 + \lambda \|\beta\|_1 \right\}$$
{{< /math >}}

Where:
- {{< math >}} $y$ {{< /math >}} is the response vector
- {{< math >}} $X$ {{< /math >}} is the design matrix
- {{< math >}} $\beta$ {{< /math >}} is the coefficient vector
- {{< math >}} $\lambda$ {{< /math >}} is the penalty parameter
- {{< math >}} $\|\beta\|_1 = \sum_{j=1}^p |\beta_j|$ {{< /math >}} is the L1 norm

The cross-validation error for a given {{< math >}} $\lambda$ {{< /math >}} is:

{{< math >}}
$$CV(\lambda) = \frac{1}{5} \sum_{k=1}^{5} L(y_k, \hat{y}_k(\lambda))$$
{{< /math >}}

Where {{< math >}} $L(y_k, \hat{y}_k(\lambda))$ {{< /math >}} is the loss function evaluated on the {{< math >}} $k$ {{< /math >}}-th fold.

# B. Prior Variance Estimation

## 🧠 Goal: Estimate Prior Variance for Gene {{< math >}} $g$ {{< /math >}}

We already estimated the prior mean {{< math >}} $\mu_{gc}$ {{< /math >}} using Poisson LASSO regression. Now, we want to estimate the prior variance of the latent expression {{< math >}} $\lambda_{gc}$ {{< /math >}}, which is modeled using a Gamma distribution:

{{< math >}} $$\lambda_{gc} \sim \text{Gamma}(\alpha_{gc}, \beta_{gc})$$ {{< /math >}}

The variance of a Gamma distribution is:

{{< math >}} $$\text{Var}(\lambda_{gc}) = \frac{\alpha_{gc}}{\beta_{gc}^2}$$ {{< /math >}}

## 🔧 Noise Models for Prior Variance

To estimate this variance, we assume a constant noise model across cells for each gene {{< math >}} $g$ {{< /math >}}, controlled by a dispersion parameter {{< math >}} $\phi_g$ {{< /math >}}. Three models are considered:

### Constant Coefficient of Variation (CV):
{{< math >}} $\phi_g^{cv}$ {{< /math >}}

Implies constant shape parameter: {{< math >}} $\alpha_{gc} = \alpha_g$ {{< /math >}}

{{< math >}} $$\text{CV} = \frac{\sqrt{\text{Var}(\lambda_{gc})}}{\text{E}[\lambda_{gc}]}$$ {{< /math >}}

### Constant Fano Factor:
{{< math >}} $\phi_g^F$ {{< /math >}}

Implies constant rate parameter: {{< math >}} $\beta_{gc} = \beta_g$ {{< /math >}}

{{< math >}} $$\text{Fano} = \frac{\text{Var}(\lambda_{gc})}{\text{E}[\lambda_{gc}]}$$ {{< /math >}}

### Constant Variance:
{{< math >}} $\phi_g^v$ {{< /math >}}

Directly assumes {{< math >}} $\text{Var}(\lambda_{gc}) = \phi_g$ {{< /math >}}

## 📈 Model Selection via Marginal Likelihood

For each model:

1. Compute the marginal likelihood of the observed data across cells.
2. Choose the model with the highest maximum likelihood.
3. The corresponding {{< math >}} $\phi_g$ {{< /math >}} is set to the maximum likelihood estimate {{< math >}} $\hat{\phi}_g$ {{< /math >}}.

## 📌 Final Step: Compute Prior Variance {{< math >}} $\hat{v}_{gc}$ {{< /math >}}

Once the best model and {{< math >}} $\hat{\phi}_g$ {{< /math >}} are selected, you can derive the prior variance {{< math >}} $\hat{v}_{gc}$ {{< /math >}} for each gene-cell pair using the Gamma distribution formula, depending on the chosen model.

# B.1 Noise Models

If {{< math >}} $\lambda_{gc} \sim \text{Gamma}(\alpha_{gc}, \beta_{gc})$ {{< /math >}}, then:

**Mean:** {{< math >}} $E[\lambda_{gc}] = \frac{\alpha_{gc}}{\beta_{gc}}$ {{< /math >}}

**Variance:** {{< math >}} $\text{Var}(\lambda_{gc}) = \frac{\alpha_{gc}}{\beta_{gc}^2}$ {{< /math >}}

## 1. 📈 Constant Coefficient of Variation (CV)

**Assumption:**
The coefficient of variation (CV) is constant across cells for gene {{< math >}} $g$ {{< /math >}}.
This implies a constant shape parameter: {{< math >}} $\alpha_{gc} = \alpha_g$ {{< /math >}}

**Derivation:**
{{< math >}}
$$\text{CV} = \frac{\sqrt{\text{Var}(\lambda_{gc})}}{E[\lambda_{gc}]} = \frac{\sqrt{\frac{\alpha_g}{\beta_{gc}^2}}}{\frac{\alpha_g}{\beta_{gc}}} = \frac{1}{\sqrt{\alpha_g}}$$
{{< /math >}}

So, under this model:

{{< math >}} $\phi_{g}^{cv} = \text{CV}^2 = \frac{1}{\alpha_g}$ {{< /math >}}

The rate parameter {{< math >}} $\beta_{gc}$ {{< /math >}} varies with {{< math >}} $\mu_{gc}$ {{< /math >}}, since {{< math >}} $\mu_{gc} = \frac{\alpha_g}{\beta_{gc}} \Rightarrow \beta_{gc} = \frac{\alpha_g}{\mu_{gc}}$ {{< /math >}}

## 2. 📊 Constant Fano Factor

**Assumption:**
The Fano factor is constant across cells for gene {{< math >}} $g$ {{< /math >}}.
This implies a constant rate parameter: {{< math >}} $\beta_{gc} = \beta_g$ {{< /math >}}

**Derivation:**
{{< math >}}
$$\text{Fano} = \frac{\text{Var}(\lambda_{gc})}{E[\lambda_{gc}]} = \frac{\frac{\alpha_{gc}}{\beta_g^2}}{\frac{\alpha_{gc}}{\beta_g}} = \frac{1}{\beta_g}$$
{{< /math >}}

So, under this model:

{{< math >}} $\phi_g^F = \text{Fano} = \frac{1}{\beta_g}$ {{< /math >}}

The shape parameter {{< math >}} $\alpha_{gc}$ {{< /math >}} varies with {{< math >}} $\mu_{gc}$ {{< /math >}}, since {{< math >}} $\mu_{gc} = \frac{\alpha_{gc}}{\beta_g} \Rightarrow \alpha_{gc} = \mu_{gc} \beta_g$ {{< /math >}}

## 3. 📐 Constant Variance

**Assumption:**
The variance of {{< math >}} $\lambda_{gc}$ {{< /math >}} is constant across cells for gene {{< math >}} $g$ {{< /math >}}.

{{< math >}} $\text{Var}(\lambda_{gc}) = \phi_g^v$ {{< /math >}}

**Derivation:**
From the Gamma variance formula:

{{< math >}} $\frac{\alpha_{gc}}{\beta_{gc}^2} = \phi_g^v$ {{< /math >}}

This model allows both {{< math >}} $\alpha_{gc}$ {{< /math >}} and {{< math >}} $\beta_{gc}$ {{< /math >}} to vary with {{< math >}} $\mu_{gc}$ {{< /math >}}, but constrains their relationship to maintain constant variance.

## ✅ Summary Table

| Model | Assumption | Fixed Parameter | Varying Parameter | Derived Relationship |
|-------|------------|----------------|-------------------|---------------------|
| **Constant CV** | {{< math >}} $\text{CV} = \frac{1}{\sqrt{\alpha_g}}$ {{< /math >}} | {{< math >}} $\alpha_{gc} = \alpha_g$ {{< /math >}} | {{< math >}} $\beta_{gc} = \frac{\alpha_g}{\mu_{gc}}$ {{< /math >}} | {{< math >}} $\text{Var} = \mu_{gc}^2/\alpha_g$ {{< /math >}} |
| **Constant Fano Factor** | {{< math >}} $\text{Fano} = \frac{1}{\beta_g}$ {{< /math >}} | {{< math >}} $\beta_{gc} = \beta_g$ {{< /math >}} | {{< math >}} $\alpha_{gc} = \mu_{gc} \beta_g$ {{< /math >}} | {{< math >}} $\text{Var} = \mu_{gc}/\beta_g$ {{< /math >}} |
| **Constant Variance** | {{< math >}} $\text{Var} = \phi_g^v$ {{< /math >}} | None | Both vary | {{< math >}} $\frac{\alpha_{gc}}{\beta_{gc}^2} = \phi_g^v$ {{< /math >}} |

## 🔁 How {{< math >}} $\hat{v}_{gc}$ {{< /math >}} is Computed

Depending on the model:

| Model | Assumption | Formula for {{< math >}} $\hat{v}_{gc} = \text{Var}(\lambda_{gc})$ {{< /math >}} |
|-------|------------|--------------------------------------------------------------------------|
| **Constant CV** | {{< math >}} $\alpha_{gc} = \alpha_g$ {{< /math >}} | {{< math >}} $\hat{v}_{gc} = \hat{\phi}_g \cdot \mu_{gc}^2$ {{< /math >}} |
| **Constant Fano Factor** | {{< math >}} $\beta_{gc} = \beta_g$ {{< /math >}} | {{< math >}} $\hat{v}_{gc} = \hat{\phi}_g \cdot \mu_{gc}$ {{< /math >}} |
| **Constant Variance** | {{< math >}} $\text{Var}(\lambda_{gc}) = \phi_g$ {{< /math >}} | {{< math >}} $\hat{v}_{gc} = \hat{\phi}_g$ {{< /math >}} |

So in all cases, **{{< math >}} $\hat{v}_{gc}$ {{< /math >}} is an estimate of the variance** of the latent expression {{< math >}} $\lambda_{gc}$ {{< /math >}}, derived from the Gamma prior and the selected noise model.

# B.2 Marginal Likelihood for Model Selection
## Poisson-Gamma to Negative Binomial: Step-by-Step Derivation

### 1. Observation Model
{{< math >}}$$Y_{gc} \sim \text{Poisson}(s_c \lambda_{gc})$${{< /math >}}

- {{< math >}}$Y_{gc}${{< /math >}}: Observed count for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$c${{< /math >}}
- {{< math >}}$s_c${{< /math >}}: Size factor for cell {{< math >}}$c${{< /math >}}
- {{< math >}}$\lambda_{gc}${{< /math >}}: True (latent) expression level

### 2. Prior on Expression
{{< math >}}$$\lambda_{gc} \sim \text{Gamma}(\alpha_{gc}, \beta_{gc})$${{< /math >}}

This introduces a **Gamma prior** on the Poisson rate, forming a **Poisson-Gamma compound model**, which is equivalent to a **Negative Binomial** model for marginal counts.

To estimate the prior mean {{< math >}}$\mu_{gc}${{< /math >}}, the model uses a **Poisson generalized linear model (GLM)** with a **log link** and **LASSO penalty**.


### 3. Poisson Likelihood

{{< math >}}
$$P(Y_{gc} \mid \lambda_{gc}) = \frac{(s_c \lambda_{gc})^{Y_{gc}} e^{-s_c \lambda_{gc}}}{Y_{gc}!}$$
{{< /math >}}

### 4. Gamma Prior

{{< math >}}
$$P(\lambda_{gc}) = \frac{\beta_{gc}^{\alpha_{gc}}}{\Gamma(\alpha_{gc})} \lambda_{gc}^{\alpha_{gc}-1} e^{-\beta_{gc} \lambda_{gc}}$$
{{< /math >}}

### 5. Marginal Distribution

Multiply and integrate:

{{< math >}}
$$P(Y_{gc}) = \int_0^{\infty} P(Y_{gc} \mid \lambda_{gc}) \cdot P(\lambda_{gc}) \, d\lambda_{gc}$$
{{< /math >}}

{{< math >}}
$$= \frac{\beta_{gc}^{\alpha_{gc}} s_c^{Y_{gc}}}{\Gamma(\alpha_{gc}) Y_{gc}!} \int_0^{\infty} \lambda_{gc}^{Y_{gc} + \alpha_{gc} - 1} e^{-(\beta_{gc} + s_c) \lambda_{gc}} \, d\lambda_{gc}$$
{{< /math >}}

This integral is the definition of the Gamma function:

{{< math >}}
$$\int_0^{\infty} x^{k-1} e^{-\theta x} dx = \frac{\Gamma(k)}{\theta^k}$$
{{< /math >}}

So we get:

{{< math >}}
$$P(Y_{gc}) = \frac{\beta_{gc}^{\alpha_{gc}} s_c^{Y_{gc}}}{\Gamma(\alpha_{gc}) Y_{gc}!} \cdot \frac{\Gamma(Y_{gc} + \alpha_{gc})}{(\beta_{gc} + s_c)^{Y_{gc} + \alpha_{gc}}}$$
{{< /math >}}

### 6. Final Simplified Form

Rewriting:

{{< math >}}
$$P(Y_{gc}) = \frac{\Gamma(Y_{gc} + \alpha_{gc})}{\Gamma(\alpha_{gc}) Y_{gc}!} \left(\frac{s_c}{\beta_{gc} + s_c}\right)^{Y_{gc}} \left(\frac{\beta_{gc}}{\beta_{gc} + s_c}\right)^{\alpha_{gc}}$$
{{< /math >}}

This is the **Negative Binomial distribution** with parameters:
- Shape parameter: {{< math >}}$\alpha_{gc}${{< /math >}}
- Success probability: {{< math >}}$p = \frac{\beta_{gc}}{\beta_{gc} + s_c}${{< /math >}}
- Mean: {{< math >}}$\mathbb{E}[Y_{gc}] = \frac{s_c \alpha_{gc}}{\beta_{gc}}${{< /math >}}
- Variance: {{< math >}}$\text{Var}[Y_{gc}] = \frac{s_c \alpha_{gc}}{\beta_{gc}} + \frac{s_c^2 \alpha_{gc}}{\beta_{gc}^2}${{< /math >}} (overdispersed relative to Poisson)

## 🧮 Result of the Integration

The symbolic integration yields:

{{< math >}}$$P(Y_{gc}) = \frac{\Gamma(Y_{gc} + \alpha_{gc})}{\Gamma(\alpha_{gc}) Y_{gc}!} \left(\frac{s_c}{\beta_{gc} + s_c}\right)^{Y_{gc}} \left(\frac{\beta_{gc}}{\beta_{gc} + s_c}\right)^{\alpha_{gc}}$${{< /math >}}

This is the probability mass function of a **Negative Binomial distribution** with:

- **Number of failures**: {{< math >}}$r = \alpha_{gc}${{< /math >}}
- **Success probability**: {{< math >}}$p = \frac{\beta_{gc}}{\beta_{gc} + s_c}${{< /math >}}

## ✅ Summary

This proves that:

{{< math >}}$$Y_{gc} \sim \text{Negative Binomial}\left(\alpha_{gc}, p = \frac{\beta_{gc}}{\beta_{gc} + s_c}\right)$${{< /math >}}

So, the **Poisson-Gamma compound model** naturally leads to a **Negative Binomial marginal distribution**, which is why it's widely used in modeling overdispersed count data like gene expression.

# B.3 Fitting to decide noise model and variance expectation

When fitting the maximum marginal likelihood to estimate the dispersion parameter {{< math >}} $\phi_g$ {{< /math >}} for a gene {{< math >}} $g$ {{< /math >}}, here's what is known and what is being estimated:

## ✅ Known Quantities

These are available from earlier steps in the model:

- **Observed counts**: {{< math >}} $Y_{gc}$ {{< /math >}} for each cell {{< math >}} $c$ {{< /math >}}
- **Size factors**: {{< math >}} $s_c$ {{< /math >}} for each cell
- **Prior means**: {{< math >}} $\mu_{gc}$ {{< /math >}}, estimated from Poisson LASSO regression
- **Noise model**: One of the three (constant CV, constant Fano factor, constant variance)

## 🧮 size factor
The size factor {{< math >}} $s_c$ {{< /math >}} is a cell-specific scaling factor used to normalize raw gene expression counts. It accounts for differences in:

- **Sequencing depth** (total number of reads per cell)
- **Capture efficiency** (how well RNA was captured from each cell)

It ensures that expression levels are comparable across cells.Without normalization, cells with more reads would appear to have higher expression for all genes, simply due to technical artifacts—not biology.

There are several methods to estimate {{< math >}} $s_c$ {{< /math >}}, but the general idea is:

### 1. Total Count Normalization (used in SAVER)

{{< math >}}
$$s_c = \frac{\text{total counts in cell } c}{\text{median total counts across all cells}}$$
{{< /math >}}

### 2. Median Ratio Method (used in DESeq2)

1. Compute the geometric mean of each gene across all cells
2. For each cell, compute the ratio of each gene's count to the gene's geometric mean
3. The median of these ratios is the size factor {{< math >}} $s_c$ {{< /math >}}

### 3. Scran's Deconvolution Method

- Pools cells to estimate size factors more robustly in sparse data

## 🔍 Parameter to Estimate

The dispersion parameter {{< math >}} $\phi_g$ {{< /math >}} for gene {{< math >}} $g$ {{< /math >}}

This is the only free parameter in the marginal likelihood optimization.

## 🧮 How the Fitting Works

1. For a given value of {{< math >}} $\phi_g$ {{< /math >}}, use the selected noise model to derive the Gamma parameters {{< math >}} $\alpha_{gc}$ {{< /math >}}, {{< math >}} $\beta_{gc}$ {{< /math >}} for each cell.

2. Use these to compute the marginal likelihood of the observed counts {{< math >}} $Y_{gc}$ {{< /math >}} under the Negative Binomial distribution.

3. Sum the log-likelihoods across all cells:

{{< math >}}
$$\log L(\phi_g) = \sum_c \log P(Y_{gc} \mid \mu_{gc}, \phi_g)$$
{{< /math >}}

4. Optimize {{< math >}} $\phi_g$ {{< /math >}} to maximize this log-likelihood.

## ✅ Summary

| Category | Items |
|----------|-------|
| **Known** | {{< math >}} $Y_{gc}$ {{< /math >}}, {{< math >}} $s_c$ {{< /math >}}, {{< math >}} $\mu_{gc}$ {{< /math >}}, noise model |
| **To Estimate** | {{< math >}} $\phi_g$ {{< /math >}} (dispersion parameter) |
| **Derived** | {{< math >}} $\alpha_{gc}$ {{< /math >}}, {{< math >}} $\beta_{gc}$ {{< /math >}} (from {{< math >}} $\mu_{gc}$ {{< /math >}} and {{< math >}} $\phi_g$ {{< /math >}}) |
| **Objective** | Maximize marginal likelihood {{< math >}} $\log L(\phi_g)$ {{< /math >}} |

# C. Full Derivation of the Bayesian Posterior of $λ_{gc}$

## 🧾 Model Setup

**Likelihood:**

{{< math >}} $$ Y_{gc} \sim \text{Poisson}(\lambda_{gc} \cdot s_c) $$ {{< /math >}}

{{< math >}} $$ p(Y_{gc} \mid \lambda_{gc}, s_c) = \frac{(s_c \lambda_{gc})^{Y_{gc}}}{Y_{gc}!} e^{-s_c \lambda_{gc}} $$ {{< /math >}}

**Prior:**

{{< math >}} $$ \lambda_{gc} \sim \text{Gamma}(\hat{\alpha}_{gc}, \hat{\beta}_{gc}) $$ {{< /math >}}

where the Gamma is in shape–rate form, i.e.,

{{< math >}} $$ p(\lambda_{gc}) = \frac{\hat{\beta}_{gc}^{\hat{\alpha}_{gc}}}{\Gamma(\hat{\alpha}_{gc})} \lambda_{gc}^{\hat{\alpha}_{gc}-1} e^{-\hat{\beta}_{gc}\lambda_{gc}} $$ {{< /math >}}

We want to compute the posterior:

{{< math >}} $$ p(\lambda_{gc} \mid Y_{gc}, s_c, \hat{\alpha}_{gc}, \hat{\beta}_{gc}) $$ {{< /math >}}

## Simplified writting
If:

**Prior:** 
{{< math >}} $\lambda \sim \text{Gamma}(\alpha, \beta)$ {{< /math >}}

**Likelihood:** 
{{< math >}} $Y \sim \text{Poisson}(s\lambda)$ {{< /math >}}

Then we will prove that the posterior is:

{{< math >}} $$\lambda | Y \sim \text{Gamma}(Y + \alpha, s + \beta)$$ {{< /math >}}

### 🔢 1. Define the Distributions

**Likelihood:**
{{< math >}} $$p(Y|\lambda) = \frac{(s\lambda)^Y e^{-s\lambda}}{Y!}$$ {{< /math >}}

**Prior (Gamma):**
{{< math >}} $$p(\lambda) = \frac{\beta^\alpha}{\Gamma(\alpha)} \lambda^{\alpha-1} e^{-\beta\lambda}$$ {{< /math >}}

**Goal:** Derive {{< math >}} $p(\lambda|Y)$ {{< /math >}}

We use Bayes' Rule:

{{< math >}} $$p(\lambda|Y) = \frac{p(Y|\lambda) \cdot p(\lambda)}{p(Y)}$$ {{< /math >}}

### Step 1: Multiply Numerator

{{< math >}} $$p(Y|\lambda) \cdot p(\lambda) = \frac{(s\lambda)^Y e^{-s\lambda}}{Y!} \cdot \frac{\beta^\alpha}{\Gamma(\alpha)} \lambda^{\alpha-1} e^{-\beta\lambda}$$ {{< /math >}}

Group terms:

{{< math >}} $$= \frac{s^Y \beta^\alpha}{Y! \Gamma(\alpha)} \lambda^{Y+\alpha-1} e^{-(s+\beta)\lambda}$$ {{< /math >}}

### Step 2: Marginal Likelihood (Evidence)

We want to compute the marginal probability of the observed data {{< math >}} $Y$ {{< /math >}}, by integrating out {{< math >}} $\lambda$ {{< /math >}}:

{{< math >}} $$p(Y) = \int_0^\infty p(Y|\lambda)p(\lambda) \, d\lambda$$ {{< /math >}}

Substitute the pdfs again:

{{< math >}} $$p(Y) = \int_0^\infty \frac{(s\lambda)^Y e^{-s\lambda}}{Y!} \cdot \frac{\beta^\alpha}{\Gamma(\alpha)} \lambda^{\alpha-1} e^{-\beta\lambda} d\lambda$$ {{< /math >}}

Simplify constants:

{{< math >}} $$p(Y) = \frac{s^Y \beta^\alpha}{Y! \Gamma(\alpha)} \int_0^\infty \lambda^{Y+\alpha-1} e^{-(s+\beta)\lambda} d\lambda$$ {{< /math >}}

Recognize the integral as the definition of the Gamma function:

{{< math >}} $$\int_0^\infty x^{k-1} e^{-\theta x} dx = \frac{\Gamma(k)}{\theta^k}$$ {{< /math >}}

So,

{{< math >}} $$\int_0^\infty \lambda^{Y+\alpha-1} e^{-(s+\beta)\lambda} d\lambda = \frac{\Gamma(Y+\alpha)}{(s+\beta)^{Y+\alpha}}$$ {{< /math >}}

Thus, the marginal likelihood is:

{{< math >}} $$p(Y) = \frac{s^Y \beta^\alpha}{Y! \Gamma(\alpha)} \cdot \frac{\Gamma(Y+\alpha)}{(s+\beta)^{Y+\alpha}}$$ {{< /math >}}


### Step 3: Combine Numerator and Denominator

{{< math >}} $$p(\lambda|Y) = \frac{\frac{s^Y \beta^\alpha}{Y! \Gamma(\alpha)} \lambda^{Y+\alpha-1} e^{-(s+\beta)\lambda}}{\frac{s^Y \beta^\alpha}{Y! \Gamma(\alpha)} \cdot \frac{\Gamma(Y+\alpha)}{(s+\beta)^{Y+\alpha}}}$$ {{< /math >}}

Cancel out constants:

{{< math >}} $$p(\lambda|Y) = \frac{(s+\beta)^{Y+\alpha}}{\Gamma(Y+\alpha)} \cdot \lambda^{Y+\alpha-1} e^{-(s+\beta)\lambda}$$ {{< /math >}}

## Final Result: Posterior is Gamma

{{< math >}} $$p(\lambda|Y) = \text{Gamma}(\lambda | Y+\alpha, \, s+\beta)$$ {{< /math >}}

Where the general Gamma pdf is:

{{< math >}} $$\text{Gamma}(\lambda|a,b) = \frac{b^a}{\Gamma(a)} \lambda^{a-1} e^{-b\lambda}$$ {{< /math >}}

# Conclusion
The above math derivation has achieved the goal of SAVER--
- derive the posterior gamma distribution for $λ_{gc}$ given the
observed counts $Y_{gc}$ 
- use the posterior mean as the normalized SAVER $ \hat{\lambda}_{gc}$
- use the variance in the posterior distribution can be thought of as a measure of uncertainty in the SAVER estimate.

