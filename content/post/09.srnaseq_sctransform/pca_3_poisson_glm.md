# Poisson GLM for Count Data and PCA Equivalence

## Poisson Deviance for Count Data

Here looks into the deviance for a Poisson Generalized Linear Model (GLM), particularly in the context of dimensionality reduction for count data (like UMI counts). It is the objective function that would be minimized instead of the squared error in PCA when modeling Poisson-distributed data.

### Understanding the Deviance Formula

The expression in Equation (1) is the Poisson deviance (or G-statistic):

{{< math >}} 
$$D=\sum_{i,j} [y_{ij} \log(\frac{y_{ij}}{\mu_{ij}})-(y_{ij}-\mu_{ij})] \quad (1)$$
{{< /math >}}

Where:
- {{< math >}} $y_{ij}$ {{< /math >}}: Represents the observed count data (e.g., UMI counts for sample i and feature j)
- {{< math >}} $\mu_{ij}$ {{< /math >}}: Represents the expected mean count for sample i and feature j under the model. This is the parameter of the Poisson distribution
- {{< math >}} $\log(\frac{y_{ij}}{\mu_{ij}})$ {{< /math >}}: Compares the observed count to the expected count
- {{< math >}} $(y_{ij}-\mu_{ij})$ {{< /math >}}: The difference between observed and expected counts

The deviance is a measure of the goodness-of-fit of a statistical model. For GLMs, minimizing the deviance is equivalent to maximizing the likelihood of the model.

### The Linear Predictor in Brief

The second part defines the linear predictor ({{< math >}} $\eta_{ij}$ {{< /math >}}) for the Poisson GLM:

{{< math >}} 
$$\eta_{ij} = \log\mu_{ij} = \log s_i + u_i' v_j \quad (2)$$
{{< /math >}}

Where:
- {{< math >}} $\eta_{ij}$ {{< /math >}}: The linear predictor. For the Poisson distribution, {{< math >}} $\eta_{ij} = \log\mu_{ij}$ {{< /math >}}. This ensures that {{< math >}} $\mu_{ij}$ {{< /math >}} (which is {{< math >}} $\exp(\eta_{ij})$ {{< /math >}}) is always positive, as required for a Poisson rate parameter.

- {{< math >}} $\log s_i$ {{< /math >}}: This term often represents an offset or a library size normalization factor. In UMI count data, {{< math >}} $s_i$ {{< /math >}} typically denotes the total number of UMIs (or reads) in sample i. Including {{< math >}} $\log s_i$ {{< /math >}} means that the model accounts for variations in sequencing depth across samples. It's often written as {{< math >}} $\log\mu_{ij} = \log s_i + \text{linear term}$ {{< /math >}}, which implies {{< math >}} $\mu_{ij} = s_i \times \exp(\text{linear term})$ {{< /math >}}. This scales the expected counts proportionally to the library size.

- {{< math >}} $u_i' v_j = v_{j1} + \sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}}: This is the low-rank approximation part, similar to what you'd see in PCA.

The expansion can be written as:

{{< math >}} 
$$u_i' v_j = v_{j1} + \sum_{l=2}^L u_{il} v_{jl} \quad (3)$$
{{< /math >}}

Where:
- {{< math >}} $u_i'$ {{< /math >}}: Represents the i-th row of the latent factor matrix (U), describing the position of sample i in the low-dimensional space
- {{< math >}} $v_j$ {{< /math >}}: Represents the j-th column of the loading matrix ({{< math >}} $V^T$ {{< /math >}} or {{< math >}} $W_k$ {{< /math >}}), describing the contribution of the j-th feature to the latent dimensions

Note that the first component {{< math >}} $v_{j1}$ {{< /math >}} might be an intercept term or a common factor, and the subsequent terms {{< math >}} $\sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}} are the actual low-rank decomposition for L latent dimensions.

### Relationship to PCA and Count Data

While PCA minimizes Euclidean distance (squared error), this Poisson GLM framework minimizes Poisson deviance.

**PCA (for Gaussian data):**
{{< math >}} 
$$D_{\text{Euclidean}} = \sum_{i,j} (y_{ij} - \mu_{ij})^2 \quad (4)$$
{{< /math >}}

where {{< math >}} $\mu_{ij}$ {{< /math >}} is the low-rank approximation. This assumes constant variance.

**Poisson GLM (for count data):**
{{< math >}} 
$$D_{\text{Poisson}} = \sum_{i,j} [y_{ij} \log(\frac{y_{ij}}{\mu_{ij}})-(y_{ij}-\mu_{ij})] \quad (5)$$
{{< /math >}}

where {{< math >}} $\mu_{ij} = s_i \exp(u_i' v_j)$ {{< /math >}}. This accounts for the mean-variance relationship of Poisson data.

The goal is still dimensionality reduction: you are modeling a complex high-dimensional count matrix ({{< math >}} $y_{ij}$ {{< /math >}}) with a lower-dimensional representation ({{< math >}} $u_i' v_j$ {{< /math >}}). However, instead of assuming Gaussian noise, it assumes Poisson noise, which is more appropriate for count data.

This type of model forms the basis of many modern single-cell RNA sequencing analysis methods (e.g., using GLM-PCA or similar factor analysis methods) because it respects the underlying statistical properties of count data, providing a more robust and biologically meaningful dimensionality reduction than standard PCA applied to transformed counts.


## Derivation of Poisson Deviance

To derive the Poisson deviance {{< math >}} $D=\sum_{i,j} [y_{ij} \log(\frac{y_{ij}}{\mu_{ij}})-(y_{ij}-\mu_{ij})]$ {{< /math >}}, we start from the log-likelihood function of the Poisson distribution. This deviance is a measure of fit in generalized linear models (GLMs) and is fundamentally linked to the likelihood.

### The Poisson Probability Mass Function (PMF)

For a single observation {{< math >}} $y$ {{< /math >}} from a Poisson distribution with mean {{< math >}} $\mu$ {{< /math >}}, the PMF is:

{{< math >}} 
$$P(Y=y|\mu) = \frac{e^{-\mu} \mu^y}{y!} \quad (1)$$
{{< /math >}}

### Log-Likelihood for a Single Observation

Taking the natural logarithm of the PMF gives the log-likelihood for a single observation:

{{< math >}} 
$$\log L(y|\mu) = \log\left(\frac{e^{-\mu} \mu^y}{y!}\right) \quad (2)$$
{{< /math >}}

Expanding this expression:

{{< math >}} 
$$\log L(y|\mu) = -\mu + y\log(\mu) - \log(y!) \quad (3)$$
{{< /math >}}

### Log-Likelihood for the Entire Dataset

For a dataset of observed counts {{< math >}} $y_{ij}$ {{< /math >}} and their corresponding model-predicted means {{< math >}} $\mu_{ij}$ {{< /math >}}, assuming independence, the total log-likelihood is the sum of individual log-likelihoods:

{{< math >}} 
$$\log L(Y|M) = \sum_{i,j} [-\mu_{ij} + y_{ij} \log(\mu_{ij}) - \log(y_{ij}!)] \quad (4)$$
{{< /math >}}

where {{< math >}} $M$ {{< /math >}} denotes the model parameters (which determine {{< math >}} $\mu_{ij}$ {{< /math >}}).

### Deviance as a Measure of Fit

In GLMs, deviance is defined as twice the difference between the log-likelihood of a saturated model and the log-likelihood of the current model.

#### Saturated Model ({{< math >}} $L_S$ {{< /math >}})

This is a hypothetical model that perfectly fits the data. For each observation {{< math >}} $y_{ij}$ {{< /math >}}, its predicted mean is exactly {{< math >}} $y_{ij}$ {{< /math >}}. So, {{< math >}} $\mu_{ij}^{\text{saturated}} = y_{ij}$ {{< /math >}}.

The log-likelihood of the saturated model is:

{{< math >}} 
$$\log L_S = \sum_{i,j} [-y_{ij} + y_{ij} \log(y_{ij}) - \log(y_{ij}!)] \quad (5)$$
{{< /math >}}

#### Current Model ({{< math >}} $L_M$ {{< /math >}})

This is the model we are evaluating, with predicted means {{< math >}} $\mu_{ij}$ {{< /math >}}.

The log-likelihood of the current model is:

{{< math >}} 
$$\log L_M = \sum_{i,j} [-\mu_{ij} + y_{ij} \log(\mu_{ij}) - \log(y_{ij}!)] \quad (6)$$
{{< /math >}}

#### Computing the Deviance

The deviance {{< math >}} $D$ {{< /math >}} is then:

{{< math >}} 
$$D = 2(\log L_S - \log L_M) \quad (7)$$
{{< /math >}}

Substituting the expressions for {{< math >}} $\log L_S$ {{< /math >}} and {{< math >}} $\log L_M$ {{< /math >}}:

{{< math >}} 
$$D = 2\sum_{i,j} [(-y_{ij} + y_{ij} \log(y_{ij}) - \log(y_{ij}!)) - (-\mu_{ij} + y_{ij} \log(\mu_{ij}) - \log(y_{ij}!))] \quad (8)$$
{{< /math >}}

Notice that the {{< math >}} $-\log(y_{ij}!)$ {{< /math >}} term cancels out:

{{< math >}} 
$$D = 2\sum_{i,j} [-y_{ij} + y_{ij} \log(y_{ij}) + \mu_{ij} - y_{ij} \log(\mu_{ij})] \quad (9)$$
{{< /math >}}

Now, we can rearrange the terms:

{{< math >}} 
$$D = 2\sum_{i,j} [y_{ij} (\log(y_{ij}) - \log(\mu_{ij})) - y_{ij} + \mu_{ij}] \quad (10)$$
{{< /math >}}

Using the logarithm property {{< math >}} $\log(a) - \log(b) = \log(\frac{a}{b})$ {{< /math >}}:

{{< math >}} 
$$D = 2\sum_{i,j} [y_{ij} \log\left(\frac{y_{ij}}{\mu_{ij}}\right) - (y_{ij} - \mu_{ij})] \quad (11)$$
{{< /math >}}

This is the standard formula for the Poisson deviance.

### Connection to the Objective Function

In GLMs, finding the best model parameters (like {{< math >}} $u_i'$ {{< /math >}} and {{< math >}} $v_j$ {{< /math >}} that define {{< math >}} $\mu_{ij}$ {{< /math >}}) is done by maximizing the log-likelihood, which is equivalent to minimizing the deviance.

So, when a model like the one for UMI counts {{< math >}} $y_{ij} \sim \text{Poi}(n_i \exp\{u_i' v_j\})$ {{< /math >}} is used, the goal is to find the {{< math >}} $u_i'$ {{< /math >}} and {{< math >}} $v_j$ {{< /math >}} that minimize this specific Poisson deviance {{< math >}} $D$ {{< /math >}}. This is the count-data equivalent of minimizing squared error in ordinary least squares or PCA for Gaussian data.

The final Poisson deviance formula becomes our optimization objective:

{{< math >}} 
$$\min_{u,v} D = \min_{u,v} \sum_{i,j} [y_{ij} \log\left(\frac{y_{ij}}{\mu_{ij}}\right) - (y_{ij} - \mu_{ij})] \quad (12)$$
{{< /math >}}

where {{< math >}} $\mu_{ij} = s_i \exp(u_i' v_j)$ {{< /math >}} in the context of dimensionality reduction for count data.

## Intercept Terms and Data Centering in GLM-PCA

Constraining the first element of each {{< math >}} $u_i'$ {{< /math >}} to equal 1 (i.e., {{< math >}} $u_{i1} = 1$ {{< /math >}} for all i) in the model {{< math >}} $\eta_{ij} = \log s_i + u_i' v_j$ {{< /math >}} effectively induces a gene-specific intercept term in the first position of each {{< math >}} $v_j$ {{< /math >}}, and this serves an analogous purpose to centering the data in standard PCA.

### Model Expansion with Fixed First Element

In our model, {{< math >}} $\eta_{ij} = \log s_i + u_i' v_j$ {{< /math >}}.

If we specify:
{{< math >}} 
$u_i' = (1, u_{i2}, u_{i3}, \ldots, u_{iL}) \quad (13)$
{{< /math >}}

And:
{{< math >}} 
$v_j = \begin{pmatrix} v_{j1} \\ v_{j2} \\ v_{j3} \\ \vdots \\ v_{jL} \end{pmatrix} \quad (14)$
{{< /math >}}

Then the term {{< math >}} $u_i' v_j$ {{< /math >}} expands to:

{{< math >}} 
$u_i' v_j = (1 \times v_{j1}) + (u_{i2} \times v_{j2}) + (u_{i3} \times v_{j3}) + \cdots + (u_{iL} \times v_{jL}) \quad (15)$
{{< /math >}}

Which simplifies to:
{{< math >}} 
$u_i' v_j = v_{j1} + \sum_{l=2}^L u_{il} v_{jl} \quad (16)$
{{< /math >}}

So the full linear predictor becomes:
{{< math >}} 
$\eta_{ij} = \log s_i + v_{j1} + \sum_{l=2}^L u_{il} v_{jl} \quad (17)$
{{< /math >}}

### The Role of {{< math >}} $v_{j1}$ {{< /math >}} as a Gene-Specific Intercept

The term {{< math >}} $v_{j1}$ {{< /math >}} is now a gene-specific intercept. For each gene j, {{< math >}} $v_{j1}$ {{< /math >}} represents a baseline (log-scale) expression level that is common across all cells, after accounting for library size {{< math >}} $s_i$ {{< /math >}}. It captures the average expression of gene j across the dataset.

The subsequent terms, {{< math >}} $\sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}}, then model the deviations from this baseline for each cell i. These are the variations in gene expression that are explained by the latent dimensions (beyond the overall average).

### Analogous to Centering

This setup is analogous to centering in standard PCA for the following reasons:

#### 1. Removing a Global Mean

In standard PCA, the first step is to center the data by subtracting the feature (gene) means from each observation. This ensures that the subsequent principal components capture the variance around that mean, rather than simply the mean itself.

**Original PCA:**
{{< math >}} 
$X_{\text{centered}} = X - \bar{X} \quad (18)$
{{< /math >}}

where {{< math >}} $\bar{X}$ {{< /math >}} contains the column means of X.

In the GLM, by having {{< math >}} $v_{j1}$ {{< /math >}} as a gene-specific intercept, the model is effectively learning and accounting for the average expression level of each gene. The remaining latent dimensions then capture the residuals or variability around these gene-specific averages.

#### 2. Modeling Deviations from Baseline

Just as centered data represents deviations from the mean, the terms {{< math >}} $\sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}} represent how the expression of gene j in cell i deviates from its expected baseline level ({{< math >}} $v_{j1}$ {{< /math >}}) after accounting for library size ({{< math >}} $\log s_i$ {{< /math >}}). These deviations are what the latent dimensions are designed to explain.

#### 3. Avoiding Capture of Mean by First Component

Without this intercept-like term, the first latent dimension (e.g., {{< math >}} $u_{i1} v_{j1}$ {{< /math >}} if {{< math >}} $u_{i1}$ {{< /math >}} wasn't fixed) would often end up capturing the largest source of variation, which for count data without a specific offset, might simply be the overall high expression of certain genes or cells, rather than biologically meaningful patterns of variability. By explicitly modeling a baseline, the subsequent latent dimensions can focus on more subtle and informative patterns.

Therefore, constraining {{< math >}} $u_{i1} = 1$ {{< /math >}} and allowing {{< math >}} $v_{j1}$ {{< /math >}} to be learned provides a flexible way to effectively "center" the gene expression data on a gene-specific baseline within the GLM framework, similar to how explicit mean subtraction centers data for traditional PCA.

## Orthogonality in GLM-PCA vs Standard PCA

In GLM-PCA (or other factor models using Poisson/Negative Binomial likelihoods), the elements of v (the loadings) won't necessarily be orthogonal.

### Why Orthogonality Is Lost

In standard PCA, the orthogonality of principal components (the v vectors or loadings) is a direct consequence of finding the eigenvectors of the covariance matrix, which is a symmetric matrix. This mathematical property guarantees orthogonal eigenvectors. The objective function in standard PCA is the minimization of squared error, which under a Gaussian assumption, aligns with this orthogonal decomposition.

However, in GLM-PCA (like the Poisson or Negative Binomial model), the objective function changes from minimizing squared error to minimizing deviance (or maximizing log-likelihood) under a non-Gaussian distribution.

Here's why this matters for orthogonality:

1. **Different Objective Function**: The optimization problem is no longer about finding eigenvectors of a covariance matrix. Instead, it's about optimizing a non-linear likelihood function (e.g., the Poisson log-likelihood involves exponentials and logarithms). The mathematical properties that guarantee orthogonality for eigenvectors of symmetric matrices simply don't apply to this new objective.

2. **No Explicit Orthogonality Constraint**: Most implementations of GLM-PCA or count-based factor analysis methods do not explicitly enforce an orthogonality constraint on the v (loading) vectors. While some methods might impose such a constraint for interpretability or uniqueness, it's not inherent to the likelihood maximization itself.

3. **Interpretation**: The "factors" or "loadings" ({{< math >}} $v_j$ {{< /math >}}) in these models represent linear combinations in the link space (e.g., log-scale for Poisson), not directly in the original count space. Their primary goal is to explain the mean of the count data through the link function, not necessarily to be uncorrelated or orthogonal axes of variance.

### What Does This Mean for Interpretation?

While the {{< math >}} $u_i'$ {{< /math >}} and {{< math >}} $v_j$ {{< /math >}} vectors still provide a useful low-dimensional representation and enable visualization similar to PCA, the lack of guaranteed orthogonality for {{< math >}} $v_j$ {{< /math >}} means:

1. **Correlation between Loadings**: The latent dimensions captured by {{< math >}} $v_j$ {{< /math >}} might be correlated with each other. This is different from standard PCA where components are uncorrelated by definition.

2. **Lack of Unique Solution**: Without an orthogonality constraint (or other identifiability constraints), there might be multiple sets of u and v that achieve the same minimum deviance, making the solution less unique. Algorithms often apply other regularization or soft constraints to aid in numerical stability and interpretability, but strict orthogonality isn't usually among them.

## Special Case of Orthogonality: Gaussian GLM-PCA with Identity Link

If GLM-PCA is specifically implemented with a Gaussian error model and an identity link function, the elements of v (the loadings) will be orthogonal. In this particular case, GLM-PCA essentially reduces to standard PCA.

### Gaussian Model + Identity Link

**Gaussian Model**: When the error distribution is Gaussian (Normal), the log-likelihood function (which GLM-PCA optimizes) is directly related to the sum of squared differences between observed and predicted values. Maximizing this likelihood is equivalent to minimizing the sum of squared errors.

**Identity Link Function**: The identity link function means that the linear predictor ({{< math >}} $\eta_{ij}$ {{< /math >}}) directly models the mean ({{< math >}} $\mu_{ij}$ {{< /math >}}), i.e., {{< math >}} $\mu_{ij} = \eta_{ij}$ {{< /math >}}. There's no transformation like a log-link.

### Equivalence to Standard PCA

When you combine these two conditions, the GLM-PCA objective function becomes:

{{< math >}} 
$\text{Minimize } \sum_{i,j} (y_{ij} - \mu_{ij})^2 \quad (19)$
{{< /math >}}

where {{< math >}} $\mu_{ij} = \log s_i + u_i' v_j$ {{< /math >}} (or just {{< math >}} $u_i' v_j$ {{< /math >}} if no offset is used, and assuming {{< math >}} $y_{ij}$ {{< /math >}} are centered).

This is precisely the sum of squared reconstruction errors that standard PCA minimizes. Standard PCA finds the optimal low-rank approximation by performing an eigendecomposition of the covariance matrix (or SVD of the centered data matrix). The eigenvectors obtained from this process are inherently orthogonal.

Therefore, if you're using GLM-PCA but specify a Gaussian family with an identity link, the underlying optimization problem is mathematically equivalent to that of standard PCA, and thus the resulting loadings (v) will be orthogonal. The "GLM" part of "GLM-PCA" refers to its generality to other exponential family distributions and link functions; if you choose the Gaussian-identity case, it specializes back to PCA.

In summary, while GLM-PCA shares the goal of dimensionality reduction and offers similar interpretive insights into latent factors, it achieves this through a different statistical model, which does not inherently enforce the orthogonality of the loading vectors unless specifically configured with Gaussian errors and identity link.


## Linear Predictor in Detail

The linear predictor ({{< math >}} $\eta_{ij}$ {{< /math >}}) for a Generalized Linear Model (GLM), specifically for count data often seen in single-cell RNA sequencing (scRNA-seq) UMI (Unique Molecular Identifier) counts. It shows how the mean of the observed counts ({{< math >}} $\mu_{ij}$ {{< /math >}}) is modeled using a combination of a size factor and a low-dimensional representation.

### 1. Link Function and Linear Predictor

The fundamental relationship in GLMs connects the expected mean to the linear predictor:

{{< math >}} 
$$\log \mu_{ij} = \eta_{ij} \quad (1)$$
{{< /math >}}

Where:
- {{< math >}} $\mu_{ij}$ {{< /math >}}: This is the expected mean count for observation i (e.g., cell i) and feature j (e.g., gene j). For count data, the mean must always be non-negative.

- {{< math >}} $\log \mu_{ij}$ {{< /math >}}: This is the link function. In GLMs, the link function connects the expected mean of the response variable to the linear predictor. For Poisson and Negative Binomial distributions (commonly used for count data), the log link is the canonical choice. It ensures that the predicted mean {{< math >}} $\mu_{ij} = e^{\eta_{ij}}$ {{< /math >}} will always be positive, which is necessary for count distributions.

- {{< math >}} $\eta_{ij}$ {{< /math >}}: This is the linear predictor. It's a linear combination of explanatory variables (or, in this case, latent factors) that directly models the linked mean.

### 2. Components of the Linear Predictor

The linear predictor is defined as:

{{< math >}} 
$$\eta_{ij} = \log s_i + \mathbf{u}'_i \mathbf{v}_j \quad (2)$$
{{< /math >}}

This expression is composed of two main components:

#### The Offset (Library Size Normalization): {{< math >}} $\log s_i$ {{< /math >}}

- {{< math >}} $s_i$ {{< /math >}}: This is often referred to as a size factor or library size. In scRNA-seq, {{< math >}} $s_i$ {{< /math >}} typically represents the total number of UMIs (or reads) counted in cell i. Cells often have vastly different sequencing depths (total counts). Without accounting for this, a gene might appear to have more counts in a deeply sequenced cell simply because that cell yielded more total counts, not because the gene is truly more highly expressed.

- {{< math >}} $\log s_i$ {{< /math >}}: This term is included as an offset in the GLM. An offset is a term whose coefficient is fixed at 1. By including {{< math >}} $\log s_i$ {{< /math >}} as an offset, the model directly incorporates the sequencing depth of each cell. This means that if a cell has twice as many total UMIs, the expected count for a given gene in that cell is also expected to be twice as high, before considering biological expression differences. This is a common form of library size normalization in count-based models.

#### The Low-Rank Approximation: {{< math >}} $\mathbf{u}'_i \mathbf{v}_j$ {{< /math >}}

- {{< math >}} $\mathbf{u}'_i$ {{< /math >}}: This is the i-th row vector from a matrix U (often called factors or cell embeddings). It represents the coordinates of sample (cell) i in a lower-dimensional latent space. Think of it as summarizing the "state" or "type" of cell i across L (the reduced dimension) latent features.

- {{< math >}} $\mathbf{v}_j$ {{< /math >}}: This is the j-th column vector from a matrix {{< math >}} $V^T$ {{< /math >}} (often called loadings or gene programs). It represents the contribution of gene j to each of the L latent dimensions.

- {{< math >}} $\mathbf{u}'_i \mathbf{v}_j$ {{< /math >}}: This is the inner product of the i-th cell's latent state and the j-th gene's loading. It forms a low-rank approximation of the (normalized) gene expression profile. This is the core dimensionality reduction component, similar to the {{< math >}} $\hat{X}_{\text{centered}}$ {{< /math >}} part in PCA, but applied to the log-transformed mean.

### 3. Detailed Structure of the Low-Rank Approximation

The low-rank approximation can be expanded as:

{{< math >}} 
$$\mathbf{u}'_i \mathbf{v}_j = v_{j1} + \sum_{l=2}^L u_{il} v_{jl} \quad (3)$$
{{< /math >}}

This decomposition provides a more detailed look at the structure:

#### Baseline Expression: {{< math >}} $v_{j1}$ {{< /math >}}

{{< math >}} $v_{j1}$ {{< /math >}} represents an intercept term or a baseline expression level for gene j in the linear predictor. It's the contribution of a first latent dimension (l=1) that is common across all cells or represents an average effect for that gene. Essentially, it models the gene's inherent expression level before considering cell-specific factors.

#### Cell-Specific Factors: {{< math >}} $\sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}}

This is the sum over the remaining {{< math >}} $L-1$ {{< /math >}} latent dimensions:

- {{< math >}} $u_{il}$ {{< /math >}}: The component for cell i on the l-th latent dimension.
- {{< math >}} $v_{jl}$ {{< /math >}}: The loading (weight) for gene j on the l-th latent dimension.

This sum models how cell i's specific characteristics ({{< math >}} $u_{il}$ {{< /math >}}) interact with gene j's responsiveness to those characteristics ({{< math >}} $v_{jl}$ {{< /math >}}) to influence its expression. These are the components that capture biological variability and allow for cell clustering or trajectory inference.

### Complete Model Formulation

Combining all components, the complete model can be written as:

{{< math >}} 
$$\log \mu_{ij} = \log s_i + v_{j1} + \sum_{l=2}^L u_{il} v_{jl} \quad (4)$$
{{< /math >}}

Or equivalently:

{{< math >}} 
$$\mu_{ij} = s_i \exp\left(v_{j1} + \sum_{l=2}^L u_{il} v_{jl}\right) \quad (5)$$
{{< /math >}}

### Model Interpretation

The entire expression defines a statistical model where the expected (log-transformed) count of gene j in cell i is explained by:

1. **Sequencing depth of cell i** ({{< math >}} $\log s_i$ {{< /math >}}): Accounts for technical variation in library size
2. **Baseline expression for gene j** ({{< math >}} $v_{j1}$ {{< /math >}}): Captures the inherent expression level of each gene
3. **Cell-gene interaction terms** ({{< math >}} $\sum_{l=2}^L u_{il} v_{jl}$ {{< /math >}}): A linear combination of latent factors that capture:
   - {{< math >}} $u_{il}$ {{< /math >}}: The cell's state across different biological dimensions
   - {{< math >}} $v_{jl}$ {{< /math >}}: How each gene contributes to those biological dimensions

This framework allows for robust dimensionality reduction of UMI count data, accounting for its specific statistical properties (non-negativity, integer counts, and often overdispersion) while simultaneously addressing technical variability (library size). The model enables:

- **Cell type identification** through clustering in the latent space defined by {{< math >}} $\mathbf{u}_i$ {{< /math >}}
- **Gene program discovery** through the loading patterns in {{< math >}} $\mathbf{v}_j$ {{< /math >}}
- **Trajectory inference** by analyzing smooth changes in the latent factors
- **Differential expression analysis** that accounts for both technical and biological sources of variation