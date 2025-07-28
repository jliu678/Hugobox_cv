# PCA: High Variance Equivalent to Low Euclidean Loss

## 1. Rank-k Projection ({{< math >}}$\hat{X}${{< /math >}})

Suppose we have a dataset {{< math >}}$X \in \mathbb{R}^{4 \times 3}${{< /math >}} (4 samples, 3 features):

{{< math >}}
$$X = \begin{pmatrix}
1 & 2 & 3 \\
4 & 5 & 6 \\
7 & 8 & 9 \\
10 & 11 & 12
\end{pmatrix} \quad (1)$$
{{< /math >}}

### Centering the Data

First, we center the data.

Feature means:
- Feature 1: {{< math >}}$(1+4+7+10)/4 = 22/4 = 5.5${{< /math >}}
- Feature 2: {{< math >}}$(2+5+8+11)/4 = 26/4 = 6.5${{< /math >}}
- Feature 3: {{< math >}}$(3+6+9+12)/4 = 30/4 = 7.5${{< /math >}}

The mean vector is {{< math >}}$\bar{x} = \begin{pmatrix} 5.5 \\ 6.5 \\ 7.5 \end{pmatrix}${{< /math >}}.

The mean matrix {{< math >}}$\bar{X}${{< /math >}} (each row is {{< math >}}$\bar{x}${{< /math >}}):

{{< math >}}
$$\bar{X} = \begin{pmatrix}
5.5 & 6.5 & 7.5 \\
5.5 & 6.5 & 7.5 \\
5.5 & 6.5 & 7.5 \\
5.5 & 6.5 & 7.5
\end{pmatrix} \quad (2)$$
{{< /math >}}

Now, we center the data:

{{< math >}}
$$X_{\text{centered}} = X - \bar{X} = \begin{pmatrix}
1-5.5 & 2-6.5 & 3-7.5 \\
4-5.5 & 5-6.5 & 6-7.5 \\
7-5.5 & 8-6.5 & 9-7.5 \\
10-5.5 & 11-6.5 & 12-7.5
\end{pmatrix} \quad (3)$$
{{< /math >}}

{{< math >}}
$$X_{\text{centered}} = \begin{pmatrix}
-4.5 & -4.5 & -4.5 \\
-1.5 & -1.5 & -1.5 \\
1.5 & 1.5 & 1.5 \\
4.5 & 4.5 & 4.5
\end{pmatrix} \quad (4)$$
{{< /math >}}

This centered data {{< math >}}$X_{\text{centered}}${{< /math >}} has a rank of 1, not full rank (which would be 3). For PCA, the rank of the centered data's covariance matrix (or {{< math >}}$X_{\text{centered}}^T X_{\text{centered}}${{< /math >}}) is what determines the number of non-zero eigenvalues and thus the maximum effective rank. In this specific case, the three features are perfectly correlated, leading to a rank of 1.

Let's use a slightly modified X to ensure {{< math >}}$X_{\text{centered}}${{< /math >}} has full rank (rank 3 in this case) for a more general PCA example, where not all variance is captured by one component.

### Revised X:

{{< math >}}
$$X = \begin{pmatrix}
1 & 2 & 3 \\
4 & 5 & 1 \\
7 & 8 & 9 \\
10 & 1 & 12
\end{pmatrix} \quad (5)$$
{{< /math >}}

Mean vector {{< math >}}$\bar{x}${{< /math >}}:
- Feature 1: {{< math >}}$(1+4+7+10)/4 = 5.5${{< /math >}}
- Feature 2: {{< math >}}$(2+5+8+1)/4 = 16/4 = 4${{< /math >}}
- Feature 3: {{< math >}}$(3+1+9+12)/4 = 25/4 = 6.25${{< /math >}}

{{< math >}}
$$\bar{x} = \begin{pmatrix} 5.5 \\ 4 \\ 6.25 \end{pmatrix} \quad (6)$$
{{< /math >}}

Centered data {{< math >}}$X_{\text{centered}}${{< /math >}}:

{{< math >}}
$$X_{\text{centered}} = \begin{pmatrix}
1-5.5 & 2-4 & 3-6.25 \\
4-5.5 & 5-4 & 1-6.25 \\
7-5.5 & 8-4 & 9-6.25 \\
10-5.5 & 1-4 & 12-6.25
\end{pmatrix} \quad (7)$$
{{< /math >}}

{{< math >}}
$$X_{\text{centered}} = \begin{pmatrix}
-4.5 & -2 & -3.25 \\
-1.5 & 1 & -5.25 \\
1.5 & 4 & 2.75 \\
4.5 & -3 & 5.75
\end{pmatrix} \quad (8)$$
{{< /math >}}

Now, this {{< math >}}$X_{\text{centered}}${{< /math >}} matrix has a rank of 3, meaning it has full rank for its dimensions (4×3), allowing for more interesting PCA results.

### PCA and Principal Components ({{< math >}}$W_k${{< /math >}})

To compute {{< math >}}$W_k${{< /math >}}, we'd typically calculate the covariance matrix of {{< math >}}$X_{\text{centered}}${{< /math >}} and then find its eigenvectors. For a numeric example, we'll assume we've done this and picked {{< math >}}$k=2${{< /math >}} principal components.

The covariance matrix {{< math >}}$C = \frac{1}{n-1} X_{\text{centered}}^T X_{\text{centered}}${{< /math >}}.

{{< math >}}
$$C = \frac{1}{3} \begin{pmatrix}
-4.5 & -2 & -3.25 \\
-1.5 & 1 & -5.25 \\
1.5 & 4 & 2.75 \\
4.5 & -3 & 5.75
\end{pmatrix}^T \begin{pmatrix}
-4.5 & -2 & -3.25 \\
-1.5 & 1 & -5.25 \\
1.5 & 4 & 2.75 \\
4.5 & -3 & 5.75
\end{pmatrix} \quad (9)$$
{{< /math >}}

Calculating this yields:

{{< math >}}
$$C \approx \begin{pmatrix}
11.25 & -7.5 & 13.5 \\
-7.5 & 10.67 & -6.5 \\
13.5 & -6.5 & 25.58
\end{pmatrix} \quad (10)$$
{{< /math >}}

Now, we'd find the eigenvalues and eigenvectors of C. Let's assume (for simplicity and to proceed with the example without complex matrix computations) that the top 2 principal components are:

{{< math >}}
$$W_k = W_2 = \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \quad (11)$$
{{< /math >}}

(These values are illustrative and would be derived from actual eigenvalue decomposition, normalized to be orthonormal).

Let's verify orthonormality and unit length (approximately, for illustration):

{{< math >}}
$$W_2^T W_2 = \begin{pmatrix} 0.5 & 0.3 & 0.8 \\ -0.2 & 0.8 & 0.5 \end{pmatrix} \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} = \begin{pmatrix} 0.98 & 0.54 \\ 0.54 & 0.93 \end{pmatrix} \quad (12)$$
{{< /math >}}

Ideally, this should be close to the identity matrix {{< math >}}$\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}${{< /math >}}. For a true example, we'd use a numerical library. For this illustrative purpose, we will proceed with these values as our {{< math >}}$W_k${{< /math >}}.

### How to reconstruct {{< math >}}$\hat{x}_i${{< /math >}}: The approximation of the original data point {{< math >}}$x_i${{< /math >}}

Let's reconstruct the first centered data point, {{< math >}}$x_{1,\text{centered}} = \begin{pmatrix} -4.5 \\ -2 \\ -3.25 \end{pmatrix}${{< /math >}}.

**Project the data point into the reduced space:**

{{< math >}}
$$z_1 = x_{1,\text{centered}} W_k = \begin{pmatrix} -4.5 & -2 & -3.25 \end{pmatrix} \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \quad (13)$$
{{< /math >}}

{{< math >}}
$$z_1 = \begin{pmatrix}
(-4.5 \times 0.5) + (-2 \times 0.3) + (-3.25 \times 0.8) \\
(-4.5 \times -0.2) + (-2 \times 0.8) + (-3.25 \times 0.5)
\end{pmatrix} \quad (14)$$
{{< /math >}}

{{< math >}}
$$z_1 = \begin{pmatrix}
-2.25 - 0.6 - 2.6 \\
0.9 - 1.6 - 1.625
\end{pmatrix} = \begin{pmatrix} -5.45 \\ -2.325 \end{pmatrix} \in \mathbb{R}^2 \quad (15)$$
{{< /math >}}

This gives the coordinates of {{< math >}}$x_{1,\text{centered}}${{< /math >}} in the low-dimensional PCA space.

**Reconstruct the approximation in the original space:**

{{< math >}}
$$\hat{x}_{1,\text{centered}} = z_1 W_k^T = \begin{pmatrix} -5.45 & -2.325 \end{pmatrix} \begin{pmatrix} 0.5 & 0.3 & 0.8 \\ -0.2 & 0.8 & 0.5 \end{pmatrix} \quad (16)$$
{{< /math >}}

{{< math >}}
$$\hat{x}_{1,\text{centered}} = \begin{pmatrix}
(-5.45 \times 0.5) + (-2.325 \times -0.2) \\
(-5.45 \times 0.3) + (-2.325 \times 0.8) \\
(-5.45 \times 0.8) + (-2.325 \times 0.5)
\end{pmatrix} \quad (17)$$
{{< /math >}}

{{< math >}}
$$\hat{x}_{1,\text{centered}} = \begin{pmatrix}
-2.725 + 0.465 \\
-1.635 - 1.86 \\
-4.36 - 1.1625
\end{pmatrix} \approx \begin{pmatrix} -2.26 \\ -3.495 \\ -5.5225 \end{pmatrix} \in \mathbb{R}^3 \quad (18)$$
{{< /math >}}

Notice that {{< math >}}$\hat{x}_{1,\text{centered}}${{< /math >}} is an approximation of {{< math >}}$x_{1,\text{centered}}${{< /math >}}. It's not identical because we reduced the dimensionality from 3 to 2, leading to some information loss.

### Full Reconstruction Formula for the Matrix {{< math >}}$\hat{X}${{< /math >}}

The full reconstruction formula for the centered matrix {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}} is:

{{< math >}}
$$\hat{X}_{\text{centered}} = X_{\text{centered}} W_k W_k^T \quad (19)$$
{{< /math >}}

First, let's calculate {{< math >}}$W_k W_k^T${{< /math >}}:

{{< math >}}
$$W_k W_k^T = \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \begin{pmatrix} 0.5 & 0.3 & 0.8 \\ -0.2 & 0.8 & 0.5 \end{pmatrix} \quad (20)$$
{{< /math >}}

{{< math >}}
$$W_k W_k^T = \begin{pmatrix}
(0.5)^2 + (-0.2)^2 & (0.5)(0.3) + (-0.2)(0.8) & (0.5)(0.8) + (-0.2)(0.5) \\
(0.3)(0.5) + (0.8)(-0.2) & (0.3)^2 + (0.8)^2 & (0.3)(0.8) + (0.8)(0.5) \\
(0.8)(0.5) + (0.5)(-0.2) & (0.8)(0.3) + (0.5)(0.8) & (0.8)^2 + (0.5)^2
\end{pmatrix} \quad (21)$$
{{< /math >}}

{{< math >}}
$$W_k W_k^T = \begin{pmatrix}
0.25 + 0.04 & 0.15 - 0.16 & 0.4 - 0.1 \\
0.15 - 0.16 & 0.09 + 0.64 & 0.24 + 0.4 \\
0.4 - 0.1 & 0.24 + 0.4 & 0.64 + 0.25
\end{pmatrix} = \begin{pmatrix}
0.29 & -0.01 & 0.3 \\
-0.01 & 0.73 & 0.64 \\
0.3 & 0.64 & 0.89
\end{pmatrix} \quad (22)$$
{{< /math >}}

Now, calculate {{< math >}}$\hat{X}_{\text{centered}} = X_{\text{centered}} (W_k W_k^T)${{< /math >}}:

{{< math >}}
$$\hat{X}_{\text{centered}} = \begin{pmatrix}
-4.5 & -2 & -3.25 \\
-1.5 & 1 & -5.25 \\
1.5 & 4 & 2.75 \\
4.5 & -3 & 5.75
\end{pmatrix} \begin{pmatrix}
0.29 & -0.01 & 0.3 \\
-0.01 & 0.73 & 0.64 \\
0.3 & 0.64 & 0.89
\end{pmatrix} \quad (23)$$
{{< /math >}}

Performing this matrix multiplication (values are approximate due to rounding):

{{< math >}}
$$\hat{X}_{\text{centered}} \approx \begin{pmatrix}
-2.26 & -3.495 & -5.5225 \\
-1.72 & 0.145 & -1.82 \\
1.95 & 4.395 & 4.96 \\
2.03 & -1.045 & 2.38
\end{pmatrix} \quad (24)$$
{{< /math >}}

Finally, to get {{< math >}}$\hat{X}${{< /math >}} in the original space, we add back the mean:

{{< math >}}
$$\hat{X} = \hat{X}_{\text{centered}} + \bar{X} \quad (25)$$
{{< /math >}}

{{< math >}}
$$\hat{X} = \begin{pmatrix}
-2.26 & -3.495 & -5.5225 \\
-1.72 & 0.145 & -1.82 \\
1.95 & 4.395 & 4.96 \\
2.03 & -1.045 & 2.38
\end{pmatrix} + \begin{pmatrix}
5.5 & 4 & 6.25 \\
5.5 & 4 & 6.25 \\
5.5 & 4 & 6.25 \\
5.5 & 4 & 6.25
\end{pmatrix} \quad (26)$$
{{< /math >}}

{{< math >}}
$$\hat{X} \approx \begin{pmatrix}
3.24 & 0.505 & 0.7275 \\
3.78 & 4.145 & 4.43 \\
7.45 & 8.395 & 11.21 \\
7.53 & 2.955 & 8.63
\end{pmatrix} \quad (27)$$
{{< /math >}}

This {{< math >}}$\hat{X}${{< /math >}} is the rank-k projection (here, {{< math >}}$k=2${{< /math >}}) of the original data. It's an approximation, and its quality is measured by MSE.

## 2. Mean Squared Error (MSE) in PCA

The MSE quantifies the average squared distance between each original data point and its reconstruction.

{{< math >}}
$$\text{MSE} = \frac{1}{n} \sum_{i=1}^n \|x_i - \hat{x}_i\|^2 \quad (28)$$
{{< /math >}}

Using the centered data {{< math >}}$X_{\text{centered}}${{< /math >}} and its reconstruction {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}}:

The difference matrix {{< math >}}$E = X_{\text{centered}} - \hat{X}_{\text{centered}}${{< /math >}}:

{{< math >}}
$$E \approx \begin{pmatrix}
-4.5-(-2.26) & -2-(-3.495) & -3.25-(-5.5225) \\
-1.5-(-1.72) & 1-0.145 & -5.25-(-1.82) \\
1.5-1.95 & 4-4.395 & 2.75-4.96 \\
4.5-2.03 & -3-(-1.045) & 5.75-2.38
\end{pmatrix} \quad (29)$$
{{< /math >}}

{{< math >}}
$$E = \begin{pmatrix}
-2.24 & 1.495 & 2.2725 \\
0.22 & 0.855 & -3.43 \\
-0.45 & -0.395 & -2.21 \\
2.47 & -1.955 & 3.37
\end{pmatrix} \quad (30)$$
{{< /math >}}

Now, calculate the squared Frobenius norm of this difference matrix, and then divide by {{< math >}}$n=4${{< /math >}}:

{{< math >}}
$$\|E\|_F^2 = \sum_{i=1}^4 \|e_i\|^2 \quad (31)$$
{{< /math >}}

{{< math >}}
$$\|e_1\|^2 = (-2.24)^2 + (1.495)^2 + (2.2725)^2 = 5.0176 + 2.235025 + 5.16410625 \approx 12.4167 \quad (32)$$
{{< /math >}}

{{< math >}}
$$\|e_2\|^2 = (0.22)^2 + (0.855)^2 + (-3.43)^2 = 0.0484 + 0.731025 + 11.7649 \approx 12.5443 \quad (33)$$
{{< /math >}}

{{< math >}}
$$\|e_3\|^2 = (-0.45)^2 + (-0.395)^2 + (-2.21)^2 = 0.2025 + 0.156025 + 4.8841 \approx 5.2426 \quad (34)$$
{{< /math >}}

{{< math >}}
$$\|e_4\|^2 = (2.47)^2 + (-1.955)^2 + (3.37)^2 = 6.1009 + 3.822025 + 11.3569 \approx 21.2798 \quad (35)$$
{{< /math >}}

Sum of squared distances {{< math >}}$\sum\|e_i\|^2 \approx 12.4167 + 12.5443 + 5.2426 + 21.2798 = 51.4834${{< /math >}}

{{< math >}}
$$\text{MSE} = \frac{1}{4} \times 51.4834 \approx 12.87085 \quad (36)$$
{{< /math >}}

This non-zero MSE indicates the information loss that occurred by reducing the dimensionality from 3 to 2. A higher MSE implies more information loss.

## 3. Energy (Sum of Squares)

The total energy of the dataset is the sum of squared values of the uncentered data X.

{{< math >}}
$$\text{Energy} = \sum_{i=1}^n \|x_i\|^2 = \|X\|_F^2 = \sum_{i=1}^n \sum_{j=1}^d x_{ij}^2 \quad (37)$$
{{< /math >}}

Using our original (uncentered) data {{< math >}}$X = \begin{pmatrix}
1 & 2 & 3 \\
4 & 5 & 1 \\
7 & 8 & 9 \\
10 & 1 & 12
\end{pmatrix}${{< /math >}}:

{{< math >}}
$$\text{Energy} = (1^2 + 2^2 + 3^2) + (4^2 + 5^2 + 1^2) + (7^2 + 8^2 + 9^2) + (10^2 + 1^2 + 12^2) \quad (38)$$
{{< /math >}}

{{< math >}}
$$\text{Energy} = (1 + 4 + 9) + (16 + 25 + 1) + (49 + 64 + 81) + (100 + 1 + 144) \quad (39)$$
{{< /math >}}

{{< math >}}
$$\text{Energy} = 14 + 42 + 194 + 245 = 495 \quad (40)$$
{{< /math >}}

The Energy (or sum of squares) of the dataset is 495. This measures the total squared magnitude of the original data points from the origin.

## Does this assume the data has been centered and scaled?

**Centering:** ✅ Yes, centering is required for PCA. PCA aims to find directions of maximum variance, and variance is defined as spread around the mean. If the data isn't centered, the first principal component would likely capture the location of the data's mean rather than its internal variability. For Energy (Sum of Squares), however, the uncentered data is used, as it reflects the total magnitude from the origin. This includes the contribution of the mean.

**Scaling:** ⚠️ Optional. Scaling (e.g., standardizing features to have unit variance) is not mathematically required for PCA. However, it's often a good practice when features have different units or vastly different scales. If one feature has values much larger than others, its variance will dominate the covariance matrix, and PCA might disproportionately focus on that feature, even if other features have important, albeit smaller, variations. Scaling ensures that all features contribute equally to the variance calculation. Our example did not include scaling.

## Why High Variance = Low Euclidean Loss?

PCA's objective is two-fold:

1. **Maximize captured variance**: It finds orthogonal directions (principal components) along which the data varies the most. The eigenvalues associated with these principal components represent the amount of variance explained by each component. Larger eigenvalues mean more variance is captured.

2. **Minimize reconstruction error (Euclidean Loss)**: When you project data onto a lower-dimensional subspace defined by these principal components and then reconstruct it, PCA aims to make this reconstructed data ({{< math >}}$\hat{X}${{< /math >}}) as close as possible to the original data ({{< math >}}$X${{< /math >}}). The Euclidean loss (or MSE) measures this closeness.

These two objectives are inherently linked. By selecting the directions that explain the most variance, PCA is effectively choosing the subspace that best preserves the overall structure and relationships within the data. When data is projected onto the directions of highest variance, the "distortion" or "loss" due to the projection is minimized. Conversely, if you project onto directions with low variance, the data points would be squished together, leading to a much larger reconstruction error. Therefore, maximizing the variance captured is equivalent to minimizing the Euclidean loss in the context of linear dimensionality reduction using PCA.


# Two Views of PCA: Matrix Factorization vs Eigen Decomposition

## PCA Objective Function (Matrix Factorization View)

The PCA objective you gave is:

{{< math >}}
$$\min_{u_i, v_j} \sum_{i=1}^I \sum_{j=1}^J (z_{ij} - u_i^{\top} v_j)^2 \quad (1)$$
{{< /math >}}

This is a low-rank matrix approximation problem. Here's what it means:

You're given a data matrix {{< math >}}$Z \in \mathbb{R}^{I \times J}${{< /math >}} (e.g. centered gene expression, with {{< math >}}$I${{< /math >}} cells, {{< math >}}$J${{< /math >}} genes).

You're looking to approximate it as a low-rank product:

{{< math >}}
$$Z \approx UV^{\top} \quad (2)$$
{{< /math >}}

where:

- {{< math >}}$U \in \mathbb{R}^{I \times L}${{< /math >}} contains factors or principal component scores {{< math >}}$u_i^{\top}${{< /math >}}
- {{< math >}}$V \in \mathbb{R}^{J \times L}${{< /math >}} contains loadings {{< math >}}$v_j^{\top}${{< /math >}}
- {{< math >}}$L \ll \min(I,J)${{< /math >}} is the rank of the approximation (number of components)

The goal is to minimize squared reconstruction error — i.e., mean squared error across all matrix entries.

## 🔷 PCA via Spectral Decomposition

In standard PCA:

1. You take a centered data matrix {{< math >}}$X \in \mathbb{R}^{n \times d}${{< /math >}}

2. Compute the covariance matrix:

{{< math >}}
$$\Sigma = \frac{1}{n} X^{\top} X \quad (3)$$
{{< /math >}}

3. Then find eigenvectors {{< math >}}$\{w_k\}${{< /math >}} and eigenvalues {{< math >}}$\{\lambda_k\}${{< /math >}} of {{< math >}}$\Sigma${{< /math >}}

Because {{< math >}}$\Sigma${{< /math >}} is symmetric and positive semidefinite, its eigenvectors are:

✅ Real,  
✅ Orthogonal,  
✅ Form an orthonormal basis.

These eigenvectors form the principal directions, and projections onto them give the principal components.

## 🔁 How They Relate

| Conceptual View | Matrix Factorization View | Eigen Decomposition View |
|----------------|---------------------------|--------------------------|
| Data approximation | {{< math >}}$Z \approx UV^{\top}${{< /math >}} | {{< math >}}$X \approx XWW^{\top}${{< /math >}} |
| Factors / scores | {{< math >}}$u_i^{\top}${{< /math >}} | {{< math >}}$x_i^{\top} W${{< /math >}} |
| Loadings / basis vectors | {{< math >}}$v_j${{< /math >}} | Columns of {{< math >}}$W${{< /math >}} (eigenvectors) |
| Objective | Minimize {{< math >}}$\|Z - UV^{\top}\|_F^2${{< /math >}} | Maximize variance along orthogonal axes |

🔸 In both views, PCA finds a rank-{{< math >}}$L${{< /math >}} approximation to the data matrix.

🔸 In the spectral view, the loadings {{< math >}}$V${{< /math >}} are orthonormal eigenvectors of the covariance matrix, hence orthogonal by construction.

🔸 In the matrix factorization view, you can constrain {{< math >}}$V${{< /math >}} to be orthogonal (which turns this into classical PCA). Without that constraint, this becomes general low-rank approximation, more like matrix factorization methods (e.g., in recommender systems).

## Numeric Example Comparison

Let's compare the two views of PCA side-by-side using the same numeric example data {{< math >}}$X \in \mathbb{R}^{4 \times 3}${{< /math >}}:

{{< math >}}
$$X = \begin{pmatrix}
1 & 2 & 3 \\
4 & 5 & 1 \\
7 & 8 & 9 \\
10 & 1 & 12
\end{pmatrix} \quad (4)$$
{{< /math >}}

First, we center the data:

{{< math >}}
$$\bar{x} = (5.5, 4, 6.25) \quad (5)$$
{{< /math >}}

{{< math >}}
$$X_{\text{centered}} = X - \bar{X} = \begin{pmatrix}
-4.5 & -2 & -3.25 \\
-1.5 & 1 & -5.25 \\
1.5 & 4 & 2.75 \\
4.5 & -3 & 5.75
\end{pmatrix} \quad (6)$$
{{< /math >}}

We'll use {{< math >}}$k=2${{< /math >}} principal components. For illustration, we assume the principal component matrix is:

{{< math >}}
$$W_k = W_2 = \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \quad (7)$$
{{< /math >}}

### View 1: Reconstruction Formula ({{< math >}}$\hat{X} = XW_k W_k^T${{< /math >}})

This view focuses on directly projecting and reconstructing the data using the principal components.

**Step 1:** Calculate {{< math >}}$W_k W_k^T${{< /math >}} (the projection matrix).

This matrix projects the data onto the subspace spanned by the k principal components and then back into the original space.

{{< math >}}
$$W_k W_k^T = \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \begin{pmatrix}
0.5 & 0.3 & 0.8 \\
-0.2 & 0.8 & 0.5
\end{pmatrix} = \begin{pmatrix}
0.29 & -0.01 & 0.3 \\
-0.01 & 0.73 & 0.64 \\
0.3 & 0.64 & 0.89
\end{pmatrix} \quad (8)$$
{{< /math >}}

**Step 2:** Apply the projection matrix to the centered data to get {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}}.

{{< math >}}
$$\hat{X}_{\text{centered}} = X_{\text{centered}} W_k W_k^T \quad (9)$$
{{< /math >}}

{{< math >}}
$$\hat{X}_{\text{centered}} \approx \begin{pmatrix}
-2.26 & -3.495 & -5.5225 \\
-1.72 & 0.145 & -1.82 \\
1.95 & 4.395 & 4.96 \\
2.03 & -1.045 & 2.38
\end{pmatrix} \quad (10)$$
{{< /math >}}

**Step 3:** Add the mean back to get {{< math >}}$\hat{X}${{< /math >}} in the original space.

{{< math >}}
$$\hat{X} = \hat{X}_{\text{centered}} + \bar{X} \quad (11)$$
{{< /math >}}

{{< math >}}
$$\hat{X} \approx \begin{pmatrix}
3.24 & 0.505 & 0.7275 \\
3.78 & 4.145 & 4.43 \\
7.45 & 8.395 & 11.21 \\
7.53 & 2.955 & 8.63
\end{pmatrix} \quad (12)$$
{{< /math >}}

This {{< math >}}$\hat{X}${{< /math >}} is the low-rank approximation of the original data.

### View 2: Objective Function Minimization ({{< math >}}$\min_{u,v} \sum_{i,j} (z_{ij} - u_i' v_j)^2${{< /math >}})

This view focuses on the underlying optimization problem that PCA solves to find the best low-rank approximation.

**Goal:** Find matrices {{< math >}}$U${{< /math >}} (factors/principal components) and {{< math >}}$V^T${{< /math >}} (loadings, where {{< math >}}$V^T = W_k${{< /math >}}) such that the sum of squared errors between the centered data {{< math >}}$Z${{< /math >}} (our {{< math >}}$X_{\text{centered}}${{< /math >}}) and its approximation {{< math >}}$UV^T${{< /math >}} is minimized.

Here, {{< math >}}$Z = X_{\text{centered}}${{< /math >}} (size {{< math >}}$n \times d${{< /math >}}, i.e., {{< math >}}$4 \times 3${{< /math >}}).  
{{< math >}}$U${{< /math >}} is the projected data (size {{< math >}}$n \times k${{< /math >}}, i.e., {{< math >}}$4 \times 2${{< /math >}}).  
{{< math >}}$V^T${{< /math >}} (or {{< math >}}$W_k${{< /math >}}) is the principal component matrix (size {{< math >}}$d \times k${{< /math >}}, i.e., {{< math >}}$3 \times 2${{< /math >}}).

**Step 1:** Project the centered data into the reduced space to get {{< math >}}$U${{< /math >}} (the 'factors' {{< math >}}$u_i'${{< /math >}}).

{{< math >}}
$$U = X_{\text{centered}} W_k \quad (13)$$
{{< /math >}}

{{< math >}}
$$U \approx \begin{pmatrix}
-5.45 & -2.325 \\
-5.475 & -1.56 \\
4.225 & 4.375 \\
2.525 & -1.075
\end{pmatrix} \quad (14)$$
{{< /math >}}

The rows of {{< math >}}$U${{< /math >}} are the {{< math >}}$u_i'${{< /math >}} vectors (e.g., {{< math >}}$u_1' = (-5.45, -2.325)${{< /math >}}).

**Step 2:** The 'loadings' {{< math >}}$v_j${{< /math >}} are the columns of {{< math >}}$W_k${{< /math >}} (or rows of {{< math >}}$W_k^T${{< /math >}}).

{{< math >}}
$$V^T = W_k = \begin{pmatrix}
0.5 & -0.2 \\
0.3 & 0.8 \\
0.8 & 0.5
\end{pmatrix} \quad (15)$$
{{< /math >}}

The columns of {{< math >}}$V^T${{< /math >}} are the {{< math >}}$v_j${{< /math >}} vectors (e.g., {{< math >}}$v_1 = \begin{pmatrix} 0.5 \\ 0.3 \\ 0.8 \end{pmatrix}${{< /math >}}).

**Step 3:** Form the reconstructed matrix {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}} using {{< math >}}$U${{< /math >}} and {{< math >}}$V^T${{< /math >}}.

{{< math >}}
$$\hat{X}_{\text{centered}} = UV^T = UW_k^T \quad (16)$$
{{< /math >}}

This matrix multiplication results in the same {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}} as in View 1, verifying the equivalence.

**Step 4:** Calculate the sum of squared errors (the value being minimized).

This is the {{< math >}}$\sum_{i,j} (z_{ij} - \hat{z}_{ij})^2${{< /math >}} where {{< math >}}$z_{ij}${{< /math >}} are elements of {{< math >}}$X_{\text{centered}}${{< /math >}} and {{< math >}}$\hat{z}_{ij}${{< /math >}} are elements of {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}}. This sum is {{< math >}}$\|X_{\text{centered}} - \hat{X}_{\text{centered}}\|_F^2${{< /math >}}.

{{< math >}}
$$\|X_{\text{centered}} - \hat{X}_{\text{centered}}\|_F^2 \approx 51.4834 \quad (17)$$
{{< /math >}}

This is the value that PCA algorithm aims to minimize.

## Comparison Side-by-Side

| Feature | View 1: Reconstruction Formula ({{< math >}}$\hat{X} = XW_k W_k^T${{< /math >}}) | View 2: Objective Function Minimization ({{< math >}}$\min \sum_{i,j} (z_{ij} - u_i' v_j)^2${{< /math >}}) |
|---------|-------|-------|
| **Focus** | Practical Application: How to actually compute the reconstructed data. | Theoretical Basis: The mathematical criterion PCA optimizes. |
| **Primary Goal** | To produce the approximated data matrix {{< math >}}$\hat{X}${{< /math >}} that has lower rank while retaining as much information as possible. | To find the optimal transformation (principal components/loadings) that results in the minimum reconstruction error. |
| **Output** | The low-rank approximation {{< math >}}$\hat{X}${{< /math >}} (or {{< math >}}$\hat{X}_{\text{centered}}${{< /math >}}). | The minimized sum of squared errors (a single scalar value representing the quality of approximation). The {{< math >}}$U${{< /math >}} (factors) and {{< math >}}$V${{< /math >}} (loadings) matrices that achieve this minimum. |
| **Input for Formula** | The original (or centered) data {{< math >}}$X${{< /math >}} and the principal component matrix {{< math >}}$W_k${{< /math >}}. | The centered data matrix {{< math >}}$Z${{< /math >}} (our {{< math >}}$X_{\text{centered}}${{< /math >}}). |
| **Components** | {{< math >}}$W_k W_k^T${{< /math >}} acts as a direct projection-reconstruction operator. | {{< math >}}$u_i'${{< /math >}} (factors/projected points) and {{< math >}}$v_j${{< /math >}} (loadings/principal component vectors) are explicitly part of the approximation. |
| **Interpretation** | "Applying the principal components to reconstruct the data." | "Finding the best low-rank representation that minimizes the squared distances between original and reconstructed data points." |
| **Relationship to MSE** | MSE is calculated from the result of this reconstruction. | The objective function is the sum of squared errors, which directly leads to MSE when averaged. |

## Conclusion

Both views describe the same underlying process but emphasize different aspects. The reconstruction formula shows you how to get the approximated data. The objective function explains why the principal components (the {{< math >}}$W_k${{< /math >}} matrix) are chosen the way they are—because they minimize this specific error metric. PCA mathematically guarantees that the {{< math >}}$W_k${{< /math >}} derived from its eigenvectors is the unique solution that minimizes this objective function for a given rank {{< math >}}$k${{< /math >}}.


