# A) Derivation: Transforming the Solution for s(t)

Let's derive the relationship $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$ from the provided information.

## Given Information

### System of Differential Equations:

{{< math >}}
$$
\frac{du(t)}{dt} = \alpha^{(k)} - \beta u(t)
$$
{{< /math >}}

{{< math >}}
$$
\frac{ds(t)}{dt} = \beta u(t) - \gamma s(t)
$$
{{< /math >}}

### Solutions for u(t) and s(t):

{{< math >}}
$$
u(t) = u_0 e^{-\beta\tau} + \frac{\alpha^{(k)}}{\beta} (1 - e^{-\beta\tau})
$$
{{< /math >}}

(This can be rewritten as $u(t) = \frac{\alpha^{(k)}}{\beta} + (u_0 - \frac{\alpha^{(k)}}{\beta})e^{-\beta\tau}$)

{{< math >}}
$$
s(t) = s_0 e^{-\gamma\tau} + \frac{\alpha^{(k)}}{\gamma} (1 - e^{-\gamma\tau}) + \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) (e^{-\gamma\tau} - e^{-\beta\tau})
$$
{{< /math >}}

(This can be expanded to $s(t) = \frac{\alpha^{(k)}}{\gamma} + (s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta})e^{-\gamma\tau} - (\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta})e^{-\beta\tau}$)

Where $\tau = t - t_0$.

### Steady-State Values:

From setting the derivatives to zero:

{{< math >}}
$$
u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}
$$
{{< /math >}}

{{< math >}}
$$
s_\infty = \frac{\beta u_\infty^{(k)}}{\gamma} = \frac{\alpha^{(k)}}{\gamma}
$$
{{< /math >}}

### New Variable Definitions:

{{< math >}}
$$
\tilde{\beta} := \frac{\beta}{\gamma - \beta}
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}(t) := s(t) - \tilde{\beta}u(t)
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}_\infty := s_\infty - \tilde{\beta}u_\infty^{(k)}
$$
{{< /math >}}

## Derivation Steps

Our goal is to demonstrate: $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$.

### 1. Express u(t) in a convenient form:

Let's rewrite u(t) using its steady-state value $u_\infty^{(k)}$:

{{< math >}}
$$
u(t) = u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}
$$
{{< /math >}}

### 2. Substitute s(t) and u(t) into the definition of $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \tilde{\beta}\left[u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}\right]
$$
{{< /math >}}

Now, substitute the definition of $\tilde{\beta} = \frac{\beta}{\gamma - \beta}$ and $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \frac{\beta}{\gamma - \beta}\left[\frac{\alpha^{(k)}}{\beta} + \left(u_0 - \frac{\alpha^{(k)}}{\beta}\right)e^{-\beta\tau}\right]
$$
{{< /math >}}

Let's simplify the terms multiplied by $\tilde{\beta}$:

{{< math >}}
$$
\tilde{\beta}u_\infty^{(k)} = \frac{\beta}{\gamma - \beta} \cdot \frac{\alpha^{(k)}}{\beta} = \frac{\alpha^{(k)}}{\gamma - \beta}
$$
{{< /math >}}

{{< math >}}
$$
\tilde{\beta}(u_0 - u_\infty^{(k)}) = \frac{\beta}{\gamma - \beta}\left(u_0 - \frac{\alpha^{(k)}}{\beta}\right) = \frac{\beta}{\gamma - \beta}\left(\frac{\beta u_0 - \alpha^{(k)}}{\beta}\right) = \frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}
$$
{{< /math >}}

Substitute these back into the expression for $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \left[\frac{\alpha^{(k)}}{\gamma - \beta} + \left(\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}\right)e^{-\beta\tau}\right]
$$
{{< /math >}}

### 3. Cancel the $e^{-\beta\tau}$ terms:

Observe the coefficients of $e^{-\beta\tau}$:

{{< math >}}
$$
-\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) - \left(\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}\right)
$$
{{< /math >}}

Since $\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta} = -\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}$, this simplifies to:

{{< math >}}
$$
-\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) - \left(-\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) = -\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) + \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) = 0
$$
{{< /math >}}

The terms containing $e^{-\beta\tau}$ indeed cancel out, as intended by the choice of $\tilde{\beta}$.

### 4. Simplify $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left(\frac{\alpha^{(k)}}{\gamma} - \frac{\alpha^{(k)}}{\gamma - \beta}\right) + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

Let's simplify the constant term:

{{< math >}}
$$
\frac{\alpha^{(k)}}{\gamma} - \frac{\alpha^{(k)}}{\gamma - \beta} = \alpha^{(k)}\left(\frac{1}{\gamma} - \frac{1}{\gamma - \beta}\right) = \alpha^{(k)}\left(\frac{(\gamma - \beta) - \gamma}{\gamma(\gamma - \beta)}\right) = \alpha^{(k)}\left(\frac{-\beta}{\gamma(\gamma - \beta)}\right) = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

So, we have:

{{< math >}}
$$
\tilde{s}(t) = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

### 5. Calculate $\tilde{s}_0$ and $\tilde{s}_\infty$:

Given {{< math >}} $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$ {{< /math >}}, {{< math >}} $ s_\infty = \frac{\beta u_\infty^{(k)}}{\gamma} = \frac{\alpha^{(k)}}{\gamma}
$ {{< /math >}} and the definition of $\tilde{s}_\infty$ and $\tilde{\beta}$, we derive below. You may also consider view $\tilde{s}_\infty$  as Steady-state of $\tilde{s}(t)$  ($\tau \to \infty$, $e^{-\gamma\tau} \to 0$ ).

{{< math >}}
$$
\tilde{s}_\infty = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} (= \lim_{\tau \to \infty} \tilde{s}(t) )
$$
{{< /math >}}

**$\tilde{s}_0$ (Initial value of $\tilde{s}(t)$ at $\tau = 0$):**

{{< math >}}
$$
\tilde{s}_0 = s(t_0) - \tilde{\beta}u(t_0) = s_0 - \tilde{\beta}u_0 = s_0 - \frac{\beta}{\gamma - \beta}u_0
$$
{{< /math >}}

### 6. Verify the target relationship: $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$

**Left-Hand Side (LHS):**

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}\right) - \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right)
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

**Right-Hand Side (RHS):**

First, calculate $(\tilde{s}_0 - \tilde{s}_\infty)$:

{{< math >}}
$$
(\tilde{s}_0 - \tilde{s}_\infty) = \left(s_0 - \frac{\beta u_0}{\gamma - \beta}\right) - \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right) = s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

Now, multiply by $e^{-\gamma\tau}$:

{{< math >}}
$$
(\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau} = \left(s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right)e^{-\gamma\tau}
$$
{{< /math >}}

Finally, let's compare the coefficients of $e^{-\gamma\tau}$ from both the LHS and RHS. We need to confirm that:

{{< math >}}
$$
s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta} \text{ equals } s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

Let's simplify the terms without $s_0$:

**From LHS coefficient:**

{{< math >}}
$$
-\frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta} = -\frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)}}{\gamma - \beta} - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{1}{\gamma - \beta} - \frac{1}{\gamma}\right) - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{\gamma - (\gamma - \beta)}{\gamma(\gamma - \beta)}\right) - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{\beta}{\gamma(\gamma - \beta)}\right) - \frac{\beta u_0}{\gamma - \beta} = \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} - \frac{\beta u_0}{\gamma - \beta}
$$
{{< /math >}}

**From RHS coefficient:**

{{< math >}}
$$
-\frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

As you can see, the simplified coefficients from both sides are indeed identical.

## Conclusion

By introducing the carefully chosen auxiliary variable $\tilde{\beta} = \frac{\beta}{\gamma - \beta}$, the transformation $\tilde{s}(t) = s(t) - \tilde{\beta}u(t)$ strategically eliminates the $e^{-\beta\tau}$ dependence from the solution of s(t). This simplification reveals a clean, exponential decay relationship for $\tilde{s}(t)$ towards its new steady-state value $\tilde{s}_\infty$, as expressed by:

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}
$$
{{< /math >}}

This shows how a more complex system can be transformed into a simpler, recognizable form by defining appropriate auxiliary variables.

# B) Derivation: $\tau$ in equation (10, 11) in [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3)

Above smoothly gives Inverse Relationship for $\tau$ in equation (10) in the [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3). Although equation (11) looks not obvious, it actually requires even simpler reasoning.

## Derive {{< math >}} $\tau$ {{< /math >}} in equation (11)
{{< math >}}
$$
\tau = -\frac{1}{\beta}\ln\left(\frac{u - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)
$$
{{< /math >}}

This expression makes sense given the solution for {{< math >}} $u(t)$ {{< /math >}} can be rewritten as below by substituting {{< math >}} $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$ {{< /math >}}:

{{< math >}}
$$
u(t) = u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}
$$
{{< /math >}}

If we solve this for {{< math >}} $\tau$ {{< /math >}}, we get:

* Rearrange: {{< math >}} $u(t) - u_\infty^{(k)} = (u_0 - u_\infty^{(k)})e^{-\beta\tau}$ {{< /math >}}

* Divide both sides by {{< math >}} $(u_0 - u_\infty^{(k)})$ {{< /math >}}: {{< math >}} $\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}} = e^{-\beta\tau}$ {{< /math >}}

* Take natural log: {{< math >}} $\ln\left(\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right) = -\beta\tau$ {{< /math >}}

* Solve for {{< math >}} $\tau$ {{< /math >}}: {{< math >}} $\tau = -\frac{1}{\beta}\ln\left(\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)$ {{< /math >}}

# C) Explicit Formula is Less Accurate than Exact Optimization in scVelo

The difference in accuracy between the explicit formula and exact optimization comes down to **model assumptions versus reality**. Let me explain:

## Why the Explicit Formula is Less Accurate

### 1. Single-Gene Assumption

The explicit formula:

{{< math >}}
$$\tau_i = -\frac{1}{\beta}\ln\left(\frac{u_i - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)$$
{{< /math >}}

- **Uses only unspliced RNA** {{< math >}} $(u)$ {{< /math >}} from one gene
- **Ignores spliced RNA** {{< math >}} $(s)$ {{< /math >}} information completely
- **Assumes perfect adherence** to the theoretical kinetic model for that single gene

### 2. No Cross-Gene Integration

- Each gene would give a **different time estimate** for the same cell
- The formula doesn't reconcile these conflicting estimates
- **Gene-specific noise** and measurement errors aren't averaged out

### 3. Model Violations

Real cells don't perfectly follow the kinetic equations because:

- **Transcriptional bursting** creates noise
- **Cell-to-cell variability** in kinetic parameters
- **Measurement noise** in RNA counts
- **Model approximations** (e.g., constant degradation rates)

## Why Exact Optimization is More Accurate

### 1. Multi-Gene Integration

The optimization objective:

{{< math >}}
$$\tau_i^* = \arg\min_\tau \sum_j \left[ (u_{i,j} - u_j(\tau))^2 + (s_{i,j} - s_j(\tau))^2 \right]$$
{{< /math >}}

- **Uses both unspliced AND spliced** RNA information
- **Integrates evidence across ALL genes** simultaneously
- **Finds the time that best explains the entire transcriptome**

### 2. Noise Averaging

- Gene-specific measurement errors **cancel out** when averaged
- **Outlier genes** have less impact on the final time estimate
- **Statistical power increases** with more observations

### 3. Comprehensive Model Fitting

- Accounts for **both RNA species** {{< math >}} $u$ {{< /math >}} and {{< math >}} $s$ {{< /math >}} jointly
- **Balances competing evidence** from different genes
- **Robust to individual gene model violations**

## Computational Trade-off

**Explicit Formula**: {{< math >}} $O(1)$ {{< /math >}} per cell
- Fast but uses limited information

**Exact Optimization**: {{< math >}} $O(n_{iter} \times n_{genes})$ {{< /math >}} per cell  
- Slower but uses all available information

## Mathematical Perspective

### Information Content Comparison

**Explicit Formula Information**:
{{< math >}}
$$I_{explicit} = f(u_{i,j}) \text{ for single gene } j$$
{{< /math >}}

**Exact Optimization Information**:
{{< math >}}
$$I_{exact} = f(\{u_{i,j}, s_{i,j}\}_{j=1}^{N_{genes}})$$
{{< /math >}}

Where {{< math >}} $I_{exact} >> I_{explicit}$ {{< /math >}} in terms of information content.

### Error Propagation

For the explicit formula, measurement error in a single gene directly affects the time estimate:

{{< math >}}
$$\sigma_{\tau}^2 \propto \frac{\sigma_u^2}{(u_i - u_\infty)^2}$$
{{< /math >}}

For exact optimization, errors are averaged across genes:

{{< math >}}
$$\sigma_{\tau}^2 \propto \frac{1}{N_{genes}} \sum_j \frac{\sigma_{u,j}^2 + \sigma_{s,j}^2}{(u_{i,j} - u_j(\tau))^2 + (s_{i,j} - s_j(\tau))^2}$$
{{< /math >}}

## Biological Analogy

Think of it like **estimating someone's age**:

**Explicit Formula** = Looking at just their hair color
- Quick, but hair color alone is unreliable

**Exact Optimization** = Looking at hair, skin, posture, clothing style, speech patterns, etc.
- Takes more time, but gives a much better estimate by integrating multiple sources of evidence

## Why the Two-Stage Strategy Works

The two-stage approach in scVelo is brilliant because:

1. **Stage 1 (Explicit)**: Gets "close enough" quickly when parameters are changing rapidly
2. **Stage 2 (Exact)**: Refines to high accuracy when parameters have stabilized

### Convergence Analysis

Early iterations: {{< math >}} $|\theta^{(t)} - \theta^{(t-1)}|$ {{< /math >}} is large
- Approximate times sufficient since parameters will change significantly anyway

Later iterations: {{< math >}} $|\theta^{(t)} - \theta^{(t-1)}| \approx 0$ {{< /math >}}
- Exact times needed for final parameter refinement

### Computational Efficiency

Total computational cost:

{{< math >}}
$$\text{Cost} = T_{approx} \times O(1) + T_{exact} \times O(n_{iter} \times n_{genes})$$
{{< /math >}}

Where {{< math >}} $T_{approx} >> T_{exact}$ {{< /math >}}, making the overall algorithm much faster than using exact optimization throughout.

## Key Insight

The explicit formula serves as an excellent **initialization and approximation** tool, while exact optimization provides the **precision needed for final convergence**. This hybrid approach captures the best of both worlds: computational efficiency and high accuracy.





# D) Cosine Similarity with Variance-Stabilized Transformation

### **Example Setup**
We have two cells \( i \) and \( j \), with gene expression vectors:

{{< math >}} 
$$ s_i = \begin{bmatrix} 2.0 \\ 3.5 \\ 1.0 \end{bmatrix}, \quad s_j = \begin{bmatrix} 3.5 \\ 2.5 \\ 1.5 \end{bmatrix} $$ 
{{< /math >}}

We also have the velocity vector for cell \( i \):

{{< math >}} 
$$ v_i = \begin{bmatrix} 1.2 \\ -0.5 \\ 0.3 \end{bmatrix} $$ 
{{< /math >}}

First, we compute the difference between the expression values of the two cells:

{{< math >}} 
$$ \delta_{ij} = s_j - s_i = \begin{bmatrix} 3.5 - 2.0 \\ 2.5 - 3.5 \\ 1.5 - 1.0 \end{bmatrix} = \begin{bmatrix} 1.5 \\ -1.0 \\ 0.5 \end{bmatrix} $$ 
{{< /math >}}

---

## **1. Standard Cosine Similarity**
The standard cosine similarity is given by:

{{< math >}} 
$$ \pi_{ij} = \cos(\delta_{ij}, v_i) = \frac{\delta_{ij}^T v_i}{\|\delta_{ij}\| \|v_i\|} $$ 
{{< /math >}}

**Step 1: Compute the dot product**
{{< math >}} 
$$ \delta_{ij}^T v_i = (1.5)(1.2) + (-1.0)(-0.5) + (0.5)(0.3) = 1.8 + 0.5 + 0.15 = 2.45 $$ 
{{< /math >}}

**Step 2: Compute the magnitudes**
{{< math >}} 
$$ \|\delta_{ij}\| = \sqrt{(1.5)^2 + (-1.0)^2 + (0.5)^2} = \sqrt{2.25 + 1.00 + 0.25} = \sqrt{3.5} \approx 1.87 $$ 
{{< /math >}}

{{< math >}} 
$$ \|v_i\| = \sqrt{(1.2)^2 + (-0.5)^2 + (0.3)^2} = \sqrt{1.44 + 0.25 + 0.09} = \sqrt{1.78} \approx 1.34 $$ 
{{< /math >}}

**Step 3: Compute cosine similarity**
{{< math >}} 
$$ \pi_{ij} = \frac{2.45}{(1.87)(1.34)} = \frac{2.45}{2.51} \approx 0.98 $$ 
{{< /math >}}

The vectors are strongly aligned in direction.

---

## **2. Variance-Stabilized Transformation**
Instead of using raw values, we apply:

{{< math >}} 
$$ \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|}, \quad \text{sign}(v_i) \sqrt{|v_i|} $$ 
{{< /math >}}

**Step 1: Transform \( \delta_{ij} \)**

{{< math >}} 
$$ \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|} = \begin{bmatrix} \text{sign}(1.5) \sqrt{1.5} \\ \text{sign}(-1.0) \sqrt{1.0} \\ \text{sign}(0.5) \sqrt{0.5} \end{bmatrix} = \begin{bmatrix} 1.22 \\ -1.00 \\ 0.71 \end{bmatrix} $$ 
{{< /math >}}

**Step 2: Transform \( v_i \)**

{{< math >}} 
$$ \text{sign}(v_i) \sqrt{|v_i|} = \begin{bmatrix} \text{sign}(1.2) \sqrt{1.2} \\ \text{sign}(-0.5) \sqrt{0.5} \\ \text{sign}(0.3) \sqrt{0.3} \end{bmatrix} = \begin{bmatrix} 1.10 \\ -0.71 \\ 0.55 \end{bmatrix} $$ 
{{< /math >}}

**Step 3: Compute transformed cosine similarity**
{{< math >}} 
$$ \pi_{ij}^{\text{stabilized}} = \frac{\left( 1.22, -1.00, 0.71 \right) \cdot \left( 1.10, -0.71, 0.55 \right)}{\| \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|}\| \|\text{sign}(v_i) \sqrt{|v_i|}\|} $$ 
{{< /math >}}

**Dot product**
{{< math >}} 
$$ 1.22(1.10) + (-1.00)(-0.71) + (0.71)(0.55) = 1.342 + 0.71 + 0.391 = 2.443 $$ 
{{< /math >}}

**Compute magnitudes**
{{< math >}} 
$$ \sqrt{(1.22)^2 + (-1.00)^2 + (0.71)^2} = \sqrt{1.488 + 1.000 + 0.504} = \sqrt{2.992} \approx 1.73 $$ 
{{< /math >}}

{{< math >}} 
$$ \sqrt{(1.10)^2 + (-0.71)^2 + (0.55)^2} = \sqrt{1.210 + 0.504 + 0.302} = \sqrt{2.016} \approx 1.42 $$ 
{{< /math >}}

**Final computation**
{{< math >}} 
$$ \pi_{ij}^{\text{stabilized}} = \frac{2.443}{(1.73)(1.42)} = \frac{2.443}{2.46} \approx 0.99 $$ 
{{< /math >}}

---

# E) Exponential Kernel Transformation for Transition Probabilities

### **1. Exponential Transformation of Cosine Correlation**
The raw cosine similarity \( \pi_{ij} \) measures directional alignment between velocity vectors, but to obtain meaningful **transition probabilities**, we apply an **exponential kernel**:

{{< math >}} 
$$ \tilde{\pi}_{ij} = \frac{1}{z_i} \exp \left(\frac{\pi_{ij}}{\sigma^2} \right) $$ 
{{< /math >}}

where:
- {{< math >}} $ \pi_{ij} $ {{< /math >}} is the **cosine similarity** between cell {{< math >}} $ i $ {{< /math >}} and {{< math >}} $ j $ {{< /math >}}.
- {{< math >}} $ \sigma $ {{< /math >}} is the **kernel width parameter**, which controls the spread of probabilities.
- {{< math >}} $ z_i $ {{< /math >}} is a **row normalization factor**, ensuring probabilities sum to 1:

{{< math >}} 
$$ z_i = \sum_j \exp \left(\frac{\pi_{ij}}{\sigma^2} \right) $$ 
{{< /math >}}

---

### **2. Why Exponential Transformation?**
The **exponential function** ensures:
✅ **Amplification of strong correlations**: Large values of {{< math >}} $ \pi_{ij} $ {{< /math >}} get boosted.  
✅ **Suppression of weak correlations**: Small or negative values decay rapidly, reducing their transition influence.  
✅ **Nonlinear scaling**: Instead of treating all similarity scores equally, this transformation **sharpens distinctions** between stronger and weaker connections.  

Think of it as a **softmax-like weighting** that enhances **directional flow probability** based on cosine similarity.

---

### **3. Row Normalization Ensures Probabilities Sum to 1**
The **normalization factor** {{< math >}} $ z_i $ {{< /math >}} guarantees valid probability distributions:

{{< math >}} 
$$ \sum_j \tilde{\pi}_{ij} = 1 $$ 
{{< /math >}}

Without this step, raw exponentiated values might grow arbitrarily large, leading to **improper probability distributions**.

---

### **4. Role of Kernel Width \( \sigma \)**
- **Large \( \sigma \)** → Softens probability differences, making transitions smoother.
- **Small \( \sigma \)** → Sharper transitions, favoring **strong** directional alignments.

Choosing an appropriate {{< math >}} $ \sigma $ {{< /math >}} ensures **biological interpretability** in applications like RNA velocity and cell fate prediction.

Would you like a **numerical example** demonstrating how this transformation works with actual cosine similarity values? 🚀
