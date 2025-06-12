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


# F) How scVelo Identifies the Starting Points of Differentiation Using Markov Chain Analysis

## **1. What Are Root Cells in Differentiation?**
In single-cell RNA velocity analysis, root cells are the **initial states** in a developmental trajectory. These are cells that:  
✅ Have **low latent time** (early progenitors).  
✅ Are **unlikely** to have transitioned from other states.  
✅ Serve as the **starting points** for differentiation pathways.  

Identifying these root cells is crucial for understanding **how differentiation originates**.

---

## **2. Stationary State Equation in RNA Velocity**
The root cells are found by solving the **stationary distribution equation** for the **transition probability matrix** {{< math >}} $ \tilde{\pi} $ {{< /math >}}:

{{< math >}} 
$$ \mu^* = \mu^* \tilde{\pi}^T $$ 
{{< /math >}}

where:
- {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} is the **transpose of the transition probability matrix** that encodes RNA velocity-driven transitions.
- {{< math >}} $ \mu^* $ {{< /math >}} represents the **stationary distribution**, which describes the **long-term probabilities** of cells staying in each state.
- The **left eigenvectors** of {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} corresponding to the **eigenvalue of 1** define these **stationary states**.

---

## **3. Why Does This Identify Root Cells?**

### **a. Probability Evolution in a Markov Process**
In a discrete-time **Markov chain**, the state of the system at time {{< math >}} $ t $ {{< /math >}} is represented by a **probability distribution** {{< math >}} $ \mu_t $ {{< /math >}}, which evolves according to the **transition probability matrix** {{< math >}} $ P $ {{< /math >}}:

{{< math >}} 
$$ \mu_{t+1} = \mu_t P $$ 
{{< /math >}}

where:
- {{< math >}} $ \mu_t $ {{< /math >}} is a row vector of probabilities over all possible states at time {{< math >}} $ t $ {{< /math >}}.
- {{< math >}} $ P $ {{< /math >}} is an {{< math >}} $ n \times n $ {{< /math >}} **transition matrix**, where {{< math >}} $ P_{ij} $ {{< /math >}} represents the probability of transitioning from state {{< math >}} $ i $ {{< /math >}} to state {{< math >}} $ j $ {{< /math >}}.

Each step updates the probability distribution, progressively altering the likelihood of being in each state.

---

### **b. Convergence to a Stable Distribution**
As the Markov chain progresses, repeated multiplications by {{< math >}} $ P $ {{< /math >}} lead the probability distribution toward a **stable equilibrium**:

{{< math >}} 
$$ \mu_{t+k} = \mu_t P^k, \quad \text{as } k \to \infty $$ 
{{< /math >}}

After many steps, the probabilities stabilize, meaning **the system reaches a final long-term behavior** where further transitions do not significantly change the probabilities.

---

### **c. Definition of the Stationary Distribution**
The **stationary distribution** {{< math >}} $ \mu^* $ {{< /math >}} is a probability vector that remains **unchanged** under transitions:

{{< math >}} 
$$ \mu^* = \mu^* P $$ 
{{< /math >}}

This implies that if a system starts in {{< math >}} $ \mu^* $ {{< /math >}}, applying the transition matrix does **not** alter its probabilities—it stays in equilibrium.

---

### **d. Why Is the Stationary Distribution Important?**
✅ **Predicts long-term state occupancy**—the fraction of time spent in each state.  
✅ **Defines steady-state probabilities** in applications like RNA velocity analysis.  
✅ **Identifies stable states (e.g., progenitor/root cells in differentiation)** when applied to biological modeling.


## **e. Example: Computing Stationary States Using Eigenvectors**
### **Step 1: Define the Transition Probability Matrix**
Consider a **3-state system** with transition probabilities:

{{< math >}} 
$$ \tilde{\pi} = \begin{bmatrix} 
0.8 & 0.1 & 0.1 \\ 
0.2 & 0.7 & 0.1 \\ 
0.1 & 0.2 & 0.7 
\end{bmatrix} $$ 
{{< /math >}}

This matrix describes **how probabilities shift** between states.

---

### **Step 2: Find the Left Eigenvector for Eigenvalue 1**
We solve:

{{< math >}} 
$$ v^T \tilde{\pi}^T = v^T $$ 
{{< /math >}}

This means we need the **eigenvector corresponding to eigenvalue \( \lambda = 1 \)**.

Computing the eigenvectors of {{< math >}} $ \tilde{\pi}^T $ {{< /math >}}, we obtain:

{{< math >}} 
$$ v = \begin{bmatrix} 0.45 \\ 0.35 \\ 0.20 \end{bmatrix} $$ 
{{< /math >}}

This **left eigenvector** represents the **stationary distribution**, meaning:
✅ **State 1 holds 45% of the long-term probability**.  
✅ **State 2 holds 35%**.  
✅ **State 3 holds 20%**.  

These **steady-state values** describe **how the system behaves in the long run**, identifying the states **least likely to transition further**.

---

### **4. Biological Interpretation in RNA Velocity**
In **single-cell analysis**:
✅ **Cells corresponding to the stationary distribution** act as **root cells**, initiating differentiation.  
✅ This eigenvector-driven method ensures an **objective, data-driven way** to define **progenitor states**.  
✅ **Eigenvalue 1 captures states that remain stable over differentiation**, meaning they are **starting points** in biological development.



# G) Is p in eq(16) dependent on gene or cell?
In scVelo's global time normalization, we have:

{{< math >}} $$t_{i;o} = Q^p_g \left( t_{i;g} - t_{o;g} \right)$$ {{< /math >}}

Where:
- {{< math >}} $t_{i;o}$ {{< /math >}}: normalized time for cell {{< math >}} $i$ {{< /math >}} relative to root cell {{< math >}} $o$ {{< /math >}}
- {{< math >}} $Q^p_g$ {{< /math >}}: {{< math >}} $p$ {{< /math >}}-quantile computed for gene {{< math >}} $g$ {{< /math >}}
- {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}: time shift between cell {{< math >}} $i$ {{< /math >}} and root cell {{< math >}} $o$ {{< /math >}} for gene {{< math >}} $g$ {{< /math >}}

## Is p Gene-Specific?

**❌ No** — {{< math >}} $p$ {{< /math >}} is **not gene-specific**.

### Explanation

{{< math >}} $Q^p_g$ {{< /math >}} means: for **each gene** {{< math >}} $g$ {{< /math >}}, you compute the **{{< math >}} $p$ {{< /math >}}-quantile** over the time shifts {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}, across all cells {{< math >}} $i$ {{< /math >}}.

But the value of **{{< math >}} $p$ {{< /math >}}** itself is the **same for all genes** — it is **not gene-specific**.

So even though the quantile is applied **per gene** (i.e., each gene gets its own {{< math >}} $Q^p_g$ {{< /math >}}), the **{{< math >}} $p$ {{< /math >}} used in all those calculations is shared**.

### Mathematical Illustration

For a fixed value {{< math >}} $p = 0.5$ {{< /math >}} (median):

{{< math >}} $$Q^{0.5}_{gene1} = \text{median}(\{t_{i;gene1} - t_{o;gene1}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{gene2} = \text{median}(\{t_{i;gene2} - t_{o;gene2}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$\vdots$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{geneM} = \text{median}(\{t_{i;geneM} - t_{o;geneM}\}_{i=1}^N)$$ {{< /math >}}

The **same** {{< math >}} $p = 0.5$ {{< /math >}} is used for all genes, but each gene gets its own quantile value.

## Is p Cell-Specific?

**❌ Also no** — {{< math >}} $p$ {{< /math >}} is not adjusted for each cell. It is only used to **transform the gene-wise temporal shifts** into a consistent time estimate **for each cell**, **relative to the root cell** {{< math >}} $o$ {{< /math >}}.

### Explanation

The only place where indexing over cells happens is in computing the latent time {{< math >}} $t_{i;o}$ {{< /math >}}, but the **value of {{< math >}} $p$ {{< /math >}} stays fixed** for all cells.

### Mathematical Process

1. **Fixed {{< math >}} $p$ {{< /math >}} across all computations**: {{< math >}} $p = p_{\text{fixed}}$ {{< /math >}}

2. **Gene-wise quantile computation**:
   {{< math >}} $$Q^{p_{\text{fixed}}}_g = \text{quantile}_{p_{\text{fixed}}}(\{t_{i;g} - t_{o;g}\}_{i=1}^N)$$ {{< /math >}}

3. **Cell-specific time assignment**:
   {{< math >}} $$t_{1;o} = Q^{p_{\text{fixed}}}_g \left( t_{1;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$t_{2;o} = Q^{p_{\text{fixed}}}_g \left( t_{2;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$\vdots$$ {{< /math >}}
   {{< math >}} $$t_{N;o} = Q^{p_{\text{fixed}}}_g \left( t_{N;g} - t_{o;g} \right)$$ {{< /math >}}

## Is p Root-Specific?

**✅ Yes** — scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness.

### Optimization Process

For each potential root cell {{< math >}} $o$ {{< /math >}}, scVelo optimizes:

{{< math >}} $$p^*_o = \arg\min_p \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

Where the loss function might be:

{{< math >}} $$\mathcal{L}_{\text{smoothness}}(p, o) = \sum_{i,j \in \text{neighbors}} \left| t_{i;o}(p) - t_{j;o}(p) \right|^2$$ {{< /math >}}

### Root Selection

The final root and {{< math >}} $p$ {{< /math >}} are chosen as:

{{< math >}} $$o^*, p^* = \arg\min_{o,p} \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

## Summary Table

| Aspect | Is {{< math >}} $p$ {{< /math >}} dependent? | Explanation |
|--------|-----|-------------|
| **Per gene** {{< math >}} $g$ {{< /math >}} | ❌ No | The quantile {{< math >}} $Q^p$ {{< /math >}} is computed per gene, but **{{< math >}} $p$ {{< /math >}} itself is fixed** |
| **Per cell** {{< math >}} $i$ {{< /math >}} | ❌ No | All cells use the same {{< math >}} $p$ {{< /math >}}; only their computed values differ |
| **Per root** {{< math >}} $o$ {{< /math >}} | ✅ Yes | scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness |


## Optimization Step

Once **scVelo** computes all {{< math >}} $ t_{i,o} $ {{< /math >}} for different candidate root cells {{< math >}} $ o $ {{< /math >}}, it looks for the root and quantile {{< math >}} $ p $ {{< /math >}} that maximize correlation between the **latent time series**  
{{< math >}} 
$$ \{ t_{i,o} \} $$ 
{{< /math >}}  
and its **convolution over a local neighborhood graph** (e.g., KNN graph of cells).

### **Mathematical Formulation**
The optimization step seeks:

{{< math >}} 
$$ \arg\max_{o,p} \text{corr} \left( t_{i,o}, \sum_{j \in N(i)} w_{ij} t_{j,o} \right) $$ 
{{< /math >}}

where:
- {{< math >}} $ N(i) $ {{< /math >}} is the **neighborhood of cell** {{< math >}} $ i $ {{< /math >}}.
- {{< math >}} $ w_{ij} $ {{< /math >}} represents **graph weights**, determining connectivity between cells.



# H) Project RNA Velocities into an Embedding (e.g., UMAP)

In scVelo, you have:

- **High-dimensional RNA velocities** (in gene expression space)
- **Low-dimensional embeddings** (e.g., UMAP, t-SNE)

The challenge is to project velocities into the embedding so that arrows in UMAP space reflect true dynamics.

## 📐 Key Variables

Let {{< math >}} $\vec{s}_i$ {{< /math >}} be the embedding position of cell {{< math >}} $i$ {{< /math >}} (e.g., in 2D UMAP space)

Let {{< math >}} $\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|}$ {{< /math >}} be the normalized direction vector from cell {{< math >}} $i$ {{< /math >}} to cell {{< math >}} $j$ {{< /math >}}

Let {{< math >}} $\pi_{ij}$ {{< /math >}} be the transition probability from cell {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} in the velocity-derived transition matrix

Let {{< math >}} $\vec{\nu}_i$ {{< /math >}} be the velocity vector in the embedding space for cell {{< math >}} $i$ {{< /math >}}

## 📊 Formula: Projected Velocity Vector

The projected velocity for cell {{< math >}} $i$ {{< /math >}} in the embedding space is:

{{< math >}} $$\vec{\nu}_i = E_{\vec{\pi}_i}[\vec{\delta}_{ij}] = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij} - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

Or in compact form:

{{< math >}} $$\vec{\nu}_i = \vec{\delta}_{i \cdot}^{\top} \vec{\pi}_i - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

## 🧠 Intuition

The **first term** is the expected movement direction in embedding space, weighted by the transition probabilities derived from RNA velocity.

The **second term**, {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}}, corrects for non-uniform density in the embedding. Without it, areas with more cells would bias the velocity field.

This ensures that the resulting velocity arrows:

- Follow the dynamics inferred from spliced/unspliced RNA
- Are not artifacts of cell density in UMAP
- Represent true directionality in cellular state transitions

## 🖼️ In Practice

This velocity vector {{< math >}} $\vec{\nu}_i$ {{< /math >}} is:

- Computed for each cell
- Overlaid as arrows on UMAP or t-SNE plots using `scv.pl.velocity_embedding()`

## ✅ Summary

| Term | Meaning |
|------|---------|
| {{< math >}} $\vec{s}_i$ {{< /math >}} | Position of cell {{< math >}} $i$ {{< /math >}} in embedding (e.g., UMAP) |
| {{< math >}} $\vec{\delta}_{ij}$ {{< /math >}} | Normalized direction vector from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\pi_{ij}$ {{< /math >}} | Velocity-based transition probability from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\vec{\nu}_i$ {{< /math >}} | Estimated velocity vector in embedding for cell {{< math >}} $i$ {{< /math >}} |
| {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}} | Density correction term |

## Mathematical Breakdown

### Step 1: Compute Direction Vectors
For each cell {{< math >}} $i$ {{< /math >}} and all its neighbors {{< math >}} $j$ {{< /math >}}:

{{< math >}} $$\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|_2}$$ {{< /math >}}

### Step 2: Get Transition Probabilities
From the velocity graph computed by scVelo:

{{< math >}} $$\pi_{ij} = P(\text{cell } i \rightarrow \text{cell } j | \text{RNA velocity})$$ {{< /math >}}

### Step 3: Weighted Average Direction
{{< math >}} $$\text{Expected direction} = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij}$$ {{< /math >}}

### Step 4: Density Correction
{{< math >}} $$\text{Uniform correction} = \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

### Step 5: Final Velocity
{{< math >}} $$\vec{\nu}_i = \text{Expected direction} - \text{Uniform correction}$$ {{< /math >}}

## Alternative Formulations

### Kernel-Based Approach
Some implementations use a kernel-weighted version:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j K(\vec{s}_i, \vec{s}_j) \pi_{ij} \vec{\delta}_{ij}}{\sum_j K(\vec{s}_i, \vec{s}_j)}$$ {{< /math >}}

where {{< math >}} $K(\vec{s}_i, \vec{s}_j)$ {{< /math >}} is a distance-based kernel (e.g., Gaussian).

### Confidence Weighting
Incorporate velocity confidence scores:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j c_{ij} \pi_{ij} \vec{\delta}_{ij}}{\sum_j c_{ij}}$$ {{< /math >}}

where {{< math >}} $c_{ij}$ {{< /math >}} represents the confidence in the transition from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}}.

## Implementation Notes

The projected velocities are typically stored in:
- `adata.obsm['velocity_umap']` for UMAP projections
- `adata.obsm['velocity_tsne']` for t-SNE projections

These can be accessed and visualized using scVelo's plotting functions to show the directional flow of cellular development in low-dimensional space.

# I) The Reconstructability Score (r)

## Mathematical Definition

{{< math >}} $$r = \text{median}_i \, \text{corr}(\pi_i, \pi_i')$$ {{< /math >}}

Where:

- {{< math >}} $\pi_i$ {{< /math >}}: The row vector of the full velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities to all other cells.

- {{< math >}} $\pi_i'$ {{< /math >}}: The row vector of the reduced velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities based only on the selected genes.

- {{< math >}} $\text{corr}(\pi_i, \pi_i')$ {{< /math >}}: This calculates the correlation between the full outgoing transition probabilities for cell {{< math >}} $i$ {{< /math >}} and the outgoing transition probabilities derived from the subset of genes for cell {{< math >}} $i$ {{< /math >}}. A high correlation means that the subset of genes largely reproduces the directional information provided by all genes for that particular cell.

- {{< math >}} $\text{median}_i$ {{< /math >}}: The score is then computed as the median of these per-cell correlations across all cells {{< math >}} $i$ {{< /math >}}. Taking the median makes the score robust to outliers (e.g., a few cells where the subset of genes might perform poorly due to local noise).

## What the Reconstructability Score Tells You

The reconstructability score quantifies how well a specific subset of genes can "reconstruct" or recapitulate the overall RNA velocity dynamics inferred from the full set of genes.

### High Score (close to 1)
If {{< math >}} $r$ {{< /math >}} is high (e.g., {{< math >}} $> 0.8$ {{< /math >}}), it suggests that the selected subset of genes (e.g., top likelihood genes) are indeed the primary drivers of the observed cellular transitions. Their dynamics alone are largely sufficient to explain the overall directionality and speed of differentiation. This validates their "driver" status and indicates that focusing on these genes provides a good summary of the system's dynamics.

### Low Score
If {{< math >}} $r$ {{< /math >}} is low, it means that the selected subset of genes does not adequately capture the full dynamic picture. The overall cellular transitions are likely influenced by a broader set of genes not included in your selection, or the selected genes might not be truly representative of the overall dynamics.

## Mathematical Components Breakdown

### Transition Probability Vectors

For each cell {{< math >}} $i$ {{< /math >}}, the transition probabilities are derived from the velocity graph:

{{< math >}} $$\pi_i = [\pi_{i,1}, \pi_{i,2}, \ldots, \pi_{i,N}]$$ {{< /math >}}

where {{< math >}} $\pi_{i,j}$ {{< /math >}} represents the probability of cell {{< math >}} $i$ {{< /math >}} transitioning to cell {{< math >}} $j$ {{< /math >}}, and {{< math >}} $N$ {{< /math >}} is the total number of cells.

### Correlation Calculation

The correlation between full and reduced transition probabilities:

{{< math >}} $$\text{corr}(\pi_i, \pi_i') = \frac{\text{cov}(\pi_i, \pi_i')}{\sigma_{\pi_i} \sigma_{\pi_i'}}$$ {{< /math >}}

where:
- {{< math >}} $\text{cov}(\pi_i, \pi_i')$ {{< /math >}} is the covariance
- {{< math >}} $\sigma_{\pi_i}$ {{< /math >}} and {{< math >}} $\sigma_{\pi_i'}$ {{< /math >}} are the standard deviations

### Median Aggregation

The final score uses median aggregation for robustness:

{{< math >}} $$r = \text{median}\{c_1, c_2, \ldots, c_N\}$$ {{< /math >}}

where {{< math >}} $c_i = \text{corr}(\pi_i, \pi_i')$ {{< /math >}} for each cell {{< math >}} $i$ {{< /math >}}.

## Use Cases and Significance

### Validation of Driver Gene Selection
It's commonly used to validate that the "top likelihood genes" (which scVelo identifies as dynamically strong) are indeed the key players shaping the global trajectory.

### Identifying Essential Gene Sets
You can test specific gene modules or pathways to see if they are collectively responsible for a particular fate decision or developmental progression.

### Dimensionality Reduction/Summarization
If a small subset of genes yields a high reconstructability score, it suggests that these genes form a powerful, low-dimensional representation of the system's dynamics.

### Biological Interpretation
High reconstructability from a specific gene set points to their significant biological role in driving the observed cellular changes.

