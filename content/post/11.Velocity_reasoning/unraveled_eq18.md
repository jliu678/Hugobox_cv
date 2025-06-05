# Equation 18 from [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492)

## 🧬 1. Joint state and generating function

We define the state of a cell as:

{{< math >}} 
$$x = (x_u, x_s)$$ 
{{< /math >}}

Where:
- {{< math >}} $x_u$ {{< /math >}}: unspliced mRNA count
- {{< math >}} $x_s$ {{< /math >}}: spliced mRNA count

We now write the generating function of their joint distribution:

{{< math >}} 
$$G(u_u, u_s, t) = \sum_x P(x,t)(u_u + 1)^{x_u}(u_s + 1)^{x_s}$$ 
{{< /math >}}

This is a modified bivariate probability generating function, where the "+1" shift is standard in certain moment-generating setups. It lets you cleanly extract moments via derivatives of {{< math >}} $G$ {{< /math >}}.

## ⚙️ 2. Expression for the "kernel" {{< math >}} $U_1$ {{< /math >}}

{{< math >}} 
$$U_1(u_u, u_s, s) = \frac{u_s b}{b-\gamma}(e^{-\gamma s} + u_u - u_s) - \frac{u_s b}{b-\gamma}e^{-bs}$$ 
{{< /math >}}

This expression arises from [solving a linear system of ODEs for the chemical master equation (CME) via generating functions](#Derive-G), and it essentially encodes how the state propagates in time. It's derived from how unspliced → spliced reactions occur over time.

## 🧠 3. Log-generating function {{< math >}} $f$ {{< /math >}}

{{< math >}} 
$$f(u_u, u_s, t) := \ln G(u_u, u_s, t) = \int_0^t \alpha(t-s) U_1(u_u, u_s, s) ds$$ 
{{< /math >}}

This integral form tells us how the log of the generating function evolves, driven by transcription rate {{< math >}} $\alpha(t-s)$ {{< /math >}} and the system dynamics encoded in {{< math >}} $U_1$ {{< /math >}}. Essentially, it's the cumulative effect of production and conversion over time.

Then the log-GF can be written in a linear form in {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}}:

{{< math >}} 
$$f(u_u, u_s, t) = \mu_u(t) u_u + \mu_s(t) u_s$$ 
{{< /math >}}

This suggests that the process is governed by Poisson distributions, since the log-GF is linear in the arguments.

## 📊 4. Explicit distribution — product of Poissons

Given the above, we can now recover the joint distribution {{< math >}} $P(x,t)$ {{< /math >}}:

{{< math >}} 
$$P(x,t) = \frac{\mu_u(t)^{x_u} e^{-\mu_u(t)}}{x_u!} \cdot \frac{\mu_s(t)^{x_s} e^{-\mu_s(t)}}{x_s!}$$ 
{{< /math >}}

So:
- {{< math >}} $x_u \sim \text{Poisson}(\mu_u(t))$ {{< /math >}}
- {{< math >}} $x_s \sim \text{Poisson}(\mu_s(t))$ {{< /math >}}
- Jointly independent

This model assumes that given time {{< math >}} $t$ {{< /math >}}, the spliced and unspliced counts are independent Poisson-distributed variables, whose rates {{< math >}} $\mu_u(t), \mu_s(t)$ {{< /math >}} evolve in time according to the underlying biochemistry.

## 🧮 5. Time-averaged distribution

Finally, since single-cell sequencing samples cells asynchronously in time, the observed distribution over counts is not {{< math >}} $P(x,t)$ {{< /math >}} at a fixed {{< math >}} $t$ {{< /math >}}, but a time-averaged version:

{{< math >}} 
$$P(x) = \frac{1}{T} \int_0^T P(x,t) dt$$ 
{{< /math >}}

This is a mixture of Poissons over time, reflecting the asynchrony of cells in scRNA-seq snapshots. This averaging introduces overdispersion, which is critical to explain the variance observed in real data — greater than what a single Poisson can model.

## ✅ Summary

This passage describes the full pipeline of deriving time-dependent joint distributions of unspliced and spliced mRNA counts in a stochastic birth–death model, and how to transform this into what single-cell experiments would observe via time-averaged Poisson mixtures.


# A) Intro: Multivariate Probability Generating Function (PGF)

## 🔢 1. Multivariate Probability Generating Function (PGF)

This function encodes the entire joint distribution of the random vector {{< math >}} $x = (x_u, x_s)$ {{< /math >}}, i.e., the number of unspliced ({{< math >}} $x_u$ {{< /math >}}) and spliced ({{< math >}} $x_s$ {{< /math >}}) transcripts.

It's a bivariate generating function, meaning it's a function of two complex variables {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}}.

The inclusion of {{< math >}} $+1$ {{< /math >}} makes it a shifted PGF, often done for technical convenience (especially when converting to moment-generating functions).

## 📈 2. Moment Extraction

From the properties of generating functions:

**First Moments:**

{{< math >}} 
$$ \frac{\partial G}{\partial u_u}\bigg|_{u_u = u_s = 0} = E[x_u], \quad \frac{\partial G}{\partial u_s}\bigg|_{u_u = u_s = 0} = E[x_s] $$ 
{{< /math >}}

**Second Moments / Covariances:**

{{< math >}} 
$$ \frac{\partial^2 G}{\partial u_u^2}\bigg|_{u_u = u_s = 0} = E[x_u(x_u - 1)], \quad \frac{\partial^2 G}{\partial u_u \partial u_s}\bigg|_{u_u = u_s = 0} = E[x_u x_s] $$ 
{{< /math >}}

These derivatives allow us to compute variances, covariances, and higher-order statistics of transcript counts.

## 🧬 3. Connection to Biological Reactions

In a linear RNA kinetic model:

- {{< math >}} $x_u$ {{< /math >}}: produced at rate {{< math >}} $\alpha$ {{< /math >}}, converted to {{< math >}} $x_s$ {{< /math >}} at rate {{< math >}} $\beta$ {{< /math >}}
- {{< math >}} $x_s$ {{< /math >}}: degraded at rate {{< math >}} $\gamma$ {{< /math >}}

The evolution of {{< math >}} $G(u_u, u_s, t)$ {{< /math >}} over time follows a partial differential equation that arises from the Chemical Master Equation (CME), which governs the time evolution of probability distributions in chemical kinetics.

## ⏳ 4. Time-Dependence

{{< math >}} $G$ {{< /math >}} is explicitly time-dependent, evolving as the distribution {{< math >}} $P(x,t)$ {{< /math >}} changes.

In some derivations, the log-generating function {{< math >}} $f = \log G$ {{< /math >}} is linear in {{< math >}} $u_u, u_s$ {{< /math >}}, which implies that {{< math >}} $x_u, x_s \sim \text{Poisson}(\mu_u(t)), \text{Poisson}(\mu_s(t))$ {{< /math >}} and are independent.

## 📉 5. Decay and Stationarity

As {{< math >}} $t \to \infty$ {{< /math >}}:

- The means {{< math >}} $\mu_u(t), \mu_s(t)$ {{< /math >}} stabilize.
- So does the generating function, converging to that of a product of Poisson distributions (one for each transcript species).

## 💡 Summary: Why It Matters

- **Encodes all statistical information** about {{< math >}} $x_u, x_s$ {{< /math >}}
- **Enables exact computation** of moments, cumulants, and correlations
- **Links stochastic biochemical kinetics** with observed scRNA-seq distributions
- **Supports modeling** of noise, burstiness, and cell-to-cell heterogeneity


# B) Chemical Master Equation Analysis<a id="Derive-G"></a>

## Model Setup

We consider two molecular species:

- {{< math >}} $x_u(t)$ {{< /math >}}: number of unspliced RNA molecules at time {{< math >}} $t$ {{< /math >}}
- {{< math >}} $x_s(t)$ {{< /math >}}: number of spliced RNA molecules at time {{< math >}} $t$ {{< /math >}}

With the following reactions:

| Reaction | Rate constant | Description |
|----------|---------------|-------------|
| Transcription | {{< math >}} $\alpha$ {{< /math >}} | Production of {{< math >}} $x_u$ {{< /math >}} |
| Splicing | {{< math >}} $\beta$ {{< /math >}} | {{< math >}} $x_u \to x_s$ {{< /math >}} |
| Degradation | {{< math >}} $\gamma$ {{< /math >}} | {{< math >}} $x_s \to \emptyset$ {{< /math >}} |

## Step 1: Chemical Master Equation (CME)

Let

{{< math >}} $$P(x_u, x_s, t) = \Pr[x_u(t) = x_u, x_s(t) = x_s] \qquad (1)$$ {{< /math >}}

The CME is:

{{< math >}} 
$$\begin{align}
\frac{d}{dt}P(x_u, x_s, t) &= \alpha[P(x_u - 1, x_s, t) - P(x_u, x_s, t)] \\
&\quad + \beta[(x_u + 1)P(x_u + 1, x_s - 1, t) - x_u P(x_u, x_s, t)] \\
&\quad + \gamma[(x_s + 1)P(x_u, x_s + 1, t) - x_s P(x_u, x_s, t)] \qquad (2)
\end{align}$$ 
{{< /math >}}

## Step 2: Generating Function Definition

Define the probability generating function (PGF) as:

{{< math >}} $$G(u_u, u_s, t) = \sum_{x_u=0}^{\infty} \sum_{x_s=0}^{\infty} P(x_u, x_s, t)(u_u + 1)^{x_u}(u_s + 1)^{x_s} \qquad (3)$$ {{< /math >}}

**Why {{< math >}} $u_u + 1$ {{< /math >}} and {{< math >}} $u_s + 1$ {{< /math >}} instead of just {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}}?**

This is a common choice in some treatments to simplify expansions around 0, but it's not essential. You can think of {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}} as shift variables.

## Step 3: Deriving PDE for {{< math >}} $G$ {{< /math >}}

Apply the CME to {{< math >}} $G$ {{< /math >}} (see derivation [here](#Derive-step-3-of-B)):

{{< math >}} 
$$\begin{align}
\frac{\partial G}{\partial t} &= \alpha(u_u + 1)G - \alpha G + \beta((u_s + 1)\frac{\partial G}{\partial u_u} - (u_u + 1)\frac{\partial G}{\partial u_u}) \\
&\quad + \gamma(-(u_s + 1)\frac{\partial G}{\partial u_s}) \qquad (4)
\end{align}$$ 
{{< /math >}}

More explicitly:

**The transcription term contributes:**
{{< math >}} $$\alpha[(u_u + 1)G - G] = \alpha u_u G \qquad (5)$$ {{< /math >}}

**The splicing term accounts for a change in {{< math >}} $x_u \to x_s$ {{< /math >}}:**
{{< math >}} $$\beta((u_s + 1)\frac{\partial G}{\partial u_u} - (u_u + 1)\frac{\partial G}{\partial u_u}) = \beta(u_s - u_u)\frac{\partial G}{\partial u_u} \qquad (6)$$ {{< /math >}}

**The degradation term:**
{{< math >}} $$\gamma(-(u_s + 1)\frac{\partial G}{\partial u_s}) = -\gamma(u_s + 1)\frac{\partial G}{\partial u_s} \qquad (7)$$ {{< /math >}}

So the PDE is:

{{< math >}} $$\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma(u_s + 1)\frac{\partial G}{\partial u_s} \qquad (8)$$ {{< /math >}}

## Step 4: Logarithmic Transform

Define:

{{< math >}} $$f(u_u, u_s, t) := \ln G(u_u, u_s, t) \qquad (9)$$ {{< /math >}}

Using chain rule:

{{< math >}} $$\frac{\partial f}{\partial t} = \frac{1}{G}\frac{\partial G}{\partial t} \qquad (10)$$ {{< /math >}}

{{< math >}} $$\frac{\partial f}{\partial u_u} = \frac{1}{G}\frac{\partial G}{\partial u_u}, \quad \frac{\partial f}{\partial u_s} = \frac{1}{G}\frac{\partial G}{\partial u_s} \qquad (11)$$ {{< /math >}}

Rewrite PDE in terms of {{< math >}} $f$ {{< /math >}}:

{{< math >}} $$\frac{\partial f}{\partial t} = \alpha u_u + \beta(u_s - u_u)\frac{\partial f}{\partial u_u} - \gamma(u_s + 1)\frac{\partial f}{\partial u_s} \qquad (12)$$ {{< /math >}}

## Step 5: Method of Characteristics

The PDE is first-order linear and can be solved by characteristics:

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s - u_u), \quad \frac{du_s}{ds} = -\gamma(u_s + 1), \quad \frac{dt}{ds} = 1 \qquad (13)$$ {{< /math >}}

and along characteristics,

{{< math >}} $$\frac{df}{ds} = \alpha u_u \qquad (14)$$ {{< /math >}}

**Boundary conditions:**

{{< math >}} $$t = 0 \Rightarrow f(u_u, u_s, 0) = 0 \Rightarrow G(u_u, u_s, 0) = 1 \qquad (15)$$ {{< /math >}}

## Step 6: Solve ODE for {{< math >}} $u_s$ {{< /math >}}

{{< math >}} $$\frac{du_s}{ds} = -\gamma(u_s + 1) \qquad (16)$$ {{< /math >}}

This is a linear ODE:

{{< math >}} $$\Rightarrow \frac{d}{ds}(u_s + 1) = -\gamma(u_s + 1) \qquad (17)$$ {{< /math >}}

**Solution:**

{{< math >}} $$u_s(s) + 1 = Ce^{-\gamma s} \qquad (18)$$ {{< /math >}}

At {{< math >}} $s = 0$ {{< /math >}}, let {{< math >}} $u_s(0) = u_{s0}$ {{< /math >}}, so:

{{< math >}} $$u_s(0) + 1 = u_{s0} + 1 = C \qquad (19)$$ {{< /math >}}

Thus:

{{< math >}} $$u_s(s) = (u_{s0} + 1)e^{-\gamma s} - 1 \qquad (20)$$ {{< /math >}}

## Step 7: Solve ODE for {{< math >}} $u_u$ {{< /math >}}

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s - u_u) \qquad (21)$$ {{< /math >}}

Rewrite:

{{< math >}} $$\frac{du_u}{ds} + \beta u_u = \beta u_s \qquad (22)$$ {{< /math >}}

**Integrating factor:**

{{< math >}} $$\mu(s) = e^{\beta s} \qquad (23)$$ {{< /math >}}

Multiply both sides:

{{< math >}} $$\frac{d}{ds}(u_u e^{\beta s}) = \beta u_s e^{\beta s} \qquad (24)$$ {{< /math >}}

Integrate from 0 to {{< math >}} $s$ {{< /math >}}:

{{< math >}} $$u_u(s)e^{\beta s} - u_u(0) = \beta \int_0^s u_s(\sigma)e^{\beta \sigma} d\sigma \qquad (25)$$ {{< /math >}}

So:

{{< math >}} $$u_u(s) = e^{-\beta s}\left(u_u(0) + \beta \int_0^s u_s(\sigma)e^{\beta \sigma} d\sigma\right) \qquad (26)$$ {{< /math >}}

Substitute {{< math >}} $u_s(\sigma) = (u_{s0} + 1)e^{-\gamma \sigma} - 1$ {{< /math >}}:

{{< math >}} 
$$\begin{align}
\int_0^s ((u_{s0} + 1)e^{-\gamma \sigma} - 1)e^{\beta \sigma} d\sigma &= (u_{s0} + 1)\int_0^s e^{(\beta - \gamma)\sigma} d\sigma - \int_0^s e^{\beta \sigma} d\sigma \qquad (27)
\end{align}$$ 
{{< /math >}}

Calculate integrals (assuming {{< math >}} $\beta \neq \gamma$ {{< /math >}}):

{{< math >}} $$\int_0^s e^{(\beta - \gamma)\sigma} d\sigma = \frac{e^{(\beta - \gamma)s} - 1}{\beta - \gamma} \qquad (28)$$ {{< /math >}}

{{< math >}} $$\int_0^s e^{\beta \sigma} d\sigma = \frac{e^{\beta s} - 1}{\beta} \qquad (29)$$ {{< /math >}}

Plug in:

{{< math >}} 
$$\begin{align}
u_u(s) &= e^{-\beta s}\left(u_u(0) + \beta(u_{s0} + 1)\frac{e^{(\beta - \gamma)s} - 1}{\beta - \gamma} - \beta\frac{e^{\beta s} - 1}{\beta}\right) \qquad (30)
\end{align}$$ 
{{< /math >}}

Simplify:

{{< math >}} 
$$\begin{align}
u_u(s) &= e^{-\beta s}u_u(0) + (u_{s0} + 1)\frac{\beta}{\beta - \gamma}e^{-\beta s}(e^{(\beta - \gamma)s} - 1) - (e^{-\beta s}(e^{\beta s} - 1)) \\
&= e^{-\beta s}u_u(0) + (u_{s0} + 1)\frac{\beta}{\beta - \gamma}(e^{-\gamma s} - e^{-\beta s}) - (1 - e^{-\beta s}) \qquad (31)
\end{align}$$ 
{{< /math >}}

## Step 8: Solve for {{< math >}} $f$ {{< /math >}}

Recall

{{< math >}} $$\frac{df}{ds} = \alpha u_u(s), \quad f(0) = 0 \qquad (32)$$ {{< /math >}}

Integrate:

{{< math >}} $$f(s) = \alpha \int_0^s u_u(\sigma) d\sigma \qquad (33)$$ {{< /math >}}

Using the expression for {{< math >}} $u_u(s)$ {{< /math >}}, integrate term by term:

{{< math >}} 
$$\begin{align}
f(s) &= \alpha \int_0^s \left[e^{-\beta \sigma}u_u(0) + (u_{s0} + 1)\frac{\beta}{\beta - \gamma}(e^{-\gamma \sigma} - e^{-\beta \sigma}) - (1 - e^{-\beta \sigma})\right] d\sigma \qquad (34)
\end{align}$$ 
{{< /math >}}

## Step 9: Interpretation

The function {{< math >}} $U_1(u_u, u_s, s)$ {{< /math >}} you asked about corresponds to the time-dependent kernel:

{{< math >}} $$U_1(u_u, u_s, s) := u_u(s) \qquad (35)$$ {{< /math >}}

This kernel describes the influence of the initial "dummy variables" {{< math >}} $u_u(0), u_{s0}$ {{< /math >}} on the system at time {{< math >}} $s$ {{< /math >}}.

The generating function is then:

{{< math >}} $$G(u_u, u_s, t) = \exp\left(\alpha \int_0^t U_1(u_u, u_s, s) ds\right) \qquad (36)$$ {{< /math >}}


# B.1) Derive step 3 of B)

## Step 1: Write the CME explicitly
<a id="Derive-step-3-of-B"></a>
The CME for the joint distribution {{< math >}} $P(x_u, x_s, t)$ {{< /math >}} of the unspliced {{< math >}} $x_u$ {{< /math >}} and spliced {{< math >}} $x_s$ {{< /math >}} RNA is:

{{< math >}}
$$
\frac{d}{dt}P(x_u, x_s, t) = \alpha[P(x_u - 1, x_s, t) - P(x_u, x_s, t)] + \beta[(x_u + 1)P(x_u + 1, x_s - 1, t) - x_u P(x_u, x_s, t)] + \gamma[(x_s + 1)P(x_u, x_s + 1, t) - x_s P(x_u, x_s, t)] \qquad (9)
$$
{{< /math >}}

where:
- {{< math >}} $\alpha$ {{< /math >}}: production rate of unspliced RNA
- {{< math >}} $\beta$ {{< /math >}}: splicing rate (unspliced → spliced)
- {{< math >}} $\gamma$ {{< /math >}}: degradation rate of spliced RNA

## Step 2: Define the generating function {{< math >}} $G(u_u, u_s, t)$ {{< /math >}}

{{< math >}}
$$
G(u_u, u_s, t) := \sum_{x_u=0}^{\infty} \sum_{x_s=0}^{\infty} P(x_u, x_s, t)(u_u + 1)^{x_u}(u_s + 1)^{x_s} \qquad (10)
$$
{{< /math >}}

**Note:** Using {{< math >}} $(u_u + 1)^{x_u}$ {{< /math >}} instead of {{< math >}} $u_u^{x_u}$ {{< /math >}} is a common shift to simplify derivatives later, but you can equivalently define with {{< math >}} $u_u^{x_u}$ {{< /math >}}.

## Step 3: Take time derivative of {{< math >}} $G$ {{< /math >}}

Using linearity of sums and derivatives:

{{< math >}}
$$
\frac{\partial G}{\partial t} = \sum_{x_u=0}^{\infty} \sum_{x_s=0}^{\infty} \frac{\partial P(x_u, x_s, t)}{\partial t} (u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Substitute the CME expression:

{{< math >}}
$$
= \sum_{x_u,x_s} [\alpha(P(x_u - 1, x_s, t) - P(x_u, x_s, t)) + \beta((x_u + 1)P(x_u + 1, x_s - 1, t) - x_u P(x_u, x_s, t)) + \gamma((x_s + 1)P(x_u, x_s + 1, t) - x_s P(x_u, x_s, t))](u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

## Step 4: Evaluate each term separately

### Term 1: Transcription

{{< math >}}
$$
\sum_{x_u,x_s} \alpha(P(x_u - 1, x_s) - P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Rewrite the first sum by shifting {{< math >}} $x_u \to x_u + 1$ {{< /math >}} in the first part:

{{< math >}}
$$
\sum_{x_u=0}^{\infty} P(x_u - 1, x_s)(u_u + 1)^{x_u} = \sum_{x_u'=-1}^{\infty} P(x_u', x_s)(u_u + 1)^{x_u' + 1}
$$
{{< /math >}}

Since {{< math >}} $P(x_u', x_s) = 0$ {{< /math >}} for {{< math >}} $x_u' < 0$ {{< /math >}}, this becomes:

{{< math >}}
$$
(u_u + 1) \sum_{x_u'=0}^{\infty} P(x_u', x_s)(u_u + 1)^{x_u'}
$$
{{< /math >}}

Therefore,

{{< math >}}
$$
\sum_{x_u,x_s} \alpha(P(x_u - 1, x_s) - P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s} = \alpha((u_u + 1)G - G) = \alpha u_u G
$$
{{< /math >}}

### Term 2: Splicing

{{< math >}}
$$
\sum_{x_u,x_s} \beta((x_u + 1)P(x_u + 1, x_s - 1) - x_u P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Split into two sums:

{{< math >}}
$$
S_1 = \beta \sum_{x_u,x_s} (x_u + 1)P(x_u + 1, x_s - 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

{{< math >}}
$$
S_2 = -\beta \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Change indices in {{< math >}} $S_1$ {{< /math >}}:

Let {{< math >}} $x_u' = x_u + 1$ {{< /math >}}, {{< math >}} $x_s' = x_s - 1$ {{< /math >}}, so {{< math >}} $x_u = x_u' - 1$ {{< /math >}}, {{< math >}} $x_s = x_s' + 1$ {{< /math >}}

Then,

{{< math >}}
$$
S_1 = \beta \sum_{x_u'=1}^{\infty} \sum_{x_s'=0}^{\infty} x_u' P(x_u', x_s')(u_u + 1)^{x_u' - 1}(u_s + 1)^{x_s' + 1}
$$
{{< /math >}}

Rearranged:<a id="why_sum_is_equivalence"></a>
 (see [note](#b2-why-double-sum-equivalence-holds-in-term-2-of-step4-in-b1))
 
{{< math >}}
$$
= \beta(u_s + 1) \sum_{x_u',x_s'} x_u' P(x_u', x_s')(u_u + 1)^{x_u' - 1}(u_s + 1)^{x_s'}
$$
{{< /math >}}

{{< math >}} $S_2$ {{< /math >}} is:

{{< math >}}
$$
S_2 = -\beta \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Now recognize:

{{< math >}}
$$
\frac{\partial G}{\partial u_u} = \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u - 1}(u_s + 1)^{x_s}
$$
{{< /math >}}

So:

{{< math >}}
$$
S_1 = \beta(u_s + 1) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

and

{{< math >}}
$$
S_2 = -\beta(u_u + 1) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

Putting together:

{{< math >}}
$$
\text{Splicing term} = \beta((u_s + 1) - (u_u + 1)) \frac{\partial G}{\partial u_u} = \beta(u_s - u_u) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

### Term 3: Degradation

{{< math >}}
$$
\sum_{x_u,x_s} \gamma((x_s + 1)P(x_u, x_s + 1) - x_s P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Split into:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u,x_s} (x_s + 1)P(x_u, x_s + 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

{{< math >}}
$$
S_4 = -\gamma \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Shift indices in {{< math >}} $S_3$ {{< /math >}}:

Let {{< math >}} $x_s' = x_s + 1 \Rightarrow x_s = x_s' - 1$ {{< /math >}}:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u=0}^{\infty} \sum_{x_s'=1}^{\infty} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s' - 1}
$$
{{< /math >}}

{{< math >}}
$$
= \gamma(u_s + 1)^{-1} \sum_{x_u,x_s'} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s'}
$$
{{< /math >}}

Recognize

{{< math >}}
$$
\frac{\partial G}{\partial u_s} = \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s - 1}
$$
{{< /math >}}

So,

{{< math >}}
$$
S_3 = \gamma(u_s + 1) \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

and

{{< math >}}
$$
S_4 = -\gamma(u_s + 1) \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

Wait, let me recalculate {{< math >}} $S_4$ {{< /math >}}:

{{< math >}}
$$
S_4 = -\gamma \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

This equals:

{{< math >}}
$$
S_4 = -\gamma(u_s + 1) \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s - 1} = -\gamma(u_s + 1) \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

Therefore, the degradation term is:

{{< math >}}
$$
\text{Degradation term} = S_3 + S_4 = \gamma(u_s + 1) \frac{\partial G}{\partial u_s} - \gamma(u_s + 1) \frac{\partial G}{\partial u_s} = 0
$$
{{< /math >}}

Actually, let me recalculate this more carefully. We have:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u,x_s} (x_s + 1)P(x_u, x_s + 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

With the substitution {{< math >}} $x_s' = x_s + 1$ {{< /math >}}, we get {{< math >}} $x_s = x_s' - 1$ {{< /math >}} and:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u=0}^{\infty} \sum_{x_s'=1}^{\infty} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s' - 1}
$$
{{< /math >}}

{{< math >}}
$$
= \frac{\gamma}{u_s + 1} \sum_{x_u,x_s'} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s'}
$$
{{< /math >}}

{{< math >}}
$$
= \frac{\gamma}{u_s + 1} \cdot (u_s + 1) \frac{\partial G}{\partial u_s} = \gamma \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

And:

{{< math >}}
$$
S_4 = -\gamma(u_s + 1) \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

So the degradation term is:

{{< math >}}
$$
\text{Degradation term} = \gamma \frac{\partial G}{\partial u_s} - \gamma(u_s + 1) \frac{\partial G}{\partial u_s} = -\gamma u_s \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

## Final Result

Combining all terms:

{{< math >}}
$$
\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u) \frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s} \qquad (11)
$$
{{< /math >}}

# B.2) Why Double Sum Equivalence Holds [in Term 2 of Step4 in B.1)](#why_sum_is_equivalence)

This equivalence depends on what values the indices range over and how the function being summed behaves.

## 1. Notation

When we write:

{{< math >}} $$\sum_{x_u', x_s'} f(x_u', x_s')$$ {{< /math >}}

This is shorthand for:

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty} f(x_u', x_s')$$ {{< /math >}}

That is, summing over all nonnegative integer pairs {{< math >}} $(x_u', x_s') \in \mathbb{N}_0 \times \mathbb{N}_0$ {{< /math >}}.

## 2. Why it matters in the derivation

In your original sum:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') (u_u + 1)^{x_u' - 1} (u_s + 1)^{x_s'}$$ {{< /math >}}

The lower bound of {{< math >}} $x_u'$ {{< /math >}} is 1 because we performed a change of variables from {{< math >}} $x_u = x_u' - 1$ {{< /math >}}, and in the original sum, {{< math >}} $x_u \geq 0$ {{< /math >}}, which implies {{< math >}} $x_u' \geq 1$ {{< /math >}}.

So the double sum with bounds:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

is not exactly the same as:

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

But if we define:

{{< math >}} $$\sum_{x_u', x_s'} := \sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

then in your derivation, the support of {{< math >}} $P(x_u', x_s')$ {{< /math >}} makes this safe because:

{{< math >}} $$x_u' P(x_u', x_s') = 0 \text{ when } x_u' = 0$$ {{< /math >}}

since the factor is 0.

So extending the lower limit to 0 adds no contribution to the sum.

## ✅ Conclusion

The expressions:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') \cdots$$ {{< /math >}}

and

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') \cdots$$ {{< /math >}}

are equal because when {{< math >}} $x_u' = 0$ {{< /math >}}, the term is 0.

Hence, we can write:

{{< math >}} $$\sum_{x_u', x_s'} := \sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

without affecting the value of the sum.