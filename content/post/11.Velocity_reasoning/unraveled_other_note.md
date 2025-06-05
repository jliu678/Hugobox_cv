# 🧭 What is an Occupation Measure?

An occupation measure describes the expected amount of time a stochastic process spends in a given state or set of states over a certain time interval.

For a stochastic process {{< math >}} $X_t$ {{< /math >}}, the occupation measure {{< math >}} $\mu_T$ {{< /math >}} over the interval {{< math >}} $[0,T]$ {{< /math >}} is defined as:

{{< math >}} 
$$\mu_T(A) = \int_0^T 1_{X_t \in A} \, dt \quad (1)$$
{{< /math >}}

where:

- {{< math >}} $A \subseteq \mathbb{R}^n$ {{< /math >}} is a set of states,
- {{< math >}} $1_{X_t \in A}$ {{< /math >}} is the indicator function (1 if {{< math >}} $X_t \in A$ {{< /math >}}, 0 otherwise),
- {{< math >}} $\mu_T(A)$ {{< /math >}} is the total time spent in region {{< math >}} $A$ {{< /math >}} up to time {{< math >}} $T$ {{< /math >}}.

## 🔁 Expected Occupation Measure

Taking the expectation over many realizations of the process gives:

{{< math >}} 
$$E[\mu_T(A)] = \int_0^T P(X_t \in A) \, dt \quad (2)$$
{{< /math >}}

This tells us how likely and for how long the process occupies a particular state (or set of states).

## 🔬 In the Chemical Master Equation (CME) Context

For a discrete-state Markov process governed by the CME, the occupation measure for a specific state {{< math >}} $x \in \mathbb{N}_0^n$ {{< /math >}} is:

{{< math >}} 
$$\mu_T(x) = \int_0^T P(x,t) \, dt \quad (3)$$
{{< /math >}}

This represents the expected time the system spends in state {{< math >}} $x$ {{< /math >}} over the interval {{< math >}} $[0,T]$ {{< /math >}}.

## 🧠 Why Is This Useful?

**Understanding System Behavior**: Reveals which states are most frequently visited or "important" over time.

**Single-Cell Dynamics**: In tools like scVelo, occupation measures relate to how long a cell stays in certain gene expression states.

**Model Evaluation**: Occupation measures can help identify dominant dynamics or regimes in large state spaces.

**Optimization and Control**: In Markov decision processes, they relate to how policies influence long-term behavior.

## 🧪 Simple Example: Birth-Death Process

Assume:

- Birth (transcription) at rate {{< math >}} $\alpha$ {{< /math >}}
- Death (degradation) at rate {{< math >}} $\beta x$ {{< /math >}}

Let {{< math >}} $P(x,t)$ {{< /math >}} be the solution of the CME. Then the occupation measure for state {{< math >}} $x$ {{< /math >}} is:

{{< math >}} 
$$\mu_T(x) = \int_0^T P(x,t) \, dt \quad (4)$$
{{< /math >}}

If {{< math >}} $\mu_T(x)$ {{< /math >}} is large for small {{< math >}} $x$ {{< /math >}}, it means the system spends most of its time in low expression states.


# 🔄 Desynchronization Assumption

The key biological assumption is:

> **Each cell is captured at a random point in its trajectory** — i.e., the sequencing process is not synchronized with the cell's internal dynamics or with other cells.

This reflects the asynchronous nature of scRNA-seq snapshots: cells are sampled independently and randomly in time, and we don't know when in their dynamic process each cell was observed.

## 📐 Translating This to a Mathematical Framework

If the process time interval over which the cell dynamics evolve is {{< math >}} $[0,T]$ {{< /math >}}, then the probability of observing a state {{< math >}} $x$ {{< /math >}} is assumed to be proportional to the fraction of time the system spends in that state.

This assumption leads to defining the empirical distribution over states from the observed dynamics using:

{{< math >}} 
$$ df = \frac{1}{T} dt $$
{{< /math >}}

This means the density of observations (across time) is proportional to time spent in each state.

## 📉 Limiting Behavior

Taking this to the limit as {{< math >}} $T \to \infty$ {{< /math >}}, we assume the system reaches stationarity. So, the probability of observing state {{< math >}} $x$ {{< /math >}}, denoted {{< math >}} $P(x)$ {{< /math >}}, is:

{{< math >}} 
$$ P(x) = \lim_{T \to \infty} \frac{1}{T} \int_0^T P(x,t) \, dt $$
{{< /math >}}

**In words:**
The long-term probability of observing a cell in state {{< math >}} $x$ {{< /math >}} is the time average of the transient probabilities {{< math >}} $P(x,t)$ {{< /math >}}.

If the system is ergodic and reaches a steady state, this simplifies further to:

{{< math >}} 
$$ P(x) = \lim_{t \to \infty} P(x,t) $$
{{< /math >}}

## 🔄 Transient vs. Steady-State Behavior

Let's say you have a stochastic process (e.g. a gene expression level evolving over time due to transcription and degradation). This process has a probability distribution over states {{< math >}} $x$ {{< /math >}} at any time {{< math >}} $t$ {{< /math >}}, which we call:

{{< math >}} $P(x,t)$ {{< /math >}} (the probability the system is in state {{< math >}} $x$ {{< /math >}} at time {{< math >}} $t$ {{< /math >}})

This is called the **transient probability** because it changes with time.

## ⏱️ Time-Averaging and Observation

In single-cell RNA-seq, we don't track individual cells over time. Instead, we take snapshots from a large population of cells that are each somewhere along their internal process. Because of this:

- We assume each observed cell is randomly sampled from a long-running process.
- This leads to the assumption that the probability of observing a state {{< math >}} $x$ {{< /math >}} is proportional to the amount of time the system spends in that state.

That's captured by the time-averaged probability:

{{< math >}} 
$$ P(x) = \lim_{T \to \infty} \frac{1}{T} \int_0^T P(x,t) \, dt $$
{{< /math >}}

This is called the **occupation measure** or **empirical distribution** — it tells us how often the system is in state {{< math >}} $x$ {{< /math >}} over time.

## ✅ If the System is Ergodic

A process is **ergodic** if its time averages and ensemble averages agree. This means:

- As {{< math >}} $t \to \infty$ {{< /math >}}, the distribution {{< math >}} $P(x,t)$ {{< /math >}} settles to a fixed value.
- That value is called the **stationary distribution**, denoted:

{{< math >}} 
$$ P(x) = \lim_{t \to \infty} P(x,t) $$
{{< /math >}}

So for ergodic systems, the time-average and the limiting distribution are the same:

{{< math >}} 
$$ \lim_{T \to \infty} \frac{1}{T} \int_0^T P(x,t) \, dt = \lim_{t \to \infty} P(x,t) $$
{{< /math >}}

## 🧬 In scRNA-seq Context

We only get one snapshot per cell, but thousands of cells. So we treat the whole collection of observations as a sample from the time-averaged distribution, which — under ergodicity — approximates the stationary distribution.

Thus, each cell's expression state is treated as a draw from {{< math >}} $P(x) = \lim_{t \to \infty} P(x,t)$ {{< /math >}}.

## 🔁 Summary Table

| Concept | Meaning |
|---------|---------|
| {{< math >}} $P(x,t)$ {{< /math >}} | Probability the system is in state {{< math >}} $x$ {{< /math >}} at time {{< math >}} $t$ {{< /math >}} |
| {{< math >}} $\frac{1}{T} \int_0^T P(x,t) dt$ {{< /math >}} | Time average — how long the system spends in {{< math >}} $x$ {{< /math >}} |
| {{< math >}} $\lim_{t \to \infty} P(x,t)$ {{< /math >}} | Stationary distribution (if the process stabilizes over time) |
| Ergodicity | Ensures time averages equal long-term behavior (makes interpretation valid) |