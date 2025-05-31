# Derivation: Transforming the Solution for s(t)

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