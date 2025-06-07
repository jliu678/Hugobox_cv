---
title: scVelo math derivation
math: true
date: '2025-05-28'
---
## GMM EM
```python
initialize mu_k, sigma_k, pi_k for k in 1..K

while not converged:

    # E-step
    for x_i in data:
        for k in 1..K:
            log_p_x_given_z = -((x_i - mu_k)^2) / (2 * sigma_k^2)
            log_p_z = log(pi_k)
            logliks[k] = log_p_x_given_z + log_p_z
        q_i = softmax(logliks)

    # M-step
    for k in 1..K:
        N_k = sum_i q_i[k]
        mu_k = sum_i q_i[k] * x_i / N_k
        sigma_k = sqrt(sum_i q_i[k] * (x_i - mu_k)^2 / N_k)
        pi_k = N_k / N

```

## scVelo EM
```python
initialize theta = (α, β, γ)

while not converged:

    # E-step
    mu_traj = [solve_ODE(t, theta) for t in time_grid]
    for x_i in data:
        for j, mu_t in enumerate(mu_traj):
            logliks[j] = -||x_i - mu_t||^2 / σ
        q_i = softmax(logliks)

    # M-step
    def loss_fn(theta):
        loss = 0
        for x_i, q_i in zip(data, q_all):
            for j, t_j in enumerate(time_grid):
                mu_tj = solve_ODE(t_j, theta)
                loss += q_i[j] * ||x_i - mu_tj||^2
        return loss

    theta = optimize(loss_fn)
````

## scVelo EM with Learnable Time Prior π_j 
```python
initialize theta = (α, β, γ)
initialize π = [1/T] * T  # uniform over time grid

while not converged:

    # E-step
    mu_traj = [solve_ODE(t, theta) for t in time_grid]
    for x_i in data:
        for j, mu_t in enumerate(mu_traj):
            logliks[j] = -||x_i - mu_t||^2 / σ + log(π[j])
        q_i = softmax(logliks)

    # M-step
    π[j] = sum_i q_i[j] / N  for all j  # update time prior

    def loss_fn(theta):
        loss = 0
        for x_i, q_i in zip(data, q_all):
            for j, t_j in enumerate(time_grid):
                mu_tj = solve_ODE(t_j, theta)
                loss += q_i[j] * ||x_i - mu_tj||^2
        return loss

    theta = optimize(loss_fn)
```