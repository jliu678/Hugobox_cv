# VeloVI Parameter Initialization and Learning Process

Illustration of the initiation and learning process for all parameters in VeloVI for our scenario of three genes (Gene A, Gene B, Gene C) in four cells (Cell 1, Cell 2, Cell 3, Cell 4).

This includes:

**Biophysical Parameters** ({{< math >}} $\theta$ {{< /math >}}): {{< math >}} $\alpha_{g1}, \beta_g, \gamma_g, t_g^s$ {{< /math >}}

**Variance Parameters** ({{< math >}} $\theta$ {{< /math >}}): {{< math >}} $\sigma_{gu}, \sigma_{gs}$ {{< /math >}} (gene-specific standard deviations for unspliced and spliced RNA)

**Neural Network Parameters** ({{< math >}} $\phi$ {{< /math >}}): Weights and biases of the encoder and decoder networks.

## 1. Initial Data and Preprocessing

Our simplified min-max scaled data:

| Cell | Gene | Unspliced (u) | Spliced (s) |
|------|------|---------------|-------------|
| 1    | A    | 0.8           | 0.2         |
| 1    | B    | 0.1           | 0.9         |
| 1    | C    | 0.5           | 0.5         |
| 2    | A    | 0.9           | 0.1         |
| 2    | B    | 0.2           | 0.8         |
| 2    | C    | 0.6           | 0.4         |
| 3    | A    | 0.2           | 0.7         |
| 3    | B    | 0.7           | 0.3         |
| 3    | C    | 0.1           | 0.9         |
| 4    | A    | 0.1           | 0.8         |
| 4    | B    | 0.8           | 0.2         |
| 4    | C    | 0.2           | 0.8         |

## 2. Initialization of All Parameters

### a. Biophysical Parameters ({{< math >}} $\theta$ {{< /math >}})

**{{< math >}} $\alpha_{g1}$ {{< /math >}} (Transcription Rate for Induction):**
- {{< math >}} $\alpha_{A1} = 0.85$ {{< /math >}} (median unspliced for top 99th percentile, e.g., Cells 1, 2 for Gene A)
- {{< math >}} $\alpha_{B1} = 0.75$ {{< /math >}} (median unspliced for top 99th percentile, e.g., Cells 3, 4 for Gene B)
- {{< math >}} $\alpha_{C1} = 0.55$ {{< /math >}} (median unspliced for top 99th percentile, e.g., Cells 1, 2 for Gene C)

**{{< math >}} $\beta_g$ {{< /math >}} (Splicing Rate):** Initialized to a small positive constant.
- {{< math >}} $\beta_A = \beta_B = \beta_C = 0.1$ {{< /math >}}

**{{< math >}} $\gamma_g$ {{< /math >}} (Degradation Rate):** Initialized to a small positive constant.
- {{< math >}} $\gamma_A = \gamma_B = \gamma_C = 0.05$ {{< /math >}}

**{{< math >}} $t_g^s$ {{< /math >}} (Switching Time):** Initialized to a constant value within {{< math >}} $[0, t_{\max}]$ {{< /math >}}.
- {{< math >}} $t_A^s = t_B^s = t_C^s = 10.0$ {{< /math >}} (assuming {{< math >}} $t_{\max} = 20$ {{< /math >}})

### b. Variance Parameters ({{< math >}} $\theta$ {{< /math >}})

These are also part of the {{< math >}} $\theta$ {{< /math >}} parameters learned by the model. They are gene-specific.

**{{< math >}} $\sigma_{gu}, \sigma_{gs}$ {{< /math >}} (Unspliced and Spliced Standard Deviations):** Typically initialized to a small, common positive constant, reflecting initial uncertainty or baseline noise.
- {{< math >}} $\sigma_{Au} = \sigma_{Bu} = \sigma_{Cu} = 0.1$ {{< /math >}}
- {{< math >}} $\sigma_{As} = \sigma_{Bs} = \sigma_{Cs} = 0.1$ {{< /math >}}

(Note: {{< math >}} $c_k$ {{< /math >}} is a fixed scalar, not a learnable parameter, e.g., {{< math >}} $c_4 = 0.1$ {{< /math >}}, others {{< math >}} $c_k = 1.0$ {{< /math >}})

### c. Neural Network Parameters ({{< math >}} $\phi$ {{< /math >}})

These are the weights and biases of all neural networks in the model.

**Encoder Network Weights/Biases:** This network takes {{< math >}} $(u_n, s_n)$ {{< /math >}} as input and outputs:
- Parameters for {{< math >}} $q_\phi(z_n | u_n, s_n)$ {{< /math >}} (e.g., {{< math >}} $\mu_{z_n}, \log\Sigma_{z_n}$ {{< /math >}} for a Gaussian).
- Parameters for {{< math >}} $q_\phi(\pi_{ng} | z_n)$ {{< /math >}} (e.g., concentration parameters for a Dirichlet).

*Initialization:* Typically initialized randomly, often using methods like Xavier/Glorot or Kaiming initialization, which help with training stability. Let's represent this as a large set of random numbers.

**Decoder Network** ({{< math >}} $h_{\text{ind}}, h_{\text{rep}}$ {{< /math >}}) **Weights/Biases:** These networks take {{< math >}} $z_n$ {{< /math >}} as input and output the scaled latent times {{< math >}} $\rho_{ng}^{(1)}$ {{< /math >}} and {{< math >}} $\rho_{ng}^{(3)}$ {{< /math >}}.

*Initialization:* Also randomly initialized, similar to the encoder.

## 3. Learning Process (One Iteration of Optimization)

Let's trace one optimization step focusing on a mini-batch of Cell 1 and Cell 3.

### a. Forward Pass (Calculate Loss {{< math >}} $L_{\text{velo}}$ {{< /math >}})

**For Cell 1** ({{< math >}} $u_1, s_1$ {{< /math >}}):

**Encoder Inference** ({{< math >}} $q_\phi(z_1 | u_1, s_1)$ {{< /math >}} and {{< math >}} $q_\phi(\pi_{1g} | z_1)$ {{< /math >}}):

The encoder network (with its current {{< math >}} $\phi$ {{< /math >}} weights/biases) takes {{< math >}} $(u_{1A}, s_{1A}, u_{1B}, s_{1B}, u_{1C}, s_{1C})$ {{< /math >}} as input.
- **Output 1:** {{< math >}} $\mu_{z_1}$ {{< /math >}} and {{< math >}} $\log\Sigma_{z_1}$ {{< /math >}} (e.g., a 10-dimensional vector for {{< math >}} $\mu$ {{< /math >}} and {{< math >}} $\log\Sigma$ {{< /math >}}).
- **Output 2:** For each gene {{< math >}} $g$ {{< /math >}}, Dirichlet concentration parameters for {{< math >}} $\pi_{1g}$ {{< /math >}} (e.g., for Gene A, {{< math >}} $\text{Dirichlet}(\text{conc}_{1A,1}, \text{conc}_{1A,2}, \text{conc}_{1A,3}, \text{conc}_{1A,4})$ {{< /math >}}).

**Sampling Latent Variables** ({{< math >}} $z_1^*$ {{< /math >}} and {{< math >}} $\pi_{1g}^*$ {{< /math >}}, leading to {{< math >}} $k_{1g}$ {{< /math >}} and {{< math >}} $t_{1g}^{(k)}$ {{< /math >}}):

- **{{< math >}} $z_1^*$ {{< /math >}}:** Sample {{< math >}} $z_1^*$ {{< /math >}} using the reparameterization trick: {{< math >}} $z_1^* = \mu_{z_1} + \exp(0.5\log\Sigma_{z_1}) \odot \epsilon_z$ {{< /math >}} (where {{< math >}} $\epsilon_z \sim \mathcal{N}(0, I_d)$ {{< /math >}}).
- **{{< math >}} $\pi_{1g}^*$ {{< /math >}} and {{< math >}} $k_{1g}$ {{< /math >}}:** For each gene, sample {{< math >}} $\pi_{1g}^*$ {{< /math >}} from its Dirichlet. Then, based on {{< math >}} $\pi_{1g}^*$ {{< /math >}}, a state {{< math >}} $k_{1g}$ {{< /math >}} is determined (e.g., if {{< math >}} $\pi_{1A}^* = [0.7, 0.1, 0.1, 0.1]$ {{< /math >}}, {{< math >}} $k_{1A}$ {{< /math >}} is likely 1).

**Latent Times:** The decoder networks ({{< math >}} $h_{\text{ind}}, h_{\text{rep}}$ {{< /math >}}), using {{< math >}} $z_1^*$ {{< /math >}} and their current {{< math >}} $\phi$ {{< /math >}} weights/biases, output {{< math >}} $\rho_{1g}^{(1)}$ {{< /math >}} and {{< math >}} $\rho_{1g}^{(3)}$ {{< /math >}}.
- If {{< math >}} $k_{1g} = 1$ {{< /math >}}, then {{< math >}} $(t_1)_{1g} = \rho_{1g}^{(1)} \times t_g^s$ {{< /math >}}.
- If {{< math >}} $k_{1g} = 3$ {{< /math >}}, then {{< math >}} $(t_3)_{1g} = t_g^s + (t_{\max} - t_g^s)\rho_{1g}^{(3)}$ {{< /math >}}.

*Example for Cell 1, Gene A:* Assume {{< math >}} $z_1^*$ {{< /math >}} leads to {{< math >}} $k_{1A} = 1$ {{< /math >}} and {{< math >}} $(t_1)_{1A} = 8.0$ {{< /math >}}.

**Calculate Predicted Abundances** ({{< math >}} $\bar{u}^{(g)}, \bar{s}^{(g)}$ {{< /math >}}):

Using the sampled {{< math >}} $k_{1g}$ {{< /math >}} and {{< math >}} $t_{1g}^{(k)}$ {{< /math >}}, along with the current biophysical parameters ({{< math >}} $\alpha_{g1}, \beta_g, \gamma_g, t_g^s$ {{< /math >}}), calculate {{< math >}} $\bar{u}^{(g)}$ {{< /math >}} and {{< math >}} $\bar{s}^{(g)}$ {{< /math >}} using the kinetic equations.

*Example for Cell 1, Gene A:* Using initial {{< math >}} $\alpha_{A1} = 0.85, \beta_A = 0.1, \gamma_A = 0.05, t_A^s = 10.0$ {{< /math >}} and {{< math >}} $t_{1A}^{(1)} = 8.0$ {{< /math >}}:
- {{< math >}} $\bar{u}^{(A)}(8.0, 1) \approx 4.68$ {{< /math >}}
- {{< math >}} $\bar{s}^{(A)}(8.0, 1) \approx 1.85$ {{< /math >}}

(Observed: {{< math >}} $u_{1A} = 0.8, s_{1A} = 0.2$ {{< /math >}}. Still a large mismatch.)

**Calculate Likelihoods** (using {{< math >}} $\sigma_{gu}, \sigma_{gs}$ {{< /math >}}):

For each gene, compare observed {{< math >}} $(u_{ng}, s_{ng})$ {{< /math >}} to predicted {{< math >}} $(\bar{u}^{(g)}, \bar{s}^{(g)})$ {{< /math >}} using Normal distributions. The variance for these Normal distributions comes from the variance parameters ({{< math >}} $\sigma_{gu}, \sigma_{gs}$ {{< /math >}}) and the state-dependent scaling {{< math >}} $c_k$ {{< /math >}}.

E.g., for Cell 1, Gene A ({{< math >}} $k=1$ {{< /math >}}):
{{< math >}} $$\text{Likelihood}_{1A} = \mathcal{N}(u_{1A} | \bar{u}^{(A)}, \sigma_{Au}^2) \times \mathcal{N}(s_{1A} | \bar{s}^{(A)}, \sigma_{As}^2)$$ {{< /math >}}

Using initial {{< math >}} $\sigma_{Au} = 0.1, \sigma_{As} = 0.1$ {{< /math >}}:
{{< math >}} $\mathcal{N}(0.8 | 4.68, 0.1^2) \times \mathcal{N}(0.2 | 1.85, 0.1^2)$ {{< /math >}} will yield an extremely tiny number, leading to a very large negative log-likelihood contribution.

The total reconstruction term is the sum of these log-likelihoods over all genes and cells in the mini-batch.

**Calculate KL Divergence Terms:**

- **KL for {{< math >}} $z_n$ {{< /math >}}:** {{< math >}} $\text{KL}(q_\phi(z_1 | u_1, s_1) \| p(z))$ {{< /math >}}. This is calculated directly from the {{< math >}} $\mu_{z_1}$ {{< /math >}} and {{< math >}} $\log\Sigma_{z_1}$ {{< /math >}} outputs of the encoder network. This term penalizes deviation from a standard normal prior.
- **KL for {{< math >}} $\pi_{ng}$ {{< /math >}}:** {{< math >}} $\mathbb{E}_{q_\phi(z_1 | u_1, s_1)}[\sum_g \text{KL}(q_\phi(\pi_{1g} | z_1) \| p(\pi_{1g}))]$ {{< /math >}}. This is approximated using {{< math >}} $z_1^*$ {{< /math >}}. It penalizes deviation of state probabilities from a uniform prior.

**Calculate Switching Time Penalty** ({{< math >}} $L_{\text{switch}}$ {{< /math >}}):

This term depends on biophysical parameters ({{< math >}} $\alpha_{g1}, \beta_g, \gamma_g, t_g^s$ {{< /math >}}) as it compares the model's steady-state predictions at {{< math >}} $t_g^s$ {{< /math >}} to the empirically derived switch points.

*Example for Gene A:* Initial {{< math >}} $u_{A3}^0 = \alpha_{A1}/\beta_A = 0.85/0.1 = 8.5$ {{< /math >}}. If {{< math >}} $u_A^* = 0.85$ {{< /math >}}, then {{< math >}} $(8.5 - 0.85)^2$ {{< /math >}} is a large penalty.

**Total Loss** ({{< math >}} $L_{\text{velo}}$ {{< /math >}}): Sum all these terms (negative ELBO components + switch penalty) for Cell 1 and Cell 3. This will be a large positive number due to poor initial fit.

### b. Backward Pass (Backpropagation to Compute Gradients)

Backpropagation efficiently computes {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial \text{parameter}}$ {{< /math >}} for every single parameter in the model ({{< math >}} $\theta$ {{< /math >}} and {{< math >}} $\phi$ {{< /math >}}).

**Start from the Loss:** The process begins with {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial L_{\text{velo}}} = 1$ {{< /math >}}.

**Gradients through Likelihoods to Predicted Means and Variances:**

Gradients flow from {{< math >}} $L_{\text{velo}}$ {{< /math >}} back through the log-likelihood terms.
- {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial \bar{u}^{(g)}}, \frac{\partial L_{\text{velo}}}{\partial \bar{s}^{(g)}}$ {{< /math >}} (how much changing predictions affects loss).
- Crucially, gradients also flow to the variance parameters: {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial \sigma_{gu}}, \frac{\partial L_{\text{velo}}}{\partial \sigma_{gs}}$ {{< /math >}}. If the current {{< math >}} $\sigma_g$ {{< /math >}} values make the likelihood too sharp (if observed is far from predicted) or too flat (if observed is very close to predicted but {{< math >}} $\sigma_g$ {{< /math >}} is huge), these gradients will indicate how to adjust {{< math >}} $\sigma_g$ {{< /math >}} to better fit the data.

**Gradients to Biophysical Parameters** ({{< math >}} $\alpha_{g1}, \beta_g, \gamma_g, t_g^s$ {{< /math >}}):

From {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial \bar{u}^{(g)}}, \frac{\partial L_{\text{velo}}}{\partial \bar{s}^{(g)}}$ {{< /math >}}, gradients propagate back through the kinetic equations (Eq. 3, 4, 9, 10, etc.) to {{< math >}} $\alpha_{g1}, \beta_g, \gamma_g$ {{< /math >}}.

The {{< math >}} $L_{\text{switch}}$ {{< /math >}} term also provides direct gradients to {{< math >}} $\alpha_{g1}, \beta_g, \gamma_g, t_g^s$ {{< /math >}} as it explicitly depends on them.

*Example:* {{< math >}} $\frac{\partial L_{\text{velo}}}{\partial \alpha_{A1}}$ {{< /math >}} will be a large positive value, indicating {{< math >}} $\alpha_{A1}$ {{< /math >}} needs to decrease to reduce the high {{< math >}} $\bar{u}^{(A)}$ {{< /math >}} prediction for Cell 1, Gene A.

**Gradients to Latent Times** ({{< math >}} $t_{ng}^{(k)}$ {{< /math >}}):

Gradients from {{< math >}} $\bar{u}^{(g)}, \bar{s}^{(g)}$ {{< /math >}} also flow back to the latent times {{< math >}} $(t_1)_{ng}$ {{< /math >}} and {{< math >}} $(t_3)_{ng}$ {{< /math >}}.

**Gradients to Decoder Network Parameters** ({{< math >}} $\phi$ {{< /math >}}):

From the latent times, gradients flow back through the decoder networks ({{< math >}} $h_{\text{ind}}, h_{\text{rep}}$ {{< /math >}}) to update their internal weights and biases (part of {{< math >}} $\phi$ {{< /math >}}). This teaches the decoder how to map {{< math >}} $z_n^*$ {{< /math >}} to meaningful latent times that best fit the data.

**Gradients to Sampled Latent Variable** ({{< math >}} $z_n^*$ {{< /math >}}):

Gradients from the decoder (and from the {{< math >}} $\text{KL}(q_\phi(\pi_{ng} | z_n) \| p(\pi_{ng}))$ {{< /math >}} term) flow back to {{< math >}} $z_n^*$ {{< /math >}}.

**Gradients to Encoder Network Outputs** ({{< math >}} $\mu_{z_n}, \log\Sigma_{z_n}$ {{< /math >}}):

This is where the reparameterization trick is vital. The gradients from {{< math >}} $z_n^*$ {{< /math >}} are propagated back to {{< math >}} $\mu_{z_n}$ {{< /math >}} and {{< math >}} $\log\Sigma_{z_n}$ {{< /math >}} (the parameters defining {{< math >}} $q_\phi(z_n | u_n, s_n)$ {{< /math >}}).

{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \mu_{z_n}} = \frac{\partial L_{\text{velo}}}{\partial z_n^*} \times 1$$ {{< /math >}}

{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \log\Sigma_{z_n}} = \frac{\partial L_{\text{velo}}}{\partial z_n^*} \times \frac{\partial z_n^*}{\partial \log\Sigma_{z_n}} = \frac{\partial L_{\text{velo}}}{\partial z_n^*} \times 0.5 \times \exp(0.5\log\Sigma_{z_n}) \times \epsilon_z$$ {{< /math >}}

Additionally, the {{< math >}} $\text{KL}(q_\phi(z_n | u_n, s_n) \| p(z))$ {{< /math >}} term directly provides gradients to {{< math >}} $\mu_{z_n}$ {{< /math >}} and {{< math >}} $\log\Sigma_{z_n}$ {{< /math >}}, pulling them towards 0 and {{< math >}} $I_d$ {{< /math >}} respectively.

**Gradients to Encoder Network Parameters** ({{< math >}} $\phi$ {{< /math >}}):

Finally, the gradients for {{< math >}} $\mu_{z_n}, \log\Sigma_{z_n}$ {{< /math >}} (and the Dirichlet parameters for {{< math >}} $\pi_{ng}$ {{< /math >}}) flow back through the encoder neural network to update its internal weights and biases (part of {{< math >}} $\phi$ {{< /math >}}). This teaches the encoder how to map raw {{< math >}} $(u_n, s_n)$ {{< /math >}} data to a latent space ({{< math >}} $z_n$ {{< /math >}}) and state probabilities ({{< math >}} $\pi_{ng}$ {{< /math >}}) that are useful for modeling.

### c. Parameter Update

All computed gradients for the mini-batch are summed.
The optimizer (Adam) updates all parameters using the learning rate: {{< math >}} $\text{Parameter}_{\text{new}} = \text{Parameter}_{\text{old}} - \text{learning\_rate} \times \frac{\partial L_{\text{velo}}}{\partial \text{Parameter}}$ {{< /math >}}

## 4. Conceptual Evolution of Parameters Over Training

This single iteration is just one small step. The entire process repeats for many mini-batches and epochs (hundreds to thousands of times). Over time, all parameters will adapt:

**Biophysical Parameters** ({{< math >}} $\alpha, \beta, \gamma, t_s$ {{< /math >}}): Will converge to gene-specific values that accurately describe the transcription, splicing, degradation, and switching kinetics, fitting the observed data and satisfying the switch point penalty.
- {{< math >}} $\alpha_{A1}$ {{< /math >}} will likely decrease from 0.85, as it was overpredicting.
- {{< math >}} $\beta_A, \gamma_A$ {{< /math >}} will adjust to shape the kinetic curves.
- {{< math >}} $t_A^s$ {{< /math >}} will shift so that {{< math >}} $\bar{u}^{(A)}(t_A^s, k=2)$ {{< /math >}} and {{< math >}} $\bar{s}^{(A)}(t_A^s, k=2)$ {{< /math >}} match {{< math >}} $u_A^*$ {{< /math >}} and {{< math >}} $s_A^*$ {{< /math >}}.

**Variance Parameters** ({{< math >}} $\sigma_{gu}, \sigma_{gs}$ {{< /math >}}): Will adjust to reflect the appropriate level of noise and variability for each gene. They will likely increase initially from 0.1 because current predictions are very far off, so a larger variance would lead to a higher (less negative) likelihood. As the predictions improve, they might stabilize at values reflecting the true biological and technical noise.

**Neural Network Parameters** ({{< math >}} $\phi$ {{< /math >}} - Encoder/Decoder weights/biases):
- **Encoder:** Learns to map observed {{< math >}} $(u_n, s_n)$ {{< /math >}} data to a meaningful latent space ({{< math >}} $z_n$ {{< /math >}}) and accurate state probabilities ({{< math >}} $\pi_{ng}$ {{< /math >}}), effectively compressing the high-dimensional data into low-dimensional biological insights. It is constrained by the KL terms to avoid overfitting and ensure a regularized latent space.
- **Decoder:** Learns to map the latent variable {{< math >}} $z_n$ {{< /math >}} to specific latent times that, in combination with the biophysical parameters, accurately reconstruct the observed RNA dynamics.