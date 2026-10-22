# Gauss-Newton Approximation: Linear vs Nonlinear Drives

## Background: what Ipopt needs

The NLP has variables $z = [\tilde{\psi}_1, \ldots, \tilde{\psi}_N, u_1, \ldots, u_N]$ and equality constraints (dynamics):

$$F_k(z) = \tilde{\psi}_{k+1} - \Phi_k(u_k, u_{k+1})\, \tilde{\psi}_k = 0, \quad k = 1, \ldots, N-1$$

The **Lagrangian** is:

$$\mathcal{L}(z, \lambda) = J(z) + \sum_k \lambda_k^\top F_k(z)$$◊

Ipopt needs $\nabla^2 \mathcal{L}$. The relevant block is the control-control Hessian:

$$\frac{\partial^2 \mathcal{L}}{\partial u_i\, \partial u_j} = \underbrace{\frac{\partial^2 J}{\partial u_i\, \partial u_j}}_{\text{objective curvature}} + \underbrace{\sum_k \lambda_k^\top \frac{\partial^2 F_k}{\partial u_i\, \partial u_j}}_{\text{constraint curvature}}$$

The **constraint curvature** term expands as:

$$\lambda_k^\top \frac{\partial^2 F_k}{\partial u_i\, \partial u_j} = -\lambda_k^\top \frac{\partial^2 \Phi_k}{\partial u_i\, \partial u_j}\, \tilde{\psi}_k$$

Computing $\partial^2 \Phi_k / \partial u_i\, \partial u_j$ exactly requires solving a second-order variational ODE — expensive. **Gauss-Newton (GN) drops this entire block**, keeping only the $(u, \psi)$ cross-terms from the Jacobian product $J^\top J$.

---

## The Hamiltonian structure

The Hamiltonian is built from drive operators:

$$H(u) = H_\text{drift} + \sum_d c_d(u)\, G_d$$

where $c_d(u)$ is the **drive coefficient function** and $G_d$ is the (fixed) operator matrix.

The second derivative of the propagator satisfies a variational equation whose forcing includes:

$$\frac{\partial^2 H}{\partial u_i\, \partial u_j} = \sum_d \frac{\partial^2 c_d}{\partial u_i\, \partial u_j}\, G_d$$

This is the key quantity. Whether GN is valid depends entirely on whether the **coefficient Hessian** $\partial^2 c_d / \partial u^2$ is zero or not.

---

## Why GN is valid for LinearDrive

For a `LinearDrive`, the coefficient is affine in the controls:

$$c_d(u) = u_{d.\text{index}}$$

Therefore:

$$\frac{\partial c_d}{\partial u_j} = \delta_{j,\, d.\text{index}} \quad \text{(constant)}, \qquad \frac{\partial^2 c_d}{\partial u_i\, \partial u_j} = 0 \quad \text{(identically zero)}$$

So $\partial^2 H / \partial u^2 = 0$. The forcing in the second-order variational ODE vanishes. The only remaining second-order terms come from **compositions of first-order sensitivities** (products $s_i \cdot s_j$), which are multiplied by $\lambda_k$ and thus small when the dynamics constraints are approximately satisfied. This is the standard regime where GN is trusted.

> **Note:** GN is not *exactly* correct for linear drives — there are still $\mathcal{O}(\|\Phi - I\|^2)$ terms from the compositions — but they vanish at feasible points and are small for short timesteps and small controls. GN is a well-justified approximation in this regime.

---

## Why GN is wrong for NonlinearDrive

The three NonlinearDrive operators in `TransmonCavitySystemDispNL` have coefficients:

| Drive operator | Coefficient $c_d(u)$ | $\partial^2 c_d/\partial \alpha_I^2$ | $\partial^2 c_d/\partial \alpha_Q^2$ | $\partial^2 c_d/\partial \alpha_I\, \partial \alpha_Q$ |
|---|---|---|---|---|
| $|\alpha|^2$ | $\alpha_I^2 + \alpha_Q^2$ | **2** | **2** | 0 |
| $\alpha_I \|\alpha\|^2$ | $\alpha_I^3 + \alpha_I \alpha_Q^2$ | $6\alpha_I$ | $2\alpha_I$ | $2\alpha_Q$ |
| $\alpha_Q \|\alpha\|^2$ | $\alpha_I^2 \alpha_Q + \alpha_Q^3$ | $2\alpha_Q$ | $6\alpha_Q$ | $2\alpha_I$ |

The $|\alpha|^2$ coefficient has **constant curvature equal to 2**, regardless of the operating point. It contributes a constant forcing term to the second-order variational ODE:

$$\frac{\partial^2 H}{\partial \alpha_I^2} = 2\, G_{|\alpha|^2}$$

where $G_{|\alpha|^2}$ is the operator $K_c \cdot \hat{n}_c(\hat{n}_c - 1)/2$ (Kerr-like, from the dispersive Hamiltonian). With $K_c = -1.259\,\text{MHz}$ and $\dim_c = 10$, the operator norm is $\mathcal{O}(1)$ — not small.

**The GN approximation ignores this constant curvature throughout the entire optimization**, not just near the solution. Ipopt receives a systematically wrong picture of the landscape in the $\alpha$-control directions, causing it to overshoot and converge poorly.

---

## Concrete picture: the $\alpha$ landscape

The infidelity objective as a function of $\alpha_I$ (with everything else fixed) looks like:

$$J(\alpha_I) \approx a\, \alpha_I + \underbrace{b\, \alpha_I^2}_{\substack{\uparrow \\ \text{curvature term} \\ \text{GN drops this}}} + \text{(higher-order terms from propagation)}$$

The contribution $b \propto K_c \cdot \|G_{|\alpha|^2}\|$ is **constant**. GN tells Ipopt that the curvature in this direction is less than it actually is $\Rightarrow$ Ipopt takes steps that are too large in the $\alpha$ directions $\Rightarrow$ overshooting $\Rightarrow$ poor convergence.

---

## Why HermitianExponentialIntegrator is unaffected

HermitianExp computes $\Phi_k = \exp(-i H(u_k)\, \Delta t)$ with $u_k$ **piecewise-constant**. Its analytic Hessian implementation includes all second-order terms exactly — it does not use GN at all for the $(u, u)$ block.

Additionally, since controls are piecewise-constant, interval $k$ depends **only on $u_k$**, not on $u_{k+1}$. Adjacent intervals are decoupled: the Hessian is block-diagonal, with each $4 \times 4$ block containing the full exact second-order information including the NonlinearDrive curvature. Ipopt inverts this exactly.

SplineODE couples $u_k$ and $u_{k+1}$ through the ODE integral (linear interpolation within each interval), creating off-diagonal blocks in the Hessian. GN drops these coupled second-order terms, compounding the error from the NonlinearDrive curvature.

---

## Summary

| | LinearDrive | NonlinearDrive |
|---|---|---|
| Coefficient $c_d(u)$ | linear (affine) | quadratic / cubic |
| $\partial^2 c_d / \partial u^2$ | $0$ (exactly) | nonzero (e.g. $2$ for $|\alpha|^2$) |
| $\partial^2 H / \partial u^2$ | $0$ | $= \sum_d (\partial^2 c_d/\partial u^2)\, G_d \neq 0$ |
| GN error (SplineODE) | small near solution | **constant, $\mathcal{O}(K_c \cdot \|\hat{n}\|)$, present everywhere** |
| HermitianExp affected? | — | No (uses exact Hessian, no GN) |

---

## What the L-BFGS test tells us

Running SplineODE with `eval_hessian=false` switches Ipopt to L-BFGS: the Hessian is estimated from gradient secant pairs rather than GN. L-BFGS implicitly captures the NonlinearDrive curvature because it sees the gradient change when stepping in the $\alpha$ directions.

- If **L-BFGS $\gg$ GN** in fidelity: the inaccurate GN Hessian is the bottleneck.
  Fix: implement exact second-order terms for SplineODE (analogous to HermitianExp).
- If **L-BFGS $\approx$ GN**: something else limits convergence (landscape, ODE tolerance,
  step coupling between knots). Need to investigate further.

Results are in `notes/benchmark_displaced_frame.md` (Solver comparison section).
