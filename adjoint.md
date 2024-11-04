https://sashalsey.github.io/adjoint
# Optimization Problem

## Objective Function
\[
\min_{\rho} \; \Phi(\rho) = F^T U
\]

## Constraints
### Volume constraint:
\[
\sum_{e=1}^{N} \rho_e v_e = v^T \rho \leq V^*
\]

### General constraint:
\[
g_1(\rho, U) \leq g^*
\]

### Bounds on density:
\[
0 < \rho \leq \rho_{\text{min}} \leq 1
\]

---

## State Equation
\[
K(\rho) U = F
\]

where
\[
K(\rho) = \sum_{e=1}^{N} K_e (\rho_e) = \sum_{e=1}^{N} \rho_e^p K_e^0
\]

---

## Sensitivity Analysis
Define:
\[
\Phi = \Phi(U(\rho))
\]

Then the derivative with respect to \( \rho \) is:
\[
\Phi' = \frac{\partial \Phi}{\partial \rho} = \frac{\partial \Phi}{\partial U} \frac{\partial U}{\partial \rho}
\]

---

## Using the Adjoint Method to Solve
Define:
\[
\hat{\Phi} = \Phi + \lambda^T (KU - F)
\]

Then:
\[
\hat{\Phi}' = \frac{\partial \Phi}{\partial U} U' + \lambda^T (K'U + KU')
\]

---

## Finding \( \lambda \) such that \( U' \) terms cancel
Set:
\[
\left( \frac{\partial \Phi}{\partial U} + \lambda^T K \right) U' = 0
\]
which implies:
\[
\frac{\partial \Phi}{\partial U} + \lambda^T K = 0
\]

So we get:
\[
K^T \lambda = -\frac{\partial \Phi^T}{\partial U}
\]

---

## Back to the Original Problem
Since:
\[
\Phi = F^T U
\]
then:
\[
\frac{\partial \Phi}{\partial U} = F^T
\]

Thus:
\[
\frac{\partial \Phi^T}{\partial U} = F = -K \lambda
\]
and therefore:
\[
\lambda = -U
\]
