# Mathematical Derivation of 6DoF GPS-IMU Fusion


This document provides detailed mathematical derivations for the Error State Iterated Kalman Filter (ESIKF) used in GPS-IMU fusion.

---

## Part I. Foundations

### 1. Algorithm Overview

The filter runs in a **predict–update** loop:

```
┌─────────────────────────────────────────────────────────────────┐
│  INITIALIZATION (Section 8)                                      │
│  origin, R₀, v₀, P₀ from first GNSS                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  PREDICTION (Section 5 & 6)                                      │
│  IMU → propagate state (Section 5) + covariance (Section 6)      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  UPDATE (Section 7)                                               │
│  GNSS → residual (Section 7.1) → iterated EKF (Section 7.7)      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              └──────────► (repeat)
```

---

### 2. Notation and Prerequisites

#### 2.1 SO(3) Lie Algebra

The **hat operator** maps $\boldsymbol{\phi} \in \mathbb{R}^3$ to $\mathfrak{so}(3)$:
$$
[\boldsymbol{\phi}]_\times = \boldsymbol{\phi}^\wedge = 
\begin{bmatrix} 0 & -\phi_z & \phi_y \\ \phi_z & 0 & -\phi_x \\ -\phi_y & \phi_x & 0 \end{bmatrix}
$$

The **A-matrix** (right Jacobian of SO(3) for exp) relates $\delta\boldsymbol{\phi}$ to rotation perturbation. For $\boldsymbol{\theta} \in \mathbb{R}^3$, $\theta = \|\boldsymbol{\theta}\|$:
$$
\mathbf{A}(\boldsymbol{\theta}) = \mathbf{I} + \frac{1-\cos\theta}{\theta^2} [\boldsymbol{\theta}]_\times + \frac{1 - \frac{\sin\theta}{\theta}}{\theta^2} [\boldsymbol{\theta}]_\times^2
$$
When $\theta \to 0$: $\mathbf{A}(\boldsymbol{\theta}) \to \mathbf{I}$.

The **inverse right Jacobian** $\mathbf{J}_r^{-1}(\boldsymbol{\phi})$ maps from tangent space of $\log(R)$ to rotation perturbation. For $\phi = \|\boldsymbol{\phi}\|$:
$$
\mathbf{J}_r^{-1}(\boldsymbol{\phi}) = \mathbf{I} + \frac{1}{2}[\boldsymbol{\phi}]_\times + \left(1 - \frac{\frac{\theta}{2}(1+\cos\theta)}{\sin\theta}\right) \left(\frac{\boldsymbol{\phi}\boldsymbol{\phi}^\top}{\|\boldsymbol{\phi}\|^2} - \mathbf{I}\right)
$$
When $\|\boldsymbol{\phi}\| \le \epsilon$: $\mathbf{J}_r^{-1}(\boldsymbol{\phi}) \approx \mathbf{I} + \frac{1}{2}[\boldsymbol{\phi}]_\times$.

#### 2.2 State Vector

State $\mathbf{x}$ has DOF = 26, DIM = 27:

| Index | Variable | Type | DOF |
|-------|----------|------|-----|
| 0–2   | $\mathbf{p}$ | $\mathbb{R}^3$ | 3 |
| 3–5   | $\mathbf{R}$ (SO3) | rotation | 3 |
| 6–8   | $\mathbf{v}$ | $\mathbb{R}^3$ | 3 |
| 9–11  | $\mathbf{b}_g$ | $\mathbb{R}^3$ | 3 |
| 12–14 | $\mathbf{b}_a$ | $\mathbb{R}^3$ | 3 |
| 15–17 | $\mathbf{t}_{gnss}$ | $\mathbb{R}^3$ | 3 |
| 18–20 | $\mathbf{R}_{gnss}$ (SO3) | rotation | 3 |
| 21–23 | $\mathbf{R}_{car}$ (SO3) | rotation | 3 |
| 24–25 | $\mathbf{g}$ (S2) | gravity dir | 2 |

---

### 3. Error State Formulation

#### 3.1 Why Error State?

The **Error State Kalman Filter (ESKF)** estimates a small error $\delta\mathbf{x}$ in the tangent space rather than the full state $\mathbf{x}$ on the manifold. Advantages:

1. **Unconstrained optimization**: Error states live in $\mathbb{R}^n$, avoiding manifold constraints during the update.
2. **Better linearization**: Small errors yield more accurate first-order approximations.
3. **Consistent covariance**: Covariance stays in the tangent space; no singularities from over-parameterization (e.g., quaternion normalization).

The nominal state is propagated in closed form; the error state is reset after each update.

#### 3.2 Manifold Operations

**Boxplus** $\oplus$: inject a tangent vector into the manifold.
- $\mathbb{R}^3$: $\mathbf{p} \oplus \delta\mathbf{p} = \mathbf{p} + \delta\mathbf{p}$
- SO3: $\mathbf{R} \oplus \boldsymbol{\phi} = \mathbf{R} \cdot \exp(\boldsymbol{\phi})$
- S2: uses exponential map on the sphere

**Boxminus** $\ominus$: extract tangent vector between two manifold points.
- $\delta\mathbf{x} = \mathbf{x}_2 \ominus \mathbf{x}_1$ such that $\mathbf{x}_2 = \mathbf{x}_1 \oplus \delta\mathbf{x}$
- SO3: $\boldsymbol{\phi} = \log(\mathbf{R}_1^{-1} \mathbf{R}_2)$

---

## Part II. Prediction Step

### 4. Continuous-Time Dynamics

The velocity field $\dot{\mathbf{x}} = f(\mathbf{x}, \mathbf{u})$ with IMU input $\mathbf{u} = (\boldsymbol{\omega}, \mathbf{a})$:

$$
\begin{aligned}
\dot{\mathbf{p}} &= \mathbf{v} \\
\dot{\mathbf{R}} &= \mathbf{R} [\boldsymbol{\omega} - \mathbf{b}_g]_\times \quad \text{(or } \dot{\mathbf{R}} = \mathbf{R} \cdot (\boldsymbol{\omega} - \mathbf{b}_g) \text{ in tangent space)} \\
\dot{\mathbf{v}} &= \mathbf{R}(\mathbf{a} - \mathbf{b}_a) + \mathbf{g} \\
\dot{\mathbf{b}}_g &= \mathbf{0}, \quad \dot{\mathbf{b}}_a = \mathbf{0}, \quad \dot{\mathbf{t}}_{gnss} = \mathbf{0}, \quad \dot{\mathbf{R}}_{gnss} = \mathbf{0}, \quad \dot{\mathbf{R}}_{car} = \mathbf{0}, \quad \dot{\mathbf{g}} = \mathbf{0}
\end{aligned}
$$

Let $\tilde{\mathbf{a}} = \mathbf{a} - \mathbf{b}_a$, $\tilde{\boldsymbol{\omega}} = \boldsymbol{\omega} - \mathbf{b}_g$. The flatted output is:
$$
f(\mathbf{x}, \mathbf{u}) = \begin{bmatrix} \mathbf{v} \\ \tilde{\boldsymbol{\omega}} \\ \mathbf{R}\tilde{\mathbf{a}} + \mathbf{g} \\ \mathbf{0} \\ \ldots \\ \mathbf{0} \end{bmatrix}
$$

---

### 5. Discrete-Time State Propagation

The continuous dynamics are integrated over $\Delta t$. Three schemes are used:

**Linear (Euler):**
$$
\mathbf{R}_{k+1} = \mathbf{R}_k \exp(\Delta t \cdot \tilde{\boldsymbol{\omega}}), \quad
\mathbf{v}_{k+1} = \mathbf{v}_k + \Delta t \cdot (\mathbf{R}_k \tilde{\mathbf{a}} + \mathbf{g})
$$
$$
\mathbf{p}_{k+1} = \mathbf{p}_k + \Delta t \cdot \mathbf{v}_k + \frac{\Delta t^2}{2} (\mathbf{R}_k \tilde{\mathbf{a}} + \mathbf{g})
$$

**Midpoint:** Uses $\tilde{\boldsymbol{\omega}}_{mid} = \frac{1}{2}(\boldsymbol{\omega}_{k} + \boldsymbol{\omega}_{k+1}) - \mathbf{b}_g$ and acceleration at $t + \Delta t/2$.

**Runge-Kutta 4th Order:** Standard RK4 for higher accuracy when $\Delta t$ is not small.

---

### 6. Covariance Propagation

#### 6.1 State Jacobian $\mathbf{F} = \frac{\partial f}{\partial \mathbf{x}}$

Only non-zero blocks are shown. $\mathbf{F}$ maps from tangent space to tangent space.

| Block | Expression |
|-------|------------|
| $\frac{\partial \dot{\mathbf{p}}}{\partial \mathbf{v}}$ | $\mathbf{I}_{3\times3}$ |
| $\frac{\partial \dot{\mathbf{v}}}{\partial \boldsymbol{\phi}^R}$ | $-\mathbf{R} [\tilde{\mathbf{a}}]_\times$ |
| $\frac{\partial \dot{\mathbf{v}}}{\partial \mathbf{b}_a}$ | $-\mathbf{R}$ |
| $\frac{\partial \dot{\mathbf{v}}}{\partial \mathbf{g}}$ | $\mathbf{M}_{g}$ (3×2 S2 tangent) |
| $\frac{\partial \dot{\boldsymbol{\phi}}^R}{\partial \mathbf{b}_g}$ | $-\mathbf{I}_{3\times3}$ |

#### 6.2 Process Noise Jacobian $\mathbf{G} = \frac{\partial f}{\partial \mathbf{w}}$

Process noise $\mathbf{w} \in \mathbb{R}^{12}$: $[\mathbf{n}_g, \mathbf{n}_a, \mathbf{n}_{bg}, \mathbf{n}_{ba}]$.

$$
\mathbf{G} = \begin{bmatrix}
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0} \\
-\mathbf{I} & \mathbf{0} & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & -\mathbf{R} & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{I} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I} \\
\vdots & \vdots & \vdots & \vdots
\end{bmatrix}
$$

#### 6.3 Process Noise Model $\mathbf{Q}$

$$
\mathbf{Q} = \text{diag}(\sigma_{ng}^2 \mathbf{I}_3, \sigma_{na}^2 \mathbf{I}_3, \sigma_{nbg}^2 \mathbf{I}_3, \sigma_{nba}^2 \mathbf{I}_3)
$$

| Parameter | Symbol | Typical (phone) | Description |
|-----------|--------|-----------------|-------------|
| Gyro noise | $\sigma_{ng}$ | ~1e-2 rad/s/√Hz | White noise on angular velocity |
| Accel noise | $\sigma_{na}$ | ~1e-1 m/s²/√Hz | White noise on linear acceleration |
| Gyro bias RW | $\sigma_{nbg}$ | ~1e-5 rad/s/√Hz | Random walk of gyro bias |
| Accel bias RW | $\sigma_{nba}$ | ~1e-4 m/s²/√Hz | Random walk of accel bias |

#### 6.4 Discrete-Time Covariance Update

For small $\Delta t$: $\mathbf{F}_1 = \mathbf{I} + \Delta t \cdot \mathbf{F}$.

For **SO3** states, the transition uses $\mathbf{A}(-\Delta t \cdot \tilde{\boldsymbol{\omega}})^\top$ to map rotation error. For **S2** (gravity), the transition uses $N_x \cdot R_{exp} \cdot M_x$ between tangent spaces.

$$
\mathbf{P}_{k+1|k} = \mathbf{F}_1 \mathbf{P}_{k|k} \mathbf{F}_1^\top + (\Delta t \cdot \mathbf{G}) \mathbf{Q} (\Delta t \cdot \mathbf{G})^\top
$$

---

## Part III. Update Step

### 7. Measurement Models and Residuals

#### 7.1 Residual Formulation

For **Euclidean** measurements (position, velocity): $\mathbf{y} = \mathbf{z} - \mathbf{h}(\mathbf{x})$.

For **rotation** on SO3, the residual must respect the manifold:

**Proper form (used in code):**
$$
\mathbf{y}_{rot} = \log\left( R_{obs} \cdot R_{pred}^{-1} \right)
$$

**Pose (position + rotation):**
$$
\mathbf{y}_{pose} = \begin{bmatrix} \mathbf{z}_p - \mathbf{h}_p(\mathbf{x}) \\ \log(R_{obs} \cdot R_{pred}^{-1}) \end{bmatrix}
$$

#### 7.2 Position Measurement (GNSS)

**Model:** $\mathbf{h}_{pos}(\mathbf{x}) = \mathbf{p} + \mathbf{R} \mathbf{t}_{gnss}$

**Jacobian:** $\frac{\partial \mathbf{h}_{pos}}{\partial \mathbf{p}} = \mathbf{I}$, $\frac{\partial \mathbf{h}_{pos}}{\partial \mathbf{t}_{gnss}} = \mathbf{R}$, $\frac{\partial \mathbf{h}_{pos}}{\partial \boldsymbol{\phi}^R} = \mathbf{R} [-\mathbf{t}_{gnss}]_\times$

**Derivation:** $\mathbf{p} + \mathbf{R} \exp(\boldsymbol{\phi}) \mathbf{t}_{gnss} \approx \mathbf{p} + \mathbf{R}(\mathbf{I} + [\boldsymbol{\phi}]_\times) \mathbf{t}_{gnss}$, hence $\frac{\partial}{\partial \boldsymbol{\phi}} \mathbf{R} [\boldsymbol{\phi}]_\times \mathbf{t}_{gnss} = \mathbf{R} [-\mathbf{t}_{gnss}]_\times$.

#### 7.3 Rotation Measurement (GNSS Heading)

**Model:** $\mathbf{h}_{rot}(\mathbf{x}) = \log(\mathbf{R} \cdot \mathbf{R}_{gnss})$

**Jacobian:** Let $\boldsymbol{\phi}_{hat} = \log(\mathbf{R} \mathbf{R}_{gnss})$. Then:
$$
\frac{\partial \log(\mathbf{R} \mathbf{R}_{gnss})}{\partial \boldsymbol{\phi}^R} = \mathbf{J}_r^{-1}(\boldsymbol{\phi}_{hat}) \cdot \mathbf{R}_{gnss}^\top, \quad
\frac{\partial \log(\mathbf{R} \mathbf{R}_{gnss})}{\partial \boldsymbol{\phi}^{R_{gnss}}} = \mathbf{J}_r^{-1}(\boldsymbol{\phi}_{hat})
$$

#### 7.4 Velocity Measurement (Zero-Velocity / Chassis)

**Model:** $\mathbf{h}_{vel}(\mathbf{x}) = \mathbf{R}_{car} (\mathbf{R}^{-1} \mathbf{v})$

**Jacobian:** $\frac{\partial \mathbf{h}_{vel}}{\partial \boldsymbol{\phi}^R} = \mathbf{R}_{car} [\mathbf{R}^{-1} \mathbf{v}]_\times$, $\frac{\partial \mathbf{h}_{vel}}{\partial \mathbf{v}} = \mathbf{R}_{car} \mathbf{R}^{-1}$, $\frac{\partial \mathbf{h}_{vel}}{\partial \boldsymbol{\phi}^{R_{car}}} = -\mathbf{R}_{car} [\mathbf{R}^{-1} \mathbf{v}]_\times$

#### 7.5 Velocity Squared Norm

**Model:** $h_{vel2}(\mathbf{x}) = \|\mathbf{v}\|^2$. **Jacobian:** $\frac{\partial h_{vel2}}{\partial \mathbf{v}} = 2 \mathbf{v}^\top$

#### 7.6 Accelerometer Bias (Zero-Acceleration)

**Model:** $\mathbf{h}_{acc}(\mathbf{x}) = \mathbf{R}(\mathbf{a}_{imu} - \mathbf{b}_a) + \mathbf{g}$ (expected zero when stationary)

**Jacobian:** $\frac{\partial \mathbf{h}_{acc}}{\partial \boldsymbol{\phi}^R} = -\mathbf{R} [\mathbf{a}_{imu} - \mathbf{b}_a]_\times$, $\frac{\partial \mathbf{h}_{acc}}{\partial \mathbf{b}_a} = -\mathbf{R}$, $\frac{\partial \mathbf{h}_{acc}}{\partial \mathbf{g}} = \mathbf{M}_g$

#### 7.7 Iterated EKF Update

For measurement $\mathbf{z}$ with model $\mathbf{h}(\mathbf{x})$ and noise covariance $\mathbf{R}$:

1. **Residual:** $\mathbf{y} = \mathbf{z} - \mathbf{h}(\mathbf{x})$ (or manifold residual for rotation)
2. **Predicted covariance:** $\mathbf{S} = \mathbf{H} \mathbf{P} \mathbf{H}^\top + \mathbf{R}$
3. **Kalman gain:** $\mathbf{K} = \mathbf{P} \mathbf{H}^\top \mathbf{S}^{-1}$
4. **State update:** $\mathbf{x} \leftarrow \mathbf{x} \oplus (\mathbf{K} \mathbf{y})$
5. **Covariance update:** $\mathbf{P} \leftarrow (\mathbf{I} - \mathbf{K} \mathbf{H}) \mathbf{P}$

The **iterated** update repeats with the current $\mathbf{x}$ and $\mathbf{H}$ until $\|\delta\mathbf{x}_i\| < \epsilon_i$ (typical $\epsilon \sim 10^{-4}$) or max iterations (5–10).

**Jacobian Summary Table:**

| Measurement | $h(\mathbf{x})$ | $\frac{\partial h}{\partial \mathbf{p}}$ | $\frac{\partial h}{\partial \mathbf{R}}$ | $\frac{\partial h}{\partial \mathbf{v}}$ | $\frac{\partial h}{\partial \mathbf{b}_a}$ | $\frac{\partial h}{\partial \mathbf{g}}$ |
|-------------|-----------------|-------------------------------------------|-------------------------------------------|-------------------------------------------|---------------------------------------------|------------------------------------------|
| Position | $\mathbf{p} + \mathbf{R}\mathbf{t}_{gnss}$ | $\mathbf{I}$ | $\mathbf{R}[-\mathbf{t}_{gnss}]_\times$ | $\mathbf{0}$ | $\mathbf{0}$ | $\mathbf{0}$ |
| Rotation | $\log(\mathbf{R}\mathbf{R}_{gnss})$ | $\mathbf{0}$ | $J_r^{-1}[\phi] R_{gnss}^\top$ | $\mathbf{0}$ | $\mathbf{0}$ | $\mathbf{0}$ |
| Velocity | $\mathbf{R}_{car}\mathbf{R}^{-1}\mathbf{v}$ | $\mathbf{0}$ | $\mathbf{R}_{car}[\mathbf{R}^{-1}\mathbf{v}]_\times$ | $\mathbf{R}_{car}\mathbf{R}^{-1}$ | $\mathbf{0}$ | $\mathbf{0}$ |
| Vel²norm | $\|\mathbf{v}\|^2$ | $\mathbf{0}$ | $\mathbf{0}$ | $2\mathbf{v}^\top$ | $\mathbf{0}$ | $\mathbf{0}$ |
| Acc bias | $\mathbf{R}(\mathbf{a}-\mathbf{b}_a)+\mathbf{g}$ | $\mathbf{0}$ | $-\mathbf{R}[\tilde{\mathbf{a}}]_\times$ | $\mathbf{0}$ | $-\mathbf{R}$ | $\mathbf{M}_g$ |

---

## Part IV. Initialization and Supplementary

### 8. Initialization

1. **Origin**: First GNSS position → `origin` in UTM.
2. **Rotation**: From GNSS heading (and optionally pitch): $\mathbf{R}_0 = R_z(\psi) \cdot R_y(\theta)$.
3. **Velocity**: $\mathbf{v}_0 = \mathbf{R}_0 \cdot \mathbf{v}_{gnss}^{body}$.
4. **Covariance**: Diagonal blocks, e.g.:
   - Position: $\sim 1e-8$
   - Rotation: yaw large (1e8), roll/pitch small (1e-5)
   - Velocity: lateral large (1e4), forward small (1e-6)
   - Biases: $\sim 1e-6$
   - Extrinsics: small (1e-8)

---

### 9. Observability

| State | Observable from | Condition |
|-------|-----------------|-----------|
| Position | GNSS | — |
| Velocity | GNSS velocity / ZUPT | Moving / stationary |
| Rotation (yaw) | GNSS heading | Moving |
| Rotation (roll/pitch) | Gravity / GNSS pitch | Prior or pitch available |
| Gyro bias $\mathbf{b}_g$ | — | Requires rotational motion |
| Accel bias $\mathbf{b}_a$ | ZUPT / position/velocity | Stationary or motion |
| Gravity $\mathbf{g}$ | Accelerometer | Orientation known |

---
