# Theoretical Framework: Explicit Transport of Excited Species in Non-Equilibrium Plasma Simulations

## Table of Contents
1. [Introduction](#introduction)
2. [The Advection-Diffusion-Reaction (ADR) System](#adr-system)
3. [Reaction Mechanisms](#reaction-mechanisms)
4. [Numerical Discretization: Scharfetter-Gummel Scheme](#numerical-discretization)
5. [Boundary Conditions](#boundary-conditions)
6. [Solver Architecture](#solver-architecture)
7. [References](#references)

---

## 1. Introduction <a name="introduction"></a>

### 1.1 The Necessity of Explicit Excited Species Transport

In low-temperature, non-equilibrium plasmas (dielectric barrier discharges, atmospheric pressure plasma jets, RF capacitive discharges), traditional "effective medium" approximations that lump excited neutral species into the background gas fail to capture critical physics. This implementation addresses this by treating excited states (metastables, vibrationally excited molecules, resonant states) as **distinct transport scalars** with separate continuity equations.

**Key Physics Requiring Explicit Treatment:**
- **Stepwise ionization**: Two-step "ladder" process with lower energy threshold
- **Penning ionization**: Chemo-ionization where excited state A* ionizes target B
- **Superelastic collisions**: Excited states heat electrons in afterglows (E ≈ 0)
- **Surface quenching**: Excited neutrals undergo complex wall interactions
- **Spatial transport of ionization potential**: Metastables carry energy over macroscopic distances

### 1.2 Comparison: Charged vs. Neutral Species Transport

| Property | Electrons/Ions | Excited Neutrals | Ground Neutrals |
|----------|---------------|------------------|-----------------|
| Electric field response | **Drift** (μE dominant) | None (q=0) | None |
| Advection | Negligible | **Dominant** (Pe >> 1) | Defines frame |
| Diffusion | High (electrons) | Moderate (D ≈ D_gas) | Self-diffusion |
| Wall interaction | Sheath/Recombination | Quenching/SEE | Pressure/Flow |
| Governing equation | Drift-Diffusion | **Advection-Diffusion** | Navier-Stokes |

---

## 2. The Advection-Diffusion-Reaction (ADR) System <a name="adr-system"></a>

### 2.1 Governing Equation

For excited species *k* with number density n_k(x,t):

```
∂n_k/∂t + ∇·Γ_k = S_k
```

where the **flux** is:

```
Γ_k = n_k·u_gas - D_k·∇n_k
       ↑           ↑
   Advection   Diffusion
```

**Key Parameters:**
- `u_gas(x,t)`: Background gas velocity field [m/s]
- `D_k`: Binary diffusion coefficient of excited species in background [m²/s]
  - Typical value: D_Ar* ≈ 1.5×10⁻⁴ m²/s at STP
  - Calculated via Chapman-Enskog theory or measured

### 2.2 The Péclet Number

The **Péclet number** characterizes advection vs. diffusion:

```
Pe = (L·u_gas) / D_k
```

**Regimes:**
- Pe << 1: Diffusion-dominated (standard central differences stable)
- Pe >> 1: Advection-dominated (requires upwinding or exponential fitting)

**Example (Plasma Jet):**
- L = 1 mm, u_gas = 50 m/s, D* = 10⁻⁴ m²/s
- Pe ≈ 500 → **Highly advection-dominated**

This necessitates the **Scharfetter-Gummel scheme** (see Section 4).

### 2.3 Source Term Structure

The source term couples excited species to electrons, ions, and other neutrals:

```
S_k = Σ_r (ν"_k,r - ν'_k,r) · R_r
```

where R_r is the rate of reaction *r* [m⁻³ s⁻¹]. See Section 3 for specific mechanisms.

---

## 3. Reaction Mechanisms <a name="reaction-mechanisms"></a>

This implementation includes **9 reaction types** critical for non-equilibrium plasma physics:

### 3.1 Direct Ionization
**Process:** e + A → 2e + A⁺

**Rate:**
```
R_iz = k_iz(Te) · ne · N_gas
```

**Energy cost:** E_iz = 15.76 eV (Argon)

**Effect:** Primary electron generation, dominant at high E/N

---

### 3.2 Excitation
**Process:** e + A → e + A*

**Rate:**
```
R_exc = k_exc(Te) · ne · N_gas
```

**Energy cost:** E_exc ~ 11.55 eV (Ar metastable)

**Effect:** Populates excited state, enables stepwise processes

---

### 3.3 Stepwise Ionization
**Process:** e + A* → 2e + A⁺

**Rate:**
```
R_step = k_step(Te) · ne · n*
```

**Energy cost:** E_step = E_iz - E_exc ≈ 4.2 eV

**Significance:** This "ladder effect" dramatically lowers the effective ionization threshold. In high-density discharges, stepwise can dominate over direct ionization.

**Without explicit n*:** The model cannot calculate R_step correctly, leading to underestimation of ionization rate and overestimation of required E/N.

---

### 3.4 Penning Ionization
**Process:** A* + B → A + B⁺ + e

**For Argon self-Penning (Hornbeck-Molnar):**
```
Ar* + Ar → Ar + Ar⁺ + e
R_penning = k_penning · n_Ar* · N_Ar
```

**Typical rate coefficient:** k_penning ~ 10⁻¹⁵ m³/s

**Significance:** Non-local ionization. Ar* generated in high-E region can advect downstream (plasma jet) or to dielectric surface (DBD) and ionize there, decoupling ionization from instantaneous E(x).

---

### 3.5 Metastable Pooling
**Process:** A* + A* → A⁺ + A + e

**Rate:**
```
R_pool = k_pool · n*²
```

**Significance:** Memory effect in pulsed discharges. Metastables accumulated during pulse ON phase provide seed electrons during pulse OFF, lowering breakdown voltage of subsequent pulse.

---

### 3.6 Superelastic Collisions
**Process:** e + A* → e(fast) + A

**Rate:**
```
R_superelastic = k_se(Te) · ne · n*
```

**Energy transfer:** Potential energy of A* → kinetic energy of electron

**Significance:** Sustains electron temperature in afterglows where E ≈ 0. Critical for predicting recombination rates and chemical reactivity in post-discharge regions.

---

### 3.7 Radiative Decay
**Process:** A* → A + hν

**Rate:**
```
R_rad = A_rad · n*
```

where A_rad is the Einstein A coefficient [s⁻¹].

**Species dependence:**
- **Metastable states** (e.g., Ar_m): A_rad ≈ 0 (forbidden transition, lifetime ~ ms)
- **Resonant states** (e.g., Ar_r): A_rad ~ 10⁸ s⁻¹ (lifetime ~ 10 ns)

**Effect:** Determines spatial extent of excited species influence. Fast radiative decay confines influence locally; slow decay enables long-range transport.

---

### 3.8 Collisional Quenching
**Process:** A* + M → A + M

**Rate:**
```
R_quench = k_Q · n* · N_gas
```

**Typical:** k_Q ~ 10⁻¹⁶ m³/s

**Effect:** Energy dissipation to gas (not tracked in electron energy equation). Important at high pressure.

---

### 3.9 Three-Body Recombination
**Process:** e + A⁺ + M → A + M

**Rate:**
```
R_rec = k_rec · ne · ni · N_gas
```

**Typical:** k_rec ~ 10⁻³⁹ m⁶/s

**Effect:** Electron sink at high density, energy source (exothermic).

---

## 4. Numerical Discretization: Scharfetter-Gummel Scheme <a name="numerical-discretization"></a>

### 4.1 The Problem with Central Differences

Standard central difference schemes (∂n/∂x ≈ (n_{i+1} - n_{i-1})/(2h)) are **unstable** for advection-dominated problems (Pe > 2), leading to:
- Spurious oscillations
- Unphysical negative densities (n < 0)
- Catastrophic failure of reaction rate calculations

### 4.2 Scharfetter-Gummel Derivation

Consider 1D steady flux between grid points i and i+1:

```
J = -D·dn/dx + u·n
```

**Solution via integrating factor:**

Rearranging: `dn/dx - (u/D)·n = -J/D`

Integrating factor: `exp(-ux/D)`

After integration from x_i to x_{i+1}:

```
J = (D/h) · [B(-α)·n_i - B(α)·n_{i+1}]
```

where:
- `α = (u·h)/D` (local Péclet number)
- `B(x) = x/(exp(x) - 1)` (Bernoulli function)

### 4.3 Bernoulli Function Properties

```c++
double Bernoulli(double x) {
    if (|x| < 1e-4) {
        return 1 - x/2 + x²/12 - x⁴/720;  // Taylor series
    }
    return x / (exp(x) - 1);
}
```

**Asymptotic behavior:**
- B(0) = 1 → Central difference (diffusion limit)
- B(α→∞) ≈ α, B(-α→∞) ≈ 0 → First-order upwind (advection limit)

**Key property:** Automatically transitions between diffusion and advection regimes, ensuring:
1. **Unconditional stability** (no CFL limit)
2. **Positivity preservation** (n ≥ 0)
3. **Second-order accuracy** in diffusion limit

### 4.4 Application to Excited Neutrals

For excited species with `u_gas ≠ 0`:

```
ν = (u_gas · dx) / D*
J_{i→i+1} = (D*/dx) · [B(-ν)·n*_i - B(ν)·n*_{i+1}]
```

This **unifies** the numerical treatment of charged (drift) and neutral (advection) species within a single exponential fitting framework.

---

## 5. Boundary Conditions <a name="boundary-conditions"></a>

### 5.1 Surface Interaction Physics

Excited neutrals do not simply "vanish" at walls. They undergo:
1. **Quenching**: A* + Wall → A + Wall (energy to surface)
2. **Reflection**: Partial reflection back to gas phase
3. **Secondary Emission**: A* + Wall → A + Wall + e⁻ (Auger de-excitation)

### 5.2 Robin-Type Boundary Condition

The appropriate condition is a **flux balance**:

```
Γ* · n̂ = (γ*·v_th*/4) · n*
```

where:
- `γ*`: Wall quenching/sticking coefficient (0 < γ* < 1)
- `v_th* = √(8kT_gas/πm*)`: Mean thermal velocity
- **Thermal flux**: (1/4)·n*·v_th* (kinetic theory)

**Implementation:**
```
Flux_wall = (γ* · v_th* / 4) · n*_boundary
```

If `u_gas < 0` (flow toward wall), add advective contribution.

### 5.3 Coupling to Electron Equation

If quenching releases secondary electrons (`γ_see* > 0`):

```
Γ_e·n̂|_wall = γ_see* · Γ*·n̂
```

This couples the "separate expression" for neutrals back to the electron continuity equation, enabling **memory effects** in DBDs where metastables at the dielectric contribute to seed electrons in the next voltage cycle.

---

## 6. Solver Architecture <a name="solver-architecture"></a>

### 6.1 Fully Implicit Newton-Krylov Framework

The addition of M excited species increases DOFs per node from 5 (ne, ni, nε, φ, σ) to 5+M. Given the **stiffness** (reaction timescales ~ ns, diffusion ~ μs, flow ~ ms), explicit time-stepping is prohibitive.

**Solution:** Fully implicit DAE solver using PETSc TS:

```
F(X) = dX/dt + ∇·Γ(X) - S(X) = 0
```

Solve via Newton iteration:
```
J·δX^k = -F(X^k)
X^{k+1} = X^k + δX^k
```

### 6.2 Jacobian Structure

The Jacobian includes:
- **Diagonal blocks**: Self-coupling via transport
- **Off-diagonal blocks**: 
  - `∂F_neutral/∂ne`: Electron impact excitation
  - `∂F_electron/∂n_neutral`: Stepwise/Penning ionization

**Including off-diagonals is critical** for Newton convergence when stepwise dominates.

### 6.3 Jacobian-Free Newton-Krylov (JFNK)

For large M, storing the full Jacobian is expensive. PETSc's JFNK approximates:

```
J·v ≈ [F(X + εv) - F(X)] / ε
```

This avoids explicit Jacobian assembly while retaining Newton convergence.

### 6.4 Preconditioning: FieldSplit

PETSc's `PCFIELDSPLIT` enables block preconditioning:
- **Block 0**: Transport equations (ne, ni, nε, n*)  
  - Preconditioner: ILU or GMRES
- **Block 1**: Poisson equation (φ)  
  - Preconditioner: Direct solve (LU)

This exploits the different mathematical character of elliptic (Poisson) vs. parabolic (transport) equations.

---

## 7. References <a name="references"></a>

1. **Hagelaar & Pitchford** (2005). "Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models." *Plasma Sources Sci. Technol.* 14, 722.

2. **Sakiyama et al.** (2012). "Plasma chemistry model of surface microdischarge in humid air and dynamics of reactive neutral species." *J. Phys. D: Appl. Phys.* 45, 425201.

3. **Balay et al.** (2023). "PETSc Users Manual." Argonne National Laboratory. ANL-21/39.

4. **Scharfetter & Gummel** (1969). "Large-signal analysis of a silicon Read diode oscillator." *IEEE Trans. Electron Devices* 16, 64-77.

5. **Boeuf & Pitchford** (1995). "Two-dimensional model of a capacitively coupled RF discharge." *Phys. Rev. E* 51, 1376.

6. **Guerra & Loureiro** (1997). "Electron and heavy particle kinetics in a low-pressure nitrogen glow discharge." *Plasma Sources Sci. Technol.* 6, 361.

7. **Porazik et al.** (2012). "Fluid simulation of electron heating and transport in atmospheric pressure plasma jets." *J. Appl. Phys.* 112, 053302.

8. **Naidis** (2011). "Modelling of plasma bullet propagation along a helium jet in ambient air." *J. Phys. D: Appl. Phys.* 44, 215203.

---

## Appendix: Implementation Checklist

✅ **Physics:**
- [x] ADR equations for excited species
- [x] 9 reaction mechanisms (ionization, excitation, stepwise, Penning, superelastic, radiative, quenching, pooling, recombination)
- [x] Robin boundary conditions with surface quenching
- [x] Secondary electron emission from excited neutrals

✅ **Numerics:**
- [x] Scharfetter-Gummel flux discretization
- [x] Bernoulli function implementation with Taylor series fallback
- [x] Local Péclet number calculation

✅ **Solver:**
- [x] PETSc TS implicit DAE solver
- [x] Newton-Krylov framework
- [x] FieldSplit preconditioning
- [x] Automatic Jacobian via finite differencing with coloring

✅ **Chemistry:**
- [x] BOLSIG+ interface for rate coefficients
- [x] LookupTable for k(mean_energy)
- [x] ReactionHandler with modular reaction types

✅ **I/O:**
- [x] HDF5/OpenPMD output for multi-species visualization
- [x] Text output for debugging
- [x] Hierarchical data organization

✅ **Configuration:**
- [x] JSON-based runtime configuration
- [x] Multiple example configs (DBD, plasma jet, Penning mixture)
- [x] Comprehensive Argon chemistry (Ar_m, Ar_r, Ar2*)

---

**Last Updated:** December 30, 2025  
**Code Version:** HydroPlas v1.0 with Explicit Excited Species Transport
