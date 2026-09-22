# Theoretical and Numerical Formulation: Reacting Flows in OpenFOAM

**Author:** Tushar Sagar  
**Specialization:** Combustion CFD, Reacting Flows & Propulsion Diagnostics  

---

## 1. Governing Navier-Stokes Equations for Compressible Multi-Species Flows

In turbulent reacting flows, the conservation laws for mass, momentum, total enthalpy, and chemical species mass fractions are solved within the Finite Volume Method (FVM) framework.

### 1.1 Continuity Equation
$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{U}) = 0$$

### 1.2 Momentum Transport
$$\frac{\partial (\rho \mathbf{U})}{\partial t} + \nabla \cdot (\rho \mathbf{U} \mathbf{U}) = -\nabla p + \nabla \cdot \boldsymbol{\tau}_{\mathrm{eff}} + \rho \mathbf{g} + \mathbf{S}_p$$
where the effective viscous stress tensor combines molecular and Reynolds stress contributions:
$$\boldsymbol{\tau}_{\mathrm{eff}} = \left( \mu + \mu_t \right) \left[ \nabla \mathbf{U} + (\nabla \mathbf{U})^T - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I} \right] - \frac{2}{3} \rho k \mathbf{I}$$

### 1.3 Total Energy / Enthalpy Equation
Solving for sensible enthalpy $h_s$:
$$\frac{\partial (\rho h_s)}{\partial t} + \nabla \cdot (\rho \mathbf{U} h_s) + \frac{\partial (\rho K)}{\partial t} + \nabla \cdot (\rho \mathbf{U} K) - \frac{\partial p}{\partial t} = \nabla \cdot \left[ \left( \frac{\mu}{\mathrm{Pr}} + \frac{\mu_t}{\mathrm{Pr}_t} \right) \nabla h_s \right] + \dot{\omega}_T + S_{\mathrm{rad}}$$
where $\dot{\omega}_T = -\sum_{i=1}^{N_s} \dot{\omega}_i \Delta h_{f,i}^{\circ}$ is the heat release rate from chemical reactions.

### 1.4 Species Mass Fraction Conservation
For each chemical species $i \in \{1, \dots, N_s - 1\}$:
$$\frac{\partial (\rho Y_i)}{\partial t} + \nabla \cdot (\rho \mathbf{U} Y_i) = \nabla \cdot \left[ \left( \rho D_{i,m} + \frac{\mu_t}{\mathrm{Sc}_t} \right) \nabla Y_i \right] + \dot{\omega}_i + S_{p,i}$$

---

## 2. Turbulence-Chemistry Interaction: The PaSR Model

In turbulent combustion, molecular mixing and chemical reaction rates occur simultaneously across disparate time and length scales. The **Partially Stirred Reactor (PaSR)** model assumes that each computational cell is split into two regions:
1. **Reaction zone:** Fraction of the cell $\kappa$ where reactants are mixed at the molecular level and react.
2. **Non-reacting zone:** Fraction $(1 - \kappa)$ where fluids remain unmixed.

### Reactive Fraction Calculation:
$$\kappa = \frac{\tau_c}{\tau_c + \tau_{\mathrm{mix}}}$$

where:
- $\tau_c$ is the characteristic chemical timescale:
  $$\tau_c = \frac{Y_i^*}{\dot{\omega}_i(T^*, Y^*)}$$
- $\tau_{\mathrm{mix}}$ is the turbulent micro-mixing timescale modeled via Kolmogorov or Taylor scales:
  $$\tau_{\mathrm{mix}} = C_{\mathrm{mix}} \sqrt{\frac{\mu_{\mathrm{eff}}}{\rho \epsilon}}$$

The mean chemical source term injected into the FVM transport equation is filtered as:
$$\overline{\dot{\omega}}_i = \kappa \cdot \dot{\omega}_i(\tilde{T}, \tilde{Y}_i)$$

---

## 3. Extension to Multiphase Metal Particle Combustion (Eulerian-Lagrangian)

In the context of **metal fuels (e.g., Iron Fe, Aluminium Al, Magnesium Mg)** for zero-carbon energy storage or solid propulsion:

1. **Carrier Gas Phase (Eulerian):** Solved via standard Navier-Stokes reacting flow equations with mass, momentum, and heat transfer source terms $\mathbf{S}_p, S_{p,i}$ exchanged with particles.
2. **Dispersed Metal Particles (Lagrangian DPM):**
   - **Particle Trajectory:**
     $$\frac{d\mathbf{x}_p}{dt} = \mathbf{u}_p, \quad m_p \frac{d\mathbf{u}_p}{dt} = \mathbf{F}_D + \mathbf{F}_G + \mathbf{F}_{\mathrm{press}}$$
   - **Heterogeneous Surface Reaction (e.g., Iron oxidation):**
     $$4\mathrm{Fe}_{(s)} + 3\mathrm{O}_{2(g)} \longrightarrow 2\mathrm{Fe}_2\mathrm{O}_{3(s)} + \Delta H_{rxn}$$
   - **Shrinking Core / Diffusion-Controlled Combustion Rate:**
     $$\frac{dm_p}{dt} = - \pi d_p \mathrm{Sh} \cdot \rho D_{\mathrm{O}_2} \ln(1 + B_M)$$
     where $B_M$ is the Spalding mass transfer number and $\mathrm{Sh}$ is the particle Sherwood number.

---

## 4. Numerical Discretization in OpenFOAM

- **Pressure-Velocity Coupling:** PIMPLE algorithm (merged PISO-SIMPLE) for transient reacting flows.
- **Spatial Discretization:** 
  - Convective terms: TVD / limited schemes (`limitedLinear01 1` for bounded species mass fractions $Y_i \in [0, 1]$).
  - Diffusion: Second-order Gaussian integration with orthogonal/non-orthogonal correction.
- **Temporal Discretization:** First-order implicit Euler or second-order backward differencing.
- **Chemistry ODE Integration:** Stiff ODE solvers using semi-implicit extrapolation (`seulex`) or Rosenbrock solvers.
