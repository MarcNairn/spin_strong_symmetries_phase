# Phase-tied Strong Symmetries

This repository contains scripts to generate the data and figures used in the paper *(arXiv:??/??2025)*. The code implements mean-field equations of motion and finite-size exact-diagonalization numerics for permutation-symmetric spin systems. A short utility script is included to automatically identify strong symmetries for collective two-level systems.

**Short description.**  
Below we summarise the physical context, notation and the essential theoretical ingredients reproduced by the simulations. The README focuses on (i) open-system / strong-symmetry concepts, (ii) the two-spin-species (collective spin-1/2) model and its spin-only reduction, and (iii) the collective three-level model and its symmetry structure. The goal is a compact, copy-ready reference with definitions, key equations, diagnostics and practical notes to reproduce and extend the analyses.

---

## Table of contents
- [Quick summary of results](#quick-summary-of-results)
- [Notation and basic definitions](#notation-and-basic-definitions)
- [Open quantum systems and strong symmetries](#open-quantum-systems-and-strong-symmetries)
- [Two-spin-species model (collective spin-1/2)](#two-spin-species-model-collective-spin-12)
  - [Microscopic Hamiltonian and dissipator](#microscopic-hamiltonian-and-dissipator)
  - [Adiabatic elimination → spin-only master equation](#adiabatic-elimination--spin-only-master-equation)
  - [Phase-dependent strong symmetry and phase-Dicke sectors](#phase-dependent-strong-symmetry-and-phase-dicke-sectors)
  - [Dynamical diagnostics and main phenomena](#dynamical-diagnostics-and-main-phenomena)
- [Collective three-level (effective spin-1) model](#collective-three-level-effective-spin-1-model)
  - [Atom-only master equation and operators](#atom-only-master-equation-and-operators)
  - [Strong symmetry generator, bright/dark eigenstates](#strong-symmetry-generator-brightdark-eigenstates)
  - [State-space dependence, multistability and DFS / eDFS](#state-space-dependence-multistability-and-dfs--edfs)
- [How to reproduce main numerics (practical checklist)](#how-to-reproduce-main-numerics-practical-checklist)
- [Background reading](#background-reading)

---

## Quick summary of results
- A controllable relative phase $\varphi$ between light–matter couplings produces a **phase-dependent strong symmetry** in collective cavity-QED models.
- The phase $\varphi$ continuously rotates the Liouvillian symmetry sectors (a phase-dependent generalized Dicke decomposition), changing the initial state's overlap with bright (drive-coupled) and dark (protected) sectors.
- As a result, the critical driving strength for entering non-stationary (time-crystalline / persistent oscillatory) phases can be strongly reduced or even vanish for particular $\varphi$ and initial states.
- Two model realizations are demonstrated:
  1. Two spin species (each a collective spin-1/2 ensemble) with coupling $g_A=g$, $g_B=g e^{-i\varphi}$ — strong symmetry: $\hat{\mathbf S}_\varphi^2$.
  2. Collective three-level (V-type) atoms sharing a ground state — strong symmetry: $\hat{\mathcal T}_\varphi$, yielding explicit bright and dark many-body eigenstates.
- Finite detunings produce decoherence-free subspaces (DFS) or an emergent DFS (eDFS) in the thermodynamic limit; the phase permits restoration/control of oscillatory phases even when detuning would otherwise suppress them.

---

## Notation and basic definitions
- $\hat{\rho}$: density matrix of the full system; $\hat{\mu}_s$, $\hat{\mu}_{\mathrm{3S}}$ denote reduced atomic (spin / 3-level) density matrices after cavity elimination.
- Cavity photon annihilation / creation: $\hat a$, $\hat a^\dagger$.
- Phase-dependent collective spin lowering:
  $$
  \hat S_\varphi = \sum_{j=1}^{N/2}\left(\hat\sigma_{jA} + e^{-i\varphi}\hat\sigma_{jB}\right),
  \qquad
  \hat\sigma_{jm} = \ket{\downarrow}_{jm}\!\bra{\uparrow}.
  $$
- Dissipator:
  $$
  \mathcal{D}[\hat L]\rho = \hat L\rho\hat L^\dagger - \tfrac{1}{2}\{\hat L^\dagger\hat L,\rho\}.
  $$
- Liouvillian:
  $$
  \mathcal L(\rho) = -\tfrac{i}{\hbar}[\hat H,\rho] + \sum_k \mathcal D[\hat L_k]\rho.
  $$
- Sector weights (projector $\Pi_\lambda$ onto eigen-sector of a strong symmetry $\hat A$):
  $$
  w_\lambda = \operatorname{Tr}\!\left[\Pi_\lambda\,\rho(0)\right].
  $$
- Order parameter for nonstationarity: normalized dominant frequency of an observable’s spectrum,
  $$
  \tilde\omega = \arg\max_\omega \mathcal F(\omega), \qquad
  \mathcal F(\omega)=\int s^z(t)e^{-i\omega t}\,dt.
  $$

---

## Open quantum systems and strong symmetries
- GKSL master equation:
  $$
  \dot{\hat\rho} = \mathcal L(\hat\rho)
  = -\frac{i}{\hbar}[\hat H,\hat\rho]
  + \sum_k \mathcal D[\hat L_k]\hat\rho.
  $$
- A **strong symmetry** $\hat A$ satisfies
  $$
  [\hat A,\hat H]=0, \qquad [\hat A,\hat L_k]=0 \;\; \forall k.
  $$
- Consequences:
  - Hilbert space decomposes into invariant subspaces $\mathcal H=\bigoplus_\lambda\mathcal H_\lambda$.
  - Dynamics decouple between sectors and the sector weights $w_\lambda$ are conserved.
  - Protected subspaces or long-lived oscillatory modes may appear depending on the Liouvillian spectrum.

---

## Two-spin-species model (collective spin-1/2)

### Microscopic Hamiltonian and dissipator
- Two species $A,B$, each with $N/2$ two-level atoms.
- Couplings $g_A=g$, $g_B=g e^{-i\varphi}$; detunings $\Delta_A=\Delta_B\equiv\Delta$.
- Hamiltonian (rotating frame):
  $$
  \begin{aligned}
  \hat H_\mathrm{cav} &= \hbar\Delta_c\,\hat a^\dagger\hat a
    + i\hbar\eta\sqrt{N}\,(\hat a^\dagger-\hat a),\\
  \hat H_\mathrm{int} &= \frac{\hbar g}{\sqrt N}\big(\hat a^\dagger\hat S_\varphi + \hat a\hat S_\varphi^\dagger\big),\\
  \hat H_\mathrm{at} &= \hbar\Delta\,\hat S_\varphi^z .
  \end{aligned}
  $$
- Full master equation:
  $$
  \dot{\hat\rho} = -\frac{i}{\hbar}[\hat H_2,\hat\rho] + \kappa\,\mathcal D[\hat a]\hat\rho.
  $$

### Adiabatic elimination → spin-only master equation
- In the bad-cavity limit $\kappa \gg g,\Delta$:
  $$
  \dot{\hat\mu}_s = -\frac{i}{\hbar}[\hat H_S,\hat\mu_s]
  + \Gamma\,\mathcal D[\hat S_\varphi]\hat\mu_s,
  $$
  with
  $$
  \hat H_S = \hbar\Delta\,\hat S_\varphi^z
  + \hbar(\Omega\hat S_\varphi + \Omega^*\hat S_\varphi^\dagger)
  + \hbar\mathcal C\,\hat S_\varphi^\dagger\hat S_\varphi.
  $$

### Phase-dependent strong symmetry and phase-Dicke sectors
- $\{\hat S_\varphi,\hat S_\varphi^\dagger,\hat S_\varphi^z\}$ generate an $\mathfrak{su}(2)$ algebra.
- The Casimir
  $$
  \hat{\mathbf S}_\varphi^2
  = \tfrac{1}{2}\big(\hat S_\varphi\hat S_\varphi^\dagger
  + \hat S_\varphi^\dagger\hat S_\varphi\big)
  + (\hat S_\varphi^z)^2
  $$
  commutes with $\hat H_S$ and $\hat S_\varphi$, defining a strong symmetry.
- This induces a **phase-dependent Dicke decomposition** $\ket{S^2,m_z}_\varphi$ with conserved sector weights.

### Dynamical diagnostics and main phenomena
- **Order parameter:** dominant oscillation frequency from Fourier spectra of collective observables.
- **Mechanism:** varying $\varphi$ rotates symmetry sectors, modifying the initial overlap with the bright subspace and lowering the effective critical drive.
- **Finite detuning:** may generate DFS or suppress oscillations; appropriate $\varphi$ can restore nonstationary dynamics.

---

## Collective three-level (effective spin-1) model

### Atom-only master equation and operators
- V-type atoms with ground $\ket{0}$ and excited states $\ket{A},\ket{B}$.
- Collective lowering operator:
  $$
  \hat\Lambda_\varphi = \hat\Lambda_A + e^{-i\varphi}\hat\Lambda_B,
  \qquad
  \hat\Lambda_k = \sum_{j=1}^N \ket{0}_j\bra{k}_j.
  $$
- Master equation:
  $$
  \dot{\hat\mu}_{\mathrm{3S}}
  = -\frac{i}{\hbar}[\hat H_{\mathrm{3S}},\hat\mu_{\mathrm{3S}}]
  + \Gamma\,\mathcal D[\hat\Lambda_\varphi]\hat\mu_{\mathrm{3S}}.
  $$

### Strong symmetry generator, bright/dark eigenstates
- Strong symmetry:
  $$
  \hat{\mathcal T}_\varphi
  = \hat N_0
  + \left(e^{-i\varphi}\hat\Lambda_A^\dagger\hat\Lambda_B
  + e^{i\varphi}\hat\Lambda_B^\dagger\hat\Lambda_A\right).
  $$
- Single-particle bright/dark states:
  $$
  \ket{\pm_{AB}} = \frac{1}{\sqrt{2}}\left(\ket{B}\pm e^{-i\varphi}\ket{A}\right).
  $$

### State-space dependence, multistability and DFS / eDFS
- Initial-state overlaps with bright and dark manifolds determine long-time dynamics.
- Phase $\varphi$ rotates the bright/dark regions in state space, leading to dynamical multistability.
- Finite detuning can produce a DFS or, in the thermodynamic limit, an emergent DFS (eDFS).

---

## How to reproduce main numerics (practical checklist)
1. Choose collective scalings consistent with the thermodynamic limit.
2. Solve either the spin-only master equation (small $N$) or the mean-field equations.
3. Scan $\varphi$ and drive strength $\eta$.
4. Compute Fourier spectra of collective observables.
5. Analyze Liouvillian eigenvalues to identify DFS / oscillatory modes.
6. Perform finite-size scaling to confirm thermodynamic behavior.

---

## Background reading
- H. P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems* (2002).
- G. Lindblad, *Commun. Math. Phys.* (1976).
