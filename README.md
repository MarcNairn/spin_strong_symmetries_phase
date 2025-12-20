# Phase-engineered Strong Symmetries

This repository contains scripts to generate the data and figures used in the paper *(arXiv:??/??2025)*. The code implements the mean-field equations of motion and finite-size exact-diagonalization numerics for permutation-symmetric spin systems. A short utility script is included to automatically identify strong symmetries in collective two- and three-level systems.

---

## Quick Summary of Main Results

- A *controllable relative phase* $\varphi$ between light–matter couplings produces a phase-dependent **strong symmetry** in collective cavity-QED models.
- The eigenstates of the combined system are determined by the choice of $\varphi$.
- By *continuously tuning* $\varphi$, the *critical driving strength* $\eta_c$ to trigger the dissipative phase transition and produce **nonstationary states** can be **arbitrarily lowered**. 
- The results are general and apply to any class of spin systems (here we study 2 and 3 level systems) that admit a collective representation, so long as the initial state does not live in a trivial/uncorrelated eigenstate.

## Mean field analysis

The workflow of the program is as follows.

We first consider solving the steady-state dynamics of the mean-field equations of motion. From this, we extract the long-time window observables $O(t)$ and their Fourier transform

$$
\mathcal{F}_O(\omega) = \int_{-\infty}^\infty dt\; e^{i\omega t} O(t).
$$

For our interests, we choose to define the **order parameter**

$$
\tilde{\omega}_O = \arg\max_\omega\left[\mathcal{F}_O(\omega)\right] \in [0,1],
$$

which, in practice, determines the position of the frequency peak in the normalized Fourier spectrum after the DC component is removed. This is the object we study to discriminate between stationary states and those with persistent oscillations into their steady state (nonstationary states).

We also need to consider the choice of initial mean field variables for both spin species and the cavity field. For this, we consider the $x$-polarized initial state for both spins, corresponding to the fully symmetric superposition $$\left(|\uparrow\downarrow\rangle + |\downarrow\uparrow\rangle\right)/\sqrt{2}$$ and an empty cavity, $\langle a^\dagger a\rangle = 0$. The real valued mean field variable vector takes the form:

$$
u_0 = (a_x, a_p, m_x^A, m_y^A, m_z^A, m_x^B, m_y^B, m_z^B)
= (0, 0, 1, 0, 0, 1, 0, 0).
$$

The scripts to implement the mean field analysis are stored in the `simul` folder, and contain a leading "`MF`" in their filename. 


## Dimensionality Reduction via Permutation Invariance

At the same time, we also consider performing finite-size dynamics. The key to simulating large system sizes (up to $N = 100$ spins) lies in exploiting the full permutation symmetry of the Hamiltonian and jump operators. For a system of $N$ identical two-level systems (spin-1/2 particles), the Hilbert space dimension naively scales as $2^N$. However, under global permutation symmetry, the state is completely characterized by its total angular momentum $j$ and the projection $m_j$, along with the multiplicity (Dicke manifold).

For permutation-invariant states, the relevant subspace is the **Dicke manifold** (or symmetric subspace), where all spins are indistinguishable. The dimension of this symmetric subspace is only $N + 1$, growing linearly with $N$ instead of exponentially.

For example:
- $N$ qubits: The symmetric subspace dimension is $N + 1$.
- $N$ qutrits (three-level systems): The symmetric subspace dimension is $\frac{(N+1)(N+2)}{2}$, growing quadratically instead of as $3^N$.

### Qutip Implementation

The `PermutationalInvariant` quantum solver (`PIQS`) in QuTiP automatically restricts the dynamics to this symmetric subspace by:
1.  Constructing the Liouvillian directly in the basis of permutation-invariant states.
2.  Using the angular momentum algebra to represent collective spin operators ($J_x, J_y, J_z, J_\pm$) that act on the Dicke states $|j, m_j\rangle$.
3.  Handling the non-collective dissipative terms (if present) by properly symmetrizing them.

This reduction allows the simulation of system sizes that are otherwise computationally intractable for exact density-matrix evolution, enabling the study of collective phenomena in large spin ensembles.

The implementation can be found within the `simul` folder, under `finite_size_spins.ipynb`.


## Strong symmetries and where to find them

In the context of open quantum systems described by a Lindblad master equation, a **strong symmetry** is defined as an operator $X$ that simultaneously commutes with both the Hamiltonian $H$ and *every* jump operator $L_k$: $[X, H] = 0$ and $[X, L_k] = 0$ for all $k$. This stringent condition ensures that the corresponding observable is conserved under the full dynamics, both coherent and dissipative, leading to a block-diagonal structure of the Liouvillian superoperator in the eigenbasis of $X$. Each block corresponds to a distinct symmetry sector, and coherences between these sectors are perfectly protected from decoherence.

Since this is a key feature of our work, it is also important to be able to find such a $X$. For this, we implement a symbolic algorithm in Python that solves for operators $X$. The script can be found in the `simul/strongsym` folder and works as follows:

1. We define the Hamiltonian $H$ and jump operator $L$ for either a two-level (spin-1/2) or three-level system, keeping the relative phase $\varphi$ as a symbolic parameter.

2. We impose $[X,H] = 0$ and $[X,L] = 0$, where $X$ is an unknown operator represented by a matrix with symbolic entries.

3. These matrix equations are expanded into a system of linear equations in the unknown matrix elements of $X$. Using sympy's `linear_eq_to_matrix`, we construct the coefficient matrix $A$ and solve for its nullspace (i.e. 0 commutation relation).

4. The nullspace basis vectors correspond to operators commuting with both $H$ and $L$. We identify the first non-identity basis element as the strong symmetry generator $X(\varphi)$. For the spin-1/2 case, this yields $\mathbf{S}^2_\varphi$; for the three-level system, we obtain $\mathcal{T}_\varphi$ as reported in the main text. If no non-trivial matrix is found it returns the identity.

The resulting operator's eigenstates reveal the invariant subspaces of the dynamics. By examining these eigenstates, we can construct the appropriate initial conditions that will remain within a single symmetry sector, enabling the stabilization of nonstationary phases at low driving strengths.

## Decoherence-free subspaces

A decoherence-free subspace (DFS) is a subspace of the total Hilbert space that is completely immune to environmental decoherence, where quantum states evolve unitarily despite the system being open. In our model, when both spin species have equal and nonzero detuning, a DFS emerges naturally. 

The presence of a DFS can be directly identified in the Liouvillian eigenvalue spectrum. The DFS can be immediately identified by a set of purely imaginary eigenvalues $\lambda_{\mathcal{L}} = \pm i \Delta$, corresponding to undamped oscillatory modes. These are traditionally independent of system size as they arise due to the symmetry of the system. 

To spot such features we devise a small script to diagonalize the full Liouvillian spectrum for a set of system sizes, for both the two level and three level systems. These are stored in `simul/alt_finite_size_spins.ipynb` and `simul/liouv_spec_spin1.ipynb`, respectively.

## Repo structure

For ease of navigation, this repository is organized as follows:

- **`simul/`**: Contains all simulation scripts for both mean-field (`MF_*.ipynb`) and finite-size exact dynamics (`finite_size_spins.ipynb`, `alt_finite_size_spins.ipynb`). 

  - **`simul/strongsym/`**: Includes helper functions and the symbolic symmetry-finding script (`find_symmetry.ipynb`).

- **`old_unused/`**: auxiliary files (not essential).

- **`finsize_data/`**: Stores data used for figures generated from the `simul` scripts

- **`figs_for_paper/`**: Contains plots and plot scripts for the figures used in the manuscript.

- **`dependencies/`**: Contains a `.yml` file for the Python dependencies and a separate `.toml` file for the julia ones.

- **`README.md`**: This file...
---

## Background Reading

- H. P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems* (2002).
- G. Lindblad, *Commun. Math. Phys.* (1976).
- V. V. Albert & L. Jiang, *Physical Review A* (2014)
- F. Minganti, A. Biella, N. Bartolo & C. Ciuti, *Physical Review A* (2018).
