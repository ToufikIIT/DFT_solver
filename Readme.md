# Simple Kohn-Sham DFT Solver for Helium

![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)

This repository contains a simple, self-contained C++ implementation of a Density Functional Theory (DFT) solver for the ground state of a Helium (He) atom. The code solves the Kohn-Sham equations within the Local Density Approximation (LDA) for the exchange-correlation functional.

The primary purpose of this code is educational: to demonstrate the core components of a DFT calculation, including the self-consistent field (SCF) procedure, numerical solution of the radial Schrödinger equation, and calculation of system energies.

## Table of Contents
1.  [Theoretical Background](#theoretical-background)
    * [Kohn-Sham DFT](#kohn-sham-dft)
    * [Self-Consistent Field (SCF) Cycle](#self-consistent-field-scf-cycle)
    * [Local Density Approximation (LDA)](#local-density-approximation-lda)
2.  [Implementation Details](#implementation-details)
3.  [How to Build and Run](#how-to-build-and-run)
    * [Prerequisites](#prerequisites)
    * [Compilation](#compilation)

---

## Theoretical Background

### Kohn-Sham DFT
Density Functional Theory is a quantum mechanical method used to investigate the electronic structure of many-body systems. It reformulates the problem from finding the complex many-electron wavefunction to finding the much simpler electron density, $n(\vec{r})$.

The Kohn-Sham approach models the interacting system as a fictitious system of non-interacting electrons moving in an effective potential, $V_{\text{eff}}(r)$. For a spherically symmetric atom like Helium, the Kohn-Sham equation for the radial part of the wavefunction $u(r) = r R(r)$ is:

$$\left[ -\frac{1}{2}\frac{d^2}{dr^2} + \frac{\ell(\ell+1)}{2r^2} + V_{\text{eff}}(r) \right] u_{nl}(r) = \epsilon_{nl} u_{nl}(r)$$

The effective potential is the sum of the external potential (from the nucleus), the Hartree potential (classical electrostatic repulsion), and the exchange-correlation potential:

$$V_{\text{eff}}(r) = V_{\text{ext}}(r) + V_{H}(r) + V_{xc}(r) = -\frac{Z}{r} + \int \frac{n(r')}{|\vec{r}-\vec{r}'|} d^3r' + V_{xc}(r)$$

### Self-Consistent Field (SCF) Cycle
Since the effective potential $V_{\text{eff}}$ depends on the density $n(r)$, which in turn depends on the wavefunctions found by solving the Kohn-Sham equation, the problem must be solved iteratively. This process is called the Self-Consistent Field (SCF) cycle.

1.  **Initial Guess**: Start with an initial guess for the electron density $n(r)$.
2.  **Calculate Potential**: Construct the effective potential $V_{\text{eff}}(r)$ from this density.
3.  **Solve KS Equations**: Solve the Kohn-Sham equations to find the new orbital energies $\epsilon$ and wavefunctions $u(r)$.
4.  **Calculate New Density**: Calculate the new electron density from the wavefunctions: $n_{\text{new}}(r) = \sum_i^{\text{occ}} \frac{|u_i(r)|^2}{4\pi r^2}$.
5.  **Check for Convergence**: Compare $n_{\text{new}}(r)$ with the previous density $n(r)$. If they are sufficiently close, the cycle has converged.
6.  **Mix and Repeat**: If not converged, mix the old and new densities ($n(r) = \alpha n_{\text{new}}(r) + (1-\alpha)n(r)$) and go back to step 2.

### Local Density Approximation (LDA)
The exact form of the exchange-correlation (XC) functional is unknown. The LDA is a simple yet effective approximation where the XC energy at a point $r$ is taken to be that of a uniform electron gas with the same density $n(r)$. This code uses the Slater exchange functional, which is the exchange part of the LDA:

$$V_x(r) = - \frac{3}{4} \left( \frac{3}{\pi} \right)^{1/3} [n(r)]^{1/3}$$

The correlation part of the potential is not included in this simple implementation.

---

## Implementation Details
* **System**: Helium atom ($Z=2$) with two electrons in the 1s orbital ($\ell=0$).
* **Numerical Grid**: A uniform radial grid is used for all functions.
* **Differential Equation Solver**: The radial Kohn-Sham equation is solved using the **Numerov method**, a highly accurate numerical technique for second-order ordinary differential equations.
* **Eigenvalue Solver**: A **shooting method** combined with a bisection search is used to find the orbital energy $\epsilon$ that satisfies the boundary condition $u(r \to \infty) = 0$.
* **Integration**: The trapezoidal rule is used for all numerical integrations.
* **Convergence**: The SCF cycle is stopped when the change in total energy and the maximum relative change in density between iterations fall below a set tolerance.

---

## How to Build and Run

### Prerequisites
* A C++ compiler that supports the C++11 standard or later (e.g., GCC, Clang).

### Compilation
Open a terminal in the project directory and compile the source file using `g++`:

```bash
g++ -std=c++17 -O2 -o dft_he main.cpp -lm
```
## Expected Output
- The output will show the iteration number, the 1s orbital energy (eps), the total energy (Etot), and the maximum relative change in density for each SCF step. The final lines will show the converged results.

- Iter   1  eps =   -0.99986591  Etot =   -2.74971217  densΔ(max rel) = 5.234693e-02
- Iter   2  eps =   -0.96347101  Etot =   -2.82527581  densΔ(max rel) = 2.158579e-02
- Iter   3  eps =   -0.94635105  Etot =   -2.84654519  densΔ(max rel) = 9.878783e-03
- ...
- ...
- Iter  17  eps =   -0.91793739  Etot =   -2.86167651  densΔ(max rel) = 1.637537e-07
- Iter  18  eps =   -0.91793725  Etot =   -2.86167675  densΔ(max rel) = 7.747514e-08

- Converged:
  * Orbital energy (1s): -0.9179372505 Ha
  * Total KS energy    : -2.861676754 Ha

## Code Structure
-The code is organized into several key functions:

### Function	Description
- main()	The main driver function. Initializes the density and runs the SCF loop.
- setup_r_grid()	Creates the uniform radial grid r.
- trapz_...()	A set of functions for numerical integration using the trapezoidal rule.
- calc_hartree()	Calculates the Hartree potential V_H(r) from the electron density n(r).
- calc_v_x()	Calculates the LDA exchange potential V_x(r) from the electron density.
- build_v_eff()	Assembles the total effective potential V_texteff(r).
- numerov_solve()	Solves the radial Schrödinger equation for a given energy E to find u(r).
- find_ground_energy()	Implements the shooting method to find the ground state orbital energy eps.
- density_from_u()	Computes the electron density n(r) from the radial wavefunction u(r).

## The Radial Schrödinger Equation

The Hamiltonian for a spherically symmetric potential is given by:
$$\hat{H} = -\frac{1}{2}\nabla^2 + V(r)$$
Because of the spherical symmetry, we can separate the variables of the wavefunction $\psi(r, \theta, \phi)$ as follows:
$$\psi(r, \theta, \phi) = \frac{u_{\ell}(r)}{r} Y_{\ell}^{m}(\theta, \phi)$$
This separation leads to the **radial equation**:
$$\frac{d^2 u_{\ell}(r)}{dr^2} = \left[ 2(E - V(r)) - \frac{\ell(\ell+1)}{r^2} \right] u_{\ell}(r)$$
Here:
* $u_{\ell}(r)$ is the radial function.
* $V(r) = -\frac{Z}{r}$ for a hydrogenic atom (with nuclear charge $Z$).
* $\ell$ is the angular momentum quantum number.

---
## Numerov Method

The radial equation is a **second-order linear ODE**, which can be solved efficiently using the Numerov algorithm. The recurrence relation is given by:
$$u_{i+1} = \frac{\left(2 - \frac{5}{6}h^2k_i^2\right)u_i - \left(1 + \frac{1}{12}h^2k_{i-1}^2\right)u_{i-1}}{1 + \frac{1}{12}h^2k_{i+1}^2}$$
where:
$$k_i^2 = 2(E - V(r_i)) - \frac{\ell(\ell+1)}{r_i^2} \quad , \quad h = \Delta r$$
This numerical scheme is **sixth-order accurate**, making it ideal for finding solutions to quantum bound states.

---
## Self-Consistent Field (SCF) Cycle

The program follows the usual **SCF procedure**:

1.  Start with a **guess density** (e.g., from a hydrogenic 1s orbital).
2.  Build the **Hartree and exchange potentials**.
3.  Construct the **effective potential** $V_{\text{eff}}(r)$.
4.  Solve the Kohn-Sham (KS) equation with the **Numerov method** to get the new wavefunction $u(r)$.
5.  Recompute the **electron density** $n(r)$.
6.  **Mix densities** to aid convergence:
    $$n^{\text{in+1}} = \alpha n^{\text{new}} + (1-\alpha)n^{\text{old}}$$
7.  Check for **convergence** in the total energy and electron density.
8.  **Repeat** the cycle until self-consistency is achieved.
