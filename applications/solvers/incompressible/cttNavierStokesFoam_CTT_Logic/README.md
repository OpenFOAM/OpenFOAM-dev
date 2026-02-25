# cttNavierStokesFoam_CTT_Logic

## Overview
`cttNavierStokesFoam_CTT` is a custom incompressible Navier-Stokes solver developed for OpenFOAM v2512. It implements a specialized PISO (Pressure-Implicit with Splitting of Operators) loop integrated with **CTT (Constant Time Transfer)** logic. 

This solver is specifically designed as a spectral reference tool to investigate the **Global Regularity** of the 3D Navier-Stokes equations using a fractal-layer regularity approach.

## Key Features
* **CTT Integration**: Custom momentum equation structure designed for spectral scale-transfer analysis.
* **Modern API**: Fully compatible with OpenFOAM-v2512 (OpenCFD).
* **Reference Handling**: Explicit pressure reference support for closed-domain simulations.
* **PISO Coupling**: Robust velocity-pressure coupling using the modern `pisoControl` class.

## Installation

### Compilation
Navigate to the solver directory and use `wmake`:

```bash
cd applications/solvers/incompressible/cttNavierStokesFoam_CTT_Logic
wmake
Theory & LogicThe solver modifies the standard momentum predictor step to account for CTT constraints:$$\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot (\phi \mathbf{U}) - \nabla \cdot (\nu \nabla \mathbf{U}) = -\nabla p$$The CTT logic ensures that the spectral transfer remains bounded, providing a numerical framework for the Global Regularity Proof.
