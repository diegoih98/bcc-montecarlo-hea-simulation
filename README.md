# bcc-montecarlo-hea-simulation

This repository contains Python code for simulating ordering/disordering in BCC-structured concentrated alloys using a Monte Carlo (Metropolis) approach.  
The implementation follows the methodology introduced in:

> **Wu, Z., Yin, F., Chen, Y., & Bei, H.**  
> *Predictive multiphase evolution in Al-containing high entropy alloys.*  
> Acta Materialia, 2017.  
> [DOI: 10.1016/j.actamat.2017.01.012](https://doi.org/10.1016/j.actamat.2017.01.012)

While the approach is inspired by the above paper, **all code in this repository was written entirely by me**.  
A key difference: in my case, **periodic boundary conditions (PBC)** are fully implemented to remove surface effects and ensure bulk-like behavior.

---

## Contents

- **`MC_BCC_CCA.py`**  
  Main script implementing the simulation pipeline:
  - Construction of a BCC lattice with two interpenetrating sublattices (α and β).
  - Randomized population of lattice sites according to target alloy composition.
  - Monte Carlo (Metropolis) swapping of nearest-neighbor pairs with simulated annealing (temperature ramp).
  - Energy calculation using nearest-neighbor enthalpy terms.
  - Output of order parameters (total and element-resolved LRO) and final configurations.

- **`corrected_vij_matrix.xlsx`**  
  Interaction parameter matrix \( V_{ij} \) derived from enthalpy of mixing data, normalized by BCC coordination number (z = 8).  
  Used in energy evaluation during Monte Carlo steps.

---

## Features

- **BCC lattice construction**: generates α and β sublattices explicitly.  
- **Initial disorder**: lattice populated randomly to match alloy composition.  
- **Monte Carlo with simulated annealing**: temperature-dependent acceptance of unfavorable swaps via Metropolis rule.  
- **Nearest-neighbor interactions only**: consistent with the enthalpy-based pair model.  
- **Periodic boundary conditions**: enforced in all three dimensions to eliminate surface artifacts and maintain correct neighbor count (z = 8 for every site).  
- **Output logging**: saves order parameter vs temperature, final configurations, and summary plots.
