# Project Overview

This project studies how solvent environments alter the conformational flexibility of small organic molecules.

Pipeline:
1. Generate conformers with RDKit
2. Optimize with ASE + xTB
3. Filter by energy
4. Compute pairwise RMSD
5. Remove duplicates
6. Compute Boltzmann populations
7. Compute entropy/flexibility metrics
8. Compare across solvents

Current goals:


Important concepts:
- Flexibility is treated as an ensemble property
- Basin analysis groups structurally similar conformers
- Entropy is used as a thermodynamic flexibility metric