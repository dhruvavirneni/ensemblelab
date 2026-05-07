conformer and stability analysis
the same molecules, not isomers, can rotate along their single bonds
* mainly, torsion and dihedral angles define conformers
computational chemistry can be used to find minimum energy conformer
* optimizes for steric clashes (when multiple atoms or groups are forced together, forcing electron cloud repulsions and less stability), bad electrostatics, and strain

Methods:
AllChem.EmbedMultipleConfs()
* creates possible 3d geometries, adjusts bond angles and torsions, tries to generate reasonable structures

AllChem.UFFOptimizeMolecule()
* relaxes geometry and lowers energy, used on rough conformers


Output:
For each conformer:
geometry
energy
maybe RMSD (structural difference)

rank them
compare distributions



Pipeline: 
SMILES
→ RDKit conformers
→ ASE optimization
→ energy calculations
→ analysis