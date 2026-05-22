# Research Log

## Entry 1

# Goal: Initialize file structure and clarify project pipeline + start basin analysis

# What I tried:
- Shannon entropy derived from Boltzmann population
# Problems:
- Issues with keeping track of conformer list from initial state to filtered, unique, etc.
- Must be streamlined 
# Decisions:
- Opted to use functional methods to separate all pipeline steps
- Pipeline found in testing.ipynb developed in notebook was adapted for this study
- Used Shannon entropy as the next obvious step from conformer Boltzmann population calculations
# Observations:

# Questions:
- Boltzmann distribution → how should it properly be visualized and used in flexibility metrics?

# Next Steps:
- Implement basin clustering analysis
- Population analysis
- Torsional analysis


## Entry 2

# Goal: Add basin clustering analysis

# What I tried:
Replacing unique_conformers initialization with basin compilation logic (based on RMSD)
Takes ‘template’ molecules and groups based on other conformers’ similarities to them, then evaluates total populations of each basin group
# Problems:
None
# Decisions:
Decided to keep unique conformer list to maintain functionality of entropy calculations and following implementations
# Observations:
Flexible molecules have more basins with relatively lower populations (not a new observation, but proven by testing current iteration of implementation) 
# Questions:

# Next Steps:
Torsion analysis
Migrate to separate file function framework
Solvent interaction additions
