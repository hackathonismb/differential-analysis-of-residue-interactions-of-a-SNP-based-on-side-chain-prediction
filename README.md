# Local Force Field Adjustment of Structure Mutations in iCn3d

### Introduction

The program iCn3D is a web application to vizualize and analyze structural data. One feature is to introduct mutations to the protein to compare the original structure and mutant variant. Currently, the structure is processed using scap, a program to predict side-chain conformations, and render the results. 

This project seeks to offer additional improvements to the mutation structure estimation to bring the total number of estimated mutation structure files a user can generate to three:

1. Scap side-chain conformation (implemented prior to this project).
2. Everything in (1) with backbone energy minimization. 
3. Everything in (2) with MD simulation to further estimate the energy minimum. 

To aid in differential analysis, information on annotated functional regions of proteins and their distance to a mutation will also be provided as a new feature.

### Methodology

PDB files for the wildtype and mutant are first analyzed to determine if any non-standard residues are present within the local vicinity of the mutation. If so, the distance of these is reported in a csv format. Water molecules and non-standard residues are then removed from the PDB files as a preprocessing step for GROMACS. *Rest of the methodology is being finalized.*

### Current Tasks and Progress

More details on the status of current tasks can be found in the wiki: https://github.com/hackathonismb/differential-analysis-of-residue-interactions-of-a-SNP-based-on-side-chain-prediction/wiki 
