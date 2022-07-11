# Local Force Field Adjustment of Structure Mutations in iCn3d

### Introduction

The program iCn3D is a web application to vizualize and analyze structural data. One feature is to introduct mutations to the protein to compare the original structure and mutant variant. Currently, the structure is processed using scap, a program to predict side-chain conformations, and render the results. This project seeks to add the following features to the iCn3D mutation analysis:

1. Backbone minimization
2. Implicit solvent molecular dynamics (MD) simulation 

This done using RESTful API calls to OpenMM, a molecular simulation program (https://openmm.org/). 

### Current Tasks

1. Implement a RESTful API of OpenMM for iCn3D to work with. 
2. Design a protocol for backbone minimization and implicit solvent MD simulation.
3. Implement and validate the protocol and API calls for OpenMM.
4. Integrate into iCn3D. 
