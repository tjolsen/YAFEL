Computational Methods for Plasticity
====================================

This is a collection of programs following the presentation in the book
Computational Methods for Plasticity, by de Souza Nieto.
It is a demonstration of the YAFEL library, along with a straightforward
implementation of the ideas presented in the text.
In general, these programs follow the theory presented in the text
rather than the suggestions on program structure, which I prefer to
discover/design on my own. Obviously, I will have read the textbook
implementation suggestions and incorporated what I feel are the good
parts.

Overall program structure:
-------------------------
- Program input data (mesh, materials, loads, boundary conditions)
- Time integration
  - For now, newton-raphson to solve nonlinear implicit quasi-static time increments
- Output (Configurable from input language/embedded lua/python?), dump all of this to HDF5 files
  - Displacements
  - Stress/strain fields
  - Mises stress

 

The organization of the programs is as follows:

- 01: Displacement-based static linear elasticity
- 02: Incremental small-strain plasticity
  - a) Isotropic von-Mises
  - b) Tresca plasticity (smoothing vs subdifferential?)
- 03: Large-strain nonlinear hyperelasticity
- 
