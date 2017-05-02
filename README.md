YAFEL
=====

Linux, OSX: [![Build Status](https://api.travis-ci.org/tjolsen/YAFEL.svg?branch=master)](https://travis-ci.org/tjolsen/YAFEL)

Yet Another Finite Element Library

YAFEL is a general-purpose finite element library written entirely in C++14.
It is being developed primarily as a learning exercise and research code.
Since this library is primarily for personal use, features are added as
the need for them arises. If you have any good ideas for something to add, please get in
contact with me to discuss them.

A quick and dirty list of features that I currently (or plan to) support:
- Eigen3 for linear algebra data structures
- ViennaCL for GPU-accelerated iterative linear system solvers
- Mixed element topologies (eg: triangles + quads)
- Arbitrary-order interpolation
- Multiple dofs per node
- Continuous-galerkin DoFManager
- (in progress) Discontinuous-galerkin DoFManager
- ASCII Output to XML-based VTU files (most simply supports single-frame data)
- HDF5 + XDMF Output for paraview (supports time-dependent data) (feature-poor, but working)
- Mesh Import from Gmsh
- (future) Mesh creation for rectilinear meshes
- (future) Mesh creation for transfinite interpolation meshes, with boundary specification with Lambda


I'd like to give credit where it's due: this library was originally inspired by the Deal.II
library (https://github.com/dealii/dealii), which is maintained by Professor Wolfgang Bangerth
at Colorado State University (formerly Texas A&M).
I had the privilege of using his library in a class, and it inspired me to really get my
hands dirty learning what goes on in these codes.


Compilation and installation
============================

(Currently in flux. Requires some reasonably recent version of cmake.)

The library and applications should be compiled using CMake. The easiest
way to incorporate new applications is to create them in a subdirectory
inside the "apps" folder. Add your directory/executable to the apps/CMakeLists.txt
file, and link against the "yafel" target.
Other libraries requiring linking (HDF5) are added automatically by CMake, if available.

Compilation is periodically tested using GCC (5.4, 6.1) and Clang (3.8, 3.9, 4.0).

At the time of writing, the library is under rapid-enough development that it is not
worth installing into a system-wide directory (eg, /usr/local/).

Warning
=======

This library is currently being re-written, and since it has no
active users (to my knowledge) the re-work is happening on the master
branch. Feel free to check out a previous version if you need something
that actually works...