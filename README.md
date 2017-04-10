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
- Mixed element topologies (eg: triangles + quads)
- Arbitrary-order interpolation
- Multiple dofs per node
- Continuous-galerkin DoFManager
- (in progress) Discontinuous-galerkin DoFManager
- (in progress) ASCII Output to XML-based VTU files (most simply supports single-frame data)
- (in progress) HDF5 + XDMF Output for paraview (supports time-dependent data)
- Mesh Import from Gmsh
- (future) Mesh creation for rectilinear meshes
- (future) Mesh creation for transfinite interpolation meshes, with boundary specification with C++11 Lambda


I'd like to give credit where it's due: this library was originally inspired by the Deal.II
library (https://github.com/dealii/dealii), which is maintained by Professor Wolfgang Bangerth
at Colorado State University (formerly Texas A&M).
I had the privilege of using his library in a class, and it inspired me to really get my
hands dirty learning what goes on in these codes.


Compilation and installation
============================

(Currently in flux. Requires some reasonably recent version of cmake.)

I recommend using the Clang compiler rather than GCC because Clang seems
to do a better job handling template-heavy code.
This library contains an implementation of arbitrary-order tensors that
makes extensive use of variadic templates and relies on compile-time computation
to generate efficient code.
Clang (3.8, 4.0) seems to do a better job than GCC (5.4, 6.0) at the time of writing
at fully inlining and simplifying the code.


At the time of writing, the library is under rapid-enough development that it is not
worth installing into a system-wide directory (eg, /usr/local/).

Using the library
=================

(Currently in flux. Link against libyafel. See Compilation section above.)

Warning
=======

This library is currently being re-written, and since it has no
active users (to my knowledge) the re-work is happening on the master
branch. Feel free to check out a previous version if you need something
that actually works...