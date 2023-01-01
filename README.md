Planetary Normal Mode Computation
==========================================================
This repository provide a highly parallel algorithm to compute the normal modes of a rotating planet. 
We apply the continous Galerkin mixed finite element method to discretize the elastic-graviational system on tetrahedral meshes and
utilize a Lanczos approach with polynomial filtering for solving this generalized eigenvalue problem. 
We use fast multipole method to deal with self-gravitation. 

Prerequisite: MPI, ParMetis, pEVSL, ExaFMM.

To compile [the forked exafmm](https://github.com/js1019/exafmm-beta), the suggested configure options:
--------------------------------------------------------- 
./configure --enable-assert --enable-mpi --enable-openmp --enable-avx CXXFLAGS=-O3  
Since our code prefers to work with double precision, **do NOT enable float!**  You can then type make. 
