# Normal Modes at Planetary Scales
<img src="https://img.shields.io/github/issues/js1019/NormalModes.svg"><img src="https://img.shields.io/github/languages/code-size/js1019/NormalModes.svg"> <img src="https://img.shields.io/github/forks/js1019/NormalModes.svg"> <img src="https://img.shields.io/github/stars/js1019/NormalModes.svg"> <img src="https://img.shields.io/github/license/js1019/NormalModes.svg">


This repository applied a combination of several highly parallel algorithms to compute the planetary interior normal modes. 
The elastic-gravitational system is discretized using the Continuous Galerkin mixed finite element method. 
A Lanczos approach with polynomial filtering is utilized for solving 
the resulting generalized eigenvalue problems. 

_Self-gravitation and rotation will be included in the future release._ 

**News: This work has been selected to serve as the benchmark for [the student cluster competition](http://www.studentclustercompetition.us/) reproducibility challenge at [SC'19](https://sc19.supercomputing.org). 
Please see [the official announcement](https://sc19.supercomputing.org/2019/04/19/from-sc-papers-to-student-cluster-competition-benchmarks-joining-forces-to-promote-reproducibility-in-hpc/) 
and [teams named](https://sc19.supercomputing.org/2019/06/10/teams-named-for-sc19-student-cluster-competition/).** 
If you look for the code for the competition, please see SCC19 version in the release. 

## Outcome of this work
Many different modes are expected. Please see below for selected normal modes computed from a standard spherically symmetric Earth model.
<p align="center">
<img src="figs/PREM2M_J2_0S2.png" width="200"/> <img src="figs/PREM2M_J2_0S3.png" width="200"/> <img src="figs/PREM2M_J2_0S5.png" width="200"/> <img src="figs/PREM2M_J2_0S7.png" width="200"/> 
</p>
<p align="center">
<img src="figs/PREM2M_J2_1S3.png" width="200"/> <img src="figs/PREM2M_J2_1S4.png" width="200"/> <img src="figs/PREM2M_J2_1S6.png" width="200"/> <img src="figs/PREM2M_J2_1S9.png" width="200"/> 
</p>
<p align="center">
<img src="figs/PREM2M_J2_2S3.png" width="200"/> <img src="figs/PREM2M_J2_3S3.png" width="200"/> <img src="figs/PREM2M_J2_4S3.png" width="200"/> <img src="figs/PREM2M_J2_5S3.png" width="200"/> 
</p>
<p align="center">
<img src="figs/PREM2M_J2_0T2.png" width="200"/> <img src="figs/PREM2M_J2_1T1.png" width="200"/> <img src="figs/PREM2M_J2_1T2.png" width="200"/> <img src="figs/PREM2M_J2_2T4.png" width="200"/> 
</p>


## A note on the design of this application 
It is not straightforward to compute the normal modes of a fully heterogeneous planet. 
We have to deal with multiphysics at large scales, 
which computationally involves many techniques in numerical linear algebra, finite element method, computer graphics, high-performance computing, etc. 
We do not expect that we can obtain everything via a single click. 
However, to solve this complicated problem, we divide the original one into several smaller subproblems. 
Indeed, via solving each subproblem separately, it actually simplifies our work significantly. 
We develop three repositories to provide a solid solution for this application:
+ [PlanetaryModels](https://github.com/js1019/PlanetaryModels): It is a planet model builder and provides scripts to generate planetary models on fully unstructured tetrahedral meshes. It provides a simple and flexible way to create a 3D body with almost arbitrary exterior and interior shapes.   
+ [pEVSL](https://github.com/eigs/pEVSL): It is a parallel eigenvalue slicing library. It provides general algebraic algorithms to solve large-scale generalized Hermitian extreme and interior eigenvalue problems. 
You may use [the forked pEVSL version](https://github.com/js1019/pEVSL) 
for this application, since it contains several modifications for your convenience. 
+ [NormalModes](https://github.com/js1019/NormalModes): It builds up the matrices from the elastic-gravitational system using finite element method and utilizes several external libraries, including [pEVSL](https://github.com/eigs/pEVSL), to solve for the normal modes of the planetary model generated by [PlanetaryModels](https://github.com/js1019/PlanetaryModels). 

The separated repositories also provide flexibility to extend our work for other purposes. 
Please let us know if you plan to utilize what we have to your work.  



## How to run this application? 
Please follow INSTALL.md to install the application. 
Please check the demos/global_conf, which shows an **extremely simple** parameter setting. 
Since the problem is deterministic, there are only a few parameters that are needed to compute the normal modes. 
You can then obtain _all_ the eigenpairs in the prescribed frequency interval. 
The values of eigenfrequencies will be shown at the end of the computation as well as their relative errors, i.e., ||Ax-&lambda;Bx||/||&lambda;||, which is typically around **10^-13**. 
The eigenfunctions will be saved in the binary format. Please check the README.md under demos/ for more details. 

**Tips**: Please always check the performance and scalability before running large-scale applications. 

**Visualization**: You can use scripts in [PlanetaryModels](https://github.com/js1019/PlanetaryModels) and [Paraview](https://www.paraview.org/) to visualize your results. Here, we show animations created by Paraview: [0S2](https://www.youtube.com/watch?v=DDfGHmqCMN0&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&index=5&t=0s), 
[0T2](https://www.youtube.com/watch?v=hxeDz8ncNH4), 
[3S9](https://www.youtube.com/watch?v=YR6N3AOTwoU&index=7&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s) and
[1T11](https://www.youtube.com/watch?v=XWY_dNAYAjE&index=6&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s). 3S9 and 1T11 are illustrated below. 

<p align="center">
<img src="figs/PREM3S9.gif" width="400"/> <img src="figs/PREM1T11.gif" width="400"/> 
</p>

You can also use the [*NMPostProcess*](https://github.com/harrymd/NMPostProcess) library to visualize modes and automatically identify them.

## Reference
+ [_**Theory, discretization and validation**_]: Jia Shi, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "A non-perturbative approach to computing seismic normal modes in rotating planets." Journal of Scientific Computing, 91:67, 2022, [the paper link](https://doi.org/10.1007/s10915-022-01836-5).
+ [_**Reproducibility, summary and discussion of the student cluster competition results**_]: Jia Shi, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "Planetary normal mode computation: Parallel algorithms, performance, and reproducibility." IEEE Transactions on Parallel and Distributed Systems, 32, no. 11 (2021): 2609-2622, [the paper link](https://ieeexplore.ieee.org/abstract/document/9319555). 
+ [_**Parallel performace and algorithms**_]: Jia Shi, Ruipeng Li, Yuanzhe Xi, Yousef Saad, and Maarten V. de Hoop. "Computing planetary interior normal modes with a highly parallel polynomial filtering eigensolver." In SC18: International Conference for High Performance Computing, Networking, Storage and Analysis, pp. 894-906. IEEE, 2018, [the paper link](https://dl.acm.org/citation.cfm?id=3291751).


## Contact 
Please report issues under this repository. Contributions are welcome. 
