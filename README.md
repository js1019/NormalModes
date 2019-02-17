# Normal Modes at Planetary Scales
This repository applied a combination of several highly parallel algorithms to compute the planetary interior normal modes. 
The elastic-gravitational system is discretized using the Continous Galerkin mixed finite element method. 
A Lanczos approach with polynomial filtering is utilized for solving 
the resulting generalized eigenvalue problems. 


_Self-gravitation and rotation will be included in the future release._ 

## Normal Modes: what will be the expected outcome?  
We show several modes that are computed from this work: [0S2](https://www.youtube.com/watch?v=DDfGHmqCMN0&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&index=5&t=0s), 
[0T2](https://www.youtube.com/watch?v=hxeDz8ncNH4), 
[3S9](https://www.youtube.com/watch?v=YR6N3AOTwoU&index=7&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s) and
[1T11](https://www.youtube.com/watch?v=XWY_dNAYAjE&index=6&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s). 3S9 and 1T11 are also illustrated below. 

<img src="figs/PREM3S9.gif" width="425"/> <img src="figs/PREM1T11.gif" width="425"/> 

Many different modes are expected. 

## Preparation
**Input**: you may use the repository [PlanetaryModels](https://github.com/js1019/PlanetaryModels) to create your planetary models. 

### How to install it? 
Prerequisite: MPI, ParMetis and pEVSL. Intel Math Kernal Library (MKL) is recommanded to use pEVSL. 

It runs on multi-CPU platforms. Currently, this code has been tested using GNU, Intel and Cray compilers and scaled up to **40k** processes. 

**Parallel Graph Partitioning** ([ParMetis](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)) is used for domain decomposition. To install ParMetis, CMake is needed and please edit metis/includemetis.h and 
change **IDXTYPEWIDTH** to **64**. You can then do 
~~~
make config; make;
~~~
Since CMake is commonly installed in modern clusters, the installation of ParMetis will be easy and simple. 

**Parallel EigenValue Slicing** ([pEVSL](https://github.com/eigs/pEVSL)) 
is used to solve the generalized eigenvalue problem. 
You may use [the forked pEVSL version](https://github.com/js1019/pEVSL) 
for this application, since it contains several modifications for your convenience. 
Please edit makefile.in for your cluster. 
We recommand users to use MKL and [Intel link line advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) 
to find the right links for you. 
If you do not have MKL, we recommand users to use Openblas or Gotoblas. 

The installation of pEVSL will be easy and simple. Please check the example to make sure that it is installed correctly. 


**This software**: once ParMetis and pEVSL are installed, you can install this software. 
Please edit makefile.in for right paths. 
We have two makefile.in examples for GNU and Intel compiler users. 
You may then go to src/ and type 
~~~
make clean; make; 
~~~
The installation will also be easy and simple. 


## How to run this application? 
Please check the demos/global_conf, which shows an **extremely simple** parameter setting. 
Since the problem is deterministic, there are only a few parameters that are needed to compute the normal modes. 
You will then be able to obtain _all_ the eigenpairs in the prescribed frequency interval. 
The values of eigenfrequencies will be shown at the end of the computation as well as their relative errors, i.e., ||Ax-&lambda;Bx||/||&lambda;||, which is typically around 10^-10. 
The eigenfunctions will be saved in the binary format. Please check the README.md under demos/ for more details. 

**Tips**: please always check the performance and scalability before you run large-scale applications. 

**Visualization**: you will need to use scripts in [PlanetaryModels](https://github.com/js1019/PlanetaryModels) and [Paraview](https://www.paraview.org/)
to visualize your results. 

## Reference
The repository and pEVSL contain codes to compute planetary normal modes for [our SuperComputing (SC'18) paper](https://dl.acm.org/citation.cfm?id=3291751), see below. 
~~~
@inproceedings{shi2018computing,
  title={Computing planetary interior normal modes with a highly parallel polynomial filtering eigensolver},
  author={Shi, Jia and Li, Ruipeng and Xi, Yuanzhe and Saad, Yousef and de Hoop, Maarten V},
  booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis, {SC}'18, Dallas, TX, November 11-16, 2018},
  pages={71:1--71:13},
  year={2018},
  organization={ACM/IEEE}
}
~~~
Parallel strong and weak scalabilities are shown in this paper. 

## Contact 
Please report issues under this repository. Please send your suggestions and comments to jia.shi@rice.edu. Contributions are welcome. 
