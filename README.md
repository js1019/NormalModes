# Normal Modes at Planetary Scales

This repository provides a highly parallel algorithm to compute the planetary normal modes. 
The elastic gravitational system is discretized using the Continous Galerkin mixed finite element method. 
A Lanczos approach with polynomial filtering is utiilzed for solving 
the resulting generalized eigenvalue problem. 

_Self-gravitation and rotation will be included in the future release._ 

## Normal Modes: what is the expected outcome?  
We show several modes that are computed from this work: [0S2](https://www.youtube.com/watch?v=DDfGHmqCMN0&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&index=5&t=0s), 
[0T2](https://www.youtube.com/watch?v=hxeDz8ncNH4), 
[3S9](https://www.youtube.com/watch?v=YR6N3AOTwoU&index=7&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s) and
[1T11](https://www.youtube.com/watch?v=XWY_dNAYAjE&index=6&list=PLUp2thaj3ruEVTLWazoRfqRK53t4hbYel&t=0s). 3S9 and 1T11 are illustrated below. 

<img src="figs/PREM3S9.gif" width="425"/> <img src="figs/PREM1T11.gif" width="425"/> 


## Preparation
**Input**: you may use the repository [PlanetaryModels](https://github.com/js1019/PlanetaryModels) to create your planetar models. 

### How to intall it? 
Prerequisite: MPI, ParMetis and pEVSL. Intel Math Kernal Library (MKL) is recommanded to use pEVSL. 

**Parallel Graph Partitioning** ([ParMetis](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)) is used for domain decomposition. To install ParMetis, you need have CMake installed and need to edit metis/includemetis.h and 
change **IDXTYPEWIDTH** to **64**. You then do 
~~~
make config; make;
~~~

**Parallel EigenValue Slicing** ([pEVSL](https://github.com/eigs/pEVSL)) is used to solve the generalized eigenvalue problem. 
You may use this [forked pEVSL version](https://github.com/js1019/pEVSL) for planetary normal mode computation. 
Please edit makefile.in for your cluster. 
We recommand users to use MKL. 
If you do not have MKL, we recommand users to use Openblas or Gotoblas. 


**This software**: once ParMetis and pEVSL are installed, you may install this software. 
Please edit makefile.in for right paths. 
We have two makefile.in examples for GNU and Intel compiler users. 
You may then go to src/ and type 
~~~
make clean; make; 
~~~

## How to run this application? 


## Reference
The repository provides codes to compute planetary normal modes for [our SuperComputing (SC'18) paper](https://dl.acm.org/citation.cfm?id=3291751), see below for details. 

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


## Contact 
Please report issues under this repository. Please send your suggestions and comments to jia.shi@rice.edu. 
