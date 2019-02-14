# Normal Modes at Planetary Scales

This repository provide a highly parallel algorithm to compute the planetary normal modes. 
The elastic gravitational system is discretized using the Continous Galerkin mixed finite element method. 
A highly parallel Lanczos approach with polynomial filtering is utiilzed for solving 
the resulting generalized eigenvalue problem. 

Self-gravitation and rotation will be included in the future release. 

## How to intall it? 
Prerequisite: MPI, ParMetis, pEVSL.


## Reference
The repository provides scripts to generate planetary models for [the SuperComputing (SC'18) paper](https://dl.acm.org/citation.cfm?id=3291751), see below for details. 

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
Please report issues under this repository. If you have comments and suggestions, please reach me at jia.shi@rice.edu. 
