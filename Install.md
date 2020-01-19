## Preparation and installation
**Input**: You need to use the repository [PlanetaryModels](https://github.com/js1019/PlanetaryModels) to create your input planetary models. Several test models are provided in the demos. 

### How to install it? 
Prerequisite: MPI, ParMetis and pEVSL. Intel Math Kernel Library (MKL) is recommended to use pEVSL. 
It runs on multi-CPU platforms. Currently, this code has been tested using GNU, Intel and Cray compilers and scaled up to **40k** processes. 

**Parallel Graph Partitioning** ([ParMetis](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)) is used for domain decomposition. To install ParMetis, CMake is needed and please edit metis/include/metis.h and 
change **IDXTYPEWIDTH** to **64**. You can then do 
~~~
make config; make;
~~~
Since CMake is commonly installed in modern clusters, the installation of ParMetis will be easy and simple. 

**Parallel EigenValue Slicing Library** ([pEVSL](https://github.com/eigs/pEVSL)) 
is used to solve the generalized eigenvalue problem. 
You may use [the forked pEVSL version](https://github.com/js1019/pEVSL) 
for this application, since it contains several modifications for your convenience. 
Please edit makefile.in for your cluster. 
We recommend users to use MKL and [Intel link line advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) 
to find the right links for you. 
If you do not have MKL, we recommend users to use OpenBLAS or GotoBLAS. 
The installation of pEVSL will be easy and simple. Please check the example to make sure that it is installed correctly. 


**This work**: Once ParMetis and pEVSL are installed, you can install this software. 
Please edit makefile.in for right paths. 
We have two makefile.in examples for GNU and Intel compiler users. 
You may then go to src/ and type 
~~~
make
~~~
The installation will also be easy and simple. 
