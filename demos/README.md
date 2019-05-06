## How to use this demo? 

Once you compile the code, you should be able to run applications in this demo. 
We have a script `Dev.sh` to run this code via using SLURM. 
You can set up the number of MPI processes and threads, etc. 

It provides four test examples: 
- a constant solid elastic ball in models/input/CONST3k/; 
- a small standard Earth model in models/input/PREM3k/; 
- a small Mood model in models/input/Mtopo100k/; 
- a small Mars model in models/input/RTMDWAK8k/. 
You can directly use the current global_conf to run your first test! 
You can run JOB = 1 or 2 and pOrder = 1 or 2 with various lowfreq and upfreq.
For the Moon model, you can only choose pOrder = 1. If you want to test pOrder = 2, 
please generate the model using the demo in  [PlanetaryModels](https://github.com/js1019/PlanetaryModels).

In general, you will need to set up a few things:
~~~ 
JOB = 1 or 2; 1 means no need for reference gravity; 2 means with reference gravity;
basename = your mesh file name; 
inputdir = your input directory; 
outputdir = your output directory;
lowfreq = the lowest eigenfrequency you need;
upfreq = the highest eigenfrequency you need;
pOrder = order of finite-element polynomial basis; you may choose 1 or 2. 
~~~

For general applications, we note that you need to generate a corresponding model with the same _pOrder_ 
and reference gravity (if needed), which needs to be precomputed using [PlanetaryModels](https://github.com/js1019/PlanetaryModels). 
You will then obtain all the eigenpairs in [lowfreq, upfreq]mHz. 
Note that the actual bounds are (2&pi;*[lowfreq, upfreq])^2*10^-6. 
Please make sure that lowfreq (smallest is 0.0) and upfreq 
are within the eigenvalue bound of the problem. 


