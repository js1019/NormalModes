## How to use this demo? 

Once you compiled the code. You should be able to run this demo. 

It provides two test examples. 
One is a constant solid elastic ball in models/input/CONST3k/. 
Another one is a small PREM in models/input/PREM3k/. 
You can directly use the current global_conf to run your first test!

In general, you will need to set up a few things:
~~~ 
JOB = 1 or 2; 1 means no need for reference gravity; means with reference gravity;
basename = your mesh file name; 
inputdir = your input directory; 
outputdir = your output directory;
CGMethod = .TRUE., which is fixed;
lowfreq = the lowest eigenfrequency you need;
upfreq = the highest eigenfrequency you need;
pOrder = order of finite-element polynomial basis; you may choose 1 or 2. 
~~~

For general applications, we note that you need to generate a corresponding model with the same _pOrder_ 
and reference gravity (if needed), which needs to be precomputed use [PlanetaryModels](https://github.com/js1019/PlanetaryModels). 
You will then obtain all the eigenpairs in [lowfreq, upfreq]. 
Note that the actual bounds are (lowfreq* 2 *\pi)^2 to (upfreq* 2 *\pi)^2. 

