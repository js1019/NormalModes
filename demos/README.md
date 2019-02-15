## How to use this demo? 

Once you compiled the code. You should be able to run this demo. 

It provides two test examples. 
One is a constant solid elastic ball in models/input/CONST3k/. 
Another one is a small PREM in models/input/PREM3k/. 
You can directly use the current global_conf to run your first test!

In general, you will need to set up a few things: 
+ JOB = 1 or 2; 1 means no need for reference gravity; 
2 means with reference gravity, which needs to be precomputed use [PlanetaryModels](https://github.com/js1019/PlanetaryModels);
+ basename = _your mesh file name_; 
+ inputdir = _your input directory_; 
+ outputdir = _your output directory_;
+ CGMethod = .TRUE. _fixed_;
+ lowfreq = _the lowest eigenfrequency you need_;
+ upfreq = _the highest eigenfrequency you need_;
+ pOrder = _order of finite-element polynomial basis_; you may choose 1 or 2. 


For general applications, we note that you need to generate a corresponding model with the same _pOrder_ and reference gravity (if needed). 
You will then obtain all the eigenpairs in [lowfreq, upfreq]. 
