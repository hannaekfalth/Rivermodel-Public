# Rivermodel Public
Public version of Rivermodel, a set of hydropower optimization models.
 
This repo contains the source code for the article *Exploring trade-offs between aggregated and turbine-level representations of hydropower in optimization models* [Manuscript submitted for publication].

The code constitutes a set of deterministic hydropower models in which the techno-physical detail can be varied, to explore trade-offs between computational time and accuracy in optimized power production. 

## USAGE
Note that the data needed to run the code is removed in this public version, due to confidentiality of the data. 

The function runmodel() is used to run the models. It takes a number of arguments, allowing to test hydropower representations with various detail. 
**runmodel(year, weeks, type, power, e, head, start)**

- **year** = any year with available data, e.g. 2019.  
- **weeks** = any unitrange of weeknumbers in a year, e.g. 1:52.
- **type** = [:LP, :NLP, :MIP]
- **power** = ["E constant head", "E taylor", "bilinear HeadE", "aggregated"]
- **e** = ["cv segments origo", "cv segments noseg", "constant eta", "ncv segments rampseg", "ncv segments zeroseg", "cv poly noseg", "ncv poly rampseg", "ncv poly zeroseg"]
- **head** = [:standardbasic, :constantmean, :constantmax, :aggregated]
- **start** = (type=:LP, power="E taylor", e="convex segments origo") or just (type=:LP,) to use same power & e as main run (but don't forget the last comma).

In the paper, the models evaluated are AMIP, A, B, B:L, C, D, and E. Here is the arguments used to run these models:  
- A_MIP: runmodel(2019, weeks=1:52, type=:NLP, power="bilinear HeadE", e="ncv poly rampseg", head=:standardbasic, start=(type=:MIP, power="bilinear HeadE", pe="ncv segments zeroseg"))
- A: runmodel(2019, weeks=1:52, type=:NLP, power="bilinear HeadE", e="ncv poly rampseg", head=:standardbasic, start=(type=:LP, power="E taylor", e="cv segments origo"))
- B: runmodel(2019, weeks=1:52, type=:NLP, power="bilinear HeadE", e="cv poly origo", head=:standardbasic, start=(type=:LP, power="E taylor", e="cv segments origo"))
- B:L: runmodel(2019, weeks=1:52, type=:LP, power="E taylor", e="cv segments origo", head=:standardbasic)
- C: runmodel(2019, weeks=1:52, type=:LP, power="E constant head", e="cv segments origo", head=:constantmean)	
- D: runmodel(2019, weeks=1:52, type=:LP, power="E constant head", e="constant eta", head=:constantmean)	
- E: runmodel(2019, weeks=1:52, type=:LP, power="aggregated")

## Content
- *Rivermodel.jl* contains the function to run the model, runmodel(), and the solver settings  
- The actual model formulations are found in the file *Jumpmodel.jl.*, except for the aggregated model, model E, which is found in *aggmodel.jl*  
- The input data are read in *rawdata.jl* and *skelleftedata.jl*  
- The model results are saved in *output.jl*  
- *compareruns.jl* and *plots.jl* are used for plotting and analyzing the model outputs









 
 
