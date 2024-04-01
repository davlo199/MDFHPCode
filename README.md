# Multidimensional Fractional Hawkes Process Manuscript, Code, Data and Results
Code, data output files and manuscript for the results in "A Multidimensional Fractional Hawkes Process for Multiple Earthquake Mainshock Aftershock Sequences". 

## Compiling the C++ code
Most scripts incorporate a C++ mex file for the "ml.m" function, all of which are contained in the directory "code/MitLef". The file "LTInversionArray.mex" should be sufficient to run the .mex file. If not follow the instructions.

Some parts of the code have been re-implemented in C++ and you may need to run
```
mex -setup
```
You will need to compile the C++ code prior to running the code. If you wish to enable OpenMP threading for fine grain parallelism,
then you will need to pass the appropriate OpenMP option
(-fopenmp below for the GNU compiler) as well as the OpenMP library to link against (-lgomp in the case below).
The full compilation step then becomes:
```
mex CXXFLAGS="$CXXFLAGS -fopenmp -fPIC" -lgomp LTInversionArray.cpp
```
Note: [You will need to load the appropriate environment on Mahuika](#environment-on-mahuika)

If you cannot compile the C++ code, in line 100 of "ml.m" set cpp_code=0;

The function "MLapp.m" incorporates a Poincare asymptotic expansion for large time, and a Power series approximation for small time.

## Data Files
Two .csv files contain the data used from Japan and the Middle America Trench within the "data" directory. These are available publicly from [USGS quake search](https://earthquake.usgs.gov/earthquakes/search/).
For convenience they are provided here. The Japan data set is contained in the file "Japan.csv" and the Middle America Trench data set is contained in the file "MAT7623.csv". Data cleaning is done at the start of the scripts when the data is read in, for example lines 5-21 in "ETASTimeFit.R".

## Slurm parallelisation

The scripts in code/slurm_parallel are only required for using the batch array parallelisation in the following MATLAB functions. If a parallel computer is unavailable simply set ``PAR=0`` for serial computing in all of the input files.

## Parameter Estimation
Scripts for parameter estimation are contained in the "code/ParamEst" directory.

Parameter estimation for the MDFHP model is performed by using the "MDFHIntensityNewSum.m" function. 
Estimating the parameters for the Japan or Middle America Trench data sets was done in two stages. First use the input files being "EstJapan.m" or "EstMAT.m" respectively with the tolerances in line 172 of "MDFHIntensityNewSum.m" set to 1e-3.
We then use the rough minima again, with the tolerances in line 172 set to 1e-6, and change the input file values 
``PAR=0; Nrand=1;`` and 
MANx0 equal to the transformed output of the first step (i.e. the unconstrained optimised parameter values which is the variable "TransEst") to save computation time. Alternatively, the code as presented in this repository will return the more accurate minima although the computational expense is greatly increased (all estimation for the MDFHP is computationally infeasible on a desktop computer).

Estimation of the ETAS model's parameter was done using "ETASTimeFit.R" and follows D. Harte closely in their cited guide. 
To run it enter in the information in the input_files for the Japan or Middle America Trench data sets verbatim (specifically, ``dataname``, ``M0``, ``A`` and ``B``). Due to instabilities in the numerical optimisation, due to poor initial values, some models are clearly wrong. We select the best model as the one with the greatest log-likelihood and such that $\tau_{N(T)}\approx N(T)$ which is suggestive of the model being fit correctly. 

Estimation of the truncated exponential parameter estimates is done by using the "TruncatedExp.m" script using the same input files as for the "MDFHIntensityNewSum.m" function. For the ETAS model 
``Nband=1; Magvec=[10,M0];``

## Residual Analysis
Scripts for residual analysis are contained in the "code/Resid&Pred" directory.

Calculation of the transformed time residual process for the ETAS model is detailed in D. Harte's cited CRAN package "PtProcess". As an example, after loading in the output file "MATETAStimefit.RData" into R use
```
library("PtProcess")
resid<-residuals(x0)
```
returns the residual process for the ETAS model fit to the Middle America Trench data set, which can also be done for the Japan data set in a similar manner.

Calculation of the transformed time residual process for the MDFHP model is done in the "MultiDimResiduals.m" script.

## Information Gain
Scripts for information gain are contained in the "code/Resid&Pred" directory.

For the ETAS model $p_i$ is computed in the script "PredCapETASMarked.m" with the input file "PredJapan.m" or "PredMat.m" for the Japan or Middle America Trench data sets respectively.

For the MDFHP model $p_i$ is computed in the script "PredCapMDFHP.m" with the input file "PredJapan.m" or "PredMat.m" for the Japan or Middle America Trench data sets respectively.

Both models are compared to the empirical Poisson process in the script "MarkedIGPT.m" which can be run by loading the output files from "PredCapETASMarked.m" or "PredCapMDFHP.m".

## Miscellaneous and Plotting Scripts

To create figures 1,2 and 3 use the script "RplotsMDFHP.R".

Computation of the expected number of offspring (e.g. Equation 9) is done using the "MiscMDFHPscripts.m" script (lines 1-13). Furthermore, it is used to calculate the mean intensity at event times for each of the subprocesses and estimate $\mathbb{P}[M\in [4.35,5.35)]$ for discussion about the predictive performance of the MDFHP model fitted to the Middle American Trench data set (approximating Equation 7).

## Output files
The "output" directory is split into three subdirectories: "Estimates","Resid" and "Pred" for the output files of the second run of the parameter estimation procedure, transformed residual processes and information gain outputs respectively.
The follow table includes the file name and a brief description of its contents. 

| File Name | Description |
| --- | --- |
|JapanETAStimefit.RData| Output of "ETASTimeFit.R" for the Japan data set (once the fitted models were selected by hand). |
| JapanMDFHP55.mat| Output of "MDFHIntensityNewSum.m" for the MDFHP5.5 model fitted to the Japan data set.|
|JapanMDFHP575.mat| Output of "MDFHIntensityNewSum.m" for the MDFHP5.75 model fitted to the Japan data set.|
|JapanMDFHP6.mat| Output of "MDFHIntensityNewSum.m" for the MDFHP6 model fitted to the Japan data set.|
| MATETAStimefit.RData  | Output of "ETASTimeFit.R" for the Middle America Trench data set (once the fitted models were selected by hand).  |
| MATMDFHP435.mat  | Output of "MDFHIntensityNewSum.m" for the MDFHP4.35 model fitted to the Middle America Trench data set.  |
| MATMDFHP455.mat  | Output of "MDFHIntensityNewSum.m" for the MDFHP4.55 model fitted to the Middle America Trench data set.  |
| MATMDFHP475.mat  | Output of "MDFHIntensityNewSum.m" for the MDFHP4.75 model fitted to the Middle America Trench data set.  |
| Japan55SP1.csv | Transformed time residual process for subprocess 1 of the MDFHP5.5 fitted to the Japan data set.|
| Japan55SP2.csv | Transformed time residual process for subprocess 2 of the MDFHP5.5 fitted to the Japan Trench data set.|
| Japan575SP1.csv | Transformed time residual process for subprocess 1 of the MDFHP5.75 fitted to the Japan data set.|
| Japan575SP2.csv | Transformed time residual process for subprocess 2 of the MDFHP5.75 fitted to the Japan Trench data set.|
| Japan6SP1.csv | Transformed time residual process for subprocess 1 of the MDFHP6 fitted to the Japan data set.|
| Japan6SP2.csv | Transformed time residual process for subprocess 2 of the MDFHP6 fitted to the Japan Trench data set.|
| MAT435SP1.csv | Transformed time residual process for subprocess 1 of the MDFHP4.35 fitted to the Middle America Trench data set.|
| MAT435SP2.csv | Transformed time residual process for subprocess 2 of the MDFHP4.35 fitted to the Middle America Trench data set.|
| MAT455SP1.csv| Transformed time residual process for subprocess 1 of the MDFHP4.55 fitted to the Middle America Trench data set.|
| MAT455SP2.csv| Transformed time residual process for subprocess 2 of the MDFHP4.55 fitted to the Middle America Trench data set.|
| MAT475SP1.csv| Transformed time residual process for subprocess 1 of the MDFHP4.75 fitted to the Middle America Trench data set.|
| MAT475SP2.csv| Transformed time residual process for subprocess 2 of the MDFHP4.75 fitted to the Middle America Trench data set.|
|JapanMDFHP55IGPT.mat |Output of "MarkedIGPT.m" for the MDFHP5.5 fitted to the Japan data set.|
|JapanETASIGPT.mat |Output of "MarkedIGPT.m" for the ETAS fitted to the Japan data set.|
|JapanMDFHP575IGPT.mat |Output of "MarkedIGPT.m" for the MDFHP5.75 fitted to the Japan data set.|
|JapanMDFHP6IGPT.mat |Output of "MarkedIGPT.m" for the MDFHP6 fitted to the Japan data set.|
|MAT435mdfhpIGPT.mat | Output of "MarkedIGPT.m" for the MDFHP4.35 fitted to the Middle America Trench data set. |
|MATetasIGPT.mat |Output of "MarkedIGPT.m" for the MDFHP4.35 fitted to the Middle America Trench data set. |
|MAT455IGPT.mat | Output of "MarkedIGPT.m" for the MDFHP4.55 fitted to the Middle America Trench data set.|
|MAT475IGPT.mat | Output of "MarkedIGPT.m" for the MDFHP4.75 fitted to the Middle America Trench data set. |


## Versions

(For estimation) R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Platform: x86_64-pc-linux-gnu (64-bit)
(For residual analysis) R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
Platform: x86_64-w64-mingw32

On Mahuika HPC (used for data estimation and predictive performance), MATLAB 2020B.
On desktop computer (used for residual analysis), MATLAB 2022B.

## Manuscript

The "manuscript" directory contains the submitted unblinded manuscript pdf, along with a .zip file as compiled by latex.  

## Comments and/or Questions
These should be directed to the corresponding authour Louis Davis (davislrs2000@gmail.com).



