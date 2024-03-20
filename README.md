# MDFHPCode
Code and data used for the results in "A Multidimensional Fractional Hawkes Process for Earthquakes". 

## Compiling the C++ code
Most scripts incorporate a C++ mex file for the "ml.m" function. The file "LTInversionArray.mex" **check correct name** should be sufficient to run the .mex file. If not follow the instructions.

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
Two .csv files contain the data used from Japan and the Middle America Trench. These are available publicly from [USGS quake search](https://earthquake.usgs.gov/earthquakes/search/).
For convenience they are provided here. The Japan data set is contained in the file "Japan.csv" and the Middle America Trench data set is contained in the file "MAT7623.csv".

## Running the scripts

**Something about directories etc**

## Parameter Estimation
Parameter estimation for the MDFHP model is performed by using the "MDFHIntensityNewSum.m" function. 
To estimate the parameters for Japan or Middle America Trench data sets use the input files being "EstJapan.m" or "EstMAT.m" respectively.

Estimation the ETAS model's parameter was done using "FitTimeETAS.R" and follows D. Harte closely in their cited guide. 
To run it, enter in the information in the input_files for the Japan or Middle America Trench data sets verbatim (specifically, "dataname", "M_0", "A" and "B"). 

Estimation of the truncated exponential parameter estimates is done by using the "TruncatedExp.m" script using the same input files as for the "MDFHIntensityNewSum.m" function.

## Residual Analysis

Calculation of the transformed time residual process for the ETAS model is detailed in D. Harte's cited CRAN package "PtProcess".

Calculation of the transformed time residual process for the MDFHP model is done in the "MultiDimResiduals.m" script.

## Information Gain

For the ETAS model $p_i$ is computed in the script "PredCapETASMarked.m" with the input file "PredJapan.m" or "PredMat.m" for the Japan or Middle America Trench data sets respectively.

For the MDFHP model $p_i$ is computed in the script "PredCapMDFHP.m" with the input file "PredJapan.m" or "PredMat.m" for the Japan or Middle America Trench data sets respectively.

Both models are compared to the empirical Poisson process in the script "MarkedIGPT.m".






