# Smoothing for Age-Period-Cohort models: a comparison between splines and random process

## Code

There are three items in this section. Two `.R` files and a folder for all the model fitting Code.

- `functions.R`: This is script for functions to fit the fitting the models, calculate model scores, etc...
- `saveData.R`: This is a script that can be run to generate the simulated data. It was used to created the file `Data/allSimulatedData.rds`.

### Model Fitting

Files for model fitting and generating the plots

- `whyIncludePenalty.R` and `pcPriorsDifference.R`: Creates the plots in Figure 2 in the manuscript.
- `fitSimulatedModels.R`: Run the simulation study and collect the results.
- `edaDeathsExample.R`: Run the analysis and model fitting for the alcohol and self harm related deaths example.
- `differentBasisExamples.R`: Creates the plots in Figure S1 in the Supplementary Material.

## Data

This folder contains two data files. The simulated data is in the `allSimulatedData.rds` file, and the alcohol and self harm related deaths data is in the `mortalityExamplesData.csv` file. 