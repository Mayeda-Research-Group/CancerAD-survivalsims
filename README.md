# CancerAD-survivalsims
Repository with scripts for cancer-AD simulations for competing risk of death and selective survival. The repository also contains a folder with the input parameterization data in CSV files.

This project is a multistate simulation model developed to examine whether competing risk of death and/or selective survival can explain the observed inverse association between cancer and dementia in the literature. The scripts are summarized below.

1.	Rate inputs.R
This file loads CSVs with data used to parameterize the model, including mortality from US lifetables, cancer and dementia incidence, cancer relative survival, and dementia mortality HR (vs. no dementia).

2. Model.R
This is the data generation script. 
    o	Part 1 inputs calibration parameters. I did this manually so you will check the calibrations by looking at the calibration plots output in Script 4 below. 
    o	Part 2 defines the rate parameters of the model, all of which are functions of age (time). This is the most annoying part of the model to check!
    o	Part 3 define the actual differential equations, initial conditions, etc. 
    o	Part 4 runs the model.
    o	Part 5 is a commented out check to make sure the function runs without errors.
    o	Note that what we call “U” in the paper is called “S” in the code.
    
3.	Calibration plot function.R
This script creates a function that checks that the model rates and rate ratios match the input parameters. outputs the plots to visually inspect the calibration of each scenario. I commented out the portion of the script at the end where I wrote a loop to inspect each parameter sequentially.

4.	Calibration plot output.Rmd
This R markdown file takes a few minutes to knit, and produces an html output that shows all the checks I did and calibration plots for the model for all the scenarios and cancer types. 

5.	Run sims.R
This script writes a function that runs the simulation model and calculates results. Then, mapply is used to run the model and save results across all possible combinations of parameter inputs. These results are saved as an R object for use in Script #6, where results and figures in the paper are generated and extracted. 

6.	CR and SB sim results.Rmd
This is the script that generates figures used in the manuscript and extracts specific results that appear in the manuscript. 

