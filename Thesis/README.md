#General description
In this folder once can find all the necessary scripts and data to replicate our findings. The accompanying figures have also been included for comparison. Below follows a description of what can be found in every folder.


#Figures
The two figures that are part of the manuscript. Figure 1 is the workflow of the package and Figure 2 is the output of our example.


#Scripts
In order to replicate the prevalence estimates found in the manuscript, two scripts are provided:
- Estimate_prevalence.R (to run on any local machine)
- Estimate_prevalence_terminal.R (developed with the intent of running on a cluster, but can also be run locally)

The main difference between these two scripts is that in the first one the default values need to be manually updated, whereas in the latter script these values are set by the user when running the script for the first time.

There is a third script (Simulations.Rmd) that allows to replicate the simulations that were run as part of the manuscript. Keep in mind that due to the stochasticity of the process changes can be found when running the simulation script multiple times.


#Supplementary material
Three files can be found in this folder:
- S1 is the supplementary figure that explains the workflow of the simulation;
- S2 is a spreadsheet that contains the result of said simulation;
- S3 is the dataset used as an example in this manuscript and that users can take to verify the working of the package.
