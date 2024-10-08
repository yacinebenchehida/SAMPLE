# SAMPLE: 

SAMPLE is an R package that can benefit many community ecologists adjust their sampling efforts when studying prevalence rates. The simplicity of SAMPLE allows researchers to run it in the field and dynamically adjust sampling schemes as needed.

## Install

``` r
library("devtools")
devtools::install_github("yacinebenchehida/SAMPLE")
```

## Dependencies

-   R (\>= 4.3.0)

SAMPLES requires: `ggplot2`, `dplyr`, `Rmisc`,`RColorBrewer`, `magrittr`.

## Example usage
### Run the full pipeline
``` r
library(SAMPLE)
data("coral_symbionts")
set.seed(812)
SAMPLE(input = coral_symbionts, output_N = "Example", replicates = 50, stability_thresh = 2, success_points = 10, diff = 1)
```

### Output
The pipeline generates two types of files:
- A pdf showing how the prevalence rate changes as the sampling increases;

<img src="Figures/Example.png" width="90%" height="90%"/>

- A text file showing for each host species and each taxa the sample size at which the prevalence rate becomes stable.

```
Host_species             Taxa                        Prevalence        thres_stability
Agaricia agaricites      Opecarcinus hypostegus      20.5              4
Agaricia lamarcki        Opecarcinus hypostegus      68.1818181818182  7
Acropora palmata         Domecia acanthophora        88                8
Acropora palmata         Spirobranchus polycerus     43                12
Millepora complanata     Domecia acanthophora        17.5555555555556  14
Millepora complanata     Acanthemblemaria spinosa    21.5384615384616  9
Millepora complanata     Megabalanus stultus         21.3846153846154  9
Orbicella faveolata_06m  Troglocarcinus corallicola  59.6923076923077  9
Orbicella faveolata_15m  Troglocarcinus corallicola  38.7142857142857  10
```

- The plot below explains what the different elements present on the plots represent:
<img src="Figures/Explanations.png" width="70%" height="70%"/>


### How to run the pipeline a step at a time

```
#) 1) Load the data
data("coral_symbionts")

# 2) Run the permutations
set.seed(812)
perm = RunPerm(input = coral_symbionts,replicates = 50)

# 3) Assess the stability based on the performed permutations
stable = stability(data = perm, stability_thresh = 2, success_points = 10, diff = 1)

# 4) Plot the results
plotstab(data = perm, info = stable, outputName = "Stability_example")
```
