# SAMPLE: 

A DESCRIPTION

## Install

``` r
library("devtools")
devtools::install_github("yacinebenchehida/SAMPLE/Package")
```

## Dependencies

-   R (\>= 4.3.0)

SAMPLES requires: `ggplot2`, `dplyr`, `Rmisc`,`RColorBrewer`, `magrittr`.

## Example usage
### Run the full pipeline
``` r
library(SAMPLE)
data("coral_symbionts")
SAMPLE(input = coral_symbionts,output_N = "Example",replicates = 50,stability_thresh = 2,sucess_points = 10,diff = 1)
```

### Output
The pipeline generates two types of files:
a pdf: 
<img src="Figures/Example.png" width="70%" height="70%"/>

EXPLANATIONS

