# I) Generate the simulated data

The R code below was used to generated the simulations:
```r
setwd("/Users/yacinebenchehida/Desktop/Other/Henrique/Crab simulations")
system("mkdir -p Simulations")
setwd("Simulations")
Pop_size = c(100,1000,10000) # Create vector of population size of 100, 1000 and 10000 individuals
Prevalence = seq(10,90,10) # create a vector of prevalence rates of 10 to 90 %

for (Pop in Pop_size){
  for (Preval in Prevalence){
    Present = Pop * Preval/100
    Absent = Pop - Present
    simdata = as.data.frame(cbind("Host1",c(rep(1,Present),rep(0,Absent))))
    data = simdata[sample(1:nrow(simdata),replace=FALSE), ]
    colnames(data) = c("Host_name","Symbiont_name")
    name_data = paste("Sim_data_Pop_size_",Pop,"_Preval_",Preval,".txt",sep="")
    write.table(x = data,file = name_data,sep="\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
  }
}
```

For each population size and prevalence rates, the R code above generates a text files that look like this: 
```
Host_name  Symbiont_name
Host1      0
Host1      0
Host1      0
Host1      0
Host1      0
Host1      1
Host1      1
Host1      0
Host1      0
...        ...

```

#  II) Assess impact of prevalence rate, number of replicates, and population size
## A) Run the pipeline on each data set generated in I)

The commands below run pipeline on each dataset generated in **I)** while varying the number of replicates from 10 to 500. The values of the number of successive points, the mean-difference and max-mean-difference parameters were kept at their default values:
```bash
for i in $(ls Sim_data*txt|tac); do
  for j in 10 20 30 40 50 60 70 80 90 100 200 500; do
    mkdir -p "$j"_replicates; Pop_size=$(echo $i|cut -d _ -f 5); Preval=$(echo $i|cut -d _ -f 7|perl -pe 's/\.txt//g')
    echo $i
    time Rscript ../Estimate_prevalence_terminal.R $i $j 2 Results_Simu_PopSize_"$Pop_size"_Preval_"$Preval" 10
    mv Res* "$j"_replicates
  done
done
```

For each set of parameters, the pipeline returns a two lines text file (and a pdf) with the results. The text file looks like this:
```
Host_species  Taxa           Prevalence  thres_stability
Host1         Symbiont_name  48.5        51
```

## B) Extract the results 

The main results of the simulations could be visualised using the following command:
```bash
(for i in Results*txt; do
  Size=$(echo $i|cut -d _ -f 4)
  Preval=$(echo $i|cut -d _ -f 6|perl -pe 's/\.txt//g')
  Obs_Preval=$(cat $i|awk 'NR > 1 {print $3}'|perl -pe 's/([0-9]+)(\.)(\d\d).+/$1$2$3/g')
  Stability=$(cat $i|awk 'NR > 1 {print $4}')
  echo -e Population size: $Size ind \| Simulated Prevalence: $Preval% \| Inferred Prevalence: $Obs_Preval% \| Stable from: $Stability Samples
done)|column -t
```

The results look like this:
```
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  50.09%  |  Stable  from:  34   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  48.36%  |  Stable  from:  13   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  50.81%  |  Stable  from:  13   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  48.98%  |  Stable  from:  93   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  49%     |  Stable  from:  35   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  48.85%  |  Stable  from:  66   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  49.26%  |  Stable  from:  29   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  51.73%  |  Stable  from:  37   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  49.36%  |  Stable  from:  48   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  49.87%  |  Stable  from:  39   Samples
Population  size:  100  ind  |  Simulated  Prevalence:  50%  |  Inferred  Prevalence:  49.71%  |  Stable  from:  12   Samples

```

## C) Gather the results for Table S2 (part A)
The command used to summarise the results in Supplementary Table S2 (part A)
```bash
(for pop in 100 1000 10000; do
  for j in 10 20 30 40 50 60 70 80 90 100 200 500; do
    cd "$j"_replicates
    for i in *_"$pop"_*.txt; do
      Samp=$(echo $i|cut -d _ -f 4|perl -pe 's/\.txt//g')
      Prev=$(echo $i|cut -d _ -f 6|perl -pe 's/\.txt//g')
      rate=$(cat $i|awk 'NR>1 {print $4}')
      prevalence=$(cat $i|awk 'NR>1 {print $3}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g')
      echo -e $Samp"\t"$Prev"\t"$rate"|"$prevalence"%"
    done|awk '{print $3}'|perl -pe 's/\n/\t/g'
    cd ..
    echo -e "\n"
  done|perl -pe 's/^\n$//g'
done) |column -t
```

#  III) Assess impact of successive points, mean-difference, max-mean-difference
## A) Run the pipeline

The commands below run pipeline on each dataset generated in **I)** while varying the number successive points (2, 10 and 50), the mean-difference (1, 2, 5, and 10) and max-mean-difference (0.5, 1, 2). The number of replicates is kept to 50 (default):
```bash
for succ in 2 10 50; do
  for meandiff in 1 2 5 10; do
    for minmax in 0.5 1 2; do
      Rscript ../Estimate_prevalence_terminal.R Sim_data_Pop_size_1000_Preval_50.txt Results_Simu_PopSize_100_Preval_50_succ_"$succ"_meandiff_"$meandiff"_minmaxdiff_"$minmax" 50 $meandiff  $succ $minmax
    done
  done
done
```

For each set of parameters, the pipeline returns a two lines text file (and a pdf) with the results. The text file looks like this:
```
Host_species	Taxa	Prevalence	thres_stability
Host1	Symbiont_name	51.1428571428571	13
```

## B) Gather the results for Table S2 (part B)

The command used to summarise the results in Supplementary Table S2 (part B):
```bash
for i in $(ls *succ_2_*.txt|sort -V); do
  rate=$(cat $i|awk 'NR>1 {print $4}')
  prevalence=$(cat $i|awk 'NR>1 {print $3}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g')
  echo -e $rate"|"$prevalence"%"
done|awk -v n=3 '{a[NR]=$0}END{ x=1; while (x<=n){ for(i=x;i<=length(a);i+=n) printf "%s",a[i]"\t"; print ""; x++; } }' |column -t
```

#  IV) Assess consistence of the results for 10 replicates using the same sets of parameters
## A) Run the pipeline:
The commands below run 10 times the pipeline on the dataset generated in **I)** with a population of a 1,000 individuals and a prevalence rate of 50%. All parameters were set to default values.
```bash
for i in seq 10; do
  time Rscript ../Estimate_prevalence_terminal.R Sim_data_Pop_size_1000_Preval_50.txt Results_Simu_PopSize_50_Preval_50_succ_10_meandiff_2_minmaxdiff_1_replicate_"$i" 50 2 10 1
done
```

## B) Gather the results for Table S2 (part C)
The command used to summarise the results in Supplementary Table S2 (part C)
```{bash}
cat *replicate_* |grep -v "Prevalence"|awk '{print $3"\t"$4}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g'
```

#  V) SAMPLE relationship Between Prevalence and Sampling Stability
## Variance in prevalence estimation using binomial sampling (Figure S5)
This simulation examines how the variance in observed prevalence changes with different true prevalence levels. We performed 1,000 simulations for each prevalence level, drawing 50 samples per simulation from a binomial distribution. The results show (check Figure S5) that when true prevalence is very low or very high, the variance in observed prevalence is small, whereas for intermediate values, the variance is larger.
```R
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Set parameters
sample_size <- 50
prevalence_levels <- c(0.01, 0.1, 0.5, 0.9, 0.99) # Testing prevalences of 1%, 10%, 50%, 90%, 99%.
num_replicates <- 1000

# Function to simulate prevalence estimation using binomial sampling
simulate_prevalence <- function(true_prevalence) {
  observed_prevalences <- numeric(num_replicates)
  
  for (i in 1:num_replicates) {
    sample <- rbinom(sample_size, 1, true_prevalence)
    observed_prevalences[i] <- mean(sample)
  }
  
  return(data.frame(TruePrevalence = true_prevalence,
                    ObservedPrevalence = observed_prevalences))
}

# Run simulations for all prevalence levels
simulated_data <- do.call(rbind, lapply(prevalence_levels, simulate_prevalence))

# Convert TruePrevalence to a factor for plotting
simulated_data$TruePrevalence <- as.factor(simulated_data$TruePrevalence)

# Plot the results using a boxplot
ggplot(simulated_data, aes(x = TruePrevalence, y = ObservedPrevalence)) +
  geom_boxplot(color = "black", fill="white",
               staplewidth = 0.15) +
  labs(x = "True Prevalence", y = "Observed Prevalence") +
  theme_bw()
```

##  Implication for SAMPLE (Figure S2)
SAMPLE evaluates prevalence stability by analyzing how prevalence estimates change as sample size increases. Because the variance of prevalence estimates depends on the true prevalence, this directly influences how quickly SAMPLE detects stability.
To illustrate this effect, we ran simulations testing known prevalence values ranging from 1% to 99% using 50 individuals per sample. Each prevalence level was simulated 1,000 times, and SAMPLE was applied using its default parameters to determine the point at which stability was reached. The script below shows how to perform a single replicate:

```R
library(SAMPLE)
args = commandArgs(trailingOnly=TRUE)
simulations_number <- args[1]

sample_size = 50 # Sample of 50 individuals
Prevalence = c(seq(1,99,1)) # test each prevalence from 1% to 99% incremented by 1%.

# Simulate sampling for low and high prevalence
simulate_data <- function(prevalence, sample_size) {
  Present <- 0
  a = 1
  
  # For each sample, update the prevalence estimate and CI
  for (i in 1:sample_size) {
    # Simulate the presence/absence of the species (1 = present, 0 = absent)
    new_sample <- rbinom(1, 1, prevalence/100)
    Present <- Present + new_sample
    estimated_prevalence <- Present /i 
    Absent <- sample_size - Present
    simdata = as.data.frame(cbind("Host1",c(rep(1,Present),rep(0,Absent))))
    simdata <- simdata[sample(1:nrow(simdata),replace=FALSE), ]

  }
  # Run SAMPLE on the simulated data 
  perm = RunPerm(input = simdata,replicates = 50)
  stable = stability(data = perm, stability_thresh = 2, success_points = 10, diff = 1)
  print(stable$Prevalence)
  print(stable$thres)
  
  return(stable$thres)
  
}

# Store the results of the replicate in a file
results <- mapply(simulate_data, MoreArgs = list(sample_size = sample_size), Prevalence)
results <- as.data.frame(cbind(Prevalence,results))
write.table(x = results,file = paste("simulations_",simulations_number,".txt",sep=""),sep="\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
```

The simulations were executed computing cluster, using a loop to submit 1,000 independent runs and speed up computations. Then the results were summarised into a single plot using the following script:
```R
# Load necessary libraries
library(dplyr)
library(ggplot2)

# List all files (assuming they are in your working directory)
file_list <- list.files(pattern = "simulations_\\d+\\.txt")

# Read all files and combine them into a single data frame
combined_data <- do.call(rbind, lapply(file_list, function(file) {
  read.table(file, header = TRUE)
}))

# Ensure the column names are correct
colnames(combined_data) <- c("Prevalence", "results")

# Convert Prevalence to a factor with levels in increasing order
combined_data$Prevalence <- factor(combined_data$Prevalence, levels = sort(unique(combined_data$Prevalence)))


# Plot the mean values with error bars
pdf("Figure_S2.pdf",15,9)
ggplot(combined_data, aes(x = Prevalence, y = results)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  scale_x_discrete(limits = sort(unique(combined_data$Prevalence)),breaks = seq(5, 95, by = 5)) +
  labs(x = "Prevalence (%)", y = "Detected Stability") +
  theme_minimal()  
dev.off()
```

The results confirm that SAMPLE consistently detects stable prevalence estimates more quickly when the true prevalence is either very low or very high (Figure S2). 
