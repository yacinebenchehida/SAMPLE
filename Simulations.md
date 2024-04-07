# I) Generate the simulated data
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
## A) Run prevalence script on each data set generated in I)
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
These commands run pipeline on each dataset generated in I while varying the number of replicates from 10 to 500. 

## B) Extract the results 
```bash

for i in Results*txt; do
  Size=$(echo $i|cut -d _ -f 4)
  Preval=$(echo $i|cut -d _ -f 6|perl -pe 's/\.txt//g')
  Obs_Preval=$(cat $i|awk 'NR > 1 {print $3}'|perl -pe 's/([0-9]+)(\.)(\d\d).+/$1$2$3/g')
  Stability=$(cat $i|awk 'NR > 1 {print $4}')
  echo -e Population size: $Size ind \| Simulated Prevalence: $Preval% \| Inferred Prevalence: $Obs_Preval% \| Stable from: $Stability Samples
done
```

## C) Gather the results (to paste in excel file):
```bash
for pop in 100 1000 10000; do for j in 10 20 30 40 50 60 70 80 90 100 200 500; do cd "$j"_replicates; for i in *_"$pop"_*.txt; do Samp=$(echo $i|cut -d _ -f 4|perl -pe 's/\.txt//g'); Prev=$(echo $i|cut -d _ -f 6|perl -pe 's/\.txt//g'); rate=$(cat $i|awk 'NR>1 {print $4}');prevalence=$(cat $i|awk 'NR>1 {print $3}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g');  echo -e $Samp"\t"$Prev"\t"$rate"|"$prevalence"%";  done|awk '{print $3}'|perl -pe 's/\n/\t/g';cd ..;  echo -e "\n"; done|perl -pe 's/^\n$//g'; done

```

#  III) Assess impact of successive points, mean-difference, max-mean-difference
## A) Run prevalence script 
```{bash}
for succ in 2 10 50; do for meandiff in 1 2 5 10; do for minmax in 0.5 1 2; do Rscript ../prevalence_script_test.R Sim_data_Pop_size_1000_Preval_50.txt Results_Simu_PopSize_100_Preval_50_succ_"$succ"_meandiff_"$meandiff"_minmaxdiff_"$minmax" 50 $meandiff  $succ $minmax; done; done; done

```

## B) Gather the results
```{bash}
for i in $(ls *succ_2_*.txt|sort -V); do rate=$(cat $i|awk 'NR>1 {print $4}');prevalence=$(cat $i|awk 'NR>1 {print $3}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g');echo -e $rate"|"$prevalence"%" ; done|awk -v n=3 '{a[NR]=$0}END{ x=1; while (x<=n){ for(i=x;i<=length(a);i+=n) printf "%s",a[i]"\t"; print ""; x++; } }' |column -t

```

#  IV) Assess consistence of the results for 10 replicates using the same sets of parameters
## A) Run prevalence script
```{bash}
for i in seq 10; do time Rscript ../prevalence_script_test.R Sim_data_Pop_size_1000_Preval_50.txt Results_Simu_PopSize_50_Preval_50_succ_10_meandiff_2_minmaxdiff_1_replicate_"$i" 50 2 10 1; done

```

## B) Gather the results
```{bash}
cat *replicate_* |grep -v "Prevalence"|awk '{print $3"\t"$4}'|perl -pe 's/(\d+)(\.)(\d)\d+/$1$2$3/g'
```
