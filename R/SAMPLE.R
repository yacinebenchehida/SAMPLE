#' Executing this function runs all the steps of the SAMPLE pipeline.
#'
#' Only the input dataframe is required to run SAMPLE, all the other arguments have default values that can be adjusted (but it is not necessary to do so). To obtain the same output as in the example file, remember to set.seed(812) after importing the dataset.
#' Alternatively, SAMPLE can be run using separate functions. To do so, please refer to RunPerm() for the first step of this process.

#' @importFrom magrittr %>%
#' @importFrom dplyr if_else
#' @import ggplot2
#' @param input Input dataframe (a dataframe object).
#' @param output_N Prefix used for the output (a character; default output_N="Results").
#' @param replicates Number of permutation replicates to perform (an integer; default replicates=50).
#' @param success_points Number of successive (mean) prevalence rates that are below a threshold (see parameter stability_thresh) used to define stability (an integer; default stability_thresh = 10).
#' @param stability_thresh Threshold used to define stability (an integer). This value will be divided by the square root of the number of replicates (a float; default stability_thresh = 2.0).
#' @param diff Difference between absolute minimum and maximum values among the all the means used to set the stability threshold (a float; default diff = 1.0).
#' @return A text file (.txt) with the output values of the analysis, and a PDF file with the generated plots from the analysis.

#'
#' @examples
#' data("coral_symbionts")
#' set.seed(812)
#' SAMPLE(input = coral_symbionts,output_N = "Example",replicates = 50,stability_thresh = 2,success_points = 10,diff = 1)

#' @export

SAMPLE <- function(input,output_N="Results",replicates=50,stability_thresh=2,success_points=10,diff=1){
  ###############
  # Upload data #
  ###############
  if(class(input) == "data.frame"){
    print("input is a data frame")
    df_name = deparse(substitute(input))
    input[is.na(input)] <- 0
    data <- input
    df.name <- deparse(substitute(data))
    print(df_name)
    colnames(data) = c("Host", colnames(data)[2:dim(data)[2]])
    # Replace data larger than 1 to 1
    if (max(unlist(data[,-c(1)])) > 1){
      print("DATA NOT ENCODED AS 0 AND 1. I AM GOING TO CONVERT VALUES LARGER THAN 1 to 1.")
      replace_larger_1 <- function(x){
        if_else(x > 1,1,x)
      }
      data <- data %>% dplyr::mutate_if(is.numeric, replace_larger_1)
    }
  } else{
    print("input is a not data frame")
    data = read.table(input, header = T, sep=c("\t",","), fill=TRUE) # Read input data file
    df_name = input
    colnames(data) = c("Host", colnames(data)[2:dim(data)[2]])

    # Replace missing data by 0
    data[is.na(data)] <- 0

    # Replace data larger than 1 to 1
    if (max(unlist(data[,-c(1)])) > 1){
      print("DATA NOT ENCODED AS 0 AND 1. I AM GOING TO CONVERT VALUES LARGER THAN 1 to 1.")
      replace_larger_1 <- function(x){
        if_else(x > 1,1,x)
      }
      data <- data %>% dplyr::mutate_if(is.numeric, replace_larger_1)
    }
  }


  #######################################
  # Set output name defined by the user #
  #######################################
  output_name = as.character(output_N)

  #######################################
  # Set the number of replicates to use #
  #######################################
  repli = replicates # Set the number of replicates to perform

  #############################################################
  # Set threshold from which to consider the system is stable #
  #############################################################
  stability = as.numeric(stability_thresh)/sqrt(as.numeric(repli))

  ###################################
  # Set number of successive points #
  ###################################
  successive_points = success_points

  #######################################
  # Set difference between mean and max #
  #######################################
  minmax = diff

  ########################
  # Print settings used #
  ########################
  cat("\n")
  cat("Running with the following settings:\n")

  cat(paste("Input name: ",df_name,sep=""),"\n")
  cat(paste("Output prefix: ",output_name,sep=""),"\n")
  cat(paste("Number of replicates: ",repli,sep=""),"\n")
  cat(paste("Mean difference tolerated/√(Number of replicates): ",stability,sep=""),"\n")
  cat(paste("Number of successive points: ",successive_points,sep=""),"\n")
  cat(paste("Difference between min and max tolerated: ",minmax,sep=""),"\n")
  cat("\n")

  ############################################################
  # Set a few variables used later by the randomization loop #
  ############################################################
  list <- list() # Create an empty list
  a = 1 # Initialize a counter what will be used to append every new element to the list
  Host = unique(data$Host) # Define unique species
  Host = Host[!is.na(Host)] # Remove the lines in the input without a defined host species
  species_counter = 0 # Initialize a counter what will count the number of species

  ########################
  # Randomization loops  #
  ########################
  print("Starting permutations")

  for (sp in Host){ # For loops over each species sp
    Host_sp = data[data$Host==sp,] # Subset the element from the raw data specific to the species sp
    colonies = dim(Host_sp)[1] # Count the total number of colonies for species sp
    species_counter = species_counter + 1
    for (taxa in 2:(dim(Host_sp)[2])){
      taxa_name = names(Host_sp[,taxa, drop = FALSE])
      taxa_size = length(Host_sp[Host_sp[taxa]==1,taxa]) # Count the total number of colonies with crabs for species sp
      print(paste("Currently running species ",species_counter,": ", sp," ","taxa: ", taxa_name, sep=""))
      if(taxa_size==0){
        print(paste(sp," has no ", taxa_name," I can't run the randomization", sep=""))
        next
      }
      if((dim(Host_sp)[1] > 5)){ # If there are at least 5 symbionts
        for (replicates in 1:repli){ # For loops to replicate
          substract_taxa = Host_sp # Make sure to reinitialize the raw input for each replicates
          Samples = 0
          for (i in 1:dim(Host_sp)[1]-1){ # Loop over each individuals in the colonies of species sp
            substract_taxa = dplyr::sample_n(substract_taxa,colonies-i,replace=FALSE) # Randomly remove 1 individual without replacement.
            Num = (length(substract_taxa[substract_taxa[taxa]==1,taxa])/(colonies-i))*100 # Looks at the number of colonies with crabs / the total number of remaining colonies (in %) for each species sp.
            list[[a]] = c(Num,taxa_name,i,replicates,sp,colonies-Samples) # Append the number of crabs counted, number of individuals removed, the replicate and the  species to the list.
            a = a + 1 # Increment by one the list so at the next iteration of the loop the information is stored in the element a + 1 of the list.
            Samples = Samples + 1
          }
        }
      }
    }
    print(paste("=====> ", sp," finished", sep=""))
  }

  ###########################################################################
  # Combine permutations results and make sure they are in the right format #
  ###########################################################################
  data = do.call(rbind, list) #  Combine the results
  data = as.data.frame(data) #  Results into a data frame format
  colnames(data) = c("Prevalence","Taxa","Substract","replicates","Host_sp","colonies") #  define column names
  data$Prevalence = as.numeric(data$Prevalence) #  Make sure the number of crabs is a numeric value
  data$Substract = as.numeric(data$Substract) #  Make sure the number of individuals removed is a numeric value

  #####################################
  # Determine threshold of stability  #
  #####################################
  #list = list() # Create an empty list for threshold based on confidence interval
  list_mean_diff = list() # Create an empty list threshold based on a mean diff between points
  list_stable = list()
  counter = 1 # Initialize a counter what will be used to append every new element to the list above
  counter_sp = 1 # Initialize a counter what will be used to append every new element to the list above
  options(warn = -1) # Prevent R to print messages and warnings not defined by this script
  meanci_previous = 0
  counter_mean_diff_low = 0
  Stability_assessement = FALSE
  suppressMessages(for (sp in unique(data$Host_sp)){ # For loops over each host species sp
    for (taxa in unique(data$Taxa)){ # Loop over each taxa (parasites or symbionts)
      Host_sp = data[data$Host_sp==sp & data$Taxa==taxa & data$colonies > 0,] # Create a data frame with the host species analysed in the loop
      Host_sp$colonies = as.numeric(Host_sp$colonies)
      Host_sp$replicates = as.numeric(Host_sp$replicates)
      tot_samples = length(unique(Host_sp$colonies)) # determine the number of host species individuals
      for (size in 1:(tot_samples)){ # Increase  sample size from 2 individuals to the total number of individuals
        occ_rate = Host_sp[Host_sp$colonies==size,1] # Extract data for the current sample size
        meanci = mean(occ_rate)  # Get the mean value of the prevalence at the current sample size
        mean_diff = abs(meanci - meanci_previous)
        meanci_previous = meanci #  Store previous mean in a different object
        if((mean_diff != "NaN") && (mean_diff < stability)){ #  If the difference of means is smaller than the stability threshold
          counter_mean_diff_low = counter_mean_diff_low + 1 # Increment the counter of mean differences below stability by 1
          list_stable[[counter_mean_diff_low]] = c(size,meanci)
        }
        else{
          counter_mean_diff_low = 0
        }
        if(counter_mean_diff_low==successive_points){
          stable_points = as.data.frame(do.call(rbind, list_stable))
          stable_point = stable_points[1,1]
          differential = max(stable_points$V2) - min(stable_points$V2)
          if(differential < minmax){
            print(paste("reached stability at ",stable_point," for ", sp," ", taxa, sep=""))
            list_mean_diff[[counter_sp]] = c(sp,stable_point,meanci,taxa)
            counter_sp = counter_sp + 1
            Stability_assessement = TRUE
            break
          }
          else{
            counter_mean_diff_low = counter_mean_diff_low - 1
            list_stable = list_stable[-1]
          }
        }
        if((size==tot_samples-1) && (counter_mean_diff_low < successive_points)){
          print(paste("Did not reach stability for ",sp," ", taxa,sep=""))
          list_mean_diff[[counter_sp]] = cbind(sp,"NA",NA,taxa) # Add NA to the list
          counter_sp = counter_sp + 1
        }
        if((size==tot_samples-1) && (Stability_assessement == FALSE) && (counter_mean_diff_low ==successive_points)){
          print(paste("Did not reach stability for ",sp," ", taxa," ", differential, " ", 1.5 * stability,sep=""))
          list_mean_diff[[counter_sp]] = cbind(sp,"NA",NA,taxa) # Add NA to the list
          counter_sp = counter_sp + 1
        }
      }
    }
  })

  info = as.data.frame(do.call(rbind,list_mean_diff)) # Transform the list with information on stability on each host species and parasites/symbionts into a single data frame.
  colnames(info) = c("Host_sp","thres","Prevalence","type_species")
  info$thres = as.numeric(info$thres)
  info$Prevalence = as.numeric(info$Prevalence)
  info$missing = "2beadded" # Create a new column stating with the stability of the system was reached or not after sampling
  info[is.na(info$thres),5] = "na" # If the system not stable in the end put "na"
  info[info$missing=="2beadded",5] = NA # If it reached stability put NA. So it is not plotted.

  #####################
  # Define color skim #
  #####################
  n <- length(unique(data$Taxa)) # Count the number of total parasites/symbionts
  list = list()
  counter = 1

  for (i in unique(info$Host_sp)){
    species_number = info[info$Host_sp==i,]
    length(species_number$type_species)
    list[[counter]] = length(species_number$type_species)
    counter = counter + 1
  }
  maximum_species = max(do.call(c,list))

  cat("\n")

  cat(paste("Number of hosts: ",length(unique(info$Host_sp)),sep=""),"\n")
  cat(paste("Number of symbionts: ", n,sep=""),"\n")
  cat(paste("Maximum number of symbionts per host: ",maximum_species,sep=""),"\n")
  cat("\n")


  #if(maximum_species < 3){ # If it smaller than 6 is you color code below
  #  couleur = viridis(begin = 0, end = .75, n)
  #}else if((maximum_species > 2) && (maximum_species < 7)){
  #couleur = c(viridis(begin = 0, end = .75, n-1),"#ff8c00")
  # print(couleur)
  #}else{ # If it large use a random color code assigned by the brewer palette
  #  print("Warning: Too many host species the plot is going to be very cluttered and hard to read")
  #  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #couleur = col_vector[1:n]
  #}

  if(n < 8){
    couleur = c("#7B3014","#D04A07","#F98C40","black","#5AA5CD","#236CA7","#26456E")
  }else{ # If it large use a random color code assigned by the brewer palette
    print("Warning: Too many host species the plot is going to be very cluttered and hard to read")
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    couleur = col_vector[1:n]
  }
  cat("Color scheme defined","\n")
  cat("Start plotting","\n")

  ###################
  # Start plotting  #
  ###################
  # The plotting works by adding one taxa (parasites/symbionts) at a time as a new layer.
  count = 1 # Counter that will by incremented for each taxa. This counter is used to choice the colors and figure out the spatial coordinate of the geom_text for the stability threshold.
  data$Taxa = gsub(pattern = "[_-]", replacement = " ", x=data$Taxa, perl = TRUE)
  data$Host_sp = gsub(pattern = "[_-]", replacement = " ", x=data$Host_sp, perl = TRUE)
  info$type_species = gsub(pattern = "[_-]", replacement = " ", x=info$type_species, perl = TRUE)
  info$Host_sp = gsub(pattern = "[_-]", replacement = " ", x=info$Host_sp, perl = TRUE)

  suppressMessages(for (taxonomy_names in sort(unique(data$Taxa))){ # Counter that will count the number of
    if(count==1){ # If count = 1 make the initial plotting with the ggplot function
      plot_data = data[data$Taxa==taxonomy_names,] #  Loop over taxa
      plot_data = plot_data %>%
        dplyr::group_by(Host_sp,Substract,colonies) %>%
        dplyr::summarise(avg = mean(Prevalence), lci = Rmisc::CI(Prevalence,ci = 0.95)[3], uci = Rmisc::CI(Prevalence,ci = 0.95)[1])
      plot_data$colonies = as.numeric(plot_data$colonies)
      plot_data$Taxa = as.factor(taxonomy_names)

      p = ggplot2::ggplot(plot_data, ggplot2::aes(x=colonies, y=avg,color=Taxa)) + ggplot2::geom_point(size=0.5) + ggplot2::scale_color_manual(values = couleur)
      p = p +  ggplot2::theme_bw() + ggplot2::facet_wrap(Host_sp~., scales = "free",  ncol = 3) +  ggplot2::xlab("Number of samples") + ggplot2::ylab("Prevalence (%)")
      p = p + ggplot2::theme(strip.background =element_rect(fill="wheat1"),strip.text.x = ggplot2::element_text(size = 14,face = "bold.italic"),legend.text=ggplot2::element_text(size=15),legend.title=ggplot2::element_text(size=16),axis.title=ggplot2::element_text(size=15),axis.text = ggplot2::element_text(size = 15))
      p = p + ggplot2::geom_ribbon(data=plot_data,aes(ymin=lci,ymax=uci),fill="grey", color="grey",alpha =0.5)
      p = p + ggplot2::geom_vline(data  = info[info$type_species==taxonomy_names,], aes(xintercept = thres),color = couleur[count] , linetype="dotted")
      p = p + ggplot2::geom_text(y = Inf, aes(x = Inf, label = ifelse(thres < 100, gsub(pattern = "(.+)", " \\1 ", thres), thres)),data = info[info$type_species==taxonomy_names,], color = couleur[count],hjust =1.2, vjust = 2,size = 4.5)
      p = p + ggplot2::geom_text(y = Inf, aes(x = Inf, label = gsub(pattern = "(.+)", " \\1 ", missing)), data = info[info$type_species==taxonomy_names,], color = couleur[count],hjust =1.5, vjust = 2,size = 4.5)
      p = p + ggplot2::geom_text(aes(x = Inf, y = Prevalence, label = round(Prevalence, digits = 1)), data =  info[info$type_species==taxonomy_names,], color = couleur[count],hjust = 1.5, vjust = 1.5,size = 4.5)
      p = p + ggplot2::guides(colour = guide_legend(override.aes = list(size=2)))
      count = count + 1
    }
    else{ #  if count > 1 add the successive layers (taxa) in each loop
      plot_data = data[data$Taxa==taxonomy_names,]
      plot_data = subset(data, data$Taxa==taxonomy_names,drop = FALSE)
      plot_data = plot_data %>%
        dplyr::group_by(Host_sp,Substract,colonies) %>%
        dplyr::summarise(avg = mean(Prevalence), lci = Rmisc::CI(Prevalence,ci = 0.95)[3], uci = Rmisc::CI(Prevalence,ci = 0.95)[1])
      plot_data$colonies = as.numeric(plot_data$colonies)
      plot_data$Taxa = as.factor(taxonomy_names)

      p = p + ggplot2::geom_point(data = plot_data, ggplot2::aes(x=colonies, y=avg, color=Taxa) ,size=0.5) +  ggplot2::scale_color_manual(values = couleur)
      p = p + ggplot2::theme_bw() + ggplot2::facet_wrap(Host_sp~., scales = "free",  ncol = 3)
      p = p + ggplot2::theme(strip.background =element_rect(fill="wheat1")) + ggplot2::theme(legend.position="bottom")
      p = p + ggplot2::theme(strip.text = element_text(face = "bold.italic",size = 14),legend.text=element_text(size=15,face="italic"),legend.title=element_text(size=16),axis.title=element_text(size=15),axis.text = element_text(size = 15))
      p = p + ggplot2::geom_ribbon(data=plot_data,aes(ymin=lci,ymax=uci),fill="grey", color="grey",alpha =0.5)
      p = p + ggplot2::geom_vline(data  = info[info$type_species==taxonomy_names,], aes(xintercept = thres),color = couleur[count] , linetype="dotted")
      p = p + ggplot2::geom_text(y = Inf, aes(x = Inf, label = ifelse((thres < 100) | (is.na(thres)), gsub(pattern = "(.+)", " \\1 ", thres), thres)), data = info[info$type_species==taxonomy_names,], color = couleur[count],hjust =  count * 1.2, vjust = 2,size = 4.5)
      p = p + ggplot2::geom_text(y = Inf, aes(x = Inf, label = gsub(pattern = "(.+)", " \\1 ", missing)), data = info[info$type_species==taxonomy_names,], color = couleur[count],hjust = count * 1.5, vjust = 2,size = 4.5)
      p = p + ggplot2::geom_text(aes(x = Inf, y = Prevalence, label = round(Prevalence, digits = 1)), data =  info[info$type_species==taxonomy_names,], color = couleur[count],hjust = 1, vjust = 1.5,size = 4.5)
      p = p + ggplot2::guides(colour = guide_legend(override.aes = list(size=2)))
      count = count + 1
    }
  })

  ###################
  # PLotting per se #
  ###################
  nb_sp = length(unique(data$Host_sp)) # Check the number of species. Import to know how large the pdf plot has to be.

  if(nb_sp > 1 && nb_sp < 4){ # If there are less than 4 species
    pdf(paste(output_name, ".pdf", sep=""),length(unique(data$Host_sp))*5,6) # this define the width and the length of the pdf if there are less than 4 species to plot
    plot(p)
    dev.off()
  } else if (nb_sp > 4){ # If there are more than 4 species
    pdf(paste(output_name, ".pdf", sep=""),12,round(length(unique(data$Host_sp)) / 4) * 4) # this define the width and the length of the pdf if there are more than 4 species to plot
    plot(p)
    dev.off()
  } else {  # If there are exactly 4 species
    pdf(paste(output_name, ".pdf", sep=""),8,5) # this define the width and the length of the pdf if there are four species to plot
    plot(p)
    dev.off()
  }

  ###############################
  # Save results as a text file #
  ###############################
  textfile = info[,c(1,4,3,2)] # Create an object textfile coutaning the same information as info but reshape in a more meaningful order.
  colnames(textfile) = c("Host_species","Taxa","Prevalence","thres_stability")
  write.table(file = paste(output_name, ".txt", sep=""), x = textfile,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE) # Save the results as a text file.

}
