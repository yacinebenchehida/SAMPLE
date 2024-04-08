#' Run permutations
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr if_else
#' @import ggplot2
#' @param input The input dataframe (a dataframe object)
#' @param replicates Nnumber of permutation replicates to perform (an integer; default replicates=50)
#' @return A dataframe.
#'
#' @examples
#' data("coral_symbionts")
#' RunPerm(input = coral_symbionts,replicates = 50)
#'
#' @export

RunPerm <- function(input,replicates=50){
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
  # Set the number of replicates to use #
  #######################################
  repli = replicates # Set the number of replicates to perform

  ########################
  # Print settings used #
  ########################
  cat("\n")
  cat("Running with the following settings:\n")

  cat(paste("Input name: ",df_name,sep=""),"\n")
  cat(paste("Number of replicates: ",repli,sep=""),"\n")
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

  return(data)
}

