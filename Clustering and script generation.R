# Load necessary libraries
library(dplyr)
library(cluster)
library(proxy)  # For Jaccard similarity calculation
require(ape)
## 1) DATA imports
setwd("C:/Users/jzs19xhz/OneDrive - The Royal Botanic Gardens, Kew/CFP_Non-native_Trees") #Set path
scripts_path<-"./Code"
predictions_path<-"~/Outputs"
host_status_df<-read.csv("./Data/CABI_HostMatrix_Cleaned.csv") # InputHSM
megatree <-read.tree("./Data/accepted_minitree_EURO_FI.tre") #Phylo, phylogenetic tree

### 2) Run lines 15-709

hosts <- host_status_df[, 1]               
status_df <- host_status_df[, -1]     
non_zero_rows <- rowSums(status_df) > 0
non_zero_cols <- colSums(status_df) > 0
filtered_status_df <- status_df[non_zero_rows, non_zero_cols]
filtered_host_status_df <- cbind(hosts[non_zero_rows], filtered_status_df)
colnames(filtered_host_status_df)[1] <- colnames(host_status_df)[1]

A<-ape::vcv.phylo(megatree)
# Step 1: Prepare the data by transposing to have cluster as rows
threat_matrix<- as.matrix(filtered_status_df)
threat_matrix<- t(threat_matrix)

# Step 2: Calculate pairwise Jaccard similarity between cluster
similarity_matrix <- as.matrix(proxy::simil(threat_matrix, method = "Jaccard"))

# Step 3: Convert similarity to dissimilarity for clustering
dissimilarity_matrix <- 1 - similarity_matrix

# Step 4: Perform hierarchical clustering
hc <- hclust(as.dist(dissimilarity_matrix), method = "average")

### Function to get a custom cluster for an individual threat
generate_cluster <- function(threat, cutoff = 0.1, top_n = 20) {
  # Ensure threat is a valid name in the similarity matrix
  if (!threat %in% colnames(similarity_matrix)) {
    stop("Threat not found in the similarity matrix.")
  }
  
  # Get similarity values for the specified threat
  threat_similarity <- similarity_matrix[threat, ]
  
  # Remove NA values (e.g., self-comparison or undefined similarity)
  threat_similarity <- threat_similarity[!is.na(threat_similarity)]
  
  # Filter cluster above the similarity cutoff
  threshold_cluster <- names(threat_similarity[threat_similarity >= cutoff])
  
  # Or select the top_n most similar cluster if fewer than top_n exceed the cutoff
  top_n_cluster <- names(sort(threat_similarity, decreasing = TRUE)[1:min(top_n, length(threat_similarity))])
  
  # Return the smaller cluster between threshold-based and top_n-based clusters
  if (length(threshold_cluster) <= length(top_n_cluster)) {
    return(threshold_cluster)
  } else {
    return(top_n_cluster)
  }
}

### Function to generate a combined cluster for a specified host
generate_combined_cluster_for_host <- function(host, similarity_matrix = similarity_matrix, host_status_df = filtered_host_status_df, cutoff = 0.1, top_n = 20, alpha = 0.05, size_limit = 30) {
  # Ensure host exists in the matrix
  if (!host %in% host_status_df[, 1]) {
    stop("Host not found in the host status matrix.")
  }
  
  # Identify cluster hosted by the specified host
  host_row <- host_status_df[host_status_df[, 1] == host, -1]
  hosted_cluster <- colnames(host_row)[host_row == 1]
  
  # Generate clusters for each hosted threat
  all_clusters <- list()
  for (threat in hosted_cluster) {
    all_clusters[[threat]] <- generate_cluster_without_redundancy(
      threat = threat,
      similarity_matrix = similarity_matrix,
      host_status_df = filtered_host_status_df,
      cutoff = cutoff,
      top_n = top_n,
      alpha = alpha
    )
  }
  
  # Combine clusters and remove duplicates
  combined_cluster <- unique(unlist(all_clusters))
  
  # If combined cluster exceeds size limit, retain most similar cluster only
  if (length(combined_cluster) > size_limit) {
    # Get similarity scores for each threat in combined cluster relative to the primary threat
    combined_similarity <- similarity_matrix[hosted_cluster[1], combined_cluster]
    
    # Sort by similarity, retain top N, and determine the similarity cut-off for excluded cluster
    sorted_cluster <- names(sort(combined_similarity, decreasing = TRUE))
    final_cluster <- sorted_cluster[1:size_limit]
    cutoff_similarity <- combined_similarity[sorted_cluster[size_limit + 1]]
    
    return(list(cluster = final_cluster, cutoff_similarity = cutoff_similarity))
  } else {
    return(list(cluster = combined_cluster, cutoff_similarity = NULL))
  }
}

### Recursive removal of redundant cluster clustering
generate_cluster_without_redundancy <- function(threat, similarity_matrix, host_status_df, cutoff = 0.1, top_n = 20, alpha = 0.05, max_size = Inf) {
  # Generate the initial cluster
  custom_cluster <- generate_cluster(threat, cutoff = cutoff, top_n = top_n)
  
  # Focal variable
  focal <- threat
  varnames <- custom_cluster
  final_cluster <- c()
  
  # Ensure varnames are valid column names in the host status matrix
  varnames <- varnames[varnames %in% colnames(host_status_df)]
  
  if (length(varnames) == 0) {
    print(paste0("No valid variables in the cluster to assess for ",threat))
    return(NULL)
  }
  
  # Initialize a list to track excluded variables
  excluded_vars <- list()
  
  while (length(varnames) > 0) {
    # Join variable names with '+' for the formula
    all_vars <- paste(c(final_cluster, varnames), collapse = "+")
    
    # Create the formula for the linear model
    ff <- as.formula(paste(focal, "~", all_vars))
    
    # Fit the linear model
    mod1 <- tryCatch(
      lm(ff, data = host_status_df),
      error = function(e) {
        warning("Model fitting failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(mod1)) break  # Stop if model fitting fails
    
    # Use drop1 to assess variable significance
    drop1_results <- drop1(mod1, test = "F")
    
    # Extract p-values (excluding intercept) and handle NA values
    p_values <- drop1_results$"Pr(>F)"[-1]
    
    # Check if any variables have NA p-values
    na_vars <- varnames[is.na(p_values)]
    
    if (length(na_vars) > 0) {
      # Remove variables with NA p-values
      varnames <- varnames[!is.na(p_values)]
      excluded_vars <- c(excluded_vars, lapply(na_vars, function(x) list(variable = x, reason = "NA_p_value")))
      #print(paste("Removed variables with NA p-values:", paste(na_vars, collapse = ", ")))
      next
    }
    
    # Check if all remaining variables are significant
    if (all(p_values < alpha)) {
      # Add remaining significant variables to the final cluster and exit the loop
      final_cluster <- unique(c(focal, varnames))
      break
    }
    
    # Find the variable with the highest p-value (least significant)
    max_p_value <- which.max(p_values)
    
    # Remove the least significant variable
    removed_var <- varnames[max_p_value]
    varnames <- varnames[-max_p_value]
    
    # Record the exclusion
    excluded_vars <- c(excluded_vars, list(list(variable = removed_var, reason = "high_p_value")))
    
    # Output the remaining variables to monitor progress
    #print(paste("Remaining variables:", paste(varnames, collapse = ", ")))
  }
  
  # Step to enforce max_size constraint
  while (length(final_cluster) > max_size) {
    # Fit a model with the current final cluster
    all_vars <- paste(final_cluster, collapse = "+")
    ff <- as.formula(paste(focal, "~", all_vars))
    
    mod2 <- tryCatch(
      lm(ff, data = host_status_df),
      error = function(e) {
        warning("Model fitting failed during size adjustment: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(mod2)) break  # Stop if model fitting fails
    
    # Use drop1 to assess variable significance
    drop1_results <- drop1(mod2, test = "F")
    p_values <- drop1_results$"Pr(>F)"[-1]
    
    # Identify the variable with the highest p-value (least significant)
    max_p_value <- which.max(p_values)
    removed_var <- final_cluster[max_p_value]
    final_cluster <- final_cluster[-max_p_value]
    
    # Record the exclusion
    excluded_vars <- c(excluded_vars, list(list(variable = removed_var, reason = "size_limit")))
    
    #print(paste("Removed variable to meet max_size:", removed_var))
  }
  
  # Validate cluster names before returning
  valid_cluster <- final_cluster[final_cluster %in% colnames(host_status_df)]
  invalid_cluster <- setdiff(final_cluster, colnames(host_status_df))
  if (length(invalid_cluster) > 0) {
    warning("The following columns in the cluster are invalid and will be removed: ", paste(invalid_cluster, collapse = ", "))
  }
  
  return(valid_cluster)
}


identify_focal_and_closest_host <- function(focal_threat, cluster = c(), full_host_status_df = filtered_host_status_df, phylo_matrix, used_hosts = c()) {
  # Extract the host names from the first column of the host_status_df
  host_names <- full_host_status_df[, 1]
  
  host_matrix_temp <- cbind(host_status_df[,1],
                            host_status_df[, colnames(host_status_df) %in% cluster == TRUE])
  names(host_matrix_temp)[1] <- "Host"
  host_matrix <- host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
  host_matrix <- as.data.frame(rbind(host_matrix,
                                     host_matrix_temp[(length(host_matrix_temp$Host)-27):
                                                        length(host_matrix_temp$Host),]))
  
  
  # Step 1: Identify the focal host
  # Filter the host names based on the focal threat
  # Step 1: Identify the focal host
  
  # Ensure the focal threat exists in the host_matrix
  if (!(focal_threat %in% colnames(host_matrix))) {
    warning(paste("Column", focal_threat, ", the focal threat, does not exist in host_matrix. Skipping this iteration."))
    return(NULL) 
  }
  
  # Filter the host names based on the focal threat
  focal_threat_hosts <- host_names[host_matrix[, focal_threat] == 1]
  
  if (length(focal_threat_hosts) == 0) {
    stop(paste("No hosts found for focal threat:", focal_threat))
  }
  
  if (length(focal_threat_hosts) == 1) {
    # If there is only one host, select it as the focal host
    focal_host <- focal_threat_hosts
  } else {
    # Validate the cluster contains valid columns
    missing_columns <- setdiff(cluster, colnames(host_matrix))
    if (length(missing_columns) > 0) {
      print("Cluster:")
      print(cluster)
      print("host_matrix Columns:")
      print(colnames(host_matrix))
      print("Missing Columns:")
      print(missing_columns)
      stop(paste("Error: The following columns in the cluster do not exist in host_matrix:", 
                 paste(missing_columns, collapse = ", ")))
    }
    
    # Identify the host with the greatest overlap with the cluster threats
    cluster_overlap <- rowSums(host_matrix[host_names %in% focal_threat_hosts, cluster, drop = FALSE])
    focal_host <- focal_threat_hosts[which.max(cluster_overlap)]
  }
  
  
  # Validate that the focal host is found in the phylogenetic matrix
  if (!focal_host %in% rownames(phylo_matrix)) {
    stop(paste("Focal host", focal_host, "not found in phylogenetic matrix."))
  }
  
  # Add the focal host to the used hosts list
  used_hosts <- unique(c(used_hosts, focal_host))
  
  # Step 2: Identify the phylogenetically closest unused host
  # Filter full host status matrix for hosts of the focal threat
  potential_hosts <- host_names[full_host_status_df[, focal_threat] == 1]
  
  # Exclude already used hosts
  unused_hosts <- setdiff(potential_hosts, used_hosts)
  
  if (length(unused_hosts) == 0) {
    # No unused hosts available for the focal threat
    warning("No unused hosts available for the focal threat. Selecting the next closest phylogenetic host.")
    
    # Identify all hosts in the phylogenetic matrix excluding used hosts
    all_unused_hosts <- setdiff(rownames(phylo_matrix), used_hosts)
    
    # Calculate phylogenetic distances between the focal host and all unused hosts
    phylo_distances_all <- phylo_matrix[focal_host, all_unused_hosts, drop = FALSE]
    
    # Identify the closest unused host
    closest_unused_host <- colnames(phylo_distances_all)[which.max(phylo_distances_all)]
    
    # Add the closest unused host to the used hosts list
    used_hosts <- unique(c(used_hosts, closest_unused_host))
    
    # Return the results
    return(list(
      focal_host = focal_host,
      closest_unused_host = closest_unused_host,
      used_hosts = used_hosts
    ))
  }
  
  # Ensure the unused hosts exist in the phylogenetic matrix
  unused_hosts <- intersect(unused_hosts, colnames(phylo_matrix))
  if (length(unused_hosts) == 0) {
    warning("No valid unused hosts in phylogenetic matrix. Selecting the next closest phylogenetic host.")
    
    # Identify all hosts in the phylogenetic matrix excluding used hosts
    all_unused_hosts <- setdiff(rownames(phylo_matrix), used_hosts)
    
    # Calculate phylogenetic distances between the focal host and all unused hosts
    phylo_distances_all <- phylo_matrix[focal_host, all_unused_hosts, drop = FALSE]
    
    # Identify the closest unused host
    closest_unused_host <- colnames(phylo_distances_all)[which.max(phylo_distances_all)]
    
    # Add the closest unused host to the used hosts list
    used_hosts <- unique(c(used_hosts, closest_unused_host))
    
    # Return the results
    return(list(
      focal_host = focal_host,
      closest_unused_host = closest_unused_host,
      used_hosts = used_hosts
    ))
  }
  
  # Calculate phylogenetic distances between the focal host and unused hosts
  phylo_distances <- phylo_matrix[focal_host, unused_hosts, drop = FALSE]
  
  # Identify the closest unused host
  closest_unused_host <- colnames(phylo_distances)[which.max(phylo_distances)]
  
  # Step 3: Return the results
  return(list(
    focal_host = focal_host,
    closest_unused_host = closest_unused_host,
    used_hosts = unique(c(used_hosts, closest_unused_host))
  ))
}



granddesign <- function(threat, similarity_matrix=similarity_matrix, 
                        host_status_df = filtered_host_status_df, 
                        cutoff = 0.1, top_n = 20, alpha = 0.05, size_limit = 30, 
                        phylo_matrix = A){
  FINAL_CLUSTER <- c()
  max_size <- Inf
  HSM_size_flag1 <- 0
  iteration_count <- 0  # Add iteration counter
  max_iterations <- 100 # Maximum allowed iterations
  
  while(HSM_size_flag1 == 0 && iteration_count < max_iterations){
    iteration_count <- iteration_count + 1
    
    ##Generate RTRCluster for focal threat
    RTRCluster1 <- generate_cluster_without_redundancy(
      threat = threat,
      similarity_matrix = similarity_matrix,
      host_status_df = filtered_host_status_df,
      max_size = max_size
    )
    
    # Check if cluster generation failed or is too small
    if(length(RTRCluster1) < 2) {
      print(paste0("RTRCluster1 too small to proceed for ",threat))
      return(RTRCluster1)
    }
    
    ##Generate host matrix for model input to judge size
    host_matrix_temp <- cbind(host_status_df[,1],
                              host_status_df[, colnames(host_status_df) %in% RTRCluster1 == TRUE])
    names(host_matrix_temp)[1] <- "Host"
    host_matrix <- host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
    host_matrix <- as.data.frame(rbind(host_matrix,
                                       host_matrix_temp[(length(host_matrix_temp$Host)-27):
                                                          length(host_matrix_temp$Host),]))
    
    matrix_size <- length(host_matrix$Host) * length(names(host_matrix))
    
    if(matrix_size > 2700){
      max_size <- length(RTRCluster1) - 1
      if(max_size < 2) {  # Check if max_size becomes too small
        print("Cannot reduce cluster size further while maintaining requirements")
        return(RTRCluster1)
      }
    } else if(2700 < matrix_size && matrix_size <= 3300){
      HSM_size_flag1 <- 1
      FINAL_CLUSTER <- RTRCluster1
      break
    } else {
      result <- RTRCluster1
      HSM_size_flag1 <- 1
      break
    }
  }
  
  if(iteration_count >= max_iterations) {
    print("Maximum iterations reached in first loop")
    print(paste0("FINAL CLUSTER IDENTIFIED for ",threat))
    FINAL_CLUSTER <- RTRCluster1
    return(FINAL_CLUSTER)
  }
  
  if(length(FINAL_CLUSTER) == 0){ 
    HSM_size_flag2 <- 0
    used_hosts <- c()
    iteration_count <- 0  # Reset iteration counter for second loop
    
    while(HSM_size_flag2 == 0 && iteration_count < max_iterations){
      iteration_count <- iteration_count + 1
      
      missing_columns <- setdiff(result, colnames(host_status_df))
      
      if (length(missing_columns) > 0) {
        print("Cluster:")
        print(cluster)
        print("host_status_df Columns:")
        print(colnames(host_status_df))
        print("Missing Columns:")
        print(missing_columns)
        stop(paste("Error: The following columns are missing from host_status_df:", 
                   paste(missing_columns, collapse = ", ")))
      }
      ##Generate new host matrix
      host_matrix_temp <- cbind(host_status_df[,1],
                                host_status_df[, colnames(host_status_df) %in% result == TRUE])
      names(host_matrix_temp)[1] <- "Host"
      host_matrix <- host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
      host_matrix <- as.data.frame(rbind(host_matrix,
                                         host_matrix_temp[(length(host_matrix_temp$Host)-27):
                                                            length(host_matrix_temp$Host),]))
      
      fhost_uhost <- identify_focal_and_closest_host(
        focal_threat = threat,
        cluster = result,
        full_host_status_df = filtered_host_status_df,
        phylo_matrix = A,
        used_hosts = used_hosts
      )
      
      if(is.null(fhost_uhost)) next 
      
      if(is.null(fhost_uhost$closest_unused_host)) {
        print("No more unused hosts available")
        return(result)
      }
      
      
      
      input_host <- fhost_uhost$closest_unused_host
      
      # Identify cluster hosted by the specified host
      host_row <- host_status_df[host_status_df[, 1] == input_host, -1]
      hosted_threat_list <- colnames(host_row)[host_row == 1]
      
      if(length(hosted_threat_list) == 0) {
        print("No cluster found for selected host")
        next
      }
      
      # Generate clusters for each hosted threat
      all_clusters <- list()
      for(hosted_threat in hosted_threat_list) {
        all_clusters[[hosted_threat]] <- generate_cluster_without_redundancy(
          threat = hosted_threat,
          similarity_matrix = similarity_matrix,
          host_status_df = filtered_host_status_df,
          cutoff = cutoff,
          top_n = top_n,
          alpha = alpha,
          max_size = max_size
        )
      }
      
      # Combine clusters and remove duplicates
      combined_cluster <- unique(unlist(all_clusters))
      
      if(length(combined_cluster) < 2) {
        print(paste0("Combined cluster too small to proceed for ",threat))
        next
      }
      
      host_matrix_temp <- cbind(host_status_df[,1],
                                host_status_df[, colnames(host_status_df) %in% combined_cluster == TRUE])
      names(host_matrix_temp)[1] <- "Host"
      host_matrix <- host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
      host_matrix <- as.data.frame(rbind(host_matrix,
                                         host_matrix_temp[(length(host_matrix_temp$Host)-27):
                                                            length(host_matrix_temp$Host),]))
      
      matrix_size <- length(host_matrix$Host) * length(names(host_matrix))
      
      if(matrix_size <= 2700){
        result <- combined_cluster
        used_hosts <- fhost_uhost$used_hosts
      } else if(2700 < matrix_size && matrix_size <= 3300){
        HSM_size_flag2 <- 1
        FINAL_CLUSTER <- combined_cluster
        break
      } else {
        cluster_lengths <- sapply(all_clusters, length)
        max_size <- max(cluster_lengths) - 1
        if(max_size == 1) {
          print("Cannot reduce cluster size further")
          FINAL_CLUSTER<-result
          host_matrix_temp <- cbind(host_status_df[,1],
                                    host_status_df[, colnames(host_status_df) %in% FINAL_CLUSTER == TRUE])
          names(host_matrix_temp)[1] <- "Host"
          host_matrix <- host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
          host_matrix <- as.data.frame(rbind(host_matrix,
                                             host_matrix_temp[(length(host_matrix_temp$Host)-27):
                                                                length(host_matrix_temp$Host),]))
          
          matrix_size <- length(host_matrix$Host) * length(names(host_matrix))
          
          print(paste("No. Hosts:", length(host_matrix$Host)))
          print(paste("Matrix size:", matrix_size))
          print(paste0("FINAL CLUSTER IDENTIFIED for ",threat))
          return(FINAL_CLUSTER)
        }
      }
    }
    
    if(iteration_count >= max_iterations) {
      print("Maximum iterations reached in second loop")
      print(paste0("FINAL CLUSTER IDENTIFIED for ",threat))
      return(FINAL_CLUSTER)
    }
  }
  print(paste0("FINAL CLUSTER IDENTIFIED for ",threat))
  return(FINAL_CLUSTER)
}

y_function <- function(x) {
  mean(plogis(rt(10^5, 2) - x))
}
# Target function to find the root for
target_function <- function(x, y_target) {
  y_function(x) - y_target
}
# Function to find a value
find_prior_mean <- function(y_target, interval = c(-20, 20)) {
  # Debugging
  print(paste("Testing lower bound:", interval[1]))
  print(paste("Target function at lower bound:", target_function(interval[1], y_target)))
  # 
  # Ensure the function is well-defined at the bounds
  lower_val <- target_function(interval[1], y_target)
  upper_val <- target_function(interval[2], y_target)
  
  if (is.na(lower_val) || is.na(upper_val)) {
    stop("The target function is not well-defined at the bounds of the interval.")
  }
  
  # Find the root
  uniroot(target_function, interval = interval, y_target = y_target)$root
}


generate_model_script <- function(cluster, threat) {

  if (length(cluster) <= 1) {
    # Write a CSV with the specified message
    final_script<-paste0("message_df <- data.frame(Message = \"No significantly related, non-redundant threats\") \n
                         write.csv(message_df, file = paste0(",predictions_path,"\"/predictions_", threat, ".csv\"), row.names = FALSE)")
    writeLines(final_script, paste0(scripts_path,"/Min_2_hosts_inputs/model_", threat, ".R"))
  }
  else{
    print(paste0("Threat cluster for ", threat, " length: ",length(cluster)))
  
  # Filter host matrix
  host_matrix_temp<-cbind(host_matrix_df[,1],host_matrix_df[, colnames(host_matrix_df) %in% cluster ==TRUE])
  names(host_matrix_temp)[1] <- "Host"
  host_matrix_gen<-host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]
  host_matrix_gen<-rbind(host_matrix_gen,host_matrix_temp[(length(host_matrix_temp$Host)-27):length(host_matrix_temp$Host),])
  
  #####Generate model input
  modelcommand<-paste0("require(dplyr)\n 
    require(ape)\n
    require(brms)\n
    require(tidyr)\n
    require(ggplot2)\n
    require(distributional)\n
    remotes::install_github(\"stan-dev/cmdstanr\")\n
    require(cmdstanr)\n
    install_cmdstan(cores=8,overwrite = TRUE)\n
    set_cmdstan_path(\"/data/home/mpx605/.cmdstan/cmdstan-2.33.1\")\n
    \n
    ## Import data\n
    # Phloygenetic tree\n
    megatree <-read.tree(\"./Data/accepted_minitree_EURO_FI.tre\")\n
    # Host status matrix (HSM) dataframe\n
    host_matrix_full<-read.csv(\"./Data/CABI_HostMatrix.csv\")\n
    ## Specify the cluster of interest\n
    cluster<-c(")
  
  for(i in 1:length(cluster)){
    if (i < length(cluster))
      modelcommand<-paste0(modelcommand,"\"",cluster[i],"\",\n")
    else
      modelcommand<-paste0(modelcommand,"\"",cluster[i],"\")\n")}
  
  modelcommand<-paste0(modelcommand,"## Cutting down the HSM to only include rows for plants with host status \n
    host_matrix_temp<-cbind(host_matrix_full[,1],host_matrix_full[, colnames(host_matrix_full) %in% cluster ==TRUE])\n
    names(host_matrix_temp)[1] <- \"Host\"\n
    host_matrix<-host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]\n
    host_matrix<-rbind(host_matrix,host_matrix_temp[(length(host_matrix_temp$Host)-27):length(host_matrix_temp$Host),])\n
    \n
    ## Generate phylogenetic covariance matrix\n
    A<-ape::vcv.phylo(megatree)\n
    A[1:length(megatree$tip.label),1:length(megatree$tip.label)]<-\n
      sqrt(A[1:length(megatree$tip.label),1:length(megatree$tip.label)])\n
    ## Set model form, cluster names without marks\n
    fullform <- bf(mvbind(\n")
  
  for(i in 1:length(cluster)){
    if (i < length(cluster))
      modelcommand<-paste0(modelcommand,cluster[i],",\n")
    else
      modelcommand<-paste0(modelcommand,cluster[i],")\n ~ 1 + (1|p|gr(Host, cov = A)),
    family = bernoulli())\n")
  }
  
  
  #####Generate priors input
  
  ## Function for finding an *approximate* value of student-t mean for a given random effect
  
  priorsvalues<-list()
  for(i in 1:length(cluster)){
    priorsvalues[i]<-round(find_prior_mean(sum(as.numeric(host_matrix_gen[,i+1]))/nrow(host_matrix_gen)), digits = 2)
  }
  cluster<-gsub("_", "", cluster)
  priorscommand<-paste0("priors <- c(")
  for(i in 1:length(cluster)){
    priorscommand<-paste0(priorscommand,"set_prior(\"student_t(2,",priorsvalues[i],",1)\", class = \"Intercept\", resp = \"",cluster[i],"\"),\n")
  }
  for(i in 1:length(cluster)){
    if (i < length(cluster))
      priorscommand<-paste0(priorscommand,"set_prior(\"normal(0,1)\", class = \"sd\", group = \"Host\", resp = \"",cluster[i],"\"),\n")
    else
      priorscommand<-paste0(priorscommand,"set_prior(\"normal(0,1)\", class = \"sd\", group = \"Host\", resp = \"",cluster[i],"\"))")
  }
  the_rest<-paste0("\n ## Train model, using full dataset with NAs included in HSM
       fullfit <- brm(fullform, data = host_matrix, data2 = list(A=A) , prior = priors,
                      chains = 1, init = 0,
                      iter=2000, warmup=500,
                      backend = \"cmdstanr\",
                      threads = threading(8)
       )
       
       ## Generate predictions for the new values
       logit_df <- as.data.frame(posterior_summary(posterior_linpred(fullfit, transform = FALSE)))
    predictions_NNT<- cbind(host_matrix_full$Host,logit_df[,seq(1,length(logit_df),4)])
    lowerquant_NNT<- cbind(host_matrix_full$Host,logit_df[,seq(3,length(logit_df),4)])
    
    colnames(predictions_NNT)[1]<-\"Host\"
    colnames(predictions_NNT)[-1]<- cluster")
final_script <- paste0(
  modelcommand, priorscommand,the_rest,
  "\n # Save predictions\n",
  "write.csv(predictions_NNT, ",predictions_path,"\"/predictions_", threat, ".csv\")\n"
)

writeLines(final_script, paste0(scripts_path,"/Min_2_hosts_inputs/model_", threat, ".R"))
  }
}


gd_gms_wrapper<-function(threat){
  
  cluster<-granddesign(threat=threat,
                       similarity_matrix = similarity_matrix,
                       host_status_df = filtered_host_status_df,
                       phylo_matrix = A)
  generate_model_script(cluster = cluster,threat=threat)
}

MASTER_FUNCTION <- function(threats_vector, min_hosts) {
  # Use the global host_status_df
  # Identify threats with more than the minimum number of hosts
  threats_with_min_hosts <- threats_vector[sapply(threats_vector, function(threat) {
    sum(host_status_df[, threat] == 1) >= min_hosts
  })]
  
  # Iterate over each threat and run the gd_gms_wrapper
  lapply(threats_with_min_hosts, function(threat) {
    gd_gms_wrapper(threat)
  })
}


## Run MASTER_FUNCTION specifying a vector for threat names to generate clusters
##   and scripts for model running. set minimum hosts as desired, 2 is recommended 
##   default
MASTER_FUNCTION(threats_vector=colnames(host_status_df[,2:200]),min_hosts=2)
