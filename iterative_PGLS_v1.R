## R Script for performing a PGLS iteratively across a list of genes (TPM~Status)
## Disclaimer: I utilized chatGPT 4o Mini (Dec 2024) to help troubleshoot error messages which arose in various 
## iterations of creating this script

# Load necessary libraries
library(ape)          # for phylogenetic tree manipulation
library(caper)        # for PGLS analysis
library(dplyr)        # for data manipulation
library(readr)        # for reading CSV files
library(tidyr)        # problem with gather function in dplyr
library(tidyverse)    # more problems!

# Step 1: Load the inputs

# Phylogenetic tree (Newick format or other supported formats)
tree <- read.tree("Tree_duck.newick")  # Phylogenetic tree in Newick format

# Gene TPM counts (rows = genes, columns = species)
#tpm_data <- read.csv("Duck_opsin_PGLS_eye.csv")  # Update with your TPM CSV file path
tpm_data <- read.csv("Duck_opsin_PGLS_eye_logtransformed.csv")  # Update with your TPM CSV file path

# Species Status data (species and their Status)
# diving status
status_data <- read.csv("Status_dive.csv")  # Update with your Status CSV file path
#sex morphology status
#status_data <- read.csv("Status_morph.csv")  # Update with your Status CSV file path
#water (marine, freshwater)
#status_data <- read.csv("Status_water.csv")  # Update with your Status CSV file path

# Gene list sorted IDs
gene_list <- read.csv("gene.csv")

# Step 2: Data Preprocessing

# Ensure species names match in all datasets
tree_species <- tree$tip.label

# Ensure that the species in the TPM data and Status data match those in the tree
status_data <- status_data %>%
  filter(Species %in% tree_species)  # Keep only species in the tree

# Convert Status into a factor (e.g., "status" should be the column in status_data)
status_data$Status <- as.factor(status_data$Status)

# Step 3: Merge TPM data with the Status data and gene list

# Add the gene names as a column in the TPM data
tpm_status_data <- tpm_data %>%
  gather(key = "Species", value = "TPM", -Genes) %>%
  left_join(status_data, by = "Species")

# Step 4: Run PGLS for each gene - if changed model, fix names in steps below
# Brownian Motion
#bm.Duck <- corPagel(0.5, phy = tree, fixed = FALSE, form = ~Species)
# Ornstein–Uhlenbeck model
ou.Duck <- corMartins(1, phy = tree, fixed = FALSE, form = ~Species)

# Define a function to run PGLS for each gene
run_pgls <- function(tpm_status_data, tree, status_data, gene) {
  # Subset the data for the current gene
  gene_tpm <- tpm_status_data %>%
    filter(Genes == gene) %>%
    select(Species, TPM) %>%
    left_join(status_data, by = "Species") %>%
    na.omit()  # Remove missing values
  
  # Ensure 'TPM' and 'Status_numeric' are numeric
  gene_tpm$TPM <- as.numeric(gene_tpm$TPM)  # Make sure TPM is numeric
  gene_tpm$Status_numeric <- as.numeric(gene_tpm$Status)  # Convert Status to numeric
  
  # Ensure 'Species' is in the final data frame before creating comparative data
  gene_tpm <- gene_tpm %>% 
    mutate(Species = as.character(Species))  # Make sure 'Species' is a character
  
  # Check for constant variables (no variation in TPM or Status_numeric)
  if (length(unique(gene_tpm$TPM)) == 1) {
    message("TPM has no variation (constant values). Skipping PGLS for gene: ", gene)
    return(NULL)  # Skip this gene
  }
  if (length(unique(gene_tpm$Status_numeric)) == 1) {
    message("Status_numeric has no variation (constant values). Skipping PGLS for gene: ", gene)
    return(NULL)  # Skip this gene
  }
  
  # Check correlation between TPM and Status_numeric
  correlation <- tryCatch({
    cor(gene_tpm$TPM, gene_tpm$Status_numeric, use = "complete.obs")
  }, error = function(e) {
    return(NA)  # Return NA if there is an error with the correlation
  })
  
  # Handle the case where correlation is NA (indicating no variation or other issues)
  if (is.na(correlation)) {
    message("Correlation is NA for gene: ", gene, ". Skipping PGLS.")
    return(NULL)  # Skip this gene
  }
  
  # Check for high correlation (greater than 0.9) and warn if found
  if (abs(correlation) > 0.9) {
    warning("High correlation between TPM and Status_numeric detected for gene: ", gene, ". This may cause issues with the model.")
  }
  
  # Create a comparative data object
  comp_data <- comparative.data(tree, gene_tpm, names.col = "Species", vcv = TRUE)
  
  # Extract the data from the comparative.data object
  comp_data_df <- comp_data$data
  
  # Ensure 'Species' is present in the dataframe for gls model
  if(!"Species" %in% names(comp_data_df)) {
    comp_data_df$Species <- gene_tpm$Species  # Add the 'Species' column if it's missing
  }
  
  # Fit the PGLS model
  tryCatch({
    pgls_model <- gls(TPM ~ Status_numeric, correlation = ou.Duck, data = comp_data_df, method = "ML", na.action = na.exclude)
    return(summary(pgls_model))
  }, error = function(e) {
    message("PGLS model fitting failed for gene: ", gene, " Error: ", e$message)
    return(NULL)
  })
}

# Step 5: Run the analysis for each gene and store results
pgls_results <- lapply(unique(tpm_status_data$Genes), function(gene) {
  result <- run_pgls(tpm_status_data, tree, status_data, gene)
  
  # Debugging: Print if the result is NULL
  if (is.null(result)) {
    message("No result for gene: ", gene)
  } else {
    message("Successfully ran PGLS for gene: ", gene)
  }
  
  # Check if the result is not NULL and add gene ID to it
  if (!is.null(result)) {
    result$Gene <- gene  # Add the gene ID to the result
  }
  
  return(result)
})

# Filter out any NULL results (if no data available for certain genes)
pgls_results <- pgls_results[!sapply(pgls_results, is.null)]

# Debugging: Check how many results are valid
message("Number of valid results: ", length(pgls_results))

# Step 6: Extract useful information from the PGLS results (e.g., p-values, coefficients)
# Prepare a data frame to store gene names and corresponding PGLS results
results_df <- do.call(rbind, lapply(pgls_results, function(pgls_summary) {
  
  # Debugging: Inspect the structure of the summary object
  if (!is.null(pgls_summary)) {
    message("Structure of summary for gene ", pgls_summary$Gene, ":")
    print(str(pgls_summary))  # Inspect the structure of the summary object
  }
  
  # Ensure the summary contains necessary components before extracting values
  if (!is.null(pgls_summary) && "tTable" %in% names(pgls_summary)) {
    # Try to extract p-value and coefficient for 'Status_numeric'
    p_value <- tryCatch({
      pgls_summary$tTable["Status_numeric", "p-value"]
    }, error = function(e) NA)
    
    coefficient <- tryCatch({
      pgls_summary$tTable["Status_numeric", "Value"]
    }, error = function(e) NA)
    
    # Check for variance-covariance (if available) or other fit statistics
    r_squared <- tryCatch({
      var_corr <- varCorr(pgls_summary)
      if (!is.null(var_corr)) {
        var_corr[1, 1]  # Modify if necessary to extract R-squared or variance components
      } else {
        NA
      }
    }, error = function(e) NA)
    
    # Return a row with gene name and the extracted metrics
    return(data.frame(Gene = pgls_summary$Gene, P_Value = p_value, Coefficient = coefficient, R_Squared = r_squared))
  } else {
    # If the summary is missing essential components, return NULL
    message("Missing essential components for gene: ", pgls_summary$Gene)
    return(NULL)
  }
}))

# Step 7: Remove NULL results and combine valid rows into a final data frame
if (is.list(results_df)) {
  # Remove NULL results
  results_df <- results_df[!sapply(results_df, is.null)]
  
  # Check if all elements are data frames
  if (all(sapply(results_df, is.data.frame))) {
    # Combine all the rows from the list into a single data frame
    results_df <- do.call(rbind, results_df)
  } else {
    stop("Some elements in results_df are not data frames. Cannot combine.")
  }
} else {
  warning("results_df is not a list at this point. Expected a list of data frames.")
}

# Check the structure of the final results_df before proceeding
str(results_df)

# Step 8: Transpose the results to have thousands of rows and a few columns
# Transpose results_df
results_df_trans <- t(results_df)

# Convert the transposed matrix into a data frame
results_df_trans <- as.data.frame(results_df_trans)

# Assign column names to the transposed data frame
colnames(results_df_trans) <- c("Gene", "P_Value", "Coefficient", "R_Squared")

# Print the final structure of results_df
str(results_df_trans)

# Step 9: Perform a FDR for pvalues to set a threshold
# Extract p-values for each gene
#results_df_trans$Adjusted_P_Value <- p.adjust(results_df_trans$P_Value, method = "BH")
results_df$Adjusted_P_Value <- p.adjust(results_df$P_Value, method = "BH")

# Write the transposed results to a CSV file
write.csv(results_df, "PGLS_results_eye_dive_transformed.csv", row.names = TRUE)

