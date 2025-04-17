
# --- 0. Setup ---
# install.packages(c("readxl", "dplyr", "tidyr", "rjags", "coda", "MCMCvis", "openxlsx")) # Run once if needed

library(readxl)
library(dplyr)
library(tidyr)
library(rjags)
library(coda)
library(MCMCvis)
library(openxlsx) # For writing to Excel

 source("functions.R", encoding = "UTF-8") # Keep if you use custom functions

# --- 1. Load Data (User's Original Structure) ---

# Set working directory (adjust path as needed)
# setwd("C:/Path/To/Your/Directory")

# Load raw data
raw_data_file <- "for 脱落率分析 20250417.xlsx" # Make sure this path is correct
cNMA_raw_data <- readxl::read_xlsx(raw_data_file)

# Get number of studies (rows) and ensure 'na' column exists and is numeric
Ns <- nrow(cNMA_raw_data)
if (!"na" %in% names(cNMA_raw_data)) {
  stop("Column 'na' (number of arms per study) not found in the data.")
}
cNMA_raw_data$na <- as.numeric(cNMA_raw_data$na)
if(any(is.na(cNMA_raw_data$na) | cNMA_raw_data$na <= 0)) {
  warning("Some studies have missing or invalid 'na' values. They might cause errors or be excluded implicitly.")
  # Optionally filter: cNMA_raw_data <- cNMA_raw_data %>% filter(!is.na(na), na > 0)
  # Ns <- nrow(cNMA_raw_data) # Update Ns if filtering
}
na_vec <- cNMA_raw_data$na

# --- 2. Extract Components (User's Original Structure - For Nc=10) ---
# c1-c3, c5-c11 needed for main effects. c2, c4 needed for interaction.
# Ensure component column names match Excel.
tryCatch({
  c1 <- cNMA_raw_data %>% separate(IVE, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c2 <- cNMA_raw_data %>% separate(reliving, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c3 <- cNMA_raw_data %>% separate(CR, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c4 <- cNMA_raw_data %>% separate(BL, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4) # For interaction ONLY
  c5 <- cNMA_raw_data %>% separate(NF, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c6 <- cNMA_raw_data %>% separate(IMAG, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c7 <- cNMA_raw_data %>% separate(PS, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c8 <- cNMA_raw_data %>% separate(BA, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c9 <- cNMA_raw_data %>% separate(RL, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c10 <- cNMA_raw_data %>% separate(HW, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
  c11 <- cNMA_raw_data %>% separate(other, c("a1","a2","a3","a4"), sep = ",", convert = TRUE) %>% select(a1:a4)
}, error = function(e) {
  stop("Error during component extraction. Check required columns (IVE, reliving, CR, BL, NF, IMAG, PS, BA, RL, HW, other). Error: ", e)
})

# *** Set Nc to 10, as intended ***
Nc <- 10 # Number of main 'd' effects

# Convert component dataframes to matrices
c1 <- as.matrix(c1); c2 <- as.matrix(c2); c3 <- as.matrix(c3); c4 <- as.matrix(c4)
c5 <- as.matrix(c5); c6 <- as.matrix(c6); c7 <- as.matrix(c7); c8 <- as.matrix(c8)
c9 <- as.matrix(c9); c10 <- as.matrix(c10); c11 <- as.matrix(c11)

max_na_comp <- ncol(c1) # Max arms based on component separation (likely 4)

# --- 3. Extract n and Calculate r ---
# (Same code as before to get n_matrix and r_matrix)
n_cols <- paste0("n.", 1:max_na_comp)
if (!all(n_cols %in% names(cNMA_raw_data))) stop("Missing required sample size columns (n.1, ...)")
n_matrix <- cNMA_raw_data %>% select(all_of(n_cols)) %>% mutate(across(everything(), as.numeric)) %>% as.matrix()

rate_cols <- paste0("dropout_rate_a", 1:max_na_comp)
if (!all(rate_cols %in% names(cNMA_raw_data))) stop("Missing required dropout rate columns (dropout_rate_a1, ...)")

parse_rate <- function(rate_entry) { # (Same robust function as before)
  if (is.numeric(rate_entry)) { return(rate_entry) }
  if (is.character(rate_entry) && grepl("/", rate_entry, fixed = TRUE)) {
    parts <- tryCatch(strsplit(rate_entry, "/")[[1]], error=function(e) NULL)
    if(!is.null(parts) && length(parts) == 2) {
      num <- suppressWarnings(as.numeric(parts[1])); den <- suppressWarnings(as.numeric(parts[2]))
      if (!is.na(num) && !is.na(den) && den != 0) { return(num / den) } } }
  return(NA_real_)
}
rate_matrix <- cNMA_raw_data %>% select(all_of(rate_cols)) %>% mutate(across(everything(), ~sapply(., parse_rate))) %>% as.matrix()

r_matrix <- round(rate_matrix * n_matrix)
r_matrix[is.na(n_matrix) | n_matrix <= 0] <- NA # Ensure r is NA if n is invalid
if(any(r_matrix < 0, na.rm=TRUE)) {
  warning("Some calculated 'r' values were negative. Setting them to NA.")
  r_matrix[r_matrix < 0 & !is.na(r_matrix)] <- NA
}
if(any(r_matrix > n_matrix, na.rm = TRUE)) {
  warning("Some calculated 'r' values were > 'n'. Capping 'r' at 'n'.")
  r_matrix <- pmin(r_matrix, n_matrix, na.rm = FALSE)
}

# --- 4. Prepare Interactions ---
Ninter <- 1
interactions_array <- array(NA, dim = c(Ns, max_na_comp, Ninter))
# Ensure c2 and c4 are valid matrices before multiplying
if(exists("c2") && exists("c4") && is.matrix(c2) && is.matrix(c4) && all(dim(c2) == dim(c4))) {
  interactions_array[, , 1] <- c2 * c4 # c2 ('reliving') * c4 ('BL')
} else {
  stop("Matrices c2 and/or c4 are missing or invalid for interaction calculation.")
}


# --- 5. Prepare Data List for JAGS (Correct for Nc=10) ---
jags_data <- list(
  r = r_matrix, n = n_matrix, na = na_vec, Ns = Ns,
  Nc = Nc,               # Now Nc=10
  Ninter = Ninter,
  # Component matrices needed for model (c1-c3, c5-c11)
  c1 = c1, c2 = c2, c3 = c3,
  # c4 is used only via interactions_array below
  c5 = c5, c6 = c6, c7 = c7, c8 = c8, c9 = c9, c10 = c10, c11 = c11,
  interactions = interactions_array # Contains c2*c4 product
)
# Check that essential data has correct dimensions
stopifnot(dim(r_matrix) == c(Ns, max_na_comp), dim(n_matrix) == c(Ns, max_na_comp))
stopifnot(length(na_vec) == Ns)
stopifnot(dim(c1) == c(Ns, max_na_comp)) # Check one component matrix dim


# --- 6. Define JAGS Model String (Correct for Nc=10) ---
# Model exactly implementing 10 main effects (d1-d10 for c1-c3, c5-c11) and 1 interaction (gamma1 for c2*c4)
model_string <- "
model {
  for(i in 1:Ns) {
    # Basic checks for valid study data (na > 0) should ideally be done before JAGS
    w[i,1] <- 0
    delta[i,1] <- 0

    # Likelihood loop for all arms k
    for (k in 1:na[i]) {
       # JAGS handles NA r[i,k] if n[i,k] is provided and > 0
       # Add safety for p[i,k] if n[i,k] could be 0? (Not standard, assume n>=1 for included arms)
       r[i,k] ~ dbin(p[i,k], n[i,k])
       logit(p[i,k]) <- theta[i,k]
       theta[i,k] <- u[i] + delta[i,k]
    }

    # Multi-arm correlation loop for arms k >= 2
    for (k in 2:na[i]) {
      delta[i,k] ~ dnorm(md[i,k], precd[i,k])
      md[i,k] <- mean[i,k] + sw[i,k]
      # Ensure precision is positive; tau prior already ensures tau > 0
      precd[i,k] <- prect * 2 * (k-1) / k
      w[i,k] <- delta[i,k] - mean[i,k]
      # Handle k=2 edge case for sw calculation more robustly
      sw[i,k] <- ifelse(k==2, w[i,1], sum(w[i,1:(k-1)]) / (k-1)) # w[i,1] is 0, so sw[i,2] = 0
      mean[i,k] <- A1[i,k] - B1[i]

      # Log-odds for ARM K: d[1..10] maps to c1-c3, c5-c11. NO d for c4.
      A1[i,k] <- d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0)) +
                 d[3]*(1-equals(c3[i,k],0)) +
                 # No term for c4 main effect
                 d[4]*(1-equals(c5[i,k],0)) + d[5]*(1-equals(c6[i,k],0)) +
                 d[6]*(1-equals(c7[i,k],0)) + d[7]*(1-equals(c8[i,k],0)) +
                 d[8]*(1-equals(c9[i,k],0)) + d[9]*(1-equals(c10[i,k],0)) +
                 d[10]*(1-equals(c11[i,k],0)) + # d[10] is the last main effect
                 gamma[1] * interactions[i,k,1]  # c2*c4 interaction using pre-calculated product
    }

    # Log-odds for ARM 1: d[1..10] maps to c1-c3, c5-c11. NO d for c4.
    B1[i] <-   d[1]*(1-equals(c1[i,1],0)) + d[2]*(1-equals(c2[i,1],0)) +
               d[3]*(1-equals(c3[i,1],0)) +
               # No term for c4 main effect
               d[4]*(1-equals(c5[i,1],0)) + d[5]*(1-equals(c6[i,1],0)) +
               d[6]*(1-equals(c7[i,1],0)) + d[7]*(1-equals(c8[i,1],0)) +
               d[8]*(1-equals(c9[i,1],0)) + d[9]*(1-equals(c10[i,1],0)) +
               d[10]*(1-equals(c11[i,1],0)) + # d[10] is the last main effect
               gamma[1] * interactions[i,1,1]  # c2*c4 interaction

    # Prior for study-specific baseline effects
    u[i] ~ dnorm(0, 0.01) # Vague prior for log-odds baseline
  }

  # Priors for remaining parameters
  tau ~ dnorm(0, 1)I(0,)      # Heterogeneity SD on log-odds scale (Half-Normal)
  prect <- 1 / pow(tau, 2)    # Heterogeneity precision

  # Prior for interaction effect(s) - Laplace (LASSO) prior
  for(inter in 1:Ninter) {
     gamma[inter] ~ ddexp(0, lambda) # Log-odds ratio interaction effect
  }
  lambda ~ dgamma(1, 0.1)        # Prior for Laplace shrinkage parameter

  # Priors for component main effects (Nc = 10)
  for(k in 1:Nc) {
    d[k] ~ dnorm(0, 0.01) # Vague prior for log-odds ratio main effects
  }
}
"

# --- 7. Run JAGS Model ---
n_chains <- 4
n_adapt <- 10000
n_iter <- 50000 # Increased iterations slightly
n_thin <- 25    # Increased thinning slightly

model_spec <- textConnection(model_string)
jags_model <- tryCatch({
  jags.model(model_spec, data = jags_data, n.chains = n_chains, n.adapt = n_adapt, quiet=FALSE)
}, error = function(e) {
  message("Error creating JAGS model. Check data inputs (dimensions, NAs) and model syntax."); message(e); return(NULL) })

if (is.null(jags_model)) { stop("JAGS model compilation failed.") }

print("Starting MCMC sampling...")
params_to_monitor <- c("tau", "d", "gamma") # Monitor effects and heterogeneity
mcmc_samples <- tryCatch({
  coda.samples(jags_model, params_to_monitor, n.iter = n_iter, thin = n_thin)
}, error = function(e) { message("Error during MCMC sampling."); message(e); return(NULL) })
print("MCMC sampling finished.")

# --- 8. Process and Summarize Results (User's Simpler Structure + Excel) ---
if (!is.null(mcmc_samples) && inherits(mcmc_samples, "mcmc.list")) {
  print("Summarizing results...")
  
  # Get the main summary table (using user's variable name)
  summary.ptsd <- MCMCsummary(mcmc_samples, Rhat = TRUE, n.eff = TRUE) # Use this name
  
  # Define the 10 names for the d parameters
  d_param_names <- c("IVE", "reliving", "CR", "NF", "IMAG", "PS", "BA", "RL", "HW", "other")
  
  # Rename the d parameters using grep (safer)
  d_indices <- grep("^d\\[", rownames(summary.ptsd))
  if (length(d_indices) == Nc && Nc == 10) {
    # *** This is the likely point of failure if lengths mismatch ***
    # Let's add a check *before* assignment
    if(length(d_param_names) == length(d_indices)) {
      rownames(summary.ptsd)[d_indices] <- d_param_names
      print("Successfully renamed 'd' parameters.")
    } else {
      warning("Length mismatch: Found ", length(d_indices), " 'd' rows but have ", length(d_param_names), " names. Renaming failed.")
    }
  } else {
    warning(paste("Found", length(d_indices), "d parameters, but expected Nc = 10. Check model/MCMC output. Using default 'd[k]' names."))
  }
  
  # Extract and rename gamma summary (simpler approach)
  # Important: Use the potentially modified summary.ptsd object here
  gamma_indices <- grep("^gamma", rownames(summary.ptsd)) # Grep for name starting with gamma
  # (might be gamma or gamma[1])
  if(length(gamma_indices) == 1) {
    gamma_summary <- summary.ptsd[gamma_indices, , drop = FALSE]
    interaction_name <- "reliving-BL" # Use the specific name based on components
    rownames(gamma_summary) <- interaction_name # Rename the single row
    print("Successfully extracted and renamed 'gamma' parameter.")
  } else {
    warning(paste("Did not find exactly 1 'gamma' parameter in summary (found", length(gamma_indices), "). Check MCMC output/model."))
    gamma_summary <- data.frame(matrix(NA, 0, ncol(summary.ptsd), dimnames=list(NULL, colnames(summary.ptsd)))) # Empty placeholder
    interaction_name <- "NOT_FOUND"
  }
  
  # Print results as requested
  print("--- Full Summary Results (Log-Odds Scale) ---")
  print(round(summary.ptsd, digits = 3)) # Print the main table (potentially with renamed rows)
  
  print(paste("--- Interaction (", interaction_name, ") Summary (Log-Odds Scale) ---"))
  if(nrow(gamma_summary) == 1) { # Check if gamma_summary is valid
    print(round(gamma_summary, digits = 3)) # Print the isolated gamma summary
  } else {
    print("< Gamma parameter not found or multiple found >")
  }
  
  # Generate trace plot specifically for d parameters (no PDF)
  print("Generating Trace Plots for 'd' parameters...")
  if(length(grep("^d\\[", colnames(as.matrix(mcmc_samples[[1]])))) > 0) {
    tryCatch({ MCMCtrace(mcmc_samples, pdf = FALSE, params = "d", wd = getwd(), n.eff = TRUE, Rhat = TRUE)
    }, error = function(e) { message("Could not generate MCMCtrace plot for 'd'."); message(e) })
  } else { print("Skipping trace plot: No 'd' parameters found in MCMC samples.") }
  
  
  # --- Write the main summary table to EXCEL ---
  output_excel_file <- paste0("./Results_Dropout_cNMA_Nc10_Corrected_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
  tryCatch({
    # Use the final summary.ptsd object which has renamed rows (if successful)
    summary_df_for_excel <- as.data.frame(round(summary.ptsd, digits = 3))
    summary_df_for_excel <- tibble::rownames_to_column(summary_df_for_excel, "Parameter") # Add final rownames as column
    
    wb <- createWorkbook(); addWorksheet(wb, "cNMA Results (LogOdds)")
    writeData(wb, "cNMA Results (LogOdds)", summary_df_for_excel)
    addStyle(wb, sheet = 1, style = createStyle(textDecoration = "bold"), rows = 1, cols = 1:(ncol(summary_df_for_excel)))
    setColWidths(wb, sheet = 1, cols = 1:ncol(summary_df_for_excel), widths = "auto")
    saveWorkbook(wb, output_excel_file, overwrite = TRUE)
    print(paste("Full summary results saved to", output_excel_file))
  }, error = function(e) { message("Error saving results to Excel file: ", output_excel_file); message(e) })
  
  # Calculate DIC
  print("Calculating DIC (may take time)...")
  dic_results <- tryCatch({ dic.samples(jags_model, n.iter = 5000, type = "pD")
  }, error = function(e){ message("Could not calculate DIC."); message(e); return(NULL) })
  if(!is.null(dic_results)){ print(dic_results) }
  
} else { print("MCMC sampling failed or produced invalid output.") }

# Close connections at the very end
# closeAllConnections() # Usually not necessary unless specific file handles were opened

print("Script finished.")

# --- 8. Process and Summarize Results (Corrected for Nc=10) ---
# 输出结果
summary.dropout <- MCMCsummary(mcmc_samples, Rhat = TRUE, n.eff = TRUE)

# 为 d 参数重命名（10 个主效应，去掉 c4）
rownames(summary.ptsd)[1:10] <- c("IVE", "reliving", "CR", "NF", "IMAG", "PS", "BA", "RL", "HW", "other")

# 为 gamma 参数命名（仅 c2-c4）
gamma_summary <- mcmc_samples[grep("gamma", rownames(summary.dropout)), , drop = FALSE]
rownames(gamma_summary) <- "c2-c4"

# 打印结果
print(round(summary.dropout, digits = 3))
print(round(gamma_summary, digits = 3))

MCMCtrace(samps.ptsd, pdf = FALSE, params = "d", wd = getwd(), n.eff = TRUE, Rhat = TRUE)
write.csv(round(summary.ptsd, digits = 3), "./outputs_20250416_2.csv")

# 计算 DIC
DIC.ptsd <- dic.samples(jags.m.ptsd, n.iter = 5000)
print(DIC.ptsd)
closeAllConnections()
