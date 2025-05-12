
library(readxl)
library(dplyr)
library(tidyr)
library(rjags)
library(coda)
library(MCMCvis)
library(openxlsx) 

cNMA_raw_data <- readxl::read_xlsx("./cNMA data dropout rate.xlsx")

Ns <- nrow(cNMA_raw_data)
if (!"na" %in% names(cNMA_raw_data)) {
  stop("Column 'na' (number of arms per study) not found in the data.")
}
cNMA_raw_data$na <- as.numeric(cNMA_raw_data$na)
na_vec <- cNMA_raw_data$na

max_na_comp <- 3 
tryCatch({
  c1  <- cNMA_raw_data %>% separate(IVE, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp)))
  c2  <- cNMA_raw_data %>% separate(repeated_recalling, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c3  <- cNMA_raw_data %>% separate(brief_recalling, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c4  <- cNMA_raw_data %>% separate(CR, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c5  <- cNMA_raw_data %>% separate(BL, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c6  <- cNMA_raw_data %>% separate(NF, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c7  <- cNMA_raw_data %>% separate(IMAG, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c8  <- cNMA_raw_data %>% separate(PS, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c9  <- cNMA_raw_data %>% separate(BA, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c10 <- cNMA_raw_data %>% separate(RL, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c11 <- cNMA_raw_data %>% separate(HW, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp))) 
  c12 <- cNMA_raw_data %>% separate(other, paste0("a", 1:max_na_comp), sep = ",", convert = TRUE, fill = "right") %>% select(all_of(paste0("a", 1:max_na_comp)))
  
}, error = function(e) {
  stop("Error during component extraction. Check required columns (IVE, repeated_recalling, etc.). Error: ", e)
})

Nc <- 11 # 12 total components - 1 (new c5, BL) with no main effect

c1 <- as.matrix(c1); c2 <- as.matrix(c2); c3 <- as.matrix(c3); c4 <- as.matrix(c4)
c5 <- as.matrix(c5); c6 <- as.matrix(c6); c7 <- as.matrix(c7); c8 <- as.matrix(c8)
c9 <- as.matrix(c9); c10 <- as.matrix(c10); c11 <- as.matrix(c11); c12 <- as.matrix(c12)

n_cols <- paste0("n.", 1:max_na_comp) 
if (!all(n_cols %in% names(cNMA_raw_data))) stop(paste("Missing required sample size columns (", paste(n_cols, collapse=", "),")"))
n_matrix <- cNMA_raw_data %>% select(all_of(n_cols)) %>% mutate(across(everything(), as.numeric)) %>% as.matrix()

rate_cols <- paste0("dropout_rate_a", 1:max_na_comp)
if (!all(rate_cols %in% names(cNMA_raw_data))) stop(paste("Missing required dropout rate columns (", paste(rate_cols, collapse=", "),")"))

parse_rate <- function(rate_entry) {
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
r_matrix[is.na(n_matrix) | n_matrix <= 0] <- NA
if(any(r_matrix < 0, na.rm=TRUE)) {
  warning("Some calculated 'r' values were negative. Setting them to NA.")
  r_matrix[r_matrix < 0 & !is.na(r_matrix)] <- NA
}
if(any(r_matrix > n_matrix, na.rm = TRUE)) {
  warning("Some calculated 'r' values were > 'n'. Capping 'r' at 'n'.")
  r_matrix <- pmin(r_matrix, n_matrix, na.rm = FALSE) 
}

Ninter <- 1 # Only one interaction term
interactions_array <- array(NA, dim = c(Ns, max_na_comp, Ninter))

if(exists("c3") && exists("c5") && is.matrix(c3) && is.matrix(c5) && all(dim(c3) == dim(c5))) {
  interactions_array[, , 1] <- c3 * c5 } else {
  stop("Matrices c3 (brief_recalling) and/or c5 (BL) are missing or invalid for interaction calculation.")
}


jags_data <- list(
  r = r_matrix, n = n_matrix, na = na_vec, Ns = Ns,
  Nc = Nc, # Now Nc=11 MODIFIED
  Ninter = Ninter,
  
  c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, c6 = c6,
  c7 = c7, c8 = c8, c9 = c9, c10 = c10, c11 = c11, c12 = c12, 
  
  interactions = interactions_array 
)
stopifnot(dim(r_matrix) == c(Ns, max_na_comp), dim(n_matrix) == c(Ns, max_na_comp))
stopifnot(length(na_vec) == Ns)
stopifnot(dim(c1) == c(Ns, max_na_comp))


model_string <- "
model {
  for(i in 1:Ns) {
    w[i,1] <- 0
    delta[i,1] <- 0

    for (k in 1:na[i]) {
      r[i,k] ~ dbin(p[i,k], n[i,k])
      logit(p[i,k]) <- theta[i,k]
      theta[i,k] <- u[i] + delta[i,k]
    }

    for (k in 2:na[i]) {
      delta[i,k] ~ dnorm(md[i,k], precd[i,k])
      md[i,k] <- mean[i,k] + sw[i,k]
      precd[i,k] <- prect * 2 * (k-1) / k
      w[i,k] <- delta[i,k] - mean[i,k]
      sw[i,k] <- ifelse(k==2, w[i,1], sum(w[i,1:(k-1)]) / (k-1))
      mean[i,k] <- A1[i,k] - B1[i]

      # Log-odds for ARM K: d[1..11] maps to new c1-c4 and c6-c12.
      # No d term for new c5 (BL).
      A1[i,k] <- d[1]*(1-equals(c1[i,k],0)) +  # IVE
                   d[2]*(1-equals(c2[i,k],0)) +  # prolonged_recalling
                   d[3]*(1-equals(c3[i,k],0)) +  # brief_recalling
                   d[4]*(1-equals(c4[i,k],0)) +  # CR
                   # c5 (BL) has no main effect
                   d[5]*(1-equals(c6[i,k],0)) +  # NF
                   d[6]*(1-equals(c7[i,k],0)) +  # IMAG
                   d[7]*(1-equals(c8[i,k],0)) +  # PS
                   d[8]*(1-equals(c9[i,k],0)) +  # BA
                   d[9]*(1-equals(c10[i,k],0)) + # RL
                   d[10]*(1-equals(c11[i,k],0)) +# HW
                   d[11]*(1-equals(c12[i,k],0)) +# other
                   gamma[1] * interactions[i,k,1] # c3*c5 interaction
    }

    # Log-odds for ARM 1 (Reference arm for contrasts within study)
    B1[i] <-    d[1]*(1-equals(c1[i,1],0)) +
                d[2]*(1-equals(c2[i,1],0)) +
                d[3]*(1-equals(c3[i,1],0)) +
                d[4]*(1-equals(c4[i,1],0)) +
                # c5 (BL) has no main effect
                d[5]*(1-equals(c6[i,1],0)) +
                d[6]*(1-equals(c7[i,1],0)) +
                d[7]*(1-equals(c8[i,1],0)) +
                d[8]*(1-equals(c9[i,1],0)) +
                d[9]*(1-equals(c10[i,1],0)) +
                d[10]*(1-equals(c11[i,1],0)) +
                d[11]*(1-equals(c12[i,1],0)) +
                gamma[1] * interactions[i,1,1]

    u[i] ~ dnorm(0, 0.01)
  }

  tau ~ dnorm(0, 1)I(0,)
  prect <- 1 / pow(tau, 2)

  for(inter in 1:Ninter) { # Ninter is 1
    gamma[inter] ~ ddexp(0, lambda)
  }
  lambda ~ dgamma(1, 0.1)

  for(k in 1:Nc) { # Nc is 11
    d[k] ~ dnorm(0, 0.01)
  }
}
"

# Run JAGS Model
n_chains <- 4
n_adapt <- 10000 
n_iter <- 50000  
n_thin <- 25     

model_spec <- textConnection(model_string)
jags_model <- tryCatch({
  jags.model(model_spec, data = jags_data, n.chains = n_chains, n.adapt = n_adapt, quiet=FALSE)
}, error = function(e) {
  message("Error creating JAGS model. Check data inputs (dimensions, NAs) and model syntax."); message(e); return(NULL) })

if (is.null(jags_model)) { stop("JAGS model compilation failed.") }

params_to_monitor <- c("tau", "d", "gamma")
mcmc_samples <- tryCatch({
  coda.samples(jags_model, params_to_monitor, n.iter = n_iter, thin = n_thin)
}, error = function(e) { message("Error during MCMC sampling."); message(e); return(NULL) })

if (!is.null(mcmc_samples) && inherits(mcmc_samples, "mcmc.list")) {

  summary.ptsd <- MCMCsummary(mcmc_samples, Rhat = TRUE, n.eff = TRUE)
  
  d_param_names <- c(
    "IVE", "repeated_recalling", "brief_recalling", "CR", 
    "NF", "IMAG", "PS", "BA", "RL", "HW", "other"      
  )
  
  d_indices <- grep("^d\\[", rownames(summary.ptsd))
  if (length(d_indices) == Nc && Nc == 11) { # MODIFIED Nc check
    if(length(d_param_names) == length(d_indices)) {
      rownames(summary.ptsd)[d_indices] <- d_param_names
      print("Successfully renamed 'd' parameters.")
    } else {
      warning("Length mismatch: Found ", length(d_indices), " 'd' rows but have ", length(d_param_names), " names. Renaming failed.")
    }
  } else {
    warning(paste("Found", length(d_indices), "d parameters, but expected Nc =", Nc, ". Check model/MCMC output. Using default 'd[k]' names."))
  }
  
  gamma_indices <- grep("^gamma", rownames(summary.ptsd))
  if(length(gamma_indices) == Ninter && Ninter == 1) {
    gamma_summary <- summary.ptsd[gamma_indices, , drop = FALSE]
    interaction_name <- "brief_recalling-BL" 
    rownames(summary.ptsd)[gamma_indices] <- interaction_name 
    print("Successfully extracted and renamed 'gamma' parameter.")
  } else {
    warning(paste("Did not find exactly", Ninter, "'gamma' parameter(s) in summary (found", length(gamma_indices), "). Check MCMC output/model."))
    interaction_name <- "NOT_FOUND"
  }
  
  print(round(summary.ptsd, digits = 3))
  
  print("Generating Trace Plots for 'd' parameters...")
  if(any(grepl("^d\\[", colnames(as.matrix(mcmc_samples[[1]]))))) {
    tryCatch({ MCMCtrace(mcmc_samples, pdf = FALSE, params = "d", wd = getwd(), n.eff = TRUE, Rhat = TRUE)
    }, error = function(e) { message("Could not generate MCMCtrace plot for 'd'."); message(e) })
  } else { print("Skipping trace plot: No 'd' parameters found in MCMC samples.") }
  
  output_excel_file <- paste0("./Results_Dropout_cNMA_Nc", Nc, "_Corrected_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
  tryCatch({
    summary_df_for_excel <- as.data.frame(round(summary.ptsd, digits = 3))
    summary_df_for_excel <- tibble::rownames_to_column(summary_df_for_excel, "Parameter")
    
    wb <- createWorkbook(); addWorksheet(wb, paste0("cNMA Results Nc", Nc," (LogOdds)")) 
    writeData(wb, paste0("cNMA Results Nc", Nc," (LogOdds)"), summary_df_for_excel)
    addStyle(wb, sheet = 1, style = createStyle(textDecoration = "bold"), rows = 1, cols = 1:(ncol(summary_df_for_excel)))
    setColWidths(wb, sheet = 1, cols = 1:ncol(summary_df_for_excel), widths = "auto")
    saveWorkbook(wb, output_excel_file, overwrite = TRUE)
    print(paste("Full summary results saved to", output_excel_file))
  }, error = function(e) { message("Error saving results to Excel file: ", output_excel_file); message(e) })
  
  print("Calculating DIC (may take time)...")
  dic_results <- tryCatch({ dic.samples(jags_model, n.iter = 5000, type = "pD") 
  }, error = function(e){ message("Could not calculate DIC."); message(e); return(NULL) })
  if(!is.null(dic_results)){ print(dic_results) }
  
} else { print("MCMC sampling failed or produced invalid output.") }
