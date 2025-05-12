# load libraries ----------------------------------------------------------
install.packages("netmeta")
library(netmeta)
install.packages("coda")
library(coda)
install.packages("rjags")
library(rjags)
install.packages("MCMCvis")
library(MCMCvis)
install.packages("tidyverse")
library(tidyverse)
install.packages("igraph")
library(igraph)
library(dplyr)
install.packages("meta")   # If not already installed
library(meta)
install.packages("gtools")
install.packages("openxlsx")

cNMA_raw_data <- readxl::read_xlsx("./cNMA data posttreatment PTSD severity.xlsx")
cNMA_raw_data <- cNMA_raw_data %>%
  mutate(na = na) %>%
  mutate(Ns = nrow(.))

cNMA_raw_data$na <- cNMA_raw_data$na %>%
  as.numeric()

y = cNMA_raw_data %>%
  select(y.1,y.2,y.3,y.4) %>%
  char_to_num()

n = cNMA_raw_data %>%
  select(n.1,n.2,n.3,n.4) %>%
  char_to_num()

sd = cNMA_raw_data %>%
  select(sd.1,sd.2,sd.3,sd.4) %>%
  char_to_num()

na = cNMA_raw_data %>%
  select(na) %>%
  unlist() %>%
  as.vector()
c1 = cNMA_raw_data %>%
  separate(IVE, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c2 = cNMA_raw_data %>%
  separate(repeated_recalling, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c3 = cNMA_raw_data %>%
  separate(brief_recalling, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c4 = cNMA_raw_data %>%
  separate(CR, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c5 = cNMA_raw_data %>%
  separate(BL,c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c6 = cNMA_raw_data %>%
  separate(NF, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c7 = cNMA_raw_data %>%
  separate(IMAG, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c8 = cNMA_raw_data %>%
  separate(PS, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c9 = cNMA_raw_data %>%
  separate(BA, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c10 = cNMA_raw_data %>%
  separate(RL, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c11 = cNMA_raw_data %>%
  separate(HW, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c12 = cNMA_raw_data %>%
  separate(other, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c13 = cNMA_raw_data %>%
  separate(waiting, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)


s_within = cNMA_raw_data %>%
  select(s_within) %>%
  unlist() %>%
  as.vector() %>%
  as.numeric()
s_m_within = cNMA_raw_data %>%
  select(s_m_within) %>%
  unlist() %>%
  as.vector() %>%
  as.numeric()



Nc=13
Ns=nrow(y)

c1 <- as.matrix(c1)
c2 <- as.matrix(c2)
c3 <- as.matrix(c3)
c4 <- as.matrix(c4)
c5 <- as.matrix(c5)
c6 <- as.matrix(c6)
c7 <- as.matrix(c7)
c8 <- as.matrix(c8)
c9 <- as.matrix(c9)
c10 <- as.matrix(c10)
c11 <- as.matrix(c11)
c12 <- as.matrix(c12)
c13 <- as.matrix(c13)

components <- list(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)

max_na <- max(na)

existing_pairs_reduced <- matrix(c(3, 5), ncol = 1)  
Ninter_reduced <- 1  

interactions_reduced <- array(NA, dim = c(Ns, max_na, Ninter_reduced))
interactions_reduced[, , 1] <- components[[5]] * components[[3]]  

Nc <- 12  # exluded BL

data <- list(
  y = y, n = n, sd = sd,
  Nc = 12, Ns = Ns, na = na,
  c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, 
  c6 = c6, c7 = c7, c8 = c8, c9 = c9, c10 = c10,
  c11 = c11, c12 = c12,c13=c13, 
  s_within = s_within, Ninter = Ninter_reduced, interactions = interactions_reduced
)

### model specification ###
model.ptsd.string <- "
model {
  for(i in 1:Ns) {
    w[i,1] <- 0
    delta[i,1] <- 0

    for (k in 1:na[i])  {
      se[i,k] <- sd[i,k]/sqrt(n[i,k])
      prec[i,k] <- (1/(se[i,k]*se[i,k]))
      y[i,k] ~ dnorm(phi[i,k], prec[i,k])
      phi[i,k] <- theta[i,k] * (s_within[i]) 
      theta[i,k] <- u[i] + delta[i,k]
    }

    # Multi-arm section 
    for (k in 2:na[i]) {
      delta[i,k] ~ dnorm(md[i,k], precd[i,k])
      md[i,k] <- mean[i,k] + sw[i,k]
      precd[i,k] <- 2 * prect * (k-1) / k
      w[i,k] <- (delta[i,k] - mean[i,k])
      sw[i,k] <- ifelse(k==2, w[i,1], sum(w[i,1:(k-1)]) / (k-1)) 

      mean[i,k] <- A1[i,k] - B1[i]

      A1[i,k] <-
        d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0)) +
        d[3]*(1-equals(c3[i,k],0)) + d[4]*(1-equals(c4[i,k],0)) +
        d[5]*(1-equals(c6[i,k],0)) +
        # No term for c5 main effect
        d[6]*(1-equals(c7[i,k],0)) + d[7]*(1-equals(c8[i,k],0)) +
        d[8]*(1-equals(c9[i,k],0)) + d[9]*(1-equals(c10[i,k],0)) +
        d[10]*(1-equals(c11[i,k],0)) + d[11]*(1-equals(c12[i,k],0)) +
        d[12]*(1-equals(c13[i,k],0))  +
        gamma[1] * interactions[i,k,1] 
    } 
    
    B1[i] <-
      d[1]*(1-equals(c1[i,1],0)) + d[2]*(1-equals(c2[i,1],0)) +
      d[3]*(1-equals(c3[i,1],0)) + d[4]*(1-equals(c4[i,1],0)) +
      d[5]*(1-equals(c6[i,1],0)) +
      # No term for c5 main effect
      d[6]*(1-equals(c7[i,1],0)) + d[7]*(1-equals(c8[i,1],0)) +
      d[8]*(1-equals(c9[i,1],0)) + d[9]*(1-equals(c10[i,1],0)) +
      d[10]*(1-equals(c11[i,1],0)) + d[11]*(1-equals(c12[i,1],0)) +
      d[12]*(1-equals(c13[i,1],0)) + 
      gamma[1] * interactions[i,1,1] 

  } 

  # Priors 
  for (i in 1:Ns) { u[i] ~ dnorm(0, .001) }
  tau ~ dnorm(0, 0.1)I(0,)
  prect <- 1/pow(tau, 2)

  # Priors for interaction effect and its shrinkage parameter
  gamma[1] ~ ddexp(0, lambda)  # Assumes only 1 interaction term (Ninter=1)
  lambda ~ dgamma(0.01, 0.01)

  for(k in 1:Nc) {
    d[k] ~ dnorm(0, .001) # Prior for the Nc=13 main effects
  }
}
"

# run model
model.ptsd.spec <- textConnection(model.ptsd.string)
jags.m.ptsd <- jags.model(model.ptsd.spec, data = data, n.chains = 4, n.adapt = 10000)
params <- c("tau", "d", "gamma")
samps.ptsd <- coda.samples(jags.m.ptsd, params, n.iter = 20000, thin = 10)

summary.ptsd <- MCMCsummary(samps.ptsd, Rhat = TRUE, n.eff = TRUE)

rownames(summary.ptsd)[1:12] <- c( "IVE", "repeated_recalling", "brief_recalling", "CR", "NF", "IMAG", "PS", "BA", "RL", "HW", "other", "waiting")

gamma_summary <- summary.ptsd[grep("gamma", rownames(summary.ptsd)), , drop = FALSE]
rownames(gamma_summary) <- "c5-c3"

print(round(summary.ptsd, digits = 3))
print(round(gamma_summary, digits = 3))

MCMCtrace(samps.ptsd, pdf = FALSE, params = "d", wd = getwd(), n.eff = TRUE, Rhat = TRUE)

DIC.ptsd <- dic.samples(jags.m.ptsd, n.iter = 5000)
print(DIC.ptsd)

### specify the posterior distribution of meaningful combinations ###

posterior_samples_matrix <- do.call(rbind, samps.ptsd)

Nc <- 12 
d_samples <- posterior_samples_matrix[, paste0("d[", 1:Nc, "]")]

component_names <- c("IVE", "repeated_recalling", "brief_recalling", "CR", "NF", 
                     "IMAG", "PS", "BA", "RL", "HW", "other", "waiting")

excluded_components <- c("waiting", "IMAG", "repeated_recalling")

component_names_for_combinations <- setdiff(component_names, excluded_components)

colnames(d_samples) <- component_names

gamma1_samples <- posterior_samples_matrix[, "gamma", drop = FALSE]
colnames(gamma1_samples) <- "BL_x_brief_recalling_interaction" # (c5*c3)

calculate_combination_summary <- function(combined_samples, combination_name) {
  mean_val <- mean(combined_samples)
  median_val <- median(combined_samples)
  ci <- quantile(combined_samples, probs = c(0.025, 0.975))
  return(data.frame(
    Combination = combination_name,
    Mean = mean_val,
    Median = median_val,
    CI_Lower = ci[1],
    CI_Upper = ci[2],
    CI_Spans_Zero = (ci[1] < 0 & ci[2] > 0) | (ci[1] > 0 & ci[2] < 0) # More robust check
  ))
}

all_combination_results <- list()

# Combinations of 2
combinations_2 <- utils::combn(component_names_for_combinations, 2, simplify = FALSE)

for (combo in combinations_2) {
  comp1_name <- combo[1]
  comp2_name <- combo[2]
  
  # Sum posterior samples for the two main effects
  combined_main_effects_samples <- d_samples[, comp1_name] + d_samples[, comp2_name]
  
  combo_name_main <- paste(comp1_name, comp2_name, sep = " + ")
  all_combination_results[[combo_name_main]] <- calculate_combination_summary(combined_main_effects_samples, combo_name_main)
  
  if ("brief_recalling" %in% combo) {
    combined_with_interaction_samples <- combined_main_effects_samples + gamma1_samples[, "BL_x_brief_recalling_interaction"]
    combo_name_interaction <- paste0(combo_name_main, " + BL (via interaction)")
    all_combination_results[[combo_name_interaction]] <- calculate_combination_summary(combined_with_interaction_samples, combo_name_interaction)
  }
}

# Combinations of 3
combinations_3 <- combn(component_names_for_combinations, 3, simplify = FALSE)

for (combo in combinations_3) {
  comp1_name <- combo[1]
  comp2_name <- combo[2]
  comp3_name <- combo[3]
  
  # Sum posterior samples for the three main effects
  combined_main_effects_samples <- d_samples[, comp1_name] + d_samples[, comp2_name] + d_samples[, comp3_name]
  
  combo_name_main <- paste(comp1_name, comp2_name, comp3_name, sep = " + ")
  all_combination_results[[combo_name_main]] <- calculate_combination_summary(combined_main_effects_samples, combo_name_main)
  
  if ("brief_recalling" %in% combo) {
    combined_with_interaction_samples <- combined_main_effects_samples + gamma1_samples[, "BL_x_brief_recalling_interaction"]
    combo_name_interaction <- paste0(combo_name_main, " + BL (via interaction)")
    all_combination_results[[combo_name_interaction]] <- calculate_combination_summary(combined_with_interaction_samples, combo_name_interaction)
  }
}

if ("brief_recalling" %in% colnames(d_samples) && "BL_x_brief_recalling_interaction" %in% colnames(gamma1_samples)) {
  
  brief_recalling_samples <- d_samples[, "brief_recalling"]
  
  br_plus_BL_interaction_samples <- brief_recalling_samples + gamma1_samples[, "BL_x_brief_recalling_interaction"]
  
  combo_name_br_BL <- "brief_recalling + BL (via interaction)"
  
  all_combination_results[[combo_name_br_BL]] <- calculate_combination_summary(
    br_plus_BL_interaction_samples, 
    combo_name_br_BL
  )
}

final_results_df <- do.call(rbind, all_combination_results)
rownames(final_results_df) <- NULL 

significant_combinations <- final_results_df %>%
  filter(CI_Spans_Zero == FALSE)

print("All Combination Results (Mean, Median, 95% CI):")
print(final_results_df, row.names = FALSE)

print("Significant Combinations (95% CI does not span zero):")
if (nrow(significant_combinations) > 0) {
  print(significant_combinations, row.names = FALSE)
} else {
  print("No two-component or three-component combinations (including BL interaction scenarios) found with 95% CI not spanning zero.")
}
