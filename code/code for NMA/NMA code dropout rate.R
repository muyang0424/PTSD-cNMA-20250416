
library(readxl)
library(tidyverse)
library(netmeta)
library(ggplot2)
library(tidyr) 
library(openxlsx)

Armlevel_data <- readxl::read_xlsx("./data for NMA dropout rate.xlsx")

df_arm <- Armlevel_data %>%
  select(NO, matches("^Arm_[1-3]$")) %>%
  pivot_longer(
    cols = -NO,
    names_to = "arm_num_str", 
    values_to = "treatment",
    names_pattern = "Arm_([1-3])", 
    values_drop_na = TRUE 
  )

df_n <- Armlevel_data %>%
  select(NO, matches("^n\\.[1-3]$")) %>% 
  pivot_longer(
    cols = -NO,
    names_to = "arm_num_str",
    values_to = "n",
    names_pattern = "n\\.([1-3])", 
    values_drop_na = TRUE
  )

df_event <- Armlevel_data %>%
  select(NO, matches("^dropout_rate_a[1-3]$")) %>%
  pivot_longer(
    cols = -NO,
    names_to = "arm_num_str",
    values_to = "event",
    names_pattern = "dropout_rate_a([1-3])", 
    values_drop_na = TRUE
  )

df_long_dropout <- df_arm %>%
  full_join(df_n, by = c("NO", "arm_num_str")) %>%
  full_join(df_event, by = c("NO", "arm_num_str")) %>%
  filter(!is.na(treatment) & !is.na(n) & !is.na(event)) %>% 
  rename(study = NO) %>%
  mutate(
    n = as.numeric(n),
    event = as.numeric(event),
    arm_num_str = as.numeric(arm_num_str) 
  ) %>%
  mutate(
    event = ifelse(event > n, n, event) 
  ) %>%
  select(study, treatment, n, event, arm_num_str) 


p_dropout <- pairwise(
  treat = treatment,
  event = event,
  n = n,
  studlab = study,
  data = df_long_dropout,
  sm = "OR",
  allstudies = TRUE
)

m1_dropout.netmeta <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = p_dropout,
  sm = "OR",
  common = FALSE,
  random = TRUE,
  reference.group = "TAU",
  details.chkmultiarm = TRUE,
  small.values = "good",
  sep.trts = " vs ",
  tol.multiarm = 0.1
)

print(m1_dropout.netmeta)

netsplit_result_dropout <- netsplit(m1_dropout.netmeta)
print(netsplit_result_dropout)

decomp_result_dropout <- decomp.design(m1_dropout.netmeta)
print(decomp_result_dropout)

netgraph(
  m1_dropout.netmeta,
  plastic = FALSE,
  thickness = "number.of.studies",
  number.of.studies = FALSE,
  points = TRUE,
  cex = 0.7,
  cex.points = 2,
  col = "black",
  col.points = "black",
  lwd.max = 5,
  labels = m1_dropout.netmeta$trts
)
title("Network Plot for Dropout Rates")

# forest plot
forest(m1_dropout.netmeta, ref = "TAU", digits = 3, xlab = "Odds Ratio (OR) for Dropout", smlab = "OR")


# result matrix
or_matrix <- round(exp(m1_dropout.netmeta$TE.random), 3)
ci_lower_or_matrix <- round(exp(m1_dropout.netmeta$lower.random), 3)
ci_upper_or_matrix <- round(exp(m1_dropout.netmeta$upper.random), 3)

treatments_ordered_dropout <- m1_dropout.netmeta$trts

result.matrix_or_ci <- matrix("", nrow = nrow(or_matrix), ncol = ncol(or_matrix),
                              dimnames = list(treatments_ordered_dropout, treatments_ordered_dropout))
for (i in 1:nrow(or_matrix)) {
  for (j in 1:ncol(or_matrix)) {
    if (i == j) {
      result.matrix_or_ci[i, j] <- "-"
    } else {
      result.matrix_or_ci[i, j] <- paste(or_matrix[i, j], " (", ci_lower_or_matrix[i, j], ", ", ci_upper_or_matrix[i, j], ")", sep = "")
    }
  }
}

direct_or_matrix <- round(exp(m1_dropout.netmeta$TE.direct.random), 3)
ci_lower_direct_or_matrix <- round(exp(m1_dropout.netmeta$lower.direct.random), 3)
ci_upper_direct_or_matrix <- round(exp(m1_dropout.netmeta$upper.direct.random), 3)

direct_result.matrix_or_ci <- matrix("", nrow = nrow(direct_or_matrix), ncol = ncol(direct_or_matrix),
                                     dimnames = list(treatments_ordered_dropout, treatments_ordered_dropout))
for (i in 1:nrow(direct_or_matrix)) {
  for (j in 1:ncol(direct_or_matrix)) {
    if (i == j) {
      direct_result.matrix_or_ci[i, j] <- "-"
    } else if (!is.na(direct_or_matrix[i,j])) {
      direct_result.matrix_or_ci[i, j] <- paste(direct_or_matrix[i, j], " (", ci_lower_direct_or_matrix[i, j], ", ", ci_upper_direct_or_matrix[i, j], ")", sep = "")
    } else {
      direct_result.matrix_or_ci[i, j] <- "" 
    }
  }
}

result.matrix_or_ci[lower.tri(result.matrix_or_ci, diag = FALSE)] <- direct_result.matrix_or_ci[lower.tri(direct_result.matrix_or_ci, diag = FALSE)]
write.csv(result.matrix_or_ci, "dropout_result_matrix_OR_CI_20250510.csv", row.names = TRUE, fileEncoding = "UTF-8")


#  ranking (P-score)
ranks_dropout_nma <- netrank(m1_dropout.netmeta, small.values = "good")
print(ranks_dropout_nma$ranking.random)

p_score_common_dropout <- ranks_dropout_nma$ranking.common
p_score_random_dropout <- ranks_dropout_nma$ranking.random

p_score_table_dropout <- data.frame(
  Treatment = names(p_score_common_dropout),
  `P-score (common)` = round(p_score_common_dropout, 4),
  `P-score (random)` = round(p_score_random_dropout, 4),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

p_score_table_dropout <- p_score_table_dropout %>%
  arrange(desc(`P-score (random)`))
print(p_score_table_dropout, row.names = FALSE)
write.xlsx(p_score_table_dropout, "dropout_P_score_table_20250510.xlsx", rowNames = FALSE)


# funnel plot
treat_order_dropout <- m1_dropout.netmeta$trts
funnel_plot_dropout <- funnel(m1_dropout.netmeta,
                              order = treat_order_dropout,
                              pch = 19,
                              col = "black",
                              method.bias = "Egger",      
                              studlab = FALSE,
                              legend = FALSE,
                              cex.studlab = 0.7,
                              xlab = "Odds Ratio (Log Scale)" 
)
