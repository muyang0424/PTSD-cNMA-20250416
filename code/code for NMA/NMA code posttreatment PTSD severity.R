library(readxl)
library(tidyverse)
library(netmeta)
library(ggplot2)
library(tidyr)
library(openxlsx)

Armlevel_data <- readxl::read_xlsx("./data for NMA PTSD posttreatment severity.xlsx")

df_long <- Armlevel_data %>%
  pivot_longer(
    cols = -`NO`,  
    names_to = c(".value", "arm"),
    names_pattern = "([^_.]+)[_.]([1-4])",
    values_drop_na = TRUE
  ) %>%
  filter(!is.na(Arm)) %>%
  rename(study = `NO`)  # 重命名为 study


# Use pairwise to generate SMDs from df_long
p1 <- pairwise(
  treat = Arm,
  n = n,
  mean = y,
  sd = sd,
  studlab = study,
  data = df_long,
  sm = "SMD"
)

# Run network meta-analysis
m1.netmeta <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = p1,
  sm = "SMD",
  common = FALSE,
  random = TRUE,
  reference.group = "WL",
  details.chkmultiarm = TRUE,
  small.values = "good",
  sep.trts = " vs ",
  tol.multiarm = 0.1
)

print(m1.netmeta)

netsplit_result <- netsplit(m1.netmeta)
print(netsplit_result)

decomp.design(m1.netmeta)

netgraph(
  m1.netmeta,
  plastic = FALSE,                 
  thickness = "number.of.studies", 
  number.of.studies = FALSE,       
  points = TRUE,                   
  cex = 0.8,                       
  cex.points = 2,                  
  col = "black",                   
  col.points = "black",            
  lwd.max = 5                      
)


forest(m1.netmeta, ref = "WL", digits = 3, xlab = "SMD")  

###result matrix

treatments <- m1.netmeta$trts  
smd_matrix <- round(m1.netmeta$TE.random, 3)
ci_lower_matrix <- round(m1.netmeta$lower.random, 3)
ci_upper_matrix <- round(m1.netmeta$upper.random, 3)
result.matrix2 <- matrix("", nrow = nrow(smd_matrix), ncol = ncol(smd_matrix), 
                         dimnames = list(treatments, treatments))  
for (i in 1:nrow(smd_matrix)) {
  for (j in 1:ncol(smd_matrix)) {
    result.matrix2[i, j] <- paste(smd_matrix[i, j], "(", ci_lower_matrix[i, j], ",", ci_upper_matrix[i, j], ")", sep = "")
  }
}

temp_direct_matrix2 <- round(m1.netmeta$TE.direct.random, 3)
ci_lower_direct_matrix2 <- round(m1.netmeta$lower.direct.random, 3)
ci_upper_direct_matrix2 <- round(m1.netmeta$upper.direct.random, 3)
direct_result.matrix2 <- matrix("", nrow = nrow(temp_direct_matrix2), ncol = ncol(temp_direct_matrix2), 
                                dimnames = list(treatments, treatments))  
for (i in 1:nrow(temp_direct_matrix2)) {
  for (j in 1:ncol(temp_direct_matrix2)) {
    direct_result.matrix2[i, j] <- paste(temp_direct_matrix2[i, j], "(", ci_lower_direct_matrix2[i, j], ",", ci_upper_direct_matrix2[i, j], ")", sep = "")
  }
}
result.matrix2[lower.tri(result.matrix2, diag = FALSE)] <- direct_result.matrix2[lower.tri(direct_result.matrix2, diag = FALSE)]

# ranking
ranks_for_nma <- netrank(m1.netmeta, small.values = "good")
print(ranks_for_nma$ranking.random)
p_score_common <- ranks_for_nma$ranking.common
p_score_random <- ranks_for_nma$ranking.random
p_score_table <- data.frame(
  Treatment = names(p_score_common),  
  `P-score (common)` = round(p_score_common, 4), 
  `P-score (random)` = round(p_score_random, 4), 
  stringsAsFactors = FALSE,
  check.names = FALSE  
)
p_score_table <- p_score_table %>%
  arrange(desc(`P-score (random)`))
print(p_score_table, row.names = FALSE)

# funnel(netmeta)
treat_order <- m1.netmeta$trts
n_treatments <- length(treat_order)
funnel <- funnel(m1.netmeta,
                 order = treat_order,
                 pch = 19,
                 col = c("black"),
                 linreg = TRUE,
                 xlim = c(-2, 2),
                 ylim = c(0.8, 0),
                 studlab = FALSE,
                 legend = FALSE,
                 cex.studlab = 0.7)



