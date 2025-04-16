

setwd("D:/OneDrive - 北京大学/PTSD治疗的高级Meta分析/component NMA/牧阳的尝试")
setwd("C:/Users/admin/OneDrive - 北京大学/PTSD治疗的高级Meta分析/component NMA/牧阳的尝试")
setwd("C:/Users/3/OneDrive - 北京大学/PTSD治疗的高级Meta分析/component NMA")

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

source("functions.R", encoding = "UTF-8")

## ????data
cNMA_raw_data <- readxl::read_xlsx("./data for trial 20250416_2.xlsx")
cNMA_raw_data <- cNMA_raw_data %>%
  mutate(na = na) %>%
  mutate(Ns = nrow(.))

cNMA_raw_data$na <- cNMA_raw_data$na %>%
  as.numeric()


# define useful variables
age =  cNMA_raw_data %>%
  separate(age, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
sex =  cNMA_raw_data %>%
  separate(sex, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
duration =  cNMA_raw_data %>%
  separate(duration, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)

cNMA_raw_data <- cNMA_raw_data %>%
  mutate(trauma = sapply(1:nrow(.), function(i) {
    paste(rep(trauma[i], na[i]), collapse = ",")
  }))
trauma = cNMA_raw_data %>%
  separate(trauma, c("a1","a2","a3","a4"), sep = ",", convert = T) %>%
  select(a1,a2,a3,a4)

cNMA_raw_data <- cNMA_raw_data %>%
  mutate(comorbidity = sapply(1:nrow(.), function(i) {
    paste(rep(comorbidity[i], na[i]), collapse = ",")
  }))
comorbidity = cNMA_raw_data %>%
  separate(comorbidity, c("a1","a2","a3","a4"), sep = ",", convert = T) %>%
  select(a1,a2,a3,a4)


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
#dropout_rate = cNMA_raw_data %>%
separate(dropout_rate, c("dropout_rate1","dropout_rate2","dropout_rate3","dropout_rate4"),sep = ";",convert = T) %>%
  select(dropout_rate1,dropout_rate2,dropout_rate3,dropout_rate4)
c1 = cNMA_raw_data %>%
  separate(IVE, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c2 = cNMA_raw_data %>%
  separate(reliving, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c3 = cNMA_raw_data %>%
  separate(CR, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c4 = cNMA_raw_data %>%
  separate(BL,c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c5 = cNMA_raw_data %>%
  separate(NF, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c6 = cNMA_raw_data %>%
  separate(IMAG, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c7 = cNMA_raw_data %>%
  separate(PS, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c8 = cNMA_raw_data %>%
  separate(BA, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c9 = cNMA_raw_data %>%
  separate(RL, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c10 = cNMA_raw_data %>%
  separate(HW, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c11 = cNMA_raw_data %>%
  separate(other, c("a1","a2","a3","a4"),sep = ",",convert = T) %>%
  select(a1,a2,a3,a4)
c12 = cNMA_raw_data %>%
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

#dropout_integ <- bind_cols(cNMA_raw_data %>%
separate(dropout_rate_a1, c("n_dropout.a1","n_recruit.a1"),sep = "/",convert = T) %>%
  select(n_dropout.a1,n_recruit.a1),
cNMA_raw_data %>%
  separate(dropout_rate_a2, c("n_dropout.a2","n_recruit.a2"),sep = "/",convert = T) %>%
  select(n_dropout.a2,n_recruit.a2),
cNMA_raw_data %>%
  separate(dropout_rate_a3, c("n_dropout.a3","n_recruit.a3"),sep = "/",convert = T) %>%
  select(n_dropout.a3,n_recruit.a3),
cNMA_raw_data %>%
  separate(dropout_rate_a4, c("n_dropout.a4","n_recruit.a4"),sep = "/",convert = T) %>%
  select(n_dropout.a4,n_recruit.a4))

#n_recruit <- dropout_integ %>%
select(starts_with("n_recruit."))

#n_dropout <- dropout_integ %>%
select(starts_with("n_dropout."))

Nc=12
Ns=nrow(y)

save(y=y,n=n,sd=sd,Nc,Ns,na=na,
     c1=c1, c2=c2,
     c3=c3,c4=c4, c5=c5, c6=c6, c7=c7,c8=c8,
     c9=c9, c10=c10, c11=c11, age=age, sex=sex, duration=duration, trauma=trauma, comorbidity = comorbidity,
     s_within,
     s_m_within,
     file = "cNMA_data_20250416_2.RData")

## load libraries ----------------------------------------------------------
library(netmeta)
library(coda)
library(rjags)
library(MCMCvis)
library(tidyverse)
source("functions.R", encoding = "UTF-8")

load("cNMA_data_20250328_10.RData")

# define interactions -----------------------------------------------------
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

# 元素矩阵列表
components <- list(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12)

# Generate all possible pairs of components
pairs <- combn(1:12, 2)

# Check which pairs exist (i.e., both components are 1 simultaneously in at least one observation)
existing_pairs <- apply(pairs, 2, function(pair) {
  i <- pair[1]
  j <- pair[2]
  any(components[[i]] * components[[j]] == 1)  # Check if the product is 1 anywhere
})

# Filter to existing pairs
existing_pairs_idx <- which(existing_pairs)
existing_pairs <- pairs[, existing_pairs_idx, drop = FALSE]
Ninter <- ncol(existing_pairs)

# Define maximum number of arms
max_na <- max(na)

# Generate interaction array for all existing pairs
interactions <- array(NA, dim = c(Ns, max_na, Ninter))
for (k in 1:Ninter) {
  i <- existing_pairs[1, k]
  j <- existing_pairs[2, k]
  interactions[, , k] <- components[[i]] * components[[j]]
}

# Calculate frequency of each interaction
interaction_freq <- numeric(Ninter)
for (k in 1:Ninter) {
  interaction_freq[k] <- sum(interactions[, , k] == 1, na.rm = TRUE)
}

# Create a data frame with interaction names and frequencies
interaction_names <- sapply(1:Ninter, function(k) {
  i <- existing_pairs[1, k]
  j <- existing_pairs[2, k]
  paste("c", i, "-c", j, sep = "")
})

results <- data.frame(
  Interaction = interaction_names,
  Frequency = interaction_freq
)

# Print results
print(paste("Total number of interaction terms (without restrictions):", Ninter))
print("Interaction terms and their frequencies:")
print(results)

# Optionally, save results to a CSV file
write.csv(results, "interaction_terms_frequencies_20250329.csv", row.names = FALSE)

# 明确定义仅有的交互项 c2-c4
existing_pairs_reduced <- matrix(c(2, 4), ncol = 1)  # 只保留 c2-c4
Ninter_reduced <- 1  # 只有一个交互项

# 定义最大臂数
max_na <- max(na)

# 生成交互项数组（仅 c2-c4）
interactions_reduced <- array(NA, dim = c(Ns, max_na, Ninter_reduced))
interactions_reduced[, , 1] <- components[[2]] * components[[4]]  # c2 * c4

# 计算交互项频率
interaction_freq <- sum(interactions_reduced[, , 1] == 1, na.rm = TRUE)

# 打印交互项数量、内容及其频率
print(paste("Number of interaction terms included:", Ninter_reduced))
print(paste("Interaction 1: c2-c4 | Frequency:", interaction_freq))

# 定义参数
Nc <- 11  # 总组件数（去掉 c4 主效应后剩 11 个主效应）

# 更新数据列表（去掉 age）
data <- list(
  y = y, n = n, sd = sd,
  Nc = 11, Ns = Ns, na = na,
  c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, 
  c6 = c6, c7 = c7, c8 = c8, c9 = c9, c10 = c10,
  c11 = c11, c12 = c12,
  s_within = s_within, Ninter = Ninter_reduced, interactions = interactions_reduced
)

# 定义模型（去掉 c4 主效应，仅保留 c2-c4 交互项）
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

    for (k in 2:na[i]) {
      delta[i,k] ~ dnorm(md[i,k], precd[i,k])   
      md[i,k] <- mean[i,k] + sw[i,k]                 
      precd[i,k] <- 2 * prect * (k-1) / k                            
      w[i,k] <- (delta[i,k] - mean[i,k])      
      sw[i,k] <- sum(w[i,1:(k-1)]) / (k-1)        

      mean[i,k] <- A1[i,k] - B1[i]
      A1[i,k] <- 
        d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0)) +           
        d[3]*(1-equals(c3[i,k],0)) + 
        d[4]*(1-equals(c5[i,k],0)) + d[5]*(1-equals(c6[i,k],0)) + 
        d[6]*(1-equals(c7[i,k],0)) + d[7]*(1-equals(c8[i,k],0)) + 
        d[8]*(1-equals(c9[i,k],0)) + d[9]*(1-equals(c10[i,k],0)) +
        d[10]*(1-equals(c11[i,k],0)) + d[11]*(1-equals(c12[i,k],0)) +
        gamma[1] * interactions[i,k,1]  # 仅 c2-c4 交互项
    }                                        

    B1[i] <- 
      d[1]*(1-equals(c1[i,1],0)) + d[2]*(1-equals(c2[i,1],0)) +
      d[3]*(1-equals(c3[i,1],0)) + 
      d[4]*(1-equals(c5[i,1],0)) + d[5]*(1-equals(c6[i,1],0)) + 
      d[6]*(1-equals(c7[i,1],0)) + d[7]*(1-equals(c8[i,1],0)) + 
      d[8]*(1-equals(c9[i,1],0)) + d[9]*(1-equals(c10[i,1],0)) + 
      d[10]*(1-equals(c11[i,1],0)) + d[11]*(1-equals(c12[i,1],0)) +
      gamma[1] * interactions[i,1,1]  # 仅 c2-c4 交互项
  }
  for (i in 1:Ns) { u[i] ~ dnorm(0, .001) }
  tau ~ dnorm(0, 0.1)I(0,)                                      
  prect <- 1/tau.sq
  tau.sq <- pow(tau, 2)
  gamma[1] ~ ddexp(0, lambda)  # 仅一个 gamma 参数
  lambda ~ dgamma(0.01, 0.01)
  for(k in 1:Nc) { d[k] ~ dnorm(0, .001) }  # 定义 d[1:11]
}
"

# 运行模型
model.ptsd.spec <- textConnection(model.ptsd.string)
jags.m.ptsd <- jags.model(model.ptsd.spec, data = data, n.chains = 4, n.adapt = 10000)
params <- c("tau", "d", "gamma")
samps.ptsd <- coda.samples(jags.m.ptsd, params, n.iter = 20000, thin = 10)

# 输出结果
summary.ptsd <- MCMCsummary(samps.ptsd, Rhat = TRUE, n.eff = TRUE)

# 为 d 参数重命名（11 个主效应，去掉 c4）
rownames(summary.ptsd)[1:11] <- c("IVE", "reliving", "CR", "NF", "IMAG", "PS", "BA", "RL", "HW", "other", "waiting")

# 为 gamma 参数命名（仅 c2-c4）
gamma_summary <- summary.ptsd[grep("gamma", rownames(summary.ptsd)), , drop = FALSE]
rownames(gamma_summary) <- "c2-c4"

# 打印结果
print(round(summary.ptsd, digits = 3))
print(round(gamma_summary, digits = 3))

MCMCtrace(samps.ptsd, pdf = FALSE, params = "d", wd = getwd(), n.eff = TRUE, Rhat = TRUE)
write.csv(round(summary.ptsd, digits = 3), "./outputs_20250416_2.csv")

# 计算 DIC
DIC.ptsd <- dic.samples(jags.m.ptsd, n.iter = 5000)
print(DIC.ptsd)
closeAllConnections()
