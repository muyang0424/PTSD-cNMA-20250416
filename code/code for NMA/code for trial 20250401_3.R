# 加载必要的库
library(readxl)
library(tidyverse)
library(netmeta)
library(ggplot2)
library(tidyr)
library(openxlsx)

# 设置工作目录（根据需要调整）
#setwd("D:/OneDrive - 北京大学/PTSD治疗的高级Meta分析/network meta analysis")
setwd("C:/Users/3/OneDrive - 北京大学/PTSD治疗的高级Meta分析/network meta analysis")
source("functions.R", encoding = "UTF-8")


# 读取 Excel 数据
Armlevel_data <- readxl::read_xlsx("./data for trial 20250416_3 只有duration数据全的.xlsx")

# 将数据转换为长格式
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

# 输出结果
print(m1.netmeta)

# 不一致性检验 (netsplit)
netsplit_result <- netsplit(m1.netmeta)
print(netsplit_result)

# 全局不一致性检验
decomp.design(m1.netmeta)




# 网络图
netgraph(
  m1.netmeta,
  plastic = FALSE,                 # 保持布局简洁
  thickness = "number.of.studies", # 线的粗细基于研究数量
  number.of.studies = FALSE,       # 不在边上显示研究数量
  points = TRUE,                   # 显示治疗节点
  cex = 0.8,                       # 缩小节点标签字体大小 (可以调整 0.8 这个值)
  cex.points = 2,                  # 保持节点本身的大小 (可以根据需要调整)
  col = "black",                   # 边的颜色保持为黑色
  col.points = "black",            # 将点的颜色设置为黑色
  lwd.max = 5                      # 最大线宽
)


# 检查网络连接性
net_connect <- netconnection(
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = p1
)

# 获取所有治疗名称
treatments <- unique(c(p1$treat1, p1$treat2))
treatments

# 1. 检查哪些治疗有非活性对照比较
non_active_controls <- c("WL", "TAU", "SCC", "sham-neurofeedback")  # 定义非活性对照
has_non_active <- sapply(treatments, function(t) {
  any((p1$treat1 == t & p1$treat2 %in% non_active_controls) | 
        (p1$treat2 == t & p1$treat1 %in% non_active_controls))
})
no_non_active <- treatments[!has_non_active]

# 2. 计算等待名单 (WL) 的比较次数
wl_comparisons <- p1 %>%
  filter(treat1 == "WL" | treat2 == "WL") %>%
  summarise(unique_treatments = length(unique(c(treat1, treat2))) - 1)  # 减去 WL 本身

# 计算每种治疗的连接数（度数）
connections <- p1 %>%
  select(treat1, treat2) %>%
  # 为每个治疗生成它连接到的其他治疗
  group_by(treat1) %>%
  summarise(connected = list(unique(treat2))) %>%
  rename(treat = treat1) %>%
  bind_rows(
    p1 %>%
      group_by(treat2) %>%
      summarise(connected = list(unique(treat1))) %>%
      rename(treat = treat2)
  ) %>%
  group_by(treat) %>%
  summarise(n_connections = n_distinct(unlist(connected))) %>%
  filter(!is.na(treat))

# 找出只通过一条边连接的治疗（度数 = 1）
single_connection <- connections %>%
  filter(n_connections == 1) %>%
  pull(treat)

# 输出结果
if (length(single_connection) == 0) {
  cat("No psychotherapies were connected via a single condition comparison.\n\n")
} else {
  cat(length(single_connection), "psychotherapies (", 
      paste(single_connection, collapse = ", "), 
      ") were connected via a single condition comparison.\n\n")
}


# 4. 计算网络连接性：实际比较数 vs 可能比较数

# 计算可能比较数
n_treatments <- length(treatments)  # Dynamically set to the number of unique treatments
possible_comparisons <- n_treatments * (n_treatments - 1) / 2  

# 从 p1 中提取独特直接比较对数
unique_pairs <- p1 %>%
  mutate(pair = paste(pmin(treat1, treat2), pmax(treat1, treat2), sep = " vs ")) %>%
  distinct(pair) %>%
  nrow()

# 输出结果
cat("网络中共有", unique_pairs, "对独特直接比较，占可能比较总数", possible_comparisons, "对。\n")

# 输出结果
cat("All therapies, except", paste(no_non_active, collapse = ", "), 
    ", had at least one non-active control comparison.\n\n")
cat("Waiting list was the most compared control condition, being compared to", 
    wl_comparisons$unique_treatments, "psychotherapies and control conditions.\n\n")
cat(length(single_connection), "psychotherapies (", 
    paste(single_connection, collapse = ", "), 
    ") were connected via a single condition comparison.\n\n")

# 结果矩阵（随机效应）
result.matrix <- m1.netmeta$TE.random
result.matrix <- round(result.matrix, 3)
temp_direct_matrix <- round(m1.netmeta$TE.direct.random, 3)
result.matrix[lower.tri(result.matrix, diag = FALSE)] <- temp_direct_matrix[lower.tri(temp_direct_matrix, diag = FALSE)]
write.csv(result.matrix, "result_matrix_20250416.csv", row.names = TRUE)


# 森林图
forest(m1.netmeta, ref = "WL", digits = 3, xlab = "SMD")  # 使用 SMD


# 获取治疗名称的正确顺序
treatments <- m1.netmeta$trts  # 使用 netmeta 对象的治疗顺序

# 结果矩阵（带 95% CI）
smd_matrix <- round(m1.netmeta$TE.random, 3)
ci_lower_matrix <- round(m1.netmeta$lower.random, 3)
ci_upper_matrix <- round(m1.netmeta$upper.random, 3)

# 创建结果矩阵，并设置行名和列名
result.matrix2 <- matrix("", nrow = nrow(smd_matrix), ncol = ncol(smd_matrix), 
                         dimnames = list(treatments, treatments))  # 使用 m1.netmeta$trts
for (i in 1:nrow(smd_matrix)) {
  for (j in 1:ncol(smd_matrix)) {
    result.matrix2[i, j] <- paste(smd_matrix[i, j], "(", ci_lower_matrix[i, j], ",", ci_upper_matrix[i, j], ")", sep = "")
  }
}

# 直接效应矩阵
temp_direct_matrix2 <- round(m1.netmeta$TE.direct.random, 3)
ci_lower_direct_matrix2 <- round(m1.netmeta$lower.direct.random, 3)
ci_upper_direct_matrix2 <- round(m1.netmeta$upper.direct.random, 3)

# 创建直接效应矩阵，并设置行名和列名
direct_result.matrix2 <- matrix("", nrow = nrow(temp_direct_matrix2), ncol = ncol(temp_direct_matrix2), 
                                dimnames = list(treatments, treatments))  # 使用 m1.netmeta$trts
for (i in 1:nrow(temp_direct_matrix2)) {
  for (j in 1:ncol(temp_direct_matrix2)) {
    direct_result.matrix2[i, j] <- paste(temp_direct_matrix2[i, j], "(", ci_lower_direct_matrix2[i, j], ",", ci_upper_direct_matrix2[i, j], ")", sep = "")
  }
}

# 将直接效应填入下三角
result.matrix2[lower.tri(result.matrix2, diag = FALSE)] <- direct_result.matrix2[lower.tri(direct_result.matrix2, diag = FALSE)]

# 保存到 CSV 文件
write.csv(result.matrix2, "result_matrix_20250416_2.csv", row.names = TRUE, fileEncoding = "UTF-8")

# 排序
ranks_for_nma <- netrank(m1.netmeta, small.values = "good")
print(ranks_for_nma$ranking.random)

# 提取 P-score
p_score_common <- ranks_for_nma$ranking.common
p_score_random <- ranks_for_nma$ranking.random

# 创建数据框，确保列名正确
p_score_table <- data.frame(
  Treatment = names(p_score_common),  # 治疗名称
  `P-score (common)` = round(p_score_common, 4),  # 固定效应模型 P-score
  `P-score (random)` = round(p_score_random, 4),  # 随机效应模型 P-score
  stringsAsFactors = FALSE,
  check.names = FALSE  # 防止 R 自动修改列名
)

# 按 P-score (random) 降序排列
p_score_table <- p_score_table %>%
  arrange(desc(`P-score (random)`))

# 打印表格
print(p_score_table, row.names = FALSE)

# 保存为 Excel 文件
write.xlsx(p_score_table, "P_score_table_20250416.xlsx", rowNames = FALSE)

#############################
# 1. 漏斗图 - funnel(netmeta)


#############################

treat_order <- m1.netmeta$trts
n_treatments <- length(treat_order)



funnel <- funnel(m1.netmeta,
                 order = treat_order,#˳???Լ??????????԰ѱȽ??µĸ?Ԥ??????
                 pch = 19,
                 col = c("black"),
                 linreg = TRUE,
                 xlim = c(-2, 2),
                 ylim = c(0.8, 0),
                 studlab = FALSE,
                 legend = FALSE,
                 cex.studlab = 0.7)




###检查transitivity
library(dplyr)
library(tidyr)
library(ggplot2)



# 分割 age、sex、duration 列并移除原始列
Armlevel_data_split <- Armlevel_data %>%
  separate(age, into = c("age_1", "age_2", "age_3", "age_4"), sep = ",", fill = "right", remove = TRUE) %>%
  separate(sex, into = c("sex_1", "sex_2", "sex_3", "sex_4"), sep = ",", fill = "right", remove = TRUE) %>%
  separate(duration, into = c("duration_1", "duration_2", "duration_3", "duration_4"), sep = ",", fill = "right", remove = TRUE)

arm_long_corrected <- Armlevel_data_split %>%
  pivot_longer(
    cols = matches("Arm_[1-4]"),
    names_to = "arm",
    names_prefix = "Arm_",
    values_to = "treatment",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = matches("age_[1-4]"),
    names_to = "age_arm",
    names_prefix = "age_",
    values_to = "age",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = matches("sex_[1-4]"),
    names_to = "sex_arm",
    names_prefix = "sex_",
    values_to = "sex",
    values_drop_na = TRUE
  ) %>%
  pivot_longer(
    cols = matches("duration_[1-4]"),
    names_to = "duration_arm",
    names_prefix = "duration_",
    values_to = "duration",
    values_drop_na = TRUE
  ) %>%
  # 确保治疗臂与对应的 age、sex、duration 对齐
  filter(arm == age_arm, arm == sex_arm, arm == duration_arm) %>%
  mutate(
    age = as.numeric(age),
    sex = as.numeric(sex),
    duration = as.numeric(duration)
  ) %>%
  rename(study = `NO`) %>%
  select(study, treatment, age, sex, duration)


# 计算每个研究的平均协变量
study_means <- arm_long_corrected %>%
  group_by(study) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    mean_sex = mean(sex, na.rm = TRUE)
  )

# 生成比较数据
comparison_data <- arm_long_corrected %>%
  group_by(study) %>%
  summarise(
    treatments = list(treatment),
    n_arms = n()
  ) %>%
  filter(n_arms >= 2) %>%  # 只保留有至少 2 个臂的研究
  mutate(
    treat_pairs = map(treatments, ~ combn(.x, 2, simplify = FALSE))
  ) %>%
  unnest(treat_pairs) %>%
  mutate(
    treat1 = map_chr(treat_pairs, 1),
    treat2 = map_chr(treat_pairs, 2),
    comparison = paste(treat1, "vs", treat2)
  ) %>%
  select(study, treat1, treat2, comparison) %>%
  left_join(study_means, by = "study")

# 检查数据
head(comparison_data)

# 绘制森林图，调整为白底、黑字、黑点、灰线
p <- ggplot(comparison_data, aes(x = mean_age, y = comparison)) +
  geom_jitter(width = 0, height = 0.2, size = 2, color = "black") +  # 黑点
  labs(x = "Age", y = "Comparison") +
  theme(
    # 白底
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 黑字
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),  # 纵轴标签黑色，字体大小 6
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    # 灰线（网格线）
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  )

# 保存为大尺寸图片
ggsave("forest_plot_age_large.png", p, width = 10, height = 25, units = "in", dpi = 300)

# sex
p <- ggplot(comparison_data, aes(x = mean_sex, y = comparison)) +
  geom_jitter(width = 0, height = 0.2, size = 2, color = "black") +  # 黑点
  labs(x = "female percent", y = "Comparison", title = "按研究平均性别分布的森林图") +
  theme(
    # 白底
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 黑字
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),  # 纵轴标签黑色，字体大小 6
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    # 灰线（网格线）
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  )

# 保存为大尺寸图片
ggsave("forest_plot_sex_large.png", p, width = 10, height = 25, units = "in", dpi = 300)



# comorbidity
# --- 前面的数据准备代码保持不变 ---
# 1. 提取 study 和 comorbidity
study_comorbidity <- Armlevel_data %>%
  rename(study = `NO`) %>% # 确保这里的列名 `NO...2` 是正确的
  select(study, comorbidity)

# 2. 将 comorbidity 合并到 comparison_data
comparison_data_comorbidity <- comparison_data %>%
   left_join(study_comorbidity, by = "study")

# --- 修改绘图代码 ---

# 3. 定义 x 轴的 breaks (原始数值) 和 labels (对应的英文标签)
x_axis_breaks <- c(0, 1, 2, 3, 4, 5)
x_axis_labels <- c(
  "0" = "None",
  "1" = "Affect Disorders",
  "2" = "Sleep Problems",
  "3" = "Substance Use Problems",
  "4" = "Severe Mental Illness",
  "5" = "Somatic Diseases"
)

# 确保标签顺序与 breaks 对应 (如果 breaks 不是严格按 0,1,2... 顺序的话，这种索引方式更安全)
correct_labels <- x_axis_labels[as.character(x_axis_breaks)]

# 绘制图形
p_comorbidity <- ggplot(comparison_data_comorbidity, aes(x = comorbidity, y = comparison)) +
  geom_jitter(width = 0, height = 0.2, size = 2, color = "black") +
  # 添加 scale_x_continuous 来指定新的标签
  scale_x_continuous(
    breaks = x_axis_breaks,
    labels = correct_labels
  ) +
  labs(
    x = "Comorbidity",
    y = "Comparison",
    # 稍微修改标题使其更符合学术习惯
    title = "Distribution of Comparisons by Study Comorbidity"
  ) +
  theme(
    # 白底
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 黑字
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    # 确保 x 轴标签也是黑色，如果标签太长可能需要调整角度
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1), # 例如旋转45度
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5), # 标题居中
    # 灰线（网格线），注意 size 已被弃用，推荐用 linewidth
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.25)
  )

# 打印或查看图形（可选）
 print(p_comorbidity)

# 4. 保存图片
ggsave(
  filename = "forest_plot_comorbidity_large_labeled.png", # 建议改个文件名
  plot     = p_comorbidity,
  width    = 12, # 可能需要加宽图幅以容纳长标签
  height   = 25,
  units    = "in",
  dpi      = 300
)



# --- 检查并执行数据准备代码 ---
# 1. 提取 study 和 trauma
# !! 请确认 Armlevel_data 中代表 study 的列名是 'NO' 还是其他名称 !!
study_trauma <- Armlevel_data %>%
  rename(study = NO) %>%  # 确认 'NO' 是正确的列名
  select(study, trauma)

# 2. 将 trauma 合并到 comparison_data
comparison_data_trauma <- comparison_data %>%
  left_join(study_trauma, by = "study")

# --- 修改绘图代码 ---

# 3. 定义 x 轴的 breaks (原始数值) 和 labels (对应的英文标签)
trauma_axis_breaks <- c(0, 1, 2, 3, 4, 5, 6)
trauma_axis_labels <- c(
  "0" = "Mixed or Unknown",
  "1" = "Motor Vehicle Accident", # 使用 MVC 的全称，更清晰
  "2" = "Somatic Diseases",              # 将“躯体疾病”翻译为英文
  "3" = "Bereaved",
  "4" = "Childhood Abuse",
  "5" = "Interpersonal Trauma (including Sexual Assault)", # 稍微调整措辞
  "6" = "War, Terror, Political Victimization"
)

# 确保标签顺序与 breaks 对应
correct_trauma_labels <- trauma_axis_labels[as.character(trauma_axis_breaks)]

# 绘制图形
p_trauma <- ggplot(comparison_data_trauma, aes(x = trauma, y = comparison)) +
  geom_jitter(width = 0, height = 0.2, size = 2, color = "black") +
  # 添加 scale_x_continuous 来指定新的标签
  scale_x_continuous(
    breaks = trauma_axis_breaks,
    labels = correct_trauma_labels
  ) +
  labs(
    x = "Trauma Type", # 修改 x 轴标题
    y = "Comparison",
    title = "Distribution of Comparisons by Study Trauma Type" # 修改图表标题
  ) +
  theme(
    # 白底
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 黑字
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    # 旋转 x 轴标签以防重叠
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.title = element_text(color = "black"),
    # 居中标题
    plot.title = element_text(color = "black", hjust = 0.5),
    # 灰线（网格线）, 使用 linewidth 替代 size
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.25)
  )

# 打印或查看图形（可选）
# print(p_trauma)

# 4. 保存图片
ggsave(
  filename = "forest_plot_trauma_large_labeled.png", # 建议改个文件名
  plot     = p_trauma,
  width    = 12,  # 可能需要调整宽度以适应标签
  height   = 25,
  units    = "in",
  dpi      = 300
)


# 基于每条治疗臂绘制 duration 的分布图
duration_summary <- arm_long_corrected %>%
  group_by(treatment) %>%
  summarise(
    mean_duration = mean(duration, na.rm = TRUE),
    sd_duration   = sd(duration, na.rm = TRUE),
    n            = n()
  ) %>%
  # 如果想要在图中画出 mean ± SD 的区间，可以预先算好区间范围
  mutate(
    lower = mean_duration - sd_duration,
    upper = mean_duration + sd_duration
  ) %>%
  # 如果想让输出按照 mean_duration 大小排序，方便在图中美观显示
  arrange(desc(mean_duration))

p_duration_forest <- ggplot(duration_summary, 
                            aes(x = mean_duration, 
                                y = reorder(treatment, mean_duration))) +
  geom_point(size = 2, color = "black") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 height = 0.2, 
                 color = "gray40") +
  labs(
    x = "Duration (mean ± SD)",
    y = "Treatment",
    title = "各治疗臂的病程（Duration）分布"
  ) +
  theme(
    # 白底
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 黑字
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    # 灰色网格线
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  )

# 如果需要保存为较大尺寸的图片
ggsave("forest_plot_duration_by_treatment.png", 
       p_duration_forest, 
       width = 10, height = 10, 
       units = "in", dpi = 300)


