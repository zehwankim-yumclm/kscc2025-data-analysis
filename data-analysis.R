# 필요한 패키지 설치/로드 --------------------------------------------------------
packages <- c("readxl", "dplyr","ggplot2","forcats","tidyr","broom",
              "stringr","scales","ggrepel","patchwork")
new <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new)) install.packages(new)
lapply(packages, library, character.only = TRUE)

# ------------------------------------------------------------------------------
# 1) 데이터 불러오기 & 정리 (데이터 수집, 저장 불러오기)ntssssssssssssssssssskfxkfxj8 ifsntaldw;eajc 
# ------------------------------------------------------------------------------

df <- read_excel("자료-20250901.xlsx", sheet = 1)

df <- df %>% mutate(
  prevalence = as.numeric(prevalence),
  SE         = as.numeric(SE)
)

# prevalence가 0~1로 변환 (데이터의 정제)
df <- df %>% mutate(prevalence = ifelse(prevalence > 1, prevalence/100, prevalence))

# 범주형 처리  (데이터의 정제)
df$year      <- factor(df$year, levels = c("2019","2020","2021","2022","2023"))
df$sex       <- as.factor(df$sex)
df$age.group <- factor(df$age.group,
                       levels = c("A10","A30","A40","A50","A60","A70"),
                       labels = c("19-29","30-39","40-49","50-59","60-69","70+"))

disease_items   <- c("anemia","CKD","DM","hyperlipidemia","hypertension")
lifestyle_items <- c("alcohol","obese_BMI","smoking")


# 회귀 분석 위해 wide 변환(데이터 변환, 분할)
df_wide <- df %>%
  select(year, sex, age.group, item, prevalence, SE) %>%
  pivot_wider(
    names_from  = item,
    values_from = c(prevalence, SE),
    names_sep   = "_"
  )

# 공용 변환 함수
logit     <- function(p) log(p/(1-p))
inv_logit <- function(x) 1/(1+exp(-x))

# ------------------------------------------------------------------------------
# 2) 기술통계 / EDA & 기본 시각화
# ------------------------------------------------------------------------------

# EDA용: 질병만 추출
df_disease <- df %>% filter(item %in% disease_items)

# logit 기반 95% CI 계산
eps <- 1e-6 # 0 근처 미 1근처의 안정화
df_disease <- df_disease %>%
  mutate(
    p_clip        = pmin(pmax(prevalence, eps), 1 - eps),
    logit_p       = log(p_clip/(1 - p_clip)),
    SE_logit      = SE / (p_clip * (1 - p_clip)),
    ci_low_logit  = inv_logit(logit_p - 1.96 * SE_logit),
    ci_high_logit = inv_logit(logit_p + 1.96 * SE_logit)
  )

# 선 그래프: 연도×연령대×질병 (성별 facet 행, 질병 facet 열)
p_ts <- ggplot(
  df_disease,
  aes(x = year, y = p_clip, color = age.group, group = age.group)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low_logit, ymax = ci_high_logit), width = 0.2, alpha = 0.8) +
  facet_grid(rows = vars(sex), cols = vars(item), scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "연도별·연령대별·질병별 유병률 (±95% CI, logit 기반)",
    x = "조사 연도", y = "유병률",
    color = "연령대"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_ts)
ggsave("Fig 1.png", p_ts, width = 14, height = 8, dpi = 300)

#-------------------------------------------------------------------------------
# 3) WLS 회귀 적합 (질병별: 질병 ~ 생활습관 3개)
# ------------------------------------------------------------------------------

# 로그릿 변환(회귀용): df_wide에 변수 생성
for (d in disease_items) {
  p  <- df_wide[[paste0("prevalence_", d)]]
  se <- df_wide[[paste0("SE_", d)]]
  p_clip <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  df_wide[[paste0("logit_", d)]] <- logit(p_clip)
  se_ok <- ifelse(is.finite(se) & se > 0, se, NA_real_)
  df_wide[[paste0("SE_logit_", d)]] <- se_ok / (p_clip * (1 - p_clip))
}
for (l in lifestyle_items) {
  p <- df_wide[[paste0("prevalence_", l)]]
  p_clip <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  df_wide[[paste0("logit_", l)]] <- logit(p_clip)
}

fit_results <- list()

for (d in disease_items) {
  formula_str <- paste0("logit_", d,
                        " ~ logit_alcohol + logit_obese_BMI + logit_smoking")
  message("Fitting: ", formula_str)
  f <- as.formula(formula_str)

  # 가중치(분산의 역수), 0 분산 방지를 위해 작은 값 더함
  w_raw <- df_wide[[paste0("SE_logit_", d)]]
  w <- 1 / ((w_raw^2) + 1e-12)

  req_cols <- c(paste0("logit_", d),
                "logit_alcohol","logit_obese_BMI","logit_smoking")
  df_model <- df_wide[, req_cols]
  df_model$weights <- w

  keep <- complete.cases(df_model) &
    is.finite(as.matrix(df_model[, req_cols])) %*% rep(1, length(req_cols)) == length(req_cols) &
    is.finite(df_model$weights) & (df_model$weights > 0)

  df_model <- df_model[keep, , drop = FALSE]

  if (nrow(df_model) == 0) {
    warning(sprintf("'%s' 회귀에 사용 가능한 행이 없습니다. (결측/무한값만 존재)", d))
    next
  }

  fit <- lm(f, data = df_model, weights = weights, na.action = na.omit)
  fit_results[[d]] <- tidy(fit) %>% mutate(N_used = nobs(fit))
}

# 결과 테이블 합치기
all_results <- bind_rows(
  lapply(names(fit_results), function(d){
    fit_results[[d]] %>% mutate(Disease = d)
  })
) %>% select(Disease, N_used, term, estimate, std.error, statistic, p.value)

print(all_results)

# ------------------------------------------------------------------------------
# 4) 회귀 결과 시각화 & 보고용 표
# ------------------------------------------------------------------------------

# (Intercept) 제외, CI/OR, 유의성 마킹
res <- all_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = term %>%
      str_replace("^logit_", "") %>%
      factor(levels = c("alcohol","obese_BMI","smoking"),
             labels = c("Alcohol","Obese (BMI)","Smoking")),
    CI_low  = estimate - 1.96*std.error,
    CI_high = estimate + 1.96*std.error,
    OR      = exp(estimate),
    OR_low  = exp(CI_low),
    OR_high = exp(CI_high),
    sig     = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    sign_dir = case_when(
      estimate > 0 & p.value < 0.05 ~ "Positive (p<0.05)",
      estimate < 0 & p.value < 0.05 ~ "Negative (p<0.05)",
      TRUE                          ~ "NS"
    )
  )

# 질환 정렬: 유의계수 많은 순
disease_order <- res %>%
  group_by(Disease) %>% summarise(sig_n = sum(p.value < 0.05)) %>%
  arrange(desc(sig_n)) %>% pull(Disease)

res <- res %>%
  mutate(
    Disease   = factor(Disease, levels = disease_order),
    term_clean = fct_rev(term_clean)
  )

theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# 4-1) Log-odds 포레스트 플롯 (경계선 추가)
p_coef <- ggplot(res, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, alpha = 0.7) +
  geom_point(aes(shape = sign_dir), size = 3) +
  facet_wrap(~ Disease, ncol = 2, scales = "free_x") +
  labs(title = "WLS Coefficients (log-odds scale)",
       x = "Estimate (log-odds)", y = NULL, shape = NULL) +
  theme_pub +
  theme(panel.border = element_rect(color = "black", fill = NA))

# 4-2) 계수 히트맵
res_heat <- res %>%
  group_by(Disease) %>% mutate(est_z = as.numeric(scale(estimate))) %>% ungroup()

p_heat <- ggplot(res_heat, aes(x = Disease, y = term_clean, fill = est_z)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), fontface = "bold") +
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444",
                       midpoint = 0, name = "Std. effect\n(within disease)") +
  labs(title = "Coefficient Heatmap with Significance", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 미리보기 및 저장
print(p_coef); print(p_heat)
ggsave("Fig 2.png", p_coef, width = 10, height = 7, dpi = 300)
ggsave("Fig 3.png", p_heat, width = 8,  height = 5, dpi = 300)

