# ============================================================================
# DIFFERENCE-IN-DIFFERENCES
# ============================================================================
#
# Run interactively with the repository root as the working directory:
#   source("did_teaching_data/workflow_did_lecture.R")
#
# The case blocks repeat code so each can be run independently. Each block
# reads its own CSV and presents the case's diagnostics, DiD estimates, and
# relevant comparisons. Comment out blocks not used in class, or copy one
# block into a scratch file.
# ============================================================================

## Global settings
FIG_WIDTH <- 8
FIG_HEIGHT <- 6

required_packages <- c(
  "dplyr",
  "fixest",
  "forcats",
  "ggfixest",
  "ggplot2",
  "did",
  "panelView",
  "patchwork",
  "WeightIt",
  "broom",
  "cobalt",
  "bacondecomp",
  "marginaleffects"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}
# synthdid is not on CRAN; install it from GitHub when missing.
if (!requireNamespace("synthdid", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("synth-inference/synthdid")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(forcats)
  library(ggplot2)
  library(did)
  library(panelView)
  library(patchwork)
  library(WeightIt)
  library(synthdid)
  library(cobalt)
  library(bacondecomp)
  library(marginaleffects)
})

options(stringsAsFactors = FALSE, scipen = 6)

# Set the data path from the interactive working directory.
if (dir.exists(file.path("did_teaching_data", "data"))) {
  data_dir <- file.path("did_teaching_data", "data")
} else {
  data_dir <- "data"
}
stopifnot(dir.exists(data_dir))

# ----------------------------------------------------------------------------
# CASE 01. Common adoption with treated and never-treated groups
# What does the treated-control change measure here?
# Compare the group means and regression estimates, then state the assumptions
# needed for a causal interpretation.
#
# Case introduction: Loyalty offer
#   A retail chain introduces a loyalty offer in some stores during week 7.
#   The data are a store-week panel, and `outcome` is log weekly revenue.
#   `treated` records whether the store is receiving the offer in that week.
#
#   Compare treated and untreated stores over time and estimate the effect of
#   the offer.
# ----------------------------------------------------------------------------

# True effect:
#   ATT = 0.30 log points in every treated store-week.

case_file_01 <- file.path(data_dir, "case_01.csv")
d01 <- read.csv(case_file_01, stringsAsFactors = FALSE)
stopifnot(
  nrow(d01) > 0L,
  !anyDuplicated(paste(d01$unit_id, d01$period, sep = "::")),
  all(d01$treated %in% c(0, 1)),
  all(d01$ever_treated %in% c(0, 1)),
  all(d01$adoption_period[d01$ever_treated == 1] == 7)
)

first_treatment_period_01 <- min(
  d01$adoption_period[d01$ever_treated == 1],
  na.rm = TRUE
)
first_treatment_period_01

# 1. Panelview: inspect treatment timing and the observed panel.
panelview(
  outcome ~ treated,
  data = d01,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 01: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Group averages and the adoption date.
d01_plot_data <- d01 |>
  mutate(plot_group = if_else(ever_treated == 1, "Treated", "Control"))
d01_group_means <- d01_plot_data |>
  summarise(
    .by = c(period, plot_group),
    mean_outcome = mean(outcome)
  )
d01_group_means
d01_group_plot <- ggplot(
  d01_group_means,
  aes(
    x = period,
    y = mean_outcome,
    color = plot_group,
    shape = plot_group,
    linetype = plot_group
  )
) +
  geom_vline(
    xintercept = first_treatment_period_01,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 01: treated and never-treated store means",
    x = "Period",
    y = "Mean outcome",
    color = "Group",
    shape = "Group",
    linetype = "Group"
  ) +
  theme_minimal()
d01_group_plot
ggsave(
  "standard_did.pdf",
  plot = d01_group_plot,
  path = "figures",
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  units = "in"
)

# 3. Baseline TWFE DiD.
twfe_fit_01 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d01
)
summary(twfe_fit_01)
twfe_table_01 <- broom::tidy(twfe_fit_01)

# 4. TWFE event study, joint Wald test of pre-period coefficients, and plot.
event_fit_01 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_01 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d01
)
summary(event_fit_01)
wald_pre_01 <- wald(event_fit_01, keep = "period::[1-5]:", print = TRUE)
print(wald_pre_01)
iplot(
  event_fit_01,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 01: TWFE event study"
)
# OR
ggfixest::ggiplot(event_fit_01) +
  labs(
    title = "Case 01: TWFE event study",
    x = "Calendar period; week 6 is the reference"
  )
# 5. Pre-treatment placebo diagnostic: use a date three weeks before adoption.
placebo_period_01 <- first_treatment_period_01 - 3L
d01_placebo <- d01 |>
  filter(period < first_treatment_period_01) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_01
    )
  )
placebo_fit_01 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d01_placebo
)
summary(placebo_fit_01)
etable(twfe_fit_01, placebo_fit_01)

feols(outcome ~ period:ever_treated | unit_id + period, data = d01_placebo)

# 6. Compare the DiD estimates with alternative estimators.
d01_cs <- d01 |>
  mutate(
    id_num = as.integer(as.factor(unit_id)),
    cohort = if_else(is.na(adoption_period), 0, adoption_period)
  )
cs_fit_01 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~1,
  data = d01_cs,
  est_method = "dr",
  control_group = "nevertreated",
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE,
  base_period = "varying", #"universal"
)
cs_fit_01
cs_simple_01 <- aggte(cs_fit_01, type = "simple")
cs_dynamic_01 <- aggte(
  cs_fit_01,
  type = "dynamic",
  bstrap = FALSE,
  cband = FALSE
)
summary(cs_simple_01)
ggdid(cs_dynamic_01, title = "Case 01: CS dynamic aggregation") /
  ggdid(cs_fit_01, title = "Case 01: CS dynamic aggregation")

d01_synthdid_panel <- panel.matrices(
  d01,
  unit = "unit_id",
  time = "period",
  outcome = "outcome",
  treatment = "treated"
)
synthdid_fit_01 <- with(
  d01_synthdid_panel,
  synthdid_estimate(Y, N0, T0)
)
synthdid_se_01 <- sqrt(vcov(synthdid_fit_01, method = "jackknife"))
synthdid_confint_01 <- c(
  as.numeric(synthdid_fit_01) - qnorm(0.975) * synthdid_se_01,
  as.numeric(synthdid_fit_01) + qnorm(0.975) * synthdid_se_01
)
print(synthdid_confint_01)
plot(synthdid_fit_01, overlay = 1) +
  ggtitle("Case 01: Synthetic DiD")

comparison_01 <- data.frame(
  specification = c(
    "TWFE",
    "Callaway-Sant'Anna aggregate",
    "Synthetic DiD"
  ),
  estimate = c(
    twfe_table_01$estimate[twfe_table_01$term == "treated"],
    unname(cs_simple_01$overall.att),
    as.numeric(synthdid_fit_01)
  ),
  standard_error = c(
    twfe_table_01$std.error[twfe_table_01$term == "treated"],
    unname(cs_simple_01$overall.se),
    synthdid_se_01
  ),
  row.names = NULL
)
comparison_01

# ----------------------------------------------------------------------------
# CASE 02. Observed covariates and outcome trends
# How do observed covariates relate to treated and comparison outcomes?
# Compare estimates with and without adjustment; assess overlap and the
# outcome paths.
#
# Case introduction: CRM adoption
#   Regional retailers adopt a customer-relationship-management system during
#   week 7. The data are a retailer-week panel, with log weekly revenue as the
#   outcome. `growth_segment` is observed before adoption.
#
#   Assess the comparability of treated and untreated retailers before
#   estimating the CRM effect.
# ----------------------------------------------------------------------------

# True effect:
#   ATT = 0.25 log points.
#   Untreated weekly slopes differ by growth segment (low -0.02, middle 0.01,
#   high 0.04); within a segment, trends are parallel.

case_file_02 <- file.path(data_dir, "case_02.csv")
d02 <- read.csv(case_file_02, stringsAsFactors = FALSE)
stopifnot(
  nrow(d02) > 0L,
  !anyDuplicated(paste(d02$unit_id, d02$period, sep = "::")),
  all(d02$treated %in% c(0, 1)),
  all(d02$ever_treated %in% c(0, 1)),
  all(!is.na(d02$growth_segment))
)
first_treatment_period_02 <- min(
  d02$adoption_period[d02$ever_treated == 1],
  na.rm = TRUE
)

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d02,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 02: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Group averages, with the common adoption date marked.
# Color by growth segment and line type by treatment, so each segment's
# treated and control paths sit side by side.
d02_plot_data <- d02 |>
  mutate(treatment_group = if_else(ever_treated == 1, "Treated", "Control"))

d02_plot_data |>
  summarize(
    .by = treatment_group,
    n(),
    n_distinct(unit_id),
    n() / n_distinct(unit_id)
  )

d02_plot_data |>
  summarise(
    .by = c(period, treatment_group),
    mean_outcome = mean(outcome)
  ) |>
  ggplot(
    aes(
      x = period,
      y = mean_outcome,
      color = treatment_group,
      shape = treatment_group
    )
  ) +
  geom_line() +
  geom_point() +
  labs(
    x = "Period",
    y = "Mean outcome",
    color = "Treatment group",
    shape = "Treatment group"
  ) +
  theme_minimal()

d02_group_means <- d02_plot_data |>
  summarise(
    .by = c(period, growth_segment, treatment_group),
    mean_outcome = mean(outcome)
  )
d02_group_means
d02_group_plot <- ggplot(
  d02_group_means,
  aes(
    x = period,
    y = mean_outcome,
    color = fct_relevel(growth_segment, c("high", "middle", "low")),
    shape = treatment_group,
    linetype = treatment_group,
    group = interaction(growth_segment, treatment_group)
  )
) +
  geom_vline(
    xintercept = first_treatment_period_02,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 02: treatment composition and growth segment",
    x = "Period",
    y = "Mean outcome",
    color = "Growth segment",
    shape = "Group",
    linetype = "Group"
  ) +
  theme_minimal()
d02_group_plot

# 3. Baseline TWFE DiD.
twfe_fit_02 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d02
)
summary(twfe_fit_02)
twfe_table_02 <- broom::tidy(twfe_fit_02)

# 4. TWFE event study, joint Wald test of pre-period coefficients, and plot.
event_fit_02 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_02 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d02
)
summary(event_fit_02)
wald_pre_02 <- wald(event_fit_02, keep = "period::[1-5]:", print = TRUE)
wald_pre_02
iplot(
  event_fit_02,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 02: TWFE event study"
)

# 5. Pre-treatment placebo diagnostic.
placebo_period_02 <- first_treatment_period_02 - 3L
d02_placebo <- d02 |>
  filter(period < first_treatment_period_02) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_02
    )
  )
placebo_fit_02 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d02_placebo
)
summary(placebo_fit_02)

feols(outcome ~ period:ever_treated | unit_id + period, data = d02_placebo)

# 6. Compare weighted estimates with the unweighted estimates.
d02 |>
  summarize(.by = c(ever_treated, growth_segment), n())

d02_baseline <- d02 |>
  summarize(
    .by = unit_id,
    ever_treated = max(ever_treated), # treatment indicator
    growth_segment = first(growth_segment) # condition parallel trends on growth segment
  )

str(d02_baseline)
ipw_object_02 <- weightit(
  ever_treated ~ growth_segment,
  data = d02_baseline,
  method = "glm",
  estimand = "ATT"
)
summary(ipw_object_02)
# Inspect overlap, extreme weights, effective sample size, and segment balance.
d02_weight_diagnostics <- d02_baseline |>
  mutate(
    group = if_else(ever_treated == 1, "Treated", "Control"),
    propensity_score = ipw_object_02$ps,
    ipw_weight = ipw_object_02$weights
  )
d02_overlap_summary <- d02_weight_diagnostics |>
  summarise(
    .by = group,
    n_units = n(),
    propensity_score_min = min(propensity_score),
    propensity_score_max = max(propensity_score),
    weight_min = min(ipw_weight),
    weight_max = max(ipw_weight),
    effective_sample_size = sum(ipw_weight)^2 / sum(ipw_weight^2)
  )
d02_overlap_summary

d02_extreme_weights <- d02_weight_diagnostics |>
  arrange(group, desc(ipw_weight), unit_id) |>
  slice_head(n = 5L, by = group) |>
  select(group, unit_id, propensity_score, ipw_weight)
d02_extreme_weights

# Segment balance before and after weighting, with effective sample sizes.
d02_balance <- bal.tab(ipw_object_02, un = TRUE, disp = "means")
d02_balance
love.plot(ipw_object_02, thresholds = c(m = 0.1))

d02_weights <- d02_baseline |>
  transmute(unit_id, ipw_weight = ipw_object_02$weights)
d02_ipw_data <- d02 |>
  left_join(d02_weights, by = "unit_id")
stopifnot(
  nrow(d02_ipw_data) == nrow(d02),
  all(!is.na(d02_ipw_data$ipw_weight))
)
ipw_fit_02 <- feols(
  outcome ~ treated | unit_id + period,
  weights = ~ipw_weight,
  vcov = ~unit_id,
  data = d02_ipw_data
)
ipw_event_fit_02 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_02 - 1L) |
    unit_id + period,
  weights = ~ipw_weight,
  vcov = ~unit_id,
  data = d02_ipw_data
)
summary(ipw_fit_02)
summary(ipw_event_fit_02)
iplot(
  ipw_event_fit_02,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 02: IPW DiD event study"
)

# Callaway-Sant'Anna IPW aggregate and dynamic estimates.
d02_cs <- d02 |>
  mutate(
    id_num = as.integer(as.factor(unit_id)),
    cohort = if_else(is.na(adoption_period), 0, adoption_period)
  )
cs_fit_02 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~growth_segment,
  data = d02_cs,
  est_method = "ipw",
  control_group = "nevertreated",
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE
)
cs_fit_02
cs_simple_02 <- aggte(cs_fit_02, type = "simple")
cs_dynamic_02 <- aggte(
  cs_fit_02,
  type = "dynamic",
  bstrap = FALSE,
  cband = FALSE
)
summary(cs_simple_02)
ggdid(cs_dynamic_02, title = "Case 02: CS-IPW dynamic aggregation")

# Synthetic DiD without covariates.
d02_synthdid_panel <- panel.matrices(
  d02,
  unit = "unit_id",
  time = "period",
  outcome = "outcome",
  treatment = "treated"
)
synthdid_fit_02 <- with(
  d02_synthdid_panel,
  synthdid_estimate(Y, N0, T0)
)
synthdid_se_02 <- sqrt(vcov(synthdid_fit_02, method = "jackknife"))
plot(synthdid_fit_02, overlay = 1) +
  ggtitle("Case 02: Synthetic DiD")

plot(synthdid_fit_02, overlay = 0)

# Synthetic DiD never sees the segment. Compare its control weights with the
# unweighted control mix and the treated mix, by segment.
d02_synthdid_weights <- data.frame(
  unit_id = rownames(d02_synthdid_panel$Y)[seq_len(d02_synthdid_panel$N0)],
  omega = attr(synthdid_fit_02, "weights")$omega
) |>
  left_join(d02_baseline, by = "unit_id")
d02_treated_mix <- d02_baseline |>
  filter(ever_treated == 1) |>
  count(growth_segment) |>
  mutate(treated_share = n / sum(n))
d02_synthdid_segments <- d02_synthdid_weights |>
  summarise(
    .by = growth_segment,
    unweighted_control_share = n() / nrow(d02_synthdid_weights),
    synthdid_weight = sum(omega)
  ) |>
  left_join(
    select(d02_treated_mix, growth_segment, treated_share),
    by = "growth_segment"
  ) |>
  arrange(growth_segment)
d02_synthdid_segments
twfe_ipw_table_02 <- broom::tidy(ipw_fit_02)
comparison_02 <- data.frame(
  specification = c(
    "Unweighted TWFE",
    "IPW weighted TWFE (growth segment)",
    "Callaway-Sant'Anna IPW (growth segment)",
    "Synthetic DiD (no covariates)"
  ),
  estimate = c(
    twfe_table_02$estimate[twfe_table_02$term == "treated"],
    twfe_ipw_table_02$estimate[twfe_ipw_table_02$term == "treated"],
    unname(cs_simple_02$overall.att),
    as.numeric(synthdid_fit_02)
  ),
  standard_error = c(
    twfe_table_02$std.error[twfe_table_02$term == "treated"],
    twfe_ipw_table_02$std.error[twfe_ipw_table_02$term == "treated"],
    unname(cs_simple_02$overall.se),
    synthdid_se_02
  ),
  row.names = NULL
)
comparison_02

# ----------------------------------------------------------------------------
# CASE 03. Adoption decisions and pre-treatment outcome paths
# Why did the treated stores adopt when they did?
# Decide which pre-treatment weeks give a valid baseline, then compare
# estimators that use different baselines.
#
# Case introduction: Markdown optimization
#   Retail stores adopt markdown-optimization software during week 7. The data
#   are a store-week panel with log weekly revenue as the outcome.
#
#   Inspect the pre-treatment outcome paths before choosing a comparison
#   window and estimating the software effect.
# ----------------------------------------------------------------------------

# True effect:
#   ATT = 0: the software has no effect.
#   Adopting stores lose 0.30 log points in weeks 5-6 and recover from week 7
#   on their own, so a baseline that includes weeks 5-6 is biased.

case_file_03 <- file.path(data_dir, "case_03.csv")
d03 <- read.csv(case_file_03, stringsAsFactors = FALSE)
stopifnot(
  nrow(d03) > 0L,
  !anyDuplicated(paste(d03$unit_id, d03$period, sep = "::")),
  all(d03$treated %in% c(0, 1)),
  all(d03$ever_treated %in% c(0, 1)),
  any(d03$ever_treated == 1)
)
first_treatment_period_03 <- min(
  d03$adoption_period[d03$ever_treated == 1],
  na.rm = TRUE
)

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d03,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 03: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Group averages.
d03_plot_data <- d03 |>
  mutate(plot_group = if_else(ever_treated == 1, "Treated", "Control"))
d03_group_means <- d03_plot_data |>
  summarise(
    .by = c(period, plot_group),
    mean_outcome = mean(outcome)
  )
d03_group_means
d03_group_plot <- ggplot(
  d03_group_means,
  aes(
    x = period,
    y = mean_outcome,
    color = plot_group,
    shape = plot_group,
    linetype = plot_group
  )
) +
  geom_vline(xintercept = first_treatment_period_03, linetype = "dashed") +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 03: mean outcomes by treatment group",
    x = "Period",
    y = "Mean outcome",
    color = "Group",
    shape = "Group",
    linetype = "Group"
  ) +
  theme_minimal()
d03_group_plot

# 3. Baseline TWFE DiD.
twfe_fit_03 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d03
)
summary(twfe_fit_03)
twfe_table_03 <- broom::tidy(twfe_fit_03)

# 4. TWFE event study, joint Wald test of pre-period coefficients, and plot.
event_fit_03 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_03 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d03
)
summary(event_fit_03)
wald_pre_03 <- wald(event_fit_03, keep = "period::[1-5]:", print = TRUE)
wald_pre_03
iplot(
  event_fit_03,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 03: TWFE event study"
)

# 5. Pre-treatment placebo diagnostic.
placebo_period_03 <- first_treatment_period_03 - 3L
d03_placebo <- d03 |>
  filter(period < first_treatment_period_03) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_03
    )
  )
placebo_fit_03 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d03_placebo
)
summary(placebo_fit_03)

# 6. Compare estimators that use different baseline weeks.
d03_cs <- d03 |>
  mutate(
    id_num = as.integer(as.factor(unit_id)),
    cohort = if_else(is.na(adoption_period), 0, adoption_period)
  )
cs_fit_03 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~1,
  data = d03_cs,
  est_method = "dr",
  control_group = "nevertreated",
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE
)
cs_fit_03
cs_simple_03 <- aggte(cs_fit_03, type = "simple")
cs_dynamic_03 <- aggte(
  cs_fit_03,
  type = "dynamic",
  bstrap = FALSE,
  cband = FALSE
)
summary(cs_simple_03)
ggdid(cs_dynamic_03, title = "Case 03: CS dynamic aggregation")

# anticipation = 2 treats weeks 5-6 as possibly affected and moves the base
# period to week 4.  This assumes the pre-adoption dip was temporary.
cs_dip_excluded_fit_03 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~1,
  data = d03_cs,
  est_method = "dr",
  control_group = "nevertreated",
  ## Account for anticipation
  anticipation = 2,
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE
)
cs_dip_excluded_03 <- aggte(cs_dip_excluded_fit_03, type = "simple")
cs_dip_excluded_dynamic_03 <- aggte(
  cs_dip_excluded_fit_03,
  type = "dynamic",
  bstrap = FALSE,
  cband = FALSE
)
summary(cs_dip_excluded_03)
ggdid(
  cs_dip_excluded_dynamic_03,
  title = "Case 03: CS dynamic aggregation, anticipation = 2"
)


twfe_fit_03_anticipation <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d03 |>
    filter(period %notin% (first_treatment_period_03 - c(1:2)))
)
summary(twfe_fit_03_anticipation)
twfe_dip_excluded_03 <- broom::tidy(twfe_fit_03_anticipation)

# Synthetic DiD without covariates.
d03_synthdid_panel <- panel.matrices(
  d03,
  unit = "unit_id",
  time = "period",
  outcome = "outcome",
  treatment = "treated"
)
synthdid_fit_03 <- with(
  d03_synthdid_panel,
  synthdid_estimate(Y, N0, T0)
)
synthdid_se_03 <- sqrt(vcov(synthdid_fit_03, method = "jackknife"))
plot(synthdid_fit_03, overlay = 1) +
  ggtitle("Case 03: Synthetic DiD")

comparison_03 <- data.frame(
  specification = c(
    "TWFE",
    "TWFE, weeks 5-6 dropped",
    "Synthetic DiD",
    "Callaway-Sant'Anna aggregate (base week 6)",
    "Callaway-Sant'Anna aggregate, anticipation = 2 (base week 4)"
  ),
  estimate = c(
    twfe_table_03$estimate[twfe_table_03$term == "treated"],
    twfe_dip_excluded_03$estimate[twfe_dip_excluded_03$term == "treated"],
    as.numeric(synthdid_fit_03),
    unname(cs_simple_03$overall.att),
    unname(cs_dip_excluded_03$overall.att)
  ),
  standard_error = c(
    twfe_table_03$std.error[twfe_table_03$term == "treated"],
    twfe_dip_excluded_03$std.error[twfe_dip_excluded_03$term == "treated"],
    synthdid_se_03,
    unname(cs_simple_03$overall.se),
    unname(cs_dip_excluded_03$overall.se)
  ),
  row.names = NULL
)
comparison_03

# ----------------------------------------------------------------------------
# CASE 04. Staggered adoption
# How does staggered timing affect the interpretation of TWFE estimates?
# Compare TWFE and cohort-specific estimates.
#
# Case introduction: Customer-data platform
#   Regional retailers adopt a customer-data platform at different weeks. The
#   data are a retailer-week panel with log weekly revenue as the outcome.
#   `adoption_period` records the first treated week; missing values indicate
#   retailers that do not adopt during the observation window.
#
#   Estimate the platform's effect and explain which units provide the
#   comparison.
# ----------------------------------------------------------------------------

# True effect:
#   ATT = 0.25 log points for every cohort and every treated week.

case_file_04 <- file.path(data_dir, "case_04.csv")
d04 <- read.csv(case_file_04, stringsAsFactors = FALSE)
stopifnot(
  nrow(d04) > 0L,
  !anyDuplicated(paste(d04$unit_id, d04$period, sep = "::")),
  all(d04$treated %in% c(0, 1)),
  all(d04$ever_treated %in% c(0, 1)),
  any(is.na(d04$adoption_period)),
  length(unique(d04$adoption_period[d04$ever_treated == 1])) > 1L
)
first_treatment_period_04 <- min(
  d04$adoption_period[d04$ever_treated == 1],
  na.rm = TRUE
)
d04_cohort_counts <- d04 |>
  filter(ever_treated == 1) |>
  distinct(unit_id, adoption_period) |>
  count(adoption_period, name = "n_units")
d04_cohort_counts

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d04,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 04: panelview by adoption timing",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Cohort averages, with one vertical line per adoption cohort.
d04_plot_data <- d04 |>
  mutate(
    cohort_group = if_else(
      is.na(adoption_period),
      "Never treated",
      paste0("Adopted in week ", adoption_period)
    )
  )
d04_cohort_means <- d04_plot_data |>
  summarise(
    .by = c(period, cohort_group),
    mean_outcome = mean(outcome)
  )
d04_cohort_means
d04_treatment_periods <- sort(unique(
  d04$adoption_period[!is.na(d04$adoption_period)]
))
d04_cohort_plot <- ggplot(
  d04_cohort_means,
  aes(
    x = period,
    y = mean_outcome,
    color = cohort_group,
    shape = cohort_group,
    linetype = cohort_group
  )
) +
  geom_vline(
    xintercept = d04_treatment_periods,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 04: mean outcomes by adoption cohort",
    x = "Period",
    y = "Mean outcome",
    color = "Cohort",
    shape = "Cohort",
    linetype = "Cohort"
  ) +
  theme_minimal()
d04_cohort_plot

# 3. Baseline TWFE DiD.
twfe_fit_04 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d04
)
summary(twfe_fit_04)
twfe_table_04 <- broom::tidy(twfe_fit_04)

# 4. Goodman-Bacon decomposition: TWFE as a weighted average of 2x2 DiDs.
# bacon() prints the weight and average estimate of each comparison type.
bacon_04 <- bacon(
  outcome ~ treated,
  data = d04,
  id_var = "unit_id",
  time_var = "period"
)

# 5. Conventional TWFE event study using time relative to adoption.
# The omitted -1000 category keeps never-treated units in the comparison group.
d04 <- d04 |>
  mutate(
    event_time = if_else(
      ever_treated == 1,
      period - adoption_period,
      -1000
    )
  )
event_fit_04 <- feols(
  outcome ~ i(event_time, ref = c(-1, -1000)) | unit_id + period,
  vcov = ~unit_id,
  data = d04
)
summary(event_fit_04)
event_table_04 <- broom::tidy(event_fit_04)
wald_pre_04 <- wald(
  event_fit_04,
  keep = "event_time::-",
  print = TRUE
)
wald_pre_04
iplot(
  event_fit_04,
  xlab = "Periods relative to adoption; -1 is the reference",
  main = "Case 04: conventional TWFE event-time coefficients"
)
twfe_dynamic_04 <- data.frame(
  event_time = as.integer(sub(
    "^event_time::",
    "",
    event_table_04$term
  )),
  twfe_estimate = event_table_04$estimate,
  twfe_standard_error = event_table_04$std.error,
  row.names = NULL
)

# 6. Pre-treatment placebo diagnostic.
placebo_period_04 <- first_treatment_period_04 - 3L
d04_placebo <- d04 |>
  filter(period < first_treatment_period_04) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_04
    )
  )
placebo_fit_04 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d04_placebo
)
summary(placebo_fit_04)

# 7. Compare TWFE with Callaway--Sant'Anna and Sun--Abraham estimates.
d04_cs <- d04 |>
  mutate(
    id_num = as.integer(as.factor(unit_id)),
    cohort = if_else(is.na(adoption_period), 0, adoption_period)
  )
cs_fit_04 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~1,
  data = d04_cs,
  est_method = "dr",
  control_group = "nevertreated",
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE
)
cs_fit_04
cs_simple_04 <- aggte(cs_fit_04, type = "simple")
cs_dynamic_04 <- aggte(
  cs_fit_04,
  type = "dynamic",
  bstrap = FALSE,
  cband = FALSE
)
summary(cs_simple_04)
ggdid(cs_dynamic_04, title = "Case 04: CS dynamic aggregation")
ggdid(cs_fit_04, title = "Case 04: No aggregation")
cs_dynamic_table_04 <- data.frame(
  event_time = as.integer(cs_dynamic_04$egt),
  cs_estimate = unname(cs_dynamic_04$att.egt),
  cs_standard_error = unname(cs_dynamic_04$se.egt),
  row.names = NULL
)
# Keep periods available from only one estimator; TWFE omits its -1 reference.
dynamic_comparison_04 <- full_join(
  twfe_dynamic_04,
  cs_dynamic_table_04,
  by = "event_time"
) |>
  arrange(event_time)
dynamic_comparison_04

# Sun--Abraham interaction-weighted estimator.  sunab() needs a cohort for
# never-treated stores, so they get one beyond the panel.
d04_sunab <- d04 |>
  mutate(cohort = if_else(is.na(adoption_period), 10000, adoption_period))
sunab_fit_04 <- feols(
  outcome ~ sunab(cohort, period) | unit_id + period,
  vcov = ~unit_id,
  data = d04_sunab
)
summary(sunab_fit_04)
summary(sunab_fit_04, agg = "ATT")
sunab_table_04 <- broom::tidy(summary(sunab_fit_04, agg = "ATT"))
comparison_04 <- data.frame(
  specification = c(
    "TWFE",
    "Callaway-Sant'Anna aggregate",
    "Sun-Abraham aggregate"
  ),
  estimate = c(
    twfe_table_04$estimate[twfe_table_04$term == "treated"],
    unname(cs_simple_04$overall.att),
    sunab_table_04$estimate
  ),
  standard_error = c(
    twfe_table_04$std.error[twfe_table_04$term == "treated"],
    unname(cs_simple_04$overall.se),
    sunab_table_04$std.error
  ),
  row.names = NULL
)
comparison_04

# ----------------------------------------------------------------------------
# CASE 05. Staggered adoption with growing effects
# Which 2x2 comparisons make up the TWFE estimate, and which of them use
# already-treated markets as controls?
#
# Case introduction: Recommendation-engine rollout
#   Markets adopt a recommendation engine at different dates. The data are a
#   market-week panel with log weekly revenue as the outcome.
#   `adoption_period` records each market's first treated week.
#
#   Estimate the effect by cohort and over exposure time, then compare it with
#   a static model.
# ----------------------------------------------------------------------------

# True effect:
#   Effect = 0.10 per week since adoption: 0.10 in the first treated week,
#   1.00 after ten weeks.
#   ATT averaged over treated market-weeks = 0.467.

case_file_05 <- file.path(data_dir, "case_05.csv")
d05 <- read.csv(case_file_05, stringsAsFactors = FALSE)
stopifnot(
  nrow(d05) > 0L,
  !anyDuplicated(paste(d05$unit_id, d05$period, sep = "::")),
  all(d05$treated %in% c(0, 1)),
  all(d05$ever_treated %in% c(0, 1)),
  any(is.na(d05$adoption_period)),
  length(unique(d05$adoption_period[d05$ever_treated == 1])) > 1L
)
first_treatment_period_05 <- min(
  d05$adoption_period[d05$ever_treated == 1],
  na.rm = TRUE
)
d05_cohort_counts <- d05 |>
  filter(ever_treated == 1) |>
  distinct(unit_id, adoption_period) |>
  count(adoption_period, name = "n_units")
d05_cohort_counts

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d05,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 05: panelview by adoption timing",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Cohort averages, with one vertical line per adoption cohort.
d05_plot_data <- d05 |>
  mutate(
    cohort_group = if_else(
      is.na(adoption_period),
      "Never treated",
      paste0("Adopted in week ", adoption_period)
    )
  )
d05_cohort_means <- d05_plot_data |>
  summarise(
    .by = c(period, cohort_group),
    mean_outcome = mean(outcome)
  )
d05_cohort_means
d05_treatment_periods <- sort(unique(
  d05$adoption_period[!is.na(d05$adoption_period)]
))
d05_cohort_plot <- ggplot(
  d05_cohort_means,
  aes(
    x = period,
    y = mean_outcome,
    color = cohort_group,
    shape = cohort_group,
    linetype = cohort_group
  )
) +
  geom_vline(
    xintercept = d05_treatment_periods,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 05: mean outcomes by adoption cohort",
    x = "Period",
    y = "Mean outcome",
    color = "Cohort",
    shape = "Cohort",
    linetype = "Cohort"
  ) +
  theme_minimal()
d05_cohort_plot

# 3. Baseline TWFE DiD.
twfe_fit_05 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d05
)
summary(twfe_fit_05)
twfe_table_05 <- broom::tidy(twfe_fit_05)

# 4. Goodman-Bacon decomposition: TWFE as a weighted average of 2x2 DiDs.
# bacon() prints the weight and average estimate of each comparison type.
bacon_05 <- bacon(
  outcome ~ treated,
  data = d05,
  id_var = "unit_id",
  time_var = "period"
)

# 5. Pre-treatment placebo diagnostic.
placebo_period_05 <- first_treatment_period_05 - 3L
d05_placebo <- d05 |>
  filter(period < first_treatment_period_05) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_05
    )
  )
placebo_fit_05 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d05_placebo
)
summary(placebo_fit_05)

# 6. Compare TWFE with the Callaway--Sant'Anna and Sun--Abraham aggregates.
d05_cs <- d05 |>
  mutate(
    id_num = as.integer(as.factor(unit_id)),
    cohort = if_else(is.na(adoption_period), 0, adoption_period)
  )
cs_fit_05 <- att_gt(
  yname = "outcome",
  tname = "period",
  idname = "id_num",
  gname = "cohort",
  xformla = ~1,
  data = d05_cs,
  est_method = "reg",
  control_group = "notyettreated",
  bstrap = FALSE,
  panel = TRUE,
  print_details = FALSE
)
cs_fit_05
cs_simple_05 <- aggte(cs_fit_05, type = "simple")
summary(cs_simple_05)

# Sun--Abraham interaction-weighted estimator.  sunab() needs a cohort for
# never-treated stores, so they get one beyond the panel.
d05_sunab <- d05 |>
  mutate(cohort = if_else(is.na(adoption_period), 10000, adoption_period))
sunab_fit_05 <- feols(
  outcome ~ sunab(cohort, period) | unit_id + period,
  vcov = ~unit_id,
  data = d05_sunab
)
summary(sunab_fit_05, agg = "ATT")
sunab_table_05 <- broom::tidy(summary(sunab_fit_05, agg = "ATT"))
comparison_05 <- data.frame(
  specification = c(
    "Conventional TWFE",
    "Callaway-Sant'Anna aggregate",
    "Sun-Abraham aggregate"
  ),
  estimate = c(
    twfe_table_05$estimate[twfe_table_05$term == "treated"],
    unname(cs_simple_05$overall.att),
    sunab_table_05$estimate
  ),
  standard_error = c(
    twfe_table_05$std.error[twfe_table_05$term == "treated"],
    unname(cs_simple_05$overall.se),
    sunab_table_05$std.error
  ),
  row.names = NULL
)
comparison_05

# ----------------------------------------------------------------------------
# CASE 06. Geographic exposure and comparison-group choice
# Does geographic exposure matter when choosing comparison stores?
# Use the exposure information and outcome paths to assess which stores
# provide a suitable comparison group.
#
# Case introduction: Local promotion
#   One store in each local market receives a promotion during week 7. The
#   data are a store-week panel with log weekly revenue as the outcome.
#   `distance_to_treated` records the store's pre-treatment distance from the
#   treated store.
#
#   Use the market structure and exposure information to define an appropriate
#   comparison group.
# ----------------------------------------------------------------------------

# True effect:
#   Direct effect on promoted stores = 0.30 log points.
#   Nearby untreated stores gain a 0.12 spillover; distant stores gain nothing.

case_file_06 <- file.path(data_dir, "case_06.csv")
d06 <- read.csv(case_file_06, stringsAsFactors = FALSE)
stopifnot(
  nrow(d06) > 0L,
  !anyDuplicated(paste(d06$unit_id, d06$period, sep = "::")),
  all(d06$treated %in% c(0, 1)),
  all(d06$ever_treated %in% c(0, 1)),
  all(!is.na(d06$distance_to_treated))
)
first_treatment_period_06 <- min(
  d06$adoption_period[d06$ever_treated == 1],
  na.rm = TRUE
)

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d06,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 06: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Inspect outcome paths by geographic exposure.
d06 <- d06 |>
  mutate(
    exposure_group = case_when(
      ever_treated == 1 ~ "Promoted store",
      distance_to_treated == 1 ~ "Nearby untreated store",
      TRUE ~ "Distant untreated store"
    )
  )
d06_exposure_means <- d06 |>
  summarise(
    .by = c(period, exposure_group),
    mean_outcome = mean(outcome)
  )
d06_exposure_means <- d06_exposure_means |>
  mutate(
    .by = exposure_group,
    change_from_pre = mean_outcome -
      mean(mean_outcome[period < first_treatment_period_06])
  )
d06_exposure_means
d06_exposure_plot <- ggplot(
  d06_exposure_means,
  aes(
    x = period,
    y = change_from_pre,
    color = exposure_group,
    shape = exposure_group,
    linetype = exposure_group
  )
) +
  geom_vline(xintercept = first_treatment_period_06, linetype = "dashed") +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 06: outcomes by geographic exposure",
    subtitle = "Change from each group's pre-promotion mean; dashed line marks launch",
    x = "Period",
    y = "Change in mean outcome",
    color = "Exposure",
    shape = "Exposure",
    linetype = "Exposure"
  ) +
  theme_minimal()
d06_exposure_plot

# 3. Baseline TWFE DiD.
twfe_fit_06 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d06
)
summary(twfe_fit_06)
twfe_table_06 <- broom::tidy(twfe_fit_06)

# 4. TWFE event study, joint Wald test of pre-period coefficients, and plot.
event_fit_06 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_06 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d06
)
summary(event_fit_06)
wald_pre_06 <- wald(event_fit_06, keep = "period::[1-5]:", print = TRUE)
wald_pre_06
iplot(
  event_fit_06,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 06: TWFE event study"
)

# 5. Pre-treatment placebo diagnostic.
placebo_period_06 <- first_treatment_period_06 - 3L
d06_placebo <- d06 |>
  filter(period < first_treatment_period_06) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_06
    )
  )
placebo_fit_06 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d06_placebo
)
summary(placebo_fit_06)

# 6. Compare estimates across control-group definitions.
# Distance 5 or greater is one candidate definition of a far control.
far_control_data_06 <- d06 |>
  filter(ever_treated == 1 | distance_to_treated >= 5)
far_controls_fit_06 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = far_control_data_06
)
far_controls_table_06 <- broom::tidy(far_controls_fit_06)
comparison_06 <- data.frame(
  specification = c(
    "TWFE using all controls",
    "TWFE using far controls only"
  ),
  estimate = c(
    twfe_table_06$estimate[twfe_table_06$term == "treated"],
    far_controls_table_06$estimate[far_controls_table_06$term == "treated"]
  ),
  standard_error = c(
    twfe_table_06$std.error[twfe_table_06$term == "treated"],
    far_controls_table_06$std.error[far_controls_table_06$term == "treated"]
  ),
  row.names = NULL
)
comparison_06

# ----------------------------------------------------------------------------
# CASE 07. Effects that differ across stores
# Does the effect of app ordering differ between urban and rural stores?
# State the parallel-trends assumption for each location, then compare a
# regression with location-specific week effects to one without them.
#
# Case introduction: Mobile-app ordering
#   A grocery chain launches mobile-app ordering in some stores during week 7.
#   The data are a store-week panel, and `outcome` is log weekly revenue.
#   `urban` records whether the store is in an urban location; it is fixed
#   before the launch.
#
#   Estimate the effect of app ordering and whether it differs between urban
#   and rural stores.
# ----------------------------------------------------------------------------

# True effect:
#   ATT = 0.15 log points in rural stores and 0.35 in urban stores
#   (difference 0.20; ATT over all treated stores 0.29).
#   Urban stores grow 0.02 per week faster whether or not they launch;
#   within a location, trends are parallel.

case_file_07 <- file.path(data_dir, "case_07.csv")
d07 <- read.csv(case_file_07, stringsAsFactors = FALSE)
stopifnot(
  nrow(d07) > 0L,
  !anyDuplicated(paste(d07$unit_id, d07$period, sep = "::")),
  all(d07$treated %in% c(0, 1)),
  all(d07$ever_treated %in% c(0, 1)),
  all(d07$urban %in% c(0, 1)),
  length(unique(d07$adoption_period[d07$ever_treated == 1])) == 1L
)
first_treatment_period_07 <- min(
  d07$adoption_period[d07$ever_treated == 1],
  na.rm = TRUE
)

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d07,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 07: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Stores per location and group, then group averages.
# Every location needs treated and control stores (overlap).
d07_store_counts <- d07 |>
  distinct(unit_id, urban, ever_treated) |>
  count(urban, ever_treated, name = "n_stores")
d07_store_counts

d07_plot_data <- d07 |>
  mutate(
    treatment_group = if_else(ever_treated == 1, "Treated", "Control"),
    location = if_else(urban == 1, "Urban", "Rural")
  )
d07_group_means <- d07_plot_data |>
  summarise(
    .by = c(period, location, treatment_group),
    mean_outcome = mean(outcome)
  )
d07_group_means
d07_group_plot <- ggplot(
  d07_group_means,
  aes(
    x = period,
    y = mean_outcome,
    color = location,
    shape = treatment_group,
    linetype = treatment_group,
    group = interaction(location, treatment_group)
  )
) +
  geom_vline(
    xintercept = first_treatment_period_07,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  labs(
    title = "Case 07: mean outcomes by location and treatment",
    x = "Period",
    y = "Mean outcome",
    color = "Location",
    shape = "Group",
    linetype = "Group"
  ) +
  theme_minimal()
d07_group_plot

# 3. Baseline TWFE DiD, pooling both locations.
twfe_fit_07 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d07
)
summary(twfe_fit_07)
twfe_table_07 <- broom::tidy(twfe_fit_07)

# 4. One 2x2 per location, by hand.
# Each location's DiD compares treated and control stores in that location
# only, so a trend shared by all urban stores cancels.
d07_changes <- d07 |>
  mutate(post = period >= first_treatment_period_07) |>
  summarise(
    .by = c(urban, ever_treated),
    change = mean(outcome[post]) - mean(outcome[!post])
  )
d07_changes
d07_did_by_location <- d07_changes |>
  summarise(
    .by = urban,
    did = change[ever_treated == 1] - change[ever_treated == 0]
  ) |>
  arrange(urban)
d07_did_by_location
did_difference_07 <- diff(d07_did_by_location$did)
did_difference_07

# 5. The same numbers from one regression.
#   outcome_it = alpha_i + lambda_{t, urban_i} + beta treated_it
#                + gamma treated_it * urban_i + e_it
# period^urban gives each location its own week effects.  beta is the rural
# effect, beta + gamma the urban effect, and gamma their difference.
d07 <- d07 |>
  mutate(
    treated_urban = treated * urban,
    ever_treated_urban = ever_treated * urban
  )
triple_fit_07 <- feols(
  outcome ~ treated + treated_urban | unit_id + period^urban,
  vcov = ~unit_id,
  data = d07
)
summary(triple_fit_07)
# The same model with one coefficient per location, for standard errors on
# each location's effect.
location_fit_07 <- feols(
  outcome ~ i(urban, treated) | unit_id + period^urban,
  vcov = ~unit_id,
  data = d07
)
summary(location_fit_07)
sum(coef(triple_fit_07))
# Delta-method standard error of beta + gamma, the urban effect.
urban_effect_07 <- marginaleffects::hypotheses(
  triple_fit_07,
  hypothesis = "treated + treated_urban = 0"
)
urban_effect_07
# Each store belongs to one location and each location has its own week
# effects, so the pooled model equals one TWFE regression per location.  The
# estimates match exactly; standard errors differ only through the
# small-sample correction, which counts all stores in the pooled model.
separate_fits_07 <- feols(
  outcome ~ treated | unit_id + period,
  split = ~urban,
  vcov = ~unit_id,
  data = d07
)
summary(separate_fits_07)
separate_comparison_07 <- data.frame(
  location = c("Rural", "Urban"),
  separate_estimate = sapply(separate_fits_07, coef),
  separate_standard_error = sapply(separate_fits_07, se),
  pooled_estimate = unname(coef(location_fit_07)),
  pooled_standard_error = unname(se(location_fit_07)),
  row.names = NULL
)
separate_comparison_07

# 6. Event study with location-specific week effects.
# The ever_treated_urban leads test whether urban treated stores already
# diverged from urban controls, relative to the rural comparison.
event_fit_07 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_07 - 1L) +
    i(period, ever_treated_urban, ref = first_treatment_period_07 - 1L) |
    unit_id + period^urban,
  vcov = ~unit_id,
  data = d07
)
summary(event_fit_07)
wald_pre_07 <- wald(
  event_fit_07,
  keep = "period::[1-5]:ever_treated_urban",
  print = TRUE
)
wald_pre_07

# 7. Without location-specific week effects.
# Urban and rural stores now share one set of week effects, so the urban
# stores' own trend is attributed to treatment in the urban coefficient.
omit_fit_07 <- feols(
  outcome ~ treated + treated_urban | unit_id + period,
  vcov = ~unit_id,
  data = d07
)
summary(omit_fit_07)
event_omit_fit_07 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_07 - 1L) +
    i(period, ever_treated_urban, ref = first_treatment_period_07 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d07
)
summary(event_omit_fit_07)
wald_pre_omit_07 <- wald(
  event_omit_fit_07,
  keep = "period::[1-5]:ever_treated_urban",
  print = TRUE
)
wald_pre_omit_07
iplot(
  list(event_fit_07, event_omit_fit_07),
  i.select = 2,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 07: urban minus rural event-study coefficients"
)
legend(
  "topleft",
  legend = c("Location-specific week effects", "Common week effects"),
  col = 1:2,
  pch = c(20, 17),
  bty = "n"
)

# 8. Compare the estimates.
triple_table_07 <- broom::tidy(triple_fit_07)
omit_table_07 <- broom::tidy(omit_fit_07)
# Delta-method urban effect without location-specific week effects.
omit_urban_effect_07 <- marginaleffects::hypotheses(
  omit_fit_07,
  "treated + treated_urban = 0",
  vcov = vcov(omit_fit_07, cluster = "unit_id")
)
comparison_07 <- data.frame(
  quantity = c(
    "Rural effect",
    "Rural effect",
    "Urban effect (delta method)",
    "Urban effect (delta method)",
    "Urban minus rural",
    "Urban minus rural",
    "ATT, all treated stores"
  ),
  specification = c(
    "Location-specific week effects",
    "Common week effects",
    "Location-specific week effects",
    "Common week effects",
    "Location-specific week effects",
    "Common week effects",
    "Pooled TWFE"
  ),
  estimate = c(
    triple_table_07$estimate[triple_table_07$term == "treated"],
    omit_table_07$estimate[omit_table_07$term == "treated"],
    urban_effect_07$estimate,
    omit_urban_effect_07$estimate,
    triple_table_07$estimate[triple_table_07$term == "treated_urban"],
    omit_table_07$estimate[omit_table_07$term == "treated_urban"],
    twfe_table_07$estimate[twfe_table_07$term == "treated"]
  ),
  standard_error = c(
    triple_table_07$std.error[triple_table_07$term == "treated"],
    omit_table_07$std.error[omit_table_07$term == "treated"],
    urban_effect_07$std.error,
    omit_urban_effect_07$std.error,
    triple_table_07$std.error[triple_table_07$term == "treated_urban"],
    omit_table_07$std.error[omit_table_07$term == "treated_urban"],
    twfe_table_07$std.error[twfe_table_07$term == "treated"]
  ),
  row.names = NULL
)
comparison_07

# ----------------------------------------------------------------------------
# CASE 08. Controlling for a post-treatment variable
# Should engagement be included in the model?
# Check whether engagement changes when advertising starts, then compare the
# estimate with and without engagement as a control.
#
# Case introduction: Personalized advertising
#   A personalized advertising intervention is introduced in week 7. The data
#   are a store-week panel with log weekly revenue as the outcome.
#   `engagement` records product-page engagement for each store-week.
#
#   Estimate the advertising effect and explain which variables should be
#   included in the model.
# ----------------------------------------------------------------------------

# True effect:
#   Total effect = 0.35 log points: a direct effect of 0.20 plus 0.15 through
#   engagement.  Advertising raises engagement by 0.80, and revenue rises
#   0.1875 per unit of engagement (0.1875 * 0.80 = 0.15).

case_file_08 <- file.path(data_dir, "case_08.csv")
d08 <- read.csv(case_file_08, stringsAsFactors = FALSE)
stopifnot(
  nrow(d08) > 0L,
  !anyDuplicated(paste(d08$unit_id, d08$period, sep = "::")),
  all(d08$treated %in% c(0, 1)),
  all(d08$ever_treated %in% c(0, 1)),
  all(!is.na(d08$engagement))
)
first_treatment_period_08 <- min(
  d08$adoption_period[d08$ever_treated == 1],
  na.rm = TRUE
)

# 1. Panelview.
panelview(
  outcome ~ treated,
  data = d08,
  index = c("unit_id", "period"),
  pre.post = TRUE,
  by.timing = TRUE,
  main = "Case 08: panelview",
  axis.adjust = TRUE,
  axis.lab = "time"
)

# 2. Group averages of revenue and engagement.
d08_group_means <- d08 |>
  mutate(group = if_else(ever_treated == 1, "Treated", "Control")) |>
  summarise(
    .by = c(period, group),
    outcome = mean(outcome),
    engagement = mean(engagement)
  )
d08_group_means
d08_plot_data <- bind_rows(
  transmute(
    d08_group_means,
    period,
    group,
    measure = "Outcome",
    value = outcome
  ),
  transmute(
    d08_group_means,
    period,
    group,
    measure = "Engagement",
    value = engagement
  )
)
d08_group_plot <- ggplot(
  d08_plot_data,
  aes(x = period, y = value, color = group, shape = group, linetype = group)
) +
  geom_vline(xintercept = first_treatment_period_08, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~measure, scales = "free_y") +
  labs(
    title = "Case 08: revenue and engagement by treatment group",
    x = "Period",
    y = "Group mean",
    color = "Group",
    shape = "Group",
    linetype = "Group"
  ) +
  theme_minimal()
d08_group_plot

# 3. Baseline TWFE DiD without engagement.
twfe_fit_08 <- feols(
  outcome ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d08
)
summary(twfe_fit_08)
twfe_table_08 <- broom::tidy(twfe_fit_08)

# 4. TWFE event study, joint Wald test of pre-period coefficients, and plot.
event_fit_08 <- feols(
  outcome ~ i(period, ever_treated, ref = first_treatment_period_08 - 1L) |
    unit_id + period,
  vcov = ~unit_id,
  data = d08
)
summary(event_fit_08)
wald_pre_08 <- wald(event_fit_08, keep = "period::[1-5]:", print = TRUE)
wald_pre_08
iplot(
  event_fit_08,
  xlab = "Calendar period; week 6 is the reference",
  main = "Case 08: TWFE event study"
)

# 5. Pre-treatment placebo diagnostic.
placebo_period_08 <- first_treatment_period_08 - 3L
d08_placebo <- d08 |>
  filter(period < first_treatment_period_08) |>
  mutate(
    placebo_treated = as.integer(
      ever_treated == 1 & period >= placebo_period_08
    )
  )
placebo_fit_08 <- feols(
  outcome ~ placebo_treated | unit_id + period,
  vcov = ~unit_id,
  data = d08_placebo
)
summary(placebo_fit_08)

# 6. Does advertising move engagement?
# A DiD with engagement as the outcome.  A clear effect means engagement is a
# post-treatment variable on the path from advertising to revenue.
engagement_fit_08 <- feols(
  engagement ~ treated | unit_id + period,
  vcov = ~unit_id,
  data = d08
)
summary(engagement_fit_08)

# 7. TWFE DiD controlling for engagement.
# Holding engagement fixed removes the part of the effect that runs through it,
# so this regression estimates the direct effect only.
mediator_fit_08 <- feols(
  outcome ~ treated + engagement | unit_id + period,
  vcov = ~unit_id,
  data = d08
)
summary(mediator_fit_08)
mediator_table_08 <- broom::tidy(mediator_fit_08)

# 8. Compare the estimates.
comparison_08 <- data.frame(
  specification = c(
    "TWFE without engagement (total effect)",
    "TWFE controlling for engagement (direct effect)"
  ),
  estimate = c(
    twfe_table_08$estimate[twfe_table_08$term == "treated"],
    mediator_table_08$estimate[mediator_table_08$term == "treated"]
  ),
  standard_error = c(
    twfe_table_08$std.error[twfe_table_08$term == "treated"],
    mediator_table_08$std.error[mediator_table_08$term == "treated"]
  ),
  row.names = NULL
)
comparison_08
