suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

get_project_root <- function() {
  script_path <- normalizePath(commandArgs(trailingOnly = FALSE), winslash = "/", mustWork = FALSE)
  script_arg <- script_path[str_detect(script_path, "^--file=")]
  if (length(script_arg) == 0) {
    stop("Run this script with Rscript so --file is available.")
  }
  script_file <- sub("^--file=", "", script_arg[1])
  # File is expected at code/R/pet/
  normalizePath(file.path(dirname(script_file), "..", "..", ".."), winslash = "/")
}

project_root <- get_project_root()

input_csv <- file.path(project_root, "data", "processed_pet", "merged_summary_csv", "summary_suvr_mean_harmonized.csv")
output_dir <- file.path(project_root, "results", "pet", "ancova_rois")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_csv)) {
  stop(paste("Input file not found:", input_csv))
}

df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c("ID", "group", "site", "sex", "age")
missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop(paste("Missing required columns:", paste(missing_required, collapse = ", ")))
}

excluded <- c("ID", "group", "site", "sex", "age", "weight", "analysis_group")
roi_cols <- setdiff(names(df), excluded)

if (length(roi_cols) == 0) {
  stop("No ROI columns detected.")
}

# Keep LC and controls from all available sites (3 cohorts in this dataset).
pet_long <- df %>%
  filter(group %in% c("LC", "Controls")) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    sex = as.factor(sex),
    site = as.factor(site),
    group = factor(group, levels = c("Controls", "LC")),
    subject = as.factor(ID)
  ) %>%
  pivot_longer(
    cols = all_of(roi_cols),
    names_to = "ROI",
    values_to = "SUVR"
  ) %>%
  mutate(
    SUVR = suppressWarnings(as.numeric(SUVR)),
    ROI = factor(ROI, levels = roi_cols)
  ) %>%
  filter(
    !is.na(SUVR),
    !is.na(age),
    !is.na(sex),
    !is.na(site),
    !is.na(group),
    !is.na(subject)
  ) %>%
  mutate(age_c = age - mean(age, na.rm = TRUE))

message("Subjects by site/group:")
print(df %>% count(site, group))
message(paste("Rows in long dataframe:", nrow(pet_long)))
message(paste("ROIs in model:", length(roi_cols)))

# Requested model:
# lmer(SUVR ~ group * ROI + age_c + sex + site + (1 | subject), data = pet_long)
fit <- lmer(
  SUVR ~ group * ROI + age_c + sex + site + (1 | subject),
  data = pet_long,
  REML = FALSE
)

summary_path <- file.path(output_dir, "mixedlm_lc_vs_controls_all_groups_lmer_summary.txt")
capture.output(summary(fit), file = summary_path)

fixed_effects <- as.data.frame(summary(fit)$coefficients)
fixed_effects$term <- rownames(fixed_effects)
rownames(fixed_effects) <- NULL
fixed_effects <- fixed_effects %>%
  rename(
    estimate = Estimate,
    se = `Std. Error`,
    df = df,
    t = `t value`,
    p = `Pr(>|t|)`
  ) %>%
  select(term, estimate, se, df, t, p) %>%
  arrange(p)

write_csv(
  fixed_effects,
  file.path(output_dir, "mixedlm_lc_vs_controls_all_groups_lmer_fixed_effects.csv")
)

emm <- emmeans(fit, ~ group | ROI)
contr <- contrast(emm, method = list(LC_minus_Controls = c(-1, 1)), adjust = "none")

roi_contrasts <- as.data.frame(contr) %>%
  rename(
    estimate_lc_minus_control = estimate,
    se = SE,
    df = df,
    t = t.ratio,
    p = p.value
  ) %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"),
    significant_fdr = p_fdr < 0.05
  ) %>%
  arrange(p)

write_csv(
  roi_contrasts,
  file.path(output_dir, "mixedlm_lc_vs_controls_all_groups_lmer_roi_contrasts.csv")
)

write_csv(
  filter(roi_contrasts, significant_fdr),
  file.path(output_dir, "mixedlm_lc_vs_controls_all_groups_lmer_roi_contrasts_significant.csv")
)

message(paste("Saved model summary:", summary_path))
message("Saved ROI contrasts to results/pet/ancova_rois/")
message(paste("Significant ROIs (FDR < 0.05):", sum(roi_contrasts$significant_fdr)))