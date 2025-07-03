library(survival)
library(dplyr)
library(readr)
library(stringr)
library(pbapply)

# === Load expression matrix ===
expr_file <- "data/processed/proteomics/iso_log/TMT10_iso_log_PDC000116_combined.csv"
expr_df <- read_csv(expr_file, show_col_types = FALSE) %>%
  filter(!is.na(feature) & feature != "")
expr_df <- as.data.frame(expr_df)
rownames(expr_df) <- expr_df$feature
expr_df <- expr_df[, -1]

# === Extract sample metadata ===
sample_cols <- colnames(expr_df)
batch <- str_extract(sample_cols, "^\\d{2}CPTAC")
batch <- ifelse(is.na(batch), "Unknown", batch)
batch <- factor(batch)

sample_case_ids <- str_extract(sample_cols, "(?<=_matrix_)([0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+)")
keep <- !is.na(sample_case_ids)
expr_df <- expr_df[, keep]
batch <- batch[keep]
sample_case_ids <- sample_case_ids[keep]
colnames(expr_df) <- sample_case_ids

# === Match to clinical ===
clin_df <- clin_df_all %>%
  filter(case_id %in% sample_case_ids) %>%
  mutate(
    tumor_stage_clean = factor(str_extract(tumor_stage, "Stage [IVX]+"),
                               levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                               ordered = TRUE),
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA)
  )

matched_ids <- intersect(colnames(expr_df), clin_df$case_id)
expr_df <- expr_df[, matched_ids]
clin_df <- clin_df %>% filter(case_id %in% matched_ids) %>%
  arrange(match(case_id, matched_ids))
batch <- batch[match(matched_ids, sample_case_ids)]




# Clinical metadata
# === Load and clean clinical data ===
clin_df_all <- read_csv(clinical_file, show_col_types = FALSE) %>%
  rename(
    case_id = Patient_ID,
    age = `consent/age`,
    sex = `consent/sex`,
    race = `consent/race`,
    ethnicity = `consent/ethnicity`,
    tumor_size_cm = `baseline/tumor_size_cm`,
    bmi = `medical_history/bmi`,
    alcohol_consumption = `medical_history/alcohol_consumption`,
    tobacco_smoking_history = `medical_history/tobacco_smoking_history`,
    vital_status = `follow-up/vital_status_at_date_of_last_contact`,
    days_to_last_contact = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`,
    days_to_death = `follow-up/number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`,
    tumor_stage = `baseline/tumor_stage_pathological`,
    histologic_subtype = `baseline/histologic_type`,
    histologic_grade = `cptac_path/histologic_grade`,
    OS_time = `Overall survival, days`,
    OS_event = `Survival status (1, dead; 0, alive)`
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),  # suppress warning if you want
    age_log2 = ifelse(!is.na(age) & age > 0, log2(age), NA),
    age_group = cut(age,
                    breaks = c(0, 40, 50, 60, 70, 80, Inf),
                    right = FALSE,
                    labels = c("<40", "40–49", "50–59", "60–69", "70–79", "80+"))
  )

clin_df_all <- clin_df_all %>%
  mutate(tumor_stage_clean = factor(str_extract(tumor_stage, "Stage [IVX]+"),
                                    levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                                    ordered = TRUE))


clin_df_all <- clin_df_all %>%
  mutate(
    tumor_stage_clean = factor(str_extract(tumor_stage, "Stage [IVX]+"),
                               levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                               ordered = TRUE),
    
    bmi = suppressWarnings(as.numeric(bmi)),  # ensure BMI is numeric
    bmi_group = cut(bmi,
                    breaks = c(-Inf, 18.5, 25, 30, 35, 40, Inf),
                    labels = c("Underweight", "Normal", "Overweight", "Obese I", "Obese II", "Obese III"),
                    right = FALSE)
  )




# Investigating Main Effects Only
    # This helps identify which variables are individually important (via their coefficients and p-values).
fit_main <- coxph(Surv(OS_time, OS_event) ~ age_log2 + tumor_stage + sex, data = clin_df_all)
summary(fit_main)


# Test PH Assumptions
# Violations suggest need for stratification (strata()) or time-dependent terms.
zph <- cox.zph(fit_main)
print(zph)


# Test for Interactions (Synergy)
  # This will estimate:
      # - Main effects: age_log2, tumor_stage_clean, sex
      # - Interaction: age_log2:tumor_stage_clean
  # The interaction p-value tells you whether the effect of age differs across tumor stages, i.e., synergistic.
fit_interact <- coxph(Surv(OS_time, OS_event) ~ age_log2 * tumor_stage + histologic_grade, data = clin_df_all)
summary(fit_interact)



# Compare Models (with vs. without interaction)
    # A significant p-value means the interaction adds predictive value, i.e., synergy is important.
anova(fit_main, fit_interact, test = "LRT")


# Repeat for All Variable Pairs
     # This tells you which metadata interactions are statistically significant.
library(survival)
library(dplyr)
library(readr)
library(stringr)

# === Define metadata variables to test ===
vars <- c("age_log2", "tumor_stage", "bmi_group", "sex")

# === Initialize result list ===
interaction_results <- list()

# === Loop over all variable pairs ===
for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    v1 <- vars[i]
    v2 <- vars[j]
    
    # Prepare data from clinical data frame
    df <- clin_df_all %>%
      filter(!is.na(OS_time), !is.na(OS_event),
             !is.na(.data[[v1]]), !is.na(.data[[v2]]))
    
    # Drop unused levels for factor variables
    if (is.character(df[[v1]]) || is.factor(df[[v1]])) {
      df[[v1]] <- droplevels(as.factor(df[[v1]]))
    }
    if (is.character(df[[v2]]) || is.factor(df[[v2]])) {
      df[[v2]] <- droplevels(as.factor(df[[v2]]))
    }
    
    # Skip if too few samples or insufficient levels
    if (nrow(df) < 30) next
    if (is.factor(df[[v1]]) && length(levels(df[[v1]])) < 2) next
    if (is.factor(df[[v2]]) && length(levels(df[[v2]])) < 2) next
    
    # Feedback
    cat("Testing:", v1, "*", v2, "with", nrow(df), "samples\n")
    
    # Create model formulas
    formula_main <- as.formula(paste("Surv(OS_time, OS_event) ~", v1, "+", v2))
    formula_interact <- as.formula(paste("Surv(OS_time, OS_event) ~", v1, "*", v2))
    
    # Fit models safely
    fit1 <- suppressWarnings(try(coxph(formula_main, data = df), silent = TRUE))
    fit2 <- suppressWarnings(try(coxph(formula_interact, data = df), silent = TRUE))
    
    # Skip if models failed
    if (inherits(fit1, "try-error") || inherits(fit2, "try-error")) {
      warning(paste("Model failed for pair:", v1, "and", v2))
      next
    }
    
    # Likelihood ratio test
    lrt <- try(anova(fit1, fit2, test = "LRT"), silent = TRUE)
    if (inherits(lrt, "try-error")) {
      warning(paste("LRT failed for pair:", v1, "and", v2))
      next
    }
    
    # Safely extract LRT p-value from the last row
    pval <- tryCatch({
      p_col <- grep("Pr\\(>\\|Chi\\|\\)", colnames(lrt), value = TRUE)
      if (length(p_col) == 0) stop("No p-value column found in LRT result.")
      lrt[[p_col]][nrow(lrt)]
    }, error = function(e) {
      warning(paste("LRT result extraction failed for", v1, "*", v2))
      return(NA)
    })
    
    cat("→ p-value:", pval, "\n")
    interaction_results[[paste(v1, v2, sep = ":")]] <- pval
    
  }
  print(interaction_results)
  length(interaction_results)  # Should be 6 if all ran
  str(interaction_results)     # Should show named entries with p-values or NAs
  
}

# === Summarize and print results ===
if (length(interaction_results) > 0) {
  synergy_df <- data.frame(
    interaction = names(interaction_results),
    p_value = unlist(interaction_results)
  ) %>%
    arrange(p_value)
  
  print(synergy_df)
  
  # Optional: filter significant interactions
  cat("\nSignificant interactions (p < 0.05):\n")
  print(synergy_df %>% filter(p_value < 0.05))
} else {
  message("No valid interaction tests were performed.")
}

# filter for significant synergies
synergy_df %>% filter(p_value < 0.05)





# Visualize Interactions
library(survminer)

# Grouping example for age_group and tumor_stage
clin_df_all$group <- interaction(clin_df_all$age_group, drop = TRUE)

fit_km <- survfit(Surv(OS_time, OS_event) ~ group, data = clin_df_all)
ggsurvplot(fit_km, data = clin_df_all, pval = TRUE)



clin_df_all$group <- interaction(clin_df_all$tumor_stage, drop = TRUE)

fit_km <- survfit(Surv(OS_time, OS_event) ~ group, data = clin_df_all)
ggsurvplot(fit_km, data = clin_df_all, pval = TRUE)






results_df <- bind_rows(results)
results_df <- results_df %>% arrange(p)
head(results_df)

pvalues <- results_df$p 
qobj <- qvalue(p = pvalues)
results_df$adj_p <- qobj$qvalues
hist(qobj)


