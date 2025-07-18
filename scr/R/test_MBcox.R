#!/usr/bin/env Rscript
# 01_run_models_dev.R — can run proteomics‐only or transcriptomics‐only

args       <- commandArgs(trailingOnly=TRUE)
mode_arg   <- sub("^mode=",   "", args[grepl("^mode=",   args)])
cancer_arg <- sub("^cancer=", "", args[grepl("^cancer=", args)])

if (!mode_arg %in% c("prot","trans")) {
  stop("You must pass --args mode=prot or mode=trans")
}
if (!cancer_arg %in% c("LUAD","BRCA","GBM")) {
  stop("You must pass --args cancer=LUAD (or BRCA/GBM)")
}
message("→ mode = ", mode_arg, "; cancer = ", cancer_arg)

# point to your personal lib
.libPaths(c( Sys.getenv("R_LIBS_USER"), .libPaths() ))

suppressPackageStartupMessages({
  library(Coxmos); library(impute); library(dplyr)
  library(stringr); library(readr); library(data.table)
})

message("→ mode = ", mode_arg, "; cancer = ", cancer_arg)

# paths
clinical_file <- "data/processed/clin_df_all_cleaned.csv"
prot_dir      <- "data/analysis_results/probatch_diagnostics/batch_correction"
trans_dir     <- "data/processed/transcriptomics"
out_dir       <- "data/coxmos_results"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# load clinical
# load clinical
clin_df <- read_csv(clinical_file, show_col_types=FALSE)

# helper to one-hot + Y
build_clinical <- function(ids) {
  df <- clin_df %>%
    filter(case_id %in% ids, !is.na(OS_time), !is.na(OS_event)) %>%
    distinct(case_id, .keep_all=TRUE)
  # factor levels...
  df <- df %>%
    mutate(
      age_group = factor(age_group,
                         levels=c("<40","40–49","50–59","60–69","70–79","80+")),
      sex = factor(sex, levels=c("female","male")),
      tumor_stage = factor(tumor_stage_clean,
                           levels=c("stage i","stage ii","stage iii","stage iv","unknown")),
      hist_grade = factor(histologic_grade,
                          levels=c("G1","G2","G3","G4","GX","[unknown]"))
    )
  Xdf <- df %>%
    select(case_id, age_group, sex, tumor_stage, hist_grade) %>%
    column_to_rownames("case_id")
  mm  <- model.matrix(~ age_group + sex + tumor_stage + hist_grade - 1, data=Xdf)
  Ydf <- df %>% transmute(time=OS_time, event=OS_event) %>% as.data.frame()
  rownames(Ydf) <- df$case_id
  list(X=mm, Y=Ydf)
}

# —— a single preprocess_mat() for _all_ blocks ——
pattern_all <- "(?:[0-9]{2}[A-Z]{2}[0-9]{3}|C3[NL]-[0-9]{5}|NX[0-9]+|TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-[A-Za-z0-9]{2,4})?)"
cat("→ Using regex:\n", pattern_all, "\n\n"); flush.console()

preprocess_mat <- function(raw_mat, blk_name) {
  cat("   [", blk_name, "] dim:", dim(raw_mat), "\n"); flush.console()
  cn  <- colnames(raw_mat)
  ids <- str_extract(cn, pattern_all)
  if (any(is.na(ids)) && grepl("^PDC", blk_name)) {
    bad <- which(is.na(ids))
    cat("     ⚠️ fallback IDs at cols:", paste(bad,collapse=","), "\n"); flush.console()
    fallback <- sub("^.*_matrix_", "", cn[bad])
    fallback <- gsub(" ", "_", fallback)
    ids[bad] <- fallback
  }
  if (any(is.na(ids))) {
    stop("❌ Unparsed IDs in ", blk_name, ": ", paste(which(is.na(ids)), collapse=","), "\n")
  }
  colnames(raw_mat) <- ids
  m   <- t(raw_mat)
  keep<- colMeans(is.na(m)) < 0.8
  m2  <- m[, keep, drop=FALSE]
  orig<- rownames(m2)
  imp <- impute.knn(m2)$data
  rownames(imp) <- orig
  colnames(imp) <- colnames(m2)
  imp[, colnames(imp) != "NA.", drop=FALSE]
}

# —— build omics_blocks depending on mode —— 
omics_blocks <- list()
if (mode_arg == "prot") {
  # 1a) list PDC files, map to TCGA
  # — PROT blocks —
  cat("→ Listing all proteomics files...\n"); flush.console()
  pf <- list.files(prot_dir,
                           pattern = "^PDC[0-9]+_(gene|iso_log|iso_frac)_batch_corrected\\.csv$",
                           full.names = TRUE)
  cat("   Found", length(pf), "proteomics files\n"); flush.console()
  if (length(pf)>0) {
    cat("   file examples:", head(basename(pf),3), "…\n"); flush.console()
  }
  
  pdc_to_tcga <- c(
    PDC000110="OV",PDC000116="COAD",PDC000120="BRCA",PDC000125="UCEC",
    PDC000127="KIRC",PDC000153="LUAD",PDC000204="GBM",PDC000221="HNSC",
    PDC000234="LUSC",PDC000270="PAAD"
  )
  prot_blocks <- setNames(
    lapply(pf, function(f) {
      df <- fread(f, data.table=FALSE)
      rownames(df) <- df[[1]]
      mat <- as.matrix(df[,-1,drop=FALSE])
      preprocess_mat(mat, basename(f))
    }),
    sapply(basename(pf), function(f) {
      pdc  <- str_extract(f, "^PDC[0-9]+")
      type <- str_extract(f, "(gene|iso_log|iso_frac)")
      tcga <- pdc_to_tcga[[pdc]]
      paste0("PROT_", tcga, "_", type)
    })
  )
  omics_blocks <- prot_blocks
  
} else {
  # mode == "trans"
  tf <- list.files(trans_dir,
                   pattern=paste0("^RNA_",cancer_arg,"_(gene|iso_log|iso_frac)\\.csv$"),
                   full.names=TRUE)
  trans_blocks <- setNames(
    lapply(tf, function(f) {
      df <- fread(f, data.table = FALSE)
      # first column is your feature‐ID, so pull it off:
      rownames(df) <- df[[1]]
      df <- df[, -1, drop = FALSE]
      mat <- t(as.matrix(df))
    }),
    # name each block “TRANS_<cancer>_<type>”
    str_replace_all(basename(tf),
                    c("^RNA_" = "TRANS_", "\\.csv$" = ""))
  )
  trans_blocks <- lapply(trans_blocks, function(mat) {
    clean_ids    <- sub("-[0-9]{2}$", "", rownames(mat))
    rownames(mat) <- clean_ids
    mat
  })
  
  omics_blocks <- trans_blocks
}
message("→ Found ---- ", paste0(length(pf)), " ---- proteomics files", "\n", "File Names:\n   ", paste(basename(pf), collapse=", "))
#message("→ Found ---- ", paste0(length(tf)), " ---- transcriptomics files", "\n", "File Names:\n   ", paste(basename(tf), collapse=", "))


# now restrict to this cancer
blocks_for <- omics_blocks[grep(paste0("^(PROT|TRANS)_",cancer_arg,"_"), names(omics_blocks))]
message("→ Using ", length(blocks_for), " blocks for cancer ", cancer_arg, ", ", mode_arg, " data")


# Filter each molecular block to only those cases that actually appear in the clinical metadata
blocks_for <- lapply(blocks_for, function(mat) {
  # keep only rows whose IDs are in the global clinical df:
  keep_ids <- intersect(rownames(mat), clin_df$case_id)
  mat[keep_ids, , drop=FALSE]
})

# now 1) get union of all samples in your processed (and pre-filtered) blocks
all_block_ids <- unique(unlist(lapply(blocks_for, rownames)))

# 2) restrict global clinical metadata to those samples
clin_df_sub <- clin_df %>%
  filter(case_id %in% all_block_ids)

# 3) build samp_anno for modelling…
samp_anno <- clin_df_sub %>%
  filter(!is.na(OS_time), !is.na(OS_event)) %>%
  distinct(case_id, .keep_all = TRUE) %>%
  mutate(
    age_group   = factor(age_group,   levels = c("<40","40–49","50–59","60–69","70–79","80+")),
    sex         = factor(sex,         levels = c("female","male")),
    tumor_stage = factor(tumor_stage_clean,
                         levels = c("stage i","stage ii","stage iii","stage iv","unknown")),
    hist_grade  = factor(histologic_grade,
                         levels = c("G1","G2","G3","G4","GX","[unknown]"))
  )

message("→ Total unique samples across blocks_for: ", length(all_block_ids))
message("→ Clinical matrix dims (X) = ", nrow(samp_anno), " × ", ncol(samp_anno))

# 4) extract the final, filtered list of case IDs:
common_ids <- samp_anno$case_id

message("→ Samples in clinical: ", length(common_ids),
        "; samples in each block will also be restricted to these.")

# 5) subset each block again to exactly those common IDs (and in the same order)
blocks_for <- lapply(blocks_for, function(mat) {
  mat2 <- mat[common_ids, , drop=FALSE]
  # sanity check:
  if (!all(rownames(mat2) == common_ids))
    stop("Mismatch in ordering after subsetting molecular block")
  mat2
})

# 6) now build your one-hot clinical matrix (exactly in the same order)
X_clin_df <- samp_anno %>%
  select(case_id, age_group, sex, tumor_stage, hist_grade) %>%
  tibble::column_to_rownames("case_id") %>%
  # ensure the rows are in the same order
  .[common_ids, ]

X_clin_onehot <- factorToBinary(X = X_clin_df, all = TRUE, sep = "_")
X_clin_num    <- as.matrix(mutate(as.data.frame(X_clin_onehot),
                                  across(everything(), as.integer)))

# 7) assemble your multiblock list
X_blocks <- c(blocks_for, list(clinical = X_clin_num))

# 8) build your Y in the same order
Y <- samp_anno %>%
  transmute(time = OS_time, event = OS_event) %>%
  as.data.frame()
rownames(Y) <- samp_anno$case_id
Y <- Y[common_ids, , drop=FALSE]

message("→ Total unique samples across blocks_for: ", length(common_ids))
message("→ Dimensions of each block in X_blocks:")
for (nm in names(X_blocks)) {
  d <- dim(X_blocks[[nm]])
  message(sprintf("   %-20s : %d × %d", nm, d[1], d[2]))
}

EPV <- getEPV.mb(X_blocks, Y)
cat("EPV")


# split once
sp   <- getTrainTest(X_blocks, Y, p=0.7, seed=1234)
Xtr  <- sp$X_train; Ytr <- sp$Y_train
Xte  <- sp$X_test;  Yte <- sp$Y_test

message("→ getTrainTest produced the following train/test sizes per block:")
for (blk in names(sp$X_train)) {
  n_tr <- nrow(sp$X_train[[blk]])
  n_te <- nrow(sp$X_test[[blk]])
  message(sprintf("   %-20s : train = %3d  /  test = %3d", blk, n_tr, n_te))
}
message("   Y                    : train =  ", nrow(sp$Y_train), "  /  test =  ",  nrow(sp$Y_test))


# fit one Coxmos model, for brevity
cvb <- cv.mb.coxmos("sb.splsicox", Xtr, Ytr, MIN_NVAR=5, MAX_NVAR=20,
                    max.ncomp=2, n_run=3, k_folds=5, MIN_EPV=0.01,
                    remove_zero_variance=TRUE,
                    remove_variance_at_fold_level = TRUE,
                    remove_non_significant = FALSE,
                    remove_near_zero_variance  = TRUE, 
                    PARALLEL=FALSE,)

mb  <- mb.coxmos("sb.splsicox", Xtr, Ytr,
                 n.comp=cvb$opt.comp, penalty=cvb$opt.penalty,
                 remove_non_significant=FALSE)
res <- list(cvb, mb)

# save
fn <- sprintf("%s_%s_%s.rds", mode_arg, cancer_arg, Sys.Date())
saveRDS(res, file.path(out_dir, fn))
message("→ saved to ", file.path(out_dir, fn))

message("→ Finished.  Full log is at: ", logfile)