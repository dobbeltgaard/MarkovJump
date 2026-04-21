library(here)
here::i_am("evaluation.R")

rm(list = ls()); gc() #clear memory
library(Rcpp)
library(RcppEigen)
library(dplyr); library(kableExtra)
sourceCpp("FUNCS_MJP_with_eigen.cpp")


files <- list.files("predictions_V3", pattern = "\\.csv$", full.names = TRUE)
file_names <- basename(files)
base_names <- sub("_fold_\\d+\\.csv$", "", file_names)
fold_numbers <- as.integer(sub(".*_fold_(\\d+)\\.csv$", "\\1", file_names))
file_info <- data.frame(file = files,name = file_names,base = base_names,fold = fold_numbers,stringsAsFactors = FALSE)
file_info <- file_info[order(file_info$base, file_info$fold), ]
predictor_list <- split(file_info, file_info$base)
predictor_data <- lapply(predictor_list, function(group) {
  out <- lapply(seq_len(nrow(group)), function(i) {df <- read.csv(group$file[i]); df;})
  names(out) <- paste0("fold_", group$fold)
  out
})
pred_list = predictor_data

#Make ensemble prediction
predictor_names <- names(pred_list)[ (grepl("_cov2_",names(pred_list)) & grepl("_no_warp",names(pred_list))) |  ( (grepl("_cov0_",names(pred_list)) | grepl("_cov1_",names(pred_list))) & !grepl("mixture",names(pred_list)) ) | grepl("mixture4",names(pred_list)) | grepl("ocllr", names(pred_list)) | grepl("olr", names(pred_list)) | grepl("opr", names(pred_list)) | grepl("uniform", names(pred_list)) | grepl("random_", names(pred_list)) ]
pred_list[["ensemble_all"]] <- vector("list", length = length(pred_list[[ predictor_names[1] ]]))
for(i in seq_along(pred_list[["ensemble_all"]])) {
  preds_i <- lapply(predictor_names, function(name) pred_list[[name]][[i]])
  pred_list[["ensemble_all"]][[i]] <- Reduce("+", preds_i) / length(preds_i)
}

#######################
### Forecast Scores ###
#######################
### Cross validation errors ###
m = 5; k = 5; 
log_err = matrix(NA, ncol = k, nrow = length(names(pred_list)))
RPS_err = matrix(NA, ncol = k, nrow = length(names(pred_list)))
Brier_e = matrix(NA, ncol = k, nrow = length(names(pred_list)))
for(i in 1:k){
  count = 0; 
  for(j in names(pred_list)){
    count = count + 1;
    pred = as.matrix(pred_list[[j]][[i]])
    obs = as.matrix(pred_list[["obs"]][[i]])
    log_err[count, i] = -mean(logscore_vectors(m, pred, obs))
    RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
    Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
  }
}
rownames(log_err) = names(pred_list)
rownames(RPS_err) = names(pred_list)
rownames(Brier_e) = names(pred_list)
errs = cbind(rowMeans(RPS_err),rowMeans(log_err))
errs[grepl("_all_", rownames(errs)) | !grepl("_exp", rownames(errs)),1:2]

mk_summary <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sd <- apply(mat, 1, sd, na.rm = TRUE)
  sprintf("%.4f (%.4f)", mu, sd)
}
errs_fmt <- cbind(
  RPS  = mk_summary(RPS_err),
  LogS = mk_summary(log_err)
)
rownames(errs_fmt) <- rownames(RPS_err)
errs_fmt[grepl("_all_", rownames(errs_fmt)) | !grepl("_exp", rownames(errs_fmt)), , drop = FALSE]

#namfoo = c("exp_exp_all", "uniform", "olr", "opr", "ocllr", "random_", "ensemble")
#idx = rowSums(sapply(namfoo, FUN = grepl, x = rownames(errs_fmt))) > 0
#foo = errs_fmt[idx, ]
# nams =
#   c("uniform", 
#     "olr", "olr_cov", "opr", "opr_cov", "ocllr", "ocllr_cov",
#     "random_forest", "random_forest_cov",
#     "gerlang_cov0_exp_exp_all_no_warp", "gerlang_cov1_exp_exp_all_no_warp", 
#     "gerlang_relax_cov0_exp_exp_all_no_warp", "gerlang_relax_cov1_exp_exp_all_no_warp", 
#     "bidiagonal_cov0_exp_exp_all_no_warp", "bidiagonal_cov1_exp_exp_all_no_warp", 
#     "tridiagonal_cov0_exp_exp_all_no_warp", "tridiagonal_cov1_exp_exp_all_no_warp", 
#     "free_upper_tri_cov0_exp_exp_all_no_warp", "free_upper_tri_cov1_exp_exp_all_no_warp",
#     "gerlang_cov0_exp_exp_all_warp", "gerlang_cov1_exp_exp_all_warp", "gerlang_cov2_exp_exp_all_no_warp",
#     "gerlang_relax_cov0_exp_exp_all_warp", "gerlang_relax_cov1_exp_exp_all_warp", "gerlang_relax_cov2_exp_exp_all_no_warp",
#     "bidiagonal_cov0_exp_exp_all_warp", "bidiagonal_cov1_exp_exp_all_warp", "bidiagonal_cov2_exp_exp_all_no_warp",
#     "tridiagonal_cov0_exp_exp_all_warp", "tridiagonal_cov1_exp_exp_all_warp", "tridiagonal_cov2_exp_exp_all_no_warp",
#     "free_upper_tri_cov0_exp_exp_all_warp", "free_upper_tri_cov1_exp_exp_all_warp", "free_upper_tri_cov2_exp_exp_all_no_warp",
#     "ensemble_all"
#     )

nams_ref = c("uniform","olr", "olr_cov", "opr", "opr_cov", "ocllr", "ocllr_cov", "DTMC")#, "nnet", "nnet_cov")#,"random_forest", "random_forest_cov")
nams_MJP = 
  c("gerlang_cov0_exp_exp_all_no_warp", "gerlang_cov1_exp_exp_all_no_warp","gerlang_cov2_exp_exp_all_no_warp", 
             "gerlang_relax_cov0_exp_exp_all_no_warp", "gerlang_relax_cov1_exp_exp_all_no_warp","gerlang_relax_cov2_exp_exp_all_no_warp",
             "bidiagonal_cov0_exp_exp_all_no_warp", "bidiagonal_cov1_exp_exp_all_no_warp","bidiagonal_cov2_exp_exp_all_no_warp", 
             "tridiagonal_cov0_exp_exp_all_no_warp", "tridiagonal_cov1_exp_exp_all_no_warp","tridiagonal_cov2_exp_exp_all_no_warp", 
             "free_upper_tri_cov0_exp_exp_all_no_warp", "free_upper_tri_cov1_exp_exp_all_no_warp", "free_upper_tri_cov2_exp_exp_all_no_warp")
nams_MJP_warp = 
  c("gerlang_cov0_exp_exp_all_warp", "gerlang_cov1_exp_exp_all_warp",
             "gerlang_relax_cov0_exp_exp_all_warp", "gerlang_relax_cov1_exp_exp_all_warp",
             "bidiagonal_cov0_exp_exp_all_warp", "bidiagonal_cov1_exp_exp_all_warp",
             "tridiagonal_cov0_exp_exp_all_warp", "tridiagonal_cov1_exp_exp_all_warp",
             "free_upper_tri_cov0_exp_exp_all_warp", "free_upper_tri_cov1_exp_exp_all_warp")
nams_MJP_mixture = 
  c("gerlang_cov0_exp_exp_all_no_warp_mixture4", "gerlang_cov1_exp_exp_all_no_warp_mixture4", 
    "gerlang_relax_cov0_exp_exp_all_no_warp_mixture4", "gerlang_relax_cov1_exp_exp_all_no_warp_mixture4", 
    "bidiagonal_cov0_exp_exp_all_no_warp_mixture4", "bidiagonal_cov1_exp_exp_all_no_warp_mixture4", 
    "tridiagonal_cov0_exp_exp_all_no_warp_mixture4", "tridiagonal_cov1_exp_exp_all_no_warp_mixture4", 
    "free_upper_tri_cov0_exp_exp_all_no_warp_mixture4", "free_upper_tri_cov1_exp_exp_all_no_warp_mixture4")
nams_ensemble = c("ensemble_all")

#nams_idx = rep(NA, length(nams));for(i in 1:length(nams)){nams_idx[i] = which(nams[i] == rownames(foo))} 
library(kableExtra)
naive = expand.grid(
  Covariates = c("No"), #covariates
  Model = c("Naïve predictor"), #näive,olr,mjp 
  param = c("Uniform distribution")
)
regression = expand.grid(
  Covariates = c("No","Yes"), #covariates
  Model = c("Cumulative link model"), #näive,olr,mjp 
  param = c("Logistic link", "Probit link", "Cloglog link")
)
forest = expand.grid(
  Covariates = c("No"),
  Model = c("Discrete-time Markov chain"),
  param = c("")
)
mjp = expand.grid(
  Covariates = c("No","Yes (global)", "Yes (state-spec.)"), #covariates
  Model = c("Markov jump process"), #näive,olr,mjp 
  param = c("$\\bm{A}_{\\RN{1}}$: Nearest class", "$\\bm{A}_{\\RN{2}}$: Relaxed Erlang", "$\\bm{A}_{\\RN{3}}$: Up to 2 classes ahead","$\\bm{A}_{\\RN{4}}$: Up to 3 classes ahead", "$\\bm{A}_{\\RN{5}}$: Free upper-triangular")
)
mjp_warp = expand.grid(
  Covariates = c("No","Yes (global)"), #covariates
  Model = c("Model with modified sojourn time estimates"), #näive,olr,mjp 
  param = c("$\\bm{A}_{\\RN{1}}$: Nearest class", "$\\bm{A}_{\\RN{2}}$: Relaxed Erlang", "$\\bm{A}_{\\RN{3}}$: Up to 2 classes ahead","$\\bm{A}_{\\RN{4}}$: Up to 3 classes ahead", "$\\bm{A}_{\\RN{5}}$: Free upper-triangular")
)
mjp_mixture = expand.grid(
  Covariates = c("No","Yes (global)"), #covariates
  Model = c("Mixture Markov jump process ($G=4$)"), #näive,olr,mjp 
  param = c("$\\bm{A}_{\\RN{1}}$: Nearest class", "$\\bm{A}_{\\RN{2}}$: Relaxed Erlang", "$\\bm{A}_{\\RN{3}}$: Up to 2 classes ahead","$\\bm{A}_{\\RN{4}}$: Up to 3 classes ahead", "$\\bm{A}_{\\RN{5}}$: Free upper-triangular")
)
ensemble = expand.grid(
  Covariates = c("Mixed"), #covariates
  Model = c("Ensemble"), #näive,olr,mjp 
  param = c("All")
)
#collapse_rows_dt = rbind(naive, regression, forest, mjp,mjp_warp ,ensemble)
#collapse_rows_dt <- collapse_rows_dt[c("Model", "param", "Covariates")]

tab_ref = rbind(naive, regression,forest); tab_ref = tab_ref[c("Model", "param", "Covariates")]; 
tab_MJP = rbind(mjp); tab_MJP = tab_MJP[c("Model", "param", "Covariates")]; 
tab_MJP_warp = rbind(mjp_warp); tab_MJP_warp = tab_MJP_warp[c("Model", "param", "Covariates")]; 
tab_MJP_mixture = rbind(mjp_mixture); tab_MJP_mixture = tab_MJP_mixture[c("Model", "param", "Covariates")]; 
tab_ensemble = rbind(ensemble); tab_ensemble = tab_ensemble[c("Model", "param", "Covariates")];

# #par_list  = pred_list
# npars = rep(0,length(nams)); count = 0;
# for(i in head(nams,-1)){
#   count = count + 1
#     dfoo = read.csv(list.files("estimates_V3", full.names = T)[grepl(paste0(i,"_fold"), list.files("estimates_V3"))][1])
#     npars[count] = NROW(dfoo)#(NCOL(par_list[[i]])); 
# }
# npars[npars == 1 | npars == 0] = "--"

get_npars <- function(nams, dir = "estimates_V3", skip_last = FALSE) {
  if (skip_last) { nams <- head(nams, -1)}
  files <- list.files(dir, full.names = TRUE)
  file_names <- basename(files)
  out <- sapply(nams, function(nm) {
    matches <- files[grepl(paste0("^", nm, "_fold"), file_names)]
    if (length(matches) == 0) {return(NA_character_)}
    dfoo <- read.csv(matches[1])
    npar <- nrow(dfoo)
    if (npar %in% c(0, 1,NA, "NA")) "--" else as.character(npar)
  }, USE.NAMES = TRUE)
  unname(out)
}

npars_MJP = get_npars(nams = nams_MJP); npars_MJP_warp = get_npars(nams = nams_MJP_warp); npars_MJP_mixture = get_npars(nams = nams_MJP_mixture); npars_ref = get_npars(nams = nams_ref); npars_ensemble = c("--")


fill_tab <- function(tab, nams, errs_fmt, dir = "estimates_V3", npars_bin = T) {
  stopifnot(nrow(tab) == length(nams))
  npars <- get_npars(nams, dir = dir)
  idx <- match(nams, rownames(errs_fmt))
  if (anyNA(idx)) {
    warning("Some model names were not found in errs_fmt: ",
            paste(nams[is.na(idx)], collapse = ", "))
  }
  if(npars_bin) tab$npars <- npars
  tab$RPS   <- errs_fmt[idx, "RPS"]
  tab$LogS  <- errs_fmt[idx, "LogS"]
  tab
}

tab_ref         <- fill_tab(tab_ref,         nams_ref,         errs_fmt)
tab_MJP         <- fill_tab(tab_MJP,         nams_MJP,         errs_fmt)
tab_MJP_warp    <- fill_tab(tab_MJP_warp,    nams_MJP_warp,    errs_fmt)
tab_MJP_mixture <- fill_tab(tab_MJP_mixture, nams_MJP_mixture, errs_fmt)
tab_ensemble    <- fill_tab(tab_ensemble,    nams_ensemble,    errs_fmt); tab_ensemble$npars= "--"




colnams = c("Method", "Model", "Exogenous info.", "\\# Pars.", "RPS", "Log S.")
#colnams = c("Method", "Model", "Exogenous info.", "RPS", "Log S.")


# 
# colnames(tab_ref) = colnams
# row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
# n_cols = length(colnams)
# kableExtra::kbl(tab_ref,booktabs = T, align = c("l","l","l","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5) %>%
#   column_spec(1, bold=T) %>%
#   collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
#   row_spec(3,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(5,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   footnote(number = c(
#     "Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included."),
#     escape = FALSE, general_title = "", threeparttable = TRUE) %>%
#   writeLines(con = paste0("figures/forecast_error_tab_ref",".tex"))
# 
# colnames(tab_MJP) = colnams
# row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
# n_cols = length(colnams)
# kableExtra::kbl(tab_MJP,booktabs = T, align = c("l","l","l","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5) %>%
#   column_spec(1, bold=T) %>%
#   collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
#   row_spec(3,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(6,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(9, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(12, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   footnote(number = c(
#     "Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}; Yes (state-spec.) uses \\\\cref{eq:trans_rates_covariates2}."),
#     escape = FALSE, general_title = "", threeparttable = TRUE) %>%
#   writeLines(con = paste0("figures/forecast_error_tab_MJP",".tex"))
# 
# colnames(tab_MJP_warp) = colnams
# row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
# n_cols = length(colnams)
# kableExtra::kbl(tab_MJP_warp,booktabs = T, align = c("l","l","l","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5) %>%
#   column_spec(1, bold=T) %>%
#   collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
#   row_spec(2,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(4,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(6, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(8, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   footnote(number = c(
#     "Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}."),
#     escape = FALSE, general_title = "", threeparttable = TRUE) %>%
#   writeLines(con = paste0("figures/forecast_error_tab_MJP_warp",".tex"))
# 
# colnames(tab_MJP_mixture) =colnams
# row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
# n_cols = length(colnams)
# kableExtra::kbl(tab_MJP_mixture,booktabs = T, align = c("l","l","l","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5) %>%
#   column_spec(1, bold=T) %>%
#   collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
#   row_spec(2,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(4,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(6, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   row_spec(8, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
#   footnote(number = c(
#     "Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}."),
#     escape = FALSE, general_title = "", threeparttable = TRUE) %>%
#   writeLines(con = paste0("figures/forecast_error_tab_MJP_mixture",".tex"))


tab_tot = rbind(tab_ref,tab_MJP,tab_MJP_warp,tab_MJP_mixture,tab_ensemble)


colnames(tab_tot) =colnams
row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
n_cols = length(colnams)
tab_tex = kableExtra::kbl(tab_tot,booktabs = T, align = c("l","l","l","c","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5, longtable = T, label = "prediction_scores",
                          caption = "Prediction scores from $k$-fold cross-prediction procedures in terms of average ranked probability score and average log score, with standard deviations across folds in parentheses. ``Exogenous info.'' indicates whether the rail characteristics described in \\Cref{sec:data} are included: Yes (global) uses \\cref{eq:trans_rates_covariates}; Yes (state-spec.) uses \\cref{eq:trans_rates_covariates2}. ``\\# Pars.'' indicates the number of estimated model parameters. Bold numbers indicate best performance within a model class (where applicable).", 
                          ) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
  row_spec(3, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(5, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(11, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(14, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(17, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(20, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(25, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(27, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(29, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(31, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(35, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(37, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(39, extra_latex_after = "\n\\addlinespace[0.5em]")   %>%
  row_spec(41, extra_latex_after = "\n\\addlinespace[0.5em]") %>%
  kable_styling(latex_options = c("repeat_header"),
                #repeat_header_continued = "\\textit{(Continued on next page...)}", 
                repeat_header_text = "\\textit{(continued)}", 
                repeat_header_method = "replace") 
tab_tex <- as.character(tab_tex)
tab_tex <- gsub("\\\\pagebreak\\[0\\]", "", tab_tex)
tab_tex <- gsub("\\\\\\\\(\\s*\n)","\\\\\\\\*\\1",tab_tex,perl = TRUE)
writeLines(tab_tex, "figures/forecast_error_tab_tot.tex")


  # row_spec(4,  extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  # row_spec(6, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  # row_spec(8, extra_latex_after = sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  #footnote(number = c("Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}."),escape = FALSE, general_title = "", threeparttable = TRUE) %>%

# tab_tex <- as.character(tab_tex)
# lines <- strsplit(tab_tex, "\n", fixed = TRUE)[[1]]
# start <- grep("^\\\\toprule$", lines)[1] 
# end   <- tail(grep("^\\\\bottomrule$", lines), 1)
# body_lines <- lines[start:end]
# writeLines(body_lines, "figures/forecast_error_tab_tot_rows.tex", useBytes = TRUE)




# colnames(tab_ensemble) = colnams
# row_group_label_fonts <- list(list(bold = T, italic = F),list(bold = F, italic = F))
# n_cols = 6
# kableExtra::kbl(tab_ensemble,booktabs = T, align = c("l","l","l","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 5) %>%
#   column_spec(1, bold=T) %>%
#   collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
#   writeLines(con = paste0("figures/forecast_error_tab_ensemble",".tex"))





##############################
### ESTIMATION DIAGNOSTICS ###
##############################

files <- list.files("diagnostics_V3", pattern = "\\.csv$", full.names = TRUE)
#files = files[!grepl("v2",files)]
file_names <- basename(files)
base_names <- sub("_fold_\\d+\\.csv$", "", file_names)
fold_numbers <- as.integer(sub(".*_fold_(\\d+)\\.csv$", "\\1", file_names))
file_info <- data.frame(file = files,name = file_names,base = base_names,fold = fold_numbers,stringsAsFactors = FALSE)
file_info <- file_info[order(file_info$base, file_info$fold), ]
diagnostics_list <- split(file_info, file_info$base)
diagnostics_data <- lapply(diagnostics_list, function(group) {
  out <- lapply(seq_len(nrow(group)), function(i) {df <- read.csv(group$file[i]); df;})
  names(out) <- paste0("fold_", group$fold)
  out
})

diag_df <- do.call(rbind, lapply(names(diagnostics_data), function(model) {
  folds_list <- diagnostics_data[[model]]
  do.call(rbind, lapply(names(folds_list), function(fold_name) {
    df <- folds_list[[fold_name]]
    df$model <- model
    df$fold  <- as.integer(sub("fold_", "", fold_name))
    df
  }))
}))
diag_df <- diag_df[, c("model", "fold", "convergence", "grad_norm", "AIC")]
diag_df <- diag_df[order(diag_df$model, diag_df$fold), ]

#187       gerlang_cov1_exp_exp_all_no_warp_mixture4    2           1 8.958703e-04 5641.536
#260 gerlang_relax_cov1_exp_exp_all_no_warp_mixture4    5           1 1.192478e-04 5673.201
#296   tridiagonal_cov0_exp_exp_all_no_warp_mixture4    1           1 6.085386e-05 5758.775


diag_df$fold <- as.integer(diag_df$fold)
AIC_w <- with(diag_df, tapply(AIC, list(model, fold), identity))
AIC_w <- AIC_w[, order(as.integer(colnames(AIC_w))), drop = FALSE]
colnames(AIC_w) <- paste0("AICf", colnames(AIC_w))
GN_w <- with(diag_df, tapply(grad_norm, list(model, fold), identity))
GN_w <- GN_w[, order(as.integer(colnames(GN_w))), drop = FALSE]
colnames(GN_w) <- paste0("gradf", colnames(GN_w))
folds <- sort(unique(diag_df$fold))
AIC_GN_tab <- cbind(AIC_w, GN_w)
AIC_GN_tab <- AIC_GN_tab[, as.vector(rbind(paste0("AICf", folds),paste0("gradf", folds))), drop = FALSE]

nams2 = c(nams_MJP, nams_MJP_warp,nams_MJP_mixture )
missing <- setdiff(nams2, rownames(AIC_GN_tab))
if(length(missing) > 0) warning("These names in nams2 are not in AIC_GN_tab: ", paste(missing, collapse=", "))
keep <- intersect(nams2, rownames(AIC_GN_tab))
AIC_GN_tab <- AIC_GN_tab[keep, , drop = FALSE]
param_map <- c(
  gerlang = "$\\bm{A}_{\\RN{1}}$: Nearest class",
  gerlang_relax = "$\\bm{A}_{\\RN{2}}$: Relaxed Erlang",
  bidiagonal = "$\\bm{A}_{\\RN{3}}$: Up to 2 classes ahead",
  tridiagonal = "$\\bm{A}_{\\RN{4}}$: Up to 3 classes ahead",
  free_upper_tri = "$\\bm{A}_{\\RN{5}}$: Free upper-triangular"
)
cov_map <- c(
  cov0 = "No",
  cov1 = "Yes (global)",
  cov2 = "Yes (state-spec.)"
)
model_map <- c(
  no_warp = "Markov jump process",
  warp    = "\\makecell[l]{\\textbf{Model with modified}\\\\ \\textbf{sojourn time estimates}}",
  mixture = "Mixture Markov jump process ($G=4$)"
)
rn <- rownames(AIC_GN_tab)
gen <- sub("_cov\\d+_.*$", "", rn)
cov <- sub("^.*_(cov\\d)_.*$", "\\1", rn)
model_type <- ifelse(
  grepl("_mixture\\d+$", rn), "mixture",
  ifelse(grepl("_no_warp$", rn), "no_warp",
         ifelse(grepl("_warp$", rn), "warp", NA))
)
pretty_df <- data.frame(
  Model = unname(model_map[model_type]),
  param = unname(param_map[gen]),
  Covariates = unname(cov_map[cov]),
  stringsAsFactors = FALSE
)
AIC_GN_df <- cbind(pretty_df, data.frame(AIC_GN_tab, check.names = FALSE))

AIC_GN_df_fmt <- AIC_GN_df
num_cols <- names(AIC_GN_df_fmt)[sapply(AIC_GN_df_fmt, is.numeric)]
grad_cols <- grep("grad_norm|GN|grad", num_cols, value = TRUE, ignore.case = TRUE)
aic_cols  <- setdiff(num_cols, grad_cols)
fmt_sci <- function(x, digits = 2) formatC(x, format = "e", digits = digits)
fmt_fix <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

AIC_GN_df_fmt[grad_cols] <- lapply(AIC_GN_df_fmt[grad_cols], fmt_sci, digits = 2)
AIC_GN_df_fmt[aic_cols]  <- lapply(AIC_GN_df_fmt[aic_cols],  fmt_fix, digits = 1)

rownames(AIC_GN_df_fmt) = 1:NROW(AIC_GN_df_fmt)
row_group_label_fonts <- list(
  list(bold = T, italic = F),
  list(bold = F, italic = F)
)

n_cols = 13
colnames(AIC_GN_df_fmt) = c(c("Method", "Model", "Exogenous info."),rep(c("AIC", "$\\frac{1}{n}||\\nabla \\mathcal{D} (\\bm \\theta)||_2$"), 5))
kableExtra::kbl(AIC_GN_df_fmt,booktabs = T, linesep = '', format = "latex",escape = FALSE, 
                align = c("l","l","l","c","c","c","c","c","c","c","c", "c", "c")) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
  footnote(number = c("Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}; Yes (state-spec.) uses \\\\cref{eq:trans_rates_covariates2}."),escape = FALSE, general_title = "", threeparttable = TRUE) %>%
  row_spec(3,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(6,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(9,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(12,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(17,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(19,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(21,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(23,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(27,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(29,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(31,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(33,  extra_latex_after = "\n\\addlinespace[0.5em]")%>%# sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  add_header_above(c(" " = 3,"Fold 1" = 2,"Fold 2" = 2,"Fold 3" = 2,"Fold 4" = 2,"Fold 5" = 2),escape = FALSE) %>%
  writeLines(con = paste0("figures/AIC_tab",".tex"))

#row_spec(3, extra_latex_after = "\n\\addlinespace[0.5em]")

###########################
### TRANSITION PATTERNS ###
###########################
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
set.seed(1); k= 5; bin_km= 0.1
d$pos_bin=floor(d$pos / bin_km); d$group_id=interaction(d$Track0, d$pos_bin, drop = TRUE); 
groups=levels(d$group_id); G=length(groups)
fold_id=sample(rep(1:k, length.out = G)); names(fold_id)=groups; fold_size=integer(k)

N_obs_by_fold <- vector("list", k)
names(N_obs_by_fold) <- paste0("fold_", 1:k)
N_exp_by_model <- setNames(vector("list", length(pred_list)), names(pred_list))
for(j in names(pred_list)) {
  N_exp_by_model[[j]] <- vector("list", k)
  names(N_exp_by_model[[j]]) <- paste0("fold_", 1:k)
}
for(i in 1:k){
  test_grps <- names(fold_id)[fold_id == i]
  pred.idx  <- which(d$group_id %in% test_grps)
  d.test    <- d[pred.idx, , drop = FALSE]
  s1i <- d.test$s1
  N_obs_by_fold[[i]] <- with(d.test, table(factor(s1, levels = states),factor(s2, levels = states)))
  for(j in names(pred_list)){
    pred <- as.matrix(pred_list[[j]][[i]])
    Nhat <- matrix(0, nrow = m, ncol = m,dimnames = list(from = states, to = states))
    for(a in states){
      idx <- which(s1i == a)
      if(length(idx) > 0){Nhat[a, ] <- colSums(pred[idx, , drop = FALSE])}
    }
    N_exp_by_model[[j]][[i]] <- Nhat
  }
}

pool_counts <- function(N_list) Reduce(`+`, lapply(N_list, as.matrix))
N_obs_pool <- pool_counts(N_obs_by_fold)
N_exp_pool_by_model <- lapply(N_exp_by_model, pool_counts)
row_norm_diff <- function(N_obs, N_hat, eps = 1e-12) {
  row_tot <- rowSums(N_obs)
  sweep(N_obs - N_hat, 1, pmax(row_tot, eps), "/")
}
R_pool_by_model <- lapply(N_exp_pool_by_model, function(Nhat_pool) {
  row_norm_diff(N_obs_pool, Nhat_pool)
})

P_obs_pool <- sweep(N_obs_pool, 1, pmax(rowSums(N_obs_pool), 1), "/")
P_obs_pool <- round(P_obs_pool, 3)

fmt_transition_table <- function(mat, states_nams, digits = 3) {
  n <- nrow(mat)
  out <- matrix("0", nrow = n, ncol = n)
  upper_idx <- upper.tri(mat, diag = TRUE)
  out[upper_idx] <- formatC(mat[upper_idx], format = "f", digits = digits)
  out <- cbind(Class = states_nams, out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- c("Class", states_nams)
  out
}

m1 = "bidiagonal_cov1_exp_exp_all_no_warp"
m2 = "bidiagonal_cov2_exp_exp_all_no_warp"
m3 = "bidiagonal_cov1_exp_exp_all_no_warp_mixture4"
states_nams = c("3", "2B", "2A", "1", "0")


P_tbl <- fmt_transition_table(P_obs_pool, states_nams, digits = 3)
kableExtra::kbl(P_tbl, booktabs = TRUE, linesep = '', format = "latex", escape = FALSE, align = "c|ccccc") %>% writeLines(con = paste0("figures/transition_matrix",".tex"))

R_bl <- fmt_transition_table(R_pool_by_model[[m1]], states_nams, digits = 3)
kableExtra::kbl(R_bl, booktabs = TRUE, linesep = '', format = "latex", escape = FALSE, align = "c|ccccc") %>% writeLines(con = paste0("figures/", "expected_transition_", m1, ".tex"))

R_bl <- fmt_transition_table(R_pool_by_model[[m2]], states_nams, digits = 3)
kableExtra::kbl(R_bl, booktabs = TRUE, linesep = '', format = "latex", escape = FALSE, align = "c|ccccc") %>% writeLines(con = paste0("figures/", "expected_transition_", m2, ".tex"))

R_bl <- fmt_transition_table(R_pool_by_model[[m3]], states_nams, digits = 3)
kableExtra::kbl(R_bl, booktabs = TRUE, linesep = '', format = "latex", escape = FALSE, align = "c|ccccc") %>% writeLines(con = paste0("figures/", "expected_transition_", m3, ".tex"))


L1 <- sapply(R_pool_by_model, function(R) {
  R2 <- R[-m, -m, drop = FALSE]
  sum(abs(R2), na.rm = TRUE)
})

nams_all <- c(nams_ref, nams_MJP, nams_MJP_warp, nams_MJP_mixture)


#Predictive entropy (sharpness)
entropy_rows <- function(P, eps = 1e-12) {P <- pmax(P, eps); -rowSums(P * log(P))}

mean_entropy_fold <- function(P, eps = 1e-12) {  mean(entropy_rows(P, eps = eps))}
mean_entropy_fold_norm <- function(P, m, eps = 1e-12) {mean_entropy_fold(P, eps = eps) / log(m)}
fold_wts <- sapply(seq_len(k), function(i) nrow(pred_list[[names(pred_list)[1]]][[i]]))
fold_wts <- fold_wts / sum(fold_wts)

H_mean_by_model <- sapply(names(pred_list), function(mod) {
  H_folds <- sapply(seq_len(k), function(i) mean_entropy_fold(as.matrix(pred_list[[mod]][[i]])))
  sum(fold_wts * H_folds)
})
Hnorm_mean_by_model <- sapply(names(pred_list), function(mod) {
  Hn_folds <- sapply(seq_len(k), function(i) mean_entropy_fold_norm(as.matrix(pred_list[[mod]][[i]]), m = m))
  sum(fold_wts * Hn_folds)
})
Htab= H_mean_by_model[names(H_mean_by_model) %in% nams_all]
Hnormtab <- Hnorm_mean_by_model[names(Hnorm_mean_by_model) %in% nams_all]


L1tab <- L1[names(L1) %in% nams_all]
collapse_rows_dt1 <- rbind(naive, regression, forest, mjp, mjp_warp, mjp_mixture)
collapse_rows_dt1 <- collapse_rows_dt1[c("Model", "param", "Covariates")]
nams_idx2 <- match(nams_all, names(L1tab))
collapse_rows_dt1$L1 <- L1tab[nams_idx2]

nams_idx3 <- match(nams_all, names(Hnormtab))

collapse_rows_dt1$Entropy_norm <- Hnormtab[nams_idx3]

colnames(collapse_rows_dt1) <- c("Method", "Model", "Exogenous info.", "L1-norm", "Entropy")
n_cols = 5
kableExtra::kbl(collapse_rows_dt1,booktabs = T, align = c("l","l","l"), linesep = '', format = "latex",escape = FALSE, digits = 4) %>%
  column_spec(1, bold=T) %>% collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts) %>%
  row_spec(3,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(5,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(11, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(14, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(17, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(20,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(25,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(27, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(29, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(31,  extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(35, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(37, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(39, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  row_spec(41, extra_latex_after = "\n\\addlinespace[0.5em]") %>%#sprintf("\\addlinespace[0.3em]\n\\cdashline{2-%d}\n\\addlinespace[0.3em]", n_cols)) %>%
  #footnote(number = c("Exogenous information indicates whether the rail characteristics described in \\\\Cref{sec:data} are included: Yes (global) uses \\\\cref{eq:trans_rates_covariates}; Yes (state-spec.) uses \\\\cref{eq:trans_rates_covariates2}."),escape = FALSE, general_title = "", threeparttable = TRUE) %>%
  writeLines(con = paste0("figures/transition_structure_L1",".tex"))



##############################
### Estimated coefficients ###
##############################
library(kableExtra)
par_list = readRDS("results/estimated_model_pars.Rdata")
nams = c("olr_cov", "opr_cov", "ocllr_cov", "gerlang_TRUE_all", "gerlang_relax_TRUE_all", "free_upper_tri_TRUE_all")
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
count = 0; k = 10;  
for(i in nams){
  mat = par_list[[i]]
  mat = mat[, (ncol(mat)-4):ncol(mat)]
  if(count == 0){ estim = apply(mat, MARGIN = 2, FUN = quantile, probs = c(0, 0.5, 1)); NAMS = rep(i, 3);} else { estim = rbind(estim, apply(mat, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.5, 0.95))); NAMS = c(NAMS, rep(i,3))}
  count = count + 1
}
df = as.data.frame(round(estim,digits=2)); 
df$X = NAMS
df$quantiles = rep(c("0%", "50%", "100%"), length(nams))
get_quantile_strings = function(df, colname){
  quants = paste(df[df$quantiles == "0%",colname], df[df$quantiles == "100%",colname], sep = ", ")
  sep_brackets1 = rep("[", length(quants)); sep_brackets2 = rep("]", length(quants)); 
  return(#paste(df[df$quantiles == "50%",colname], 
               paste(paste(sep_brackets1, quants, sep = ""), sep_brackets2, sep = "")
               #)
  )
}
foo = sapply(exo.cols,FUN = get_quantile_strings, df = df)
row.names(foo) = nams
colnames(foo) = c("Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")
foo = as.data.frame((foo))
foo$Model = NA
foo[grepl("cov", rownames(foo) ), "Model"] = "Ordered multinomial regression"
foo[!grepl("cov", rownames(foo) ), "Model"] = "Markov jump process"
foo$param = NA
foo[grepl("cov", rownames(foo) ), "param"] = c("Logistic link", "Probit link", "Cloglog link")
foo[!grepl("cov", rownames(foo) ), "param"] = c("Generalized Erlang, $\\bm{A}'$", "Parameterized upper, $\\bm{A}''$", "Free upper, $\\bm{A}'''$")

foo = foo[, c("Model", "param","Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")]
rownames(foo) = 1:NROW(foo)
colnames(foo) = c("param", "Model","Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")
row_group_label_fonts <- list(
  list(bold = T, italic = F),
  list(bold = F, italic = F)
)
kableExtra::kbl(foo,booktabs = T, align = c("l","l","c","c","c","c","c"), linesep = '', format = "latex",escape = FALSE) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts)


###############################################
### plot of transition probability matrices ###
###############################################
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
library(Matrix)
library(ggplot2)
library(ggpubr)
library(MASS)
library(reshape2)
library(tidyr)


m = 5
states = c("3", "2B", "2A", "1", "0")
#par_list = readRDS("results/estimated_model_pars.Rdata")
d = read.csv("defect_data.csv")
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
z = as.matrix(d[,exo.cols])
idx = 1 #sample(1:NROW(d), size = 1)
text.size <- 11
ndays = 365*10


source("plot_functions.R")
w <- 2; h <- 1.5  # target PDF size in inches
#k <- size_scaler(w, h, 3, 3) 
k = 1; base_size=6;
nam = "bidiagonal_cov1_exp_exp_all_no_warp"; generator = "bidiagonal"; warp_indicator = F; mixtures = FALSE; K = 4; state_covs = F;
p1 <- trans_dist_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
p2 <- trans_prob_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
ggplot2::ggsave(paste0("figures/trans_dist_", nam,".pdf"), p1,width = w, height = h, units = "in")
ggplot2::ggsave(paste0("figures/tpm_dist_", nam,".pdf"), p2,width = w, height = h, units = "in")


nam = "bidiagonal_cov2_exp_exp_all_no_warp"; generator = "bidiagonal"; warp_indicator = F; mixtures = FALSE; K = 1; state_covs = T;
p1 <- trans_dist_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
p1
p2 <- trans_prob_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
ggplot2::ggsave(paste0("figures/trans_dist_", nam,".pdf"), p1,width = w, height = h, units = "in")
ggplot2::ggsave(paste0("figures/tpm_dist_", nam,".pdf"), p2,width = w, height = h, units = "in")

nam = "bidiagonal_cov1_exp_exp_all_no_warp_mixture4"; generator = "bidiagonal"; warp_indicator = F; mixtures = TRUE; K = 4; state_covs = F;
p1 <- trans_dist_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
p2 <- trans_prob_fig(m, states, 1, nam, ndays, generator, warp_indicator,base_size = base_size, k = k, mixture = mixtures, state_covs = state_covs, K = K)
ggplot2::ggsave(paste0("figures/trans_dist_", nam,".pdf"), p1,width = w, height = h, units = "in")
ggplot2::ggsave(paste0("figures/tpm_dist_", nam,".pdf"), p2,width = w, height = h, units = "in")


### SOME TESTING
# #test on effective time
# nam = "bidiagonal_cov2_exp_exp_all_no_warp"; generator = "bidiagonal"; source("plot_functions.R")
# tt <- (0:(ndays-1))/365
# teff_list = rep(NA, length(tt)); count = 0
# w_list = matrix(NA, nrow = length(tt), ncol = 5)
# A_list = list()
# dist = matrix(NA, nrow = (length(tt)), ncol =  5); count = 0
# for(ttt in tt){
#   count = count + 1
#   teff <- MJP_effective_time(m, s1 = c(1), u = ttt,pars = get_pars(nam), z = matrix(z[idx, , drop = FALSE], nrow = 1), generator = generator, state_covs = T, covs_bin = T, mixture = F, warping = F,transient_dist_method = "pade")
#   dist[count,] = MJP_predict(m, s1 = c(1), u = ttt,pars = get_pars(nam), z = matrix(z[idx, , drop = FALSE], nrow = 1), generator = generator, state_covs = T, covs_bin = T, mixture = F, warping = F,transient_dist_method = "pade", K=4)
#   teff_list[count] = teff$t_eff
#   w_list[count,] = teff$w
#   A_list[[count]] = teff$A
# }
# 
# 
# count = 0; init =t(as.matrix(c(1,0,0,0,0)))
# dist = matrix(NA, nrow = (length(tt)), ncol =  5); count = 0
# for(ttt in tt){
#   count = count + 1
#   A = A_list[[count]]
#   dist[count,] = init %*% as.matrix(expm(A*teff_list[count]))
# }
# 
# plot(tt, dist[,1], type = "l")
# lines(tt, dist[,2], col = 2)
# lines(tt, dist[,3], col = 3)
# lines(tt, dist[,4], col = 4)
# lines(tt, dist[,5], col = 5)
# 
# 
# plot(tt, w_list[,1], type = "l")
# lines(tt, w_list[,2], col = 2)
# lines(tt, w_list[,3], col = 3)
# lines(tt, w_list[,4], col = 4)
# lines(tt, w_list[,5], col = 5)
# 
# plot(tt, teff_list, type="l")
# 
# plot(tt, w_list[,1], type ="l")
# lines(tt, w_list[,2], col = 2)
# lines(tt, w_list[,3], col = 3)
# lines(tt, w_list[,4], col = 4)
# lines(tt, w_list[,5], col = 5)
# 
# 
# A = make_A3(m = 5, c(1,0,100,0,
#                     2,0,0,
#                      .3,0,
#                      50))
#  
# 
# dist = matrix(NA, nrow = (length(tt)), ncol =  5); count = 0
# for(ttt in tt/10){
#   count = count + 1
#   dist[count,] = t(as.matrix(c(1,0,0,0,0))) %*% as.matrix(expm(A*ttt))
# }
# 
# plot(tt, dist[,1], type = "l")
# lines(tt, dist[,2], col = 2)
# lines(tt, dist[,3], col = 3)
# lines(tt, dist[,4], col = 4)
# lines(tt, dist[,5], col = 5)

### OLR logistic ###
olr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d, method = "logistic")
d.test = d[idx, ]
d.test$s1 = 1
sol1 = matrix(NA, nrow = ndays, ncol = m)
count = 0
for(t in (1:ndays)/365){
  count = count + 1
  d.test$t = t
  sol1[count, ] = predict(olr_cov, newdata = d.test, type = "p");
}
dpp <- as.data.frame(sol1)
colnames(dpp) <- states
dpp$time <- 1:ndays
dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
p4 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75) +
  theme(
    text = element_text(size = text.size, family = "serif"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal", # Set legend direction to horizontal
    legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")

##########################################
### Multicategory reliability diagrams ###
##########################################
library(reshape2)
library(ggplot2)
library(cowplot)

forecast_category <- function(forecast, quantiles) {
  cumulative <- cumsum(forecast)
  res = sapply(quantiles, function(q) which(cumulative >= q)[1])
  return(res)
}
compute_reliability <- function(obs, z, q, qmin, qmax) {
  if(obs > z){
    p = 0
  } else if(obs < z){
    p = 1
  } else if(obs == z){
    if(qmax != qmin){
      p = (q-qmin) / (qmax - qmin) 
    } else{
      p = 0.5
    }
  }
  return(p)
}
multi.reliable = function(obs, y, quantiles){
  L = length(quantiles)
  forecast_quantile_matrix <- t(apply(y, 1, forecast_category, quantiles = quantiles))
  Cq = matrix(NA, nrow = NROW(y), ncol = L)
  for(i in 1:NROW(y)){
    Q = which(forecast_quantile_matrix[i,] == obs[i])
    if(length(Q) > 0) {
      qmin = quantiles[min(Q)]
      qmax = quantiles[max(Q)]
    } else {
      qmin = NA
      qmax = NA
    }
    for(j in 1:L){
      Cq[i,j] = compute_reliability(obs[i], forecast_quantile_matrix[i,j],quantiles[j], qmin = qmin, qmax = qmax)
    }
  }
  Cqave = colMeans(Cq)
  
  forecast_diff_matrix <- forecast_quantile_matrix - obs
  abs_error_matrix = abs(forecast_diff_matrix)
  ave_cat_err = mean(colMeans(abs_error_matrix))
  
  return(list(
    Cq = Cqave, 
    ave_cat_err = ave_cat_err
  ))
}
make_boots = function(obs, y, quantiles, iters){
  res = matrix(NA, nrow = iters, ncol = length(quantiles))
  err = matrix(NA, nrow = iters, ncol = length(quantiles))
  for(i in 1:iters){
    idx = sample(1:NROW(obs), replace = T, size = NROW(obs) )
    res[i,] = multi.reliable(obs[idx], y[idx, ], quantiles)$Cq
    forecast_quantile_matrix <- t(apply(y[idx, ], 1, forecast_category, quantiles = quantiles))
    forecast_diff_matrix <- forecast_quantile_matrix - obs[idx]
    abs_error_matrix = abs(forecast_diff_matrix)
    err[i, ] = colMeans(abs_error_matrix)
  }
  return(list(
    Cq = apply(res, 2, quantile,c(0.1,0.9)), 
    quan_cat_err = apply(err, 2, quantile,c(0.1,0.9)) 
    ) )
}


checkerboard_plot = function(obs, y, quantiles, base_size = 11){
  forecast_quantile_matrix <- t(apply(y, 1, forecast_category, quantiles = quantiles))
  forecast_diff_matrix <- forecast_quantile_matrix - obs
  diff_df <- data.frame(forecast_diff_matrix)
  colnames(diff_df) = quantiles
  diff_df$Observation <- obs
  diff_df$ID <- 1:nrow(diff_df)  # Create an ID for each observation
  diff_df_long <- melt(diff_df, id.vars = c("ID", "Observation"), variable.name = "Quantile", value.name = "Forecast_Diff")
  text.size=base_size #10

  diff_df_long$Quantile <- as.numeric(as.character(diff_df_long$Quantile))
  
  p1 = ggplot(diff_df_long, aes(x = Quantile, y = Forecast_Diff)) +
    geom_tile(aes(fill = after_stat(count)), color = "white", stat = "bin2d", binwidth =  c(length(quantiles)/100, 1)) +
    scale_fill_gradient(low = "grey90", high = "black") +
    #scale_fill_gradient(low = "white", high = "black") +
    labs(x = "Forecast quantiles", y = "Category error") +
    scale_x_continuous(breaks = seq(min(diff_df_long$Quantile), max(diff_df_long$Quantile), by = 0.4)) + # Adjust 'by' as needed
    theme(
      legend.position = "none", 
      text = element_text(size = text.size, family = "serif"),  
      panel.background = element_rect(fill = NA, color = NA),  
      plot.background = element_rect(fill = NA, color = NA),  
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank(),
      panel.border = element_rect(color = "grey", fill = NA, linewidth = .5)
    )
  return(p1)
}


reliability_plot = function(Cq, CqCI, quantiles, base_size = 11 ){
  data <- data.frame(
    Time = quantiles,
    Value = r$Cq,
    Lower = rb$Cq[1,],
    Upper = rb$Cq[2,]
  )
  text.size = base_size
  p2 = ggplot(data, aes(x = Time, y = Value)) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  
    geom_line(color = "black") +                
    geom_point(color = "black", size = 0.5) + 
    #xlim(0,1) + 
    ylim(0,1) +
    #geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "black", alpha = 0.2) +  
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.05, color = "black") +
    labs(x = "Quantile of forecast distribution", 
         y = "Observations below forecast quantile") +
    theme(
      legend.position = "none", 
      text = element_text(size = text.size, family = "serif"),  
      panel.background = element_rect(fill = "white", color = "grey"),
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank(),  
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(breaks = quantiles) + 
    annotate("text", x = 0.25, y = 1, 
             label = paste("Average category error =", format(round(r$ave_cat_err,2), nsmall = 2)), 
             size = 3.5, hjust = 0.5, color = "black", family = "serif") +  
    annotate("text", x = 0.25, y = 0.925, 
             label = paste0("CI 90% = [", format(round(rowMeans(rb$quan_cat_err)[1],2), nsmall = 2), ", ", format(round(rowMeans(rb$quan_cat_err)[2],2), nsmall = 2), "]"), 
             size = 3.5, hjust = 0.5, color = "black", family = "serif")
  return(p2)
}

pool_cv_preds <- function(pred_list, model, k) {
  OBS_all <- integer(0)
  Y_all   <- NULL
  for (fold in 1:k) {
    obs <- pred_list[["obs"]][[fold]]
    OBS <- apply(obs, 1, which.max)
    y   <- as.matrix(pred_list[[model]][[fold]])
    OBS_all <- c(OBS_all, OBS)
    Y_all   <- rbind(Y_all, y)
  }
  list(OBS = OBS_all, Y = Y_all)
}

#pred_list = readRDS("results/model_predictions_V3.Rdata")
#pred_list[["ensemble"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_softplus_softplus_all"]] )/2#+ pred_list[["empirical_dist_corr"]])/3
#pred_list[["ensemble2"]] =(pred_list[["free_upper_tri_TRUE_softplus_softplus_all"]] + pred_list[["empirical_dist_smart"]])/2
#pred_list[["ensemble3"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_softplus_softplus_all"]] + pred_list[["empirical_dist_smart"]])/3
#pred_list[["ensemble"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_all"]] )/2


nams <- c("uniform","olr_cov", "random_forest_cov",
          "bidiagonal_cov1_exp_exp_all_no_warp",
          "bidiagonal_cov2_exp_exp_all_no_warp",
          "bidiagonal_cov1_exp_exp_all_warp",
          "bidiagonal_cov1_exp_exp_all_no_warp_mixture4",
          "ensemble_all")

quantiles <- seq(0.05, 0.95, by = 0.1)
iters <- 200
text.size <- 11
k <- 5

for (mod in nams) {
  pooled <- pool_cv_preds(pred_list, mod, k = k)
  OBS <- pooled$OBS
  y   <- pooled$Y
  
  r  <- multi.reliable(OBS, y, quantiles)
  rb <- make_boots(OBS, y, quantiles, iters)
  
  p1 <- checkerboard_plot(OBS, y, quantiles)
  p2 <- reliability_plot(r, rb, quantiles)
  p3 <- ggdraw() + draw_plot(p2) + draw_plot(p1, x = 0.59, y = 0.15, width = .4, height = .5)
  
  nfil <- paste0("figures/reliability_diag_", mod, "_pooledCV.pdf")
  pdf(file = nfil, width = 4, height = 3)
  print(p3)
  dev.off()
  graphics.off()
  
  print(mod)
}


# 
# 
# nams = c("uniform","olr_cov", "random_forest_cov", "bidiagonal_cov1_exp_exp_all_no_warp", "bidiagonal_cov2_exp_exp_all_no_warp", "bidiagonal_cov1_exp_exp_all_warp", "bidiagonal_cov1_exp_exp_all_no_warp_mixture4", "ensemble_all")
# fold_numbers = 1
# for(i in 1:length(nams)){
#   obs = pred_list[["obs"]][[fold_numbers]]
#   OBS = apply(obs, 1, which.max)
#   y = pred_list[[nams[i]]][[fold_numbers]]
#   quantiles = seq(0.05, 0.95, by = 0.1)
#   iters = 200
#   text.size = 11
#   
#   r = multi.reliable(OBS, y,quantiles)
#   rb = make_boots(OBS, y, quantiles, iters)
#   
#   p1 = checkerboard_plot(OBS, y, quantiles)
#   p2 = reliability_plot(r, rb, quantiles)
#   p3 = ggdraw() + draw_plot(p2) + draw_plot(p1, x = 0.59, y = 0.15, width = .4, height = .5)
#   
#   nfil = paste0("figures/reliability_diag_", nams[i],".pdf")
#   pdf(file = nfil,width = 4, height = 3) 
#   print(p3)
#   dev.off()
#   graphics.off()
#   
#   print(nams[i])
# }
#good

####################
### THE TRIPTYCH ###
####################
#i.e. binary evaluation
#How good are the respective models to forecast 0 or 1 events?
bin_predictions = function(col1, col2, preds){ #convert multicategorial predictions into binary for event in col1 to col2
  if(col1 == col2){
    res = preds[,col1]
  } else if(col1 < col2){
    res = rowSums(preds[,col1:col2])
  }
  res[res > 1] = 1;
  res[res < 0] = 0;
  return(res)
}

zero_bin_obs = function(col1, col2, obs){
  idx = apply(obs, 1, function(row) which(row == 1))
  return(idx >= col1 & idx <= col2)
}

library(ggplot2)
m = 5; int1 = m-1; int2 = m; 

### Reliability diagrams ###
source("binary_eval_funcs/reliability_diag.R")
pdf(file = "figures/reliability_diagram_empirical_corr.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd1 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_olr_cov.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd2 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["olr_cov"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_mjp_free.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd3 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_ensemble.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd4 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["ensemble"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()

### ROC curves ###
library("ROSE")
rc0 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["empirical_dist"]]), n.thresholds = 200)
rc1 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]),add.roc = T, n.thresholds = 200)
rc2 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["olr_cov"]]), add.roc = T, n.thresholds = 200)
rc3 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), add.roc = T, n.thresholds = 200)
rc4 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["ensemble"]]), add.roc = T, n.thresholds = 200)
nams = c(paste("Empi. dist. = ", format(round(rc1$auc, 3), nsmall = 2)), 
         paste("Cumu. link = ", format(round(rc2$auc, 3), nsmall = 2)),
         paste("MJP = ", format(round(rc3$auc, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(rc4$auc, 3), nsmall = 2))
        )
RC1 = data.frame(rc1$false.positive.rate,rc1$true.positive.rate, rep(nams[1], length(rc1$true.positive.rate)) ); colnames(RC1) = c("false.positive.rate", "true.positive.rate", "group");
RC2 = data.frame(rc2$false.positive.rate,rc2$true.positive.rate, rep(nams[2], length(rc2$true.positive.rate)) ); colnames(RC2) = c("false.positive.rate", "true.positive.rate", "group");
RC3 = data.frame(rc3$false.positive.rate,rc3$true.positive.rate, rep(nams[3], length(rc3$true.positive.rate)) ); colnames(RC3) = c("false.positive.rate", "true.positive.rate", "group");
RC4 = data.frame(rc4$false.positive.rate,rc4$true.positive.rate, rep(nams[4], length(rc4$true.positive.rate)) ); colnames(RC4) = c("false.positive.rate", "true.positive.rate", "group");
combined_data <- rbind(RC1, RC2, RC3, RC4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 11
pdf(file = "figures/ROCs.pdf",width = 4, height = 3) 
ggplot(combined_data, aes(x = false.positive.rate, y = true.positive.rate, color = group)) +
  geom_line(linewidth = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Group") +
  theme_minimal() +
  theme(text = element_text(size = text.size, family = "serif"), 
                        panel.background = element_rect(fill = "white", color = "black"), 
                        panel.grid.minor = element_line(color = "lightgray"),
                        legend.position = c(0.75, 0.25),  # Moves the legend closer to the bottom
                        legend.direction = "horizontal",
                        legend.background = element_rect(fill = "transparent", color = NA),  # Set transparent background for legend
                        legend.key = element_rect(fill = "transparent", color = NA)
        ) + 
  guides(color = guide_legend(title = "", ncol = 1))
dev.off()

### Murphy diagrams ###
source("binary_eval_funcs/murphy_diag.R")
md1 = murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]), f2 = bin_predictions(int1,int2,pred_list[["olr_cov"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))
md2 = murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), f2 = bin_predictions(int1,int2,pred_list[["ensemble"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))
a1 = sum(md1$y[c(-1),1] * diff(md1$tsep))
a2 =sum(md1$y[c(-1),2] * diff(md1$tsep))
a3 =sum(md2$y[c(-1),1] * diff(md2$tsep))
a4 =sum(md2$y[c(-1),2] * diff(md2$tsep))
nams = c(paste("Empi. dist. = ", format(round(a1, 3), nsmall = 2)), 
         paste("Cumu. link = ", format(round(a2, 3), nsmall = 2)),
         paste("MJP = ", format(round(a3, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(a4, 3), nsmall = 2))
)
MD1 = data.frame(md1$tsep,md1$y[,1], rep(nams[1], length(md1$y[,1])) ); colnames(MD1) = c("par", "val", "group");
MD2 = data.frame(md1$tsep,md1$y[,2], rep(nams[2], length(md1$y[,2])) ); colnames(MD2) = c("par", "val", "group");
MD3 = data.frame(md2$tsep,md2$y[,1], rep(nams[3], length(md2$y[,1])) ); colnames(MD3) = c("par", "val", "group");
MD4 = data.frame(md2$tsep,md2$y[,2], rep(nams[4], length(md2$y[,2])) ); colnames(MD4) = c("par", "val", "group");
combined_data <- rbind(MD1, MD2, MD3, MD4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 11
pdf(file = "figures/muprhy.pdf",width = 4, height = 3) 
ggplot(combined_data, aes(x = par, y = val, color = group)) +
  geom_line(linewidth = 0.75) +
  labs(x = expression(omega), y = expression(paste("Empirical score, ", S(omega),"")), color = "Group") +
  theme_minimal() + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.position = c(0.75, 0.25),  # Moves the legend closer to the bottom
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),  # Set transparent background for legend
        legend.key = element_rect(fill = "transparent", color = NA)
  ) + 
  guides(color = guide_legend(title = "", ncol = 1))
dev.off()


################################################
### TESTING SENSITIVITY ON INITIAL CONDITION ###
################################################

rm(list = ls()) #clear memory
d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
track <- unique(d$Track)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
library(TMB)
compile("FUNCS_MJP_with_TMB.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB"))

library(lhs)
N <- 5
npar <- 25
set.seed(123)
initial_points <- randomLHS(N, npar)

gen = "free_upper_tri"; covslink = "softplus"; cov = T; log_bin = T; rps_bin = T; brier_bin = T
results <- vector("list", N)

for (i in 1:N) {
  init <- initial_points[i, ]
  time_taken <- system.time(
    opt_res <- tryCatch(
      optim( par = init, fn = MJP_score, m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = as.matrix(d[,exo.cols]), generator = gen, link_type_base = covslink, link_type_covs = covslink, covs_bin = cov, likelihood_bin = log_bin, rps_bin = rps_bin, brier_bin = brier_bin, transient_dist_method = "pade", method = "BFGS", control = list(maxit = 1000)),
      error = function(e) NULL
    )
  )
  if (!is.null(opt_res)) {results[[i]] <- list(par = opt_res$par,value = opt_res$value,convergence = opt_res$convergence,counts = opt_res$counts,time = time_taken[["elapsed"]])} 
  else {results[[i]] <- list(par = rep(NA, npar),value = NA,convergence = NA,counts = NA,time = NA)}
}
df <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    t(res$par),
    objective = res$value,
    convergence = res$convergence,
    time = res$time
  )
}))

library(dplyr)
library(ggplot2)
library(cluster)

summary(df)

# Clustering of parameter estimates
param_matrix <- as.matrix(df[, !(names(df) %in% c("objective", "convergence","time"))])
k <- 3  # Adjust number of clusters as needed
clust <- kmeans(param_matrix, centers = k)
df$cluster <- factor(clust$cluster)

ggplot(df, aes(x = objective)) +
  geom_histogram(bins = 40, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Final Objective Values")

ggplot(df, aes(x = cluster, y = objective)) +
  geom_boxplot() +
  theme_minimal()











