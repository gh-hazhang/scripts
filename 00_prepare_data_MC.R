setwd("/ghsfa/projects/pharma/projects/sirius_pharma/hazhang_projects/Treatement_Effect_RUOMRD_CSO_call_05052024/DMR")
dataDir <- "/ghsfa/projects/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/METHYLATION_DATA"
saveDir <- "/ghsfa/projects/pharma/projects/sirius_pharma/hazhang_projects/Treatement_Effect_RUOMRD_CSO_call_05052024/DMR/METHYLATION_DATA"
if(!file.exists(saveDir)) dir.create(saveDir, recursive = T)
outdir <-  "/ghsfa/projects/pharma/projects/sirius_pharma/hazhang_projects/Treatement_Effect_RUOMRD_CSO_call_05052024/DMR/METHYLATION_DATA"
if(!file.exists(outdir)) dir.create(outdir, recursive = T)

# the following can be commented out if libraries are already installed
# devtools::install("/home/mcai/workspace/ghpipeline/ghcnv")
install.packages("quantmod")
install.packages("fastglm")

# library(ghcnv)
library(ggplot2)
library(quantmod)
library(plotly)
library(GenomicRanges)
library(RColorBrewer)
library(MASS)
library(parallel)
library(readr)
library(fastglm)

############################
##  Load panels annotation
############################

repo <- "/ghsfa/projects/bioinformatics/02_DEVELOPMENT/230325_S3_APPLICATION_PDL1/results"

# epi_panel, epi_probes, genomic_panel, genomic_probes
load(file = file.path(repo, "230403_PDL1_PREDICTION/panel_annotation.RData"))

#################
#  Sample info
#################

load(file = "/ghds/bioinformatics/02_DEVELOPMENT/230325_S3_APPLICATION_PDL1/results/230403_PDL1_PREDICTION/sample_sirius.RData")

sample_ambrx_raw = read.table("/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/Ambrx_all_batches_all_col_sample_level.csv",
                          sep=",", header=T, stringsAsFactors = F, comment.char = "")

sample_ambrx = sample_ambrx_raw[c("GHSampleID", "runid", "fc_dir", "Sample_status")]
sample_ambrx["QC_call_merged"] = ifelse(sample_ambrx$Sample_status=="SUCCESS", "PASS", "FAIL")
sample_ambrx = subset(sample_ambrx, select = c("GHSampleID", "runid", "fc_dir", "QC_call_merged"))
colnames(sample_ambrx) = c("run_sample_id", "runid", "bip_path", "QC_call_merged")
sample_ambrx["cancer_type"] = "Prostate"
sample_ambrx["category"] = "prostate"
sample_ambrx["source"] = "ARX_SOW01"

# # merge with samples_sirius
# sample_sirius = dplyr::bind_rows(sample_sirius, sample_ambrx)
# save(sample_sirius, file = file.path(dataDir, "sample_sirius.RData"))

# ############################
# ##  LDT re-processed 
# ##  MSRE normalized data
# ###########################

# #sample_sai <- read.table(file = "/ghds/omni_v2/users/schen/methylation_data/all-methyl-eap-dev-training-meta.all_samples.tsv",
# #                          sep = "\t", header = T, stringsAsFactors = F, comment.char = "")

# # A. samples curated by Sai (N=4439)
# sample_sai <- read.table(file = "/ghds//omni_v2/users/schen/methylation_data/intermediate/msre-LDT-merge-result.tsv",
#                          sep = "\t", header = T, stringsAsFactors = F, comment.char = "")
# sample_sai["fc_dir"] = paste("/ghds/omni_v2/users/schen/methylation_data/msre_call", sample_sai$runid, sample_sai$sample_id, sep="/")

# sample_ambrx_for_sai = sample_ambrx[c("run_sample_id", "runid")]
# colnames(sample_ambrx_for_sai) = c("sample_id", "runid")
# # sample_ambrx_for_sai["repo"] = "/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/explore_1/ldt_methyl/msre_call"
# sample_ambrx_for_sai["fc_dir"] = paste("/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/explore_1/ldt_methyl/msre_call", 
#                                        sample_ambrx_for_sai$sample_id, sep="/")
# sample_sai = dplyr::bind_rows(sample_sai, sample_ambrx_for_sai)
# sample_sai <- sample_sai[match(sample_sirius$run_sample_id, sample_sai$sample_id),]

# # repo <- "/ghds/omni_v2/users/schen/methylation_data/msre_call"
# msre_lr_feature_files <- file.path(sample_sai$fc_dir, 
#                                    paste0(sample_sai$sample_id, ".msre_lr_feature.tsv"))

# # # B. samples from Ambrx cohort (N=109)
# # 
# # ambrx_repo = "/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/explore_1/ldt_methyl/msre_call"
# # msre_lr_feature_files_ambrx <- file.path(ambrx_repo, sample_ambrx$run_sample_id, 
# #                                        paste0(sample_ambrx$run_sample_id, ".msre_lr_feature.tsv"))
# # 
# # # 4439 + 109 = 4548
# # msre_lr_feature_files = c(msre_lr_feature_files_sai, msre_lr_feature_files_ambrx)
# l <- sapply(msre_lr_feature_files, file.exists)
# sample_sai <- sample_sai[l,]
# msre_lr_feature <- mclapply(msre_lr_feature_files[l], read.table, sep = "\t", header = T, stringsAsFactors = F, mc.cores = 12)

# # 
# msre_all_call_files <- file.path(sample_sai$fc_dir, 
#                                    paste0(sample_sai$sample_id, ".msre_all_call.hdr.tsv"))

# l <- sapply(msre_all_call_files, file.exists)
# s <- sapply(msre_all_call_files, file.size)
# msre_all_call <- mclapply(msre_all_call_files[l], function(f) read.table(f, sep = "\t", header = T, stringsAsFactors = F)[1,], mc.cores = 12)
# msre_all_call <- do.call(rbind, msre_all_call)

# save(sample_sai, msre_lr_feature, msre_all_call, file = file.path(dataDir, "msre_lr_feature.RData"))



# ##
# # A. samples curated by Sai (N=4439)
# msre_all_call_files_sai <- file.path(repo, sample_sai$runid, sample_sai$sample_id,
#                                  paste0(sample_sai$sample_id, ".msre_all_call.hdr.tsv"))
# l <- sapply(msre_all_call_files, file.exists)
# s <- sapply(msre_all_call_files, file.size)
# msre_all_call_sai <- mclapply(msre_all_call_files_sai[l], function(f) read.table(f, sep = "\t", header = T, stringsAsFactors = F)[1,], mc.cores = 12)
# msre_all_call_sai <- do.call(rbind, msre_all_call_sai)

# # B. samples from Ambrx cohort (N=109)
# msre_all_call_files_ambrx <- file.path(ambrx_repo, sample_ambrx$GHSampleID,
#                                      paste0(sample_ambrx$GHSampleID, ".msre_all_call.hdr.tsv"))
# l <- sapply(msre_all_call_files_ambrx, file.exists)
# s <- sapply(msre_all_call_files_ambrx, file.size)
# msre_all_call_ambrx <- mclapply(msre_all_call_files_ambrx[l], function(f) read.table(f, sep = "\t", header = T, stringsAsFactors = F)[1,], mc.cores = 12)
# msre_all_call_ambrx <- do.call(rbind, msre_all_call_ambrx)

# msre_all_call = rbind(msre_all_call_sai, msre_all_call_ambrx)

# save(sample_sai, sample_ambrx, msre_lr_feature, msre_all_call, file = file.path(dataDir, "msre_lr_feature.RData"))



# counts_bip_info <- data.frame(chrom = msre_lr_feature[[1]]$chr,
#                               start = msre_lr_feature[[1]]$start,
#                               end = msre_lr_feature[[1]]$end,
#                               region_id = paste(msre_lr_feature[[1]]$chr, msre_lr_feature[[1]]$start, msre_lr_feature[[1]]$end, sep = "_"))
# counts_bip_info <- counts_bip_info %>% mutate(type = epi_panel$type[match(counts_bip_info$region_id, epi_panel$region_id)],
#                                               gene = epi_panel$gene[match(counts_bip_info$region_id, epi_panel$region_id)])

# sample_sirius = sample_sirius[sample_sirius["run_sample_id"] != "B00225087", ]

# # i_match <- match(sample_sirius$run_sample_id, sample_sai$sample_id)
# # i_match <- match(sample_sirius$run_sample_id, sample_sai$sample_id)
# i_match <- match(sample_sirius$run_sample_id, sample_sai$sample_id)
# # i_match = 1:length(msre_lr_feature)
# counts_bip_norm <- lapply(i_match, function(i) {
#   out <- rep(NA, nrow(msre_lr_feature[[2]]));
#   if(!is.na(i)) {
#     if(!is.null(msre_lr_feature[[i]]) && nrow( msre_lr_feature[[i]]) > 0)
#       out <- msre_lr_feature[[i]]$log_norm_region;
#   }
#   out;})
# counts_bip_norm <- do.call(cbind, counts_bip_norm)
# colnames(counts_bip_norm) <- sample_sirius$run_sample_id
# rownames(counts_bip_norm) <- paste(msre_lr_feature[[1]]$chr, msre_lr_feature[[1]]$start, msre_lr_feature[[1]]$end, sep = "_")
# # is the next line really needed? run it anyways for now
# msre_all_call <- msre_all_call[match(sample_sirius$run_sample_id, msre_all_call$run_sample_id),]
# save(counts_bip_norm, counts_bip_info, msre_all_call, file = file.path(dataDir, "counts_ldt_bip_norm.RData"))

# counts_bip <- lapply(i_match, function(i) {
#   out <- rep(NA, nrow(msre_lr_feature[[2]]));
#   if(!is.na(i)) {
#     if(!is.null(msre_lr_feature[[i]]) &&nrow( msre_lr_feature[[i]]) > 0)
#       out <- msre_lr_feature[[i]]$num_molecules;
#   }
#   out;})
# counts_bip <- do.call(cbind, counts_bip)
# colnames(counts_bip) <- sample_sirius$run_sample_id
# rownames(counts_bip) <- paste(msre_lr_feature[[1]]$chr, msre_lr_feature[[1]]$start, msre_lr_feature[[1]]$end, sep = "_")
# save(counts_bip, counts_bip_info, msre_all_call, file = file.path(dataDir, "counts_ldt_bip.RData"))

# counts_hyper_epi_bip <- counts_bip
# counts_scaled_hyper_epi_bip <- counts_bip_norm
# counts_hyper_epi_bip_info <- counts_bip_info


# ##############################
# ##  Normalize bip scores
# ##  Adjust for batch effects
# ##############################

# i_samples <- which(apply(counts_hyper_epi_bip, 2, function(x) !any(is.na(x))))
# i_targets <- which(apply(counts_hyper_epi_bip, 1, mean, na.rm = T) > 0)

# runid <- factor(sample_sirius$runid[i_samples])
# dd <- data.frame(runid = runid)
# m <- model.matrix(~ runid, dd)

# if(0) {
#   targets_batch_coeffs <- mclapply(i_targets, function(i) {
#     y <- counts_scaled_hyper_epi_bip[i,i_samples];
#     l <- fastglm(x = m, y = y, family =  gaussian())
#     as.data.frame(summary(l)$coefficients)
#   }, mc.cores = 12)
#   #p_alk <- do.call(rbind, p_alk)
#   save(targets_batch_coeffs, file = file.path(saveDir, "targets_batch_coeffs.230424.RData"))
# }

# targets_batch_norm_score <- mclapply(i_targets, function(i) {
#   y <- counts_scaled_hyper_epi_bip[i,i_samples];
#   l <- fastglm(x = m, y = y, family =  gaussian())
#   l$residuals
# }, mc.cores = 12)
# targets_batch_norm_score <- do.call(rbind, targets_batch_norm_score)
# colnames(targets_batch_norm_score) <- sample_sirius$run_sample_id[i_samples]
# save(targets_batch_norm_score, i_samples, file = file.path(saveDir, "targets_batch_norm_score.230424.RData"))

# # the following 2 lines are not run, as boxplotc is not found
# # i <- which(epi_panel$type == "positive")[20]
# # boxplotc(targets_batch_norm_score[i,] ~ sample_sirius$runid[i_samples], las = 2, cex.axis = 0.1)

# ##########################################
# ##  Generate CpG strata molecule counts
# ##########################################
# catalin_dir = "/ghds/bioinformatics/02_DEVELOPMENT/230420_CRC_MOLECULAR_SUBTYPES/results/PROCESSED_SAMPLES"
# ambrx_dir = "/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype"
# sample_sirius["count_cg_strata_tsv1"] = paste(catalin_dir, sample_sirius$runid, sample_sirius$run_sample_id,
#                                           paste0(sample_sirius$run_sample_id, ".hyper.epi.counts.cg.strata.tsv"), sep="/")

# sample_sirius["count_cg_strata_tsv2"] = paste(ambrx_dir, sample_sirius$runid, sample_sirius$run_sample_id,
#                                            paste0(sample_sirius$run_sample_id, ".hyper.epi.counts.cg.strata.tsv"), sep="/")

# get_count_cg_strata <- function(a, b) {
#   if (file.exists(a)) {
#     return(a)
#   } else if (file.exists(b)) {
#     return(b)
#   } else {
#     return(NA)
#   }
# }

# sample_sirius["count_cg_strata_tsv"] = apply(sample_sirius[c("count_cg_strata_tsv1", "count_cg_strata_tsv2")], 1, 
#                                           function(row) get_count_cg_strata(row[1], row[2]))
# counts_cg_strata_files = sample_sirius["count_cg_strata_tsv"]
# l <- sapply(counts_cg_strata_files, file.exists)
# counts_cg_strata <- mclapply(counts_cg_strata_files[l], read.table, sep = "\t", header = T, stringsAsFactors = F,
#                              mc.cores = 12)
# save(counts_cg_strata, file = file.path(dataDir, "counts_cg_strata.RData"))

# #########################
# ##  Normalize scores
# ##  vs pos control
# ##  regions strata
# #########################

# cg_strata_lbls <- paste0("cg", c(as.character(6: 12), "12.30"))
# pseudo_count <- 1e-5
# i_pos_control <- which(epi_panel$type == "positive")

# counts_cg_strata_pos_controls <- lapply(counts_cg_strata, function(d) apply(d[i_pos_control, cg_strata_lbls], 2, sum, na.rm = T))
# counts_cg_strata_pos_controls <- do.call(rbind, counts_cg_strata_pos_controls)

# ##  Total pos cg 12:30
# M <- counts_cg_strata_pos_controls[,"cg12.30"]

# ##  cg 6:12 scores for
# ##  pos control regions
# Q <- lapply(1:nrow(sample_sirius), function(i)
#   log10(counts_cg_strata_pos_controls[i, -ncol(counts_cg_strata_pos_controls)]/M[i] + pseudo_count))
# Q <- do.call(rbind, Q)

# ##  normalize cg 6:12 scores for
# ##  pos control regions
# i_trn <- which((sample_sirius$cancer_type == "normal" | sample_sirius$source == "TND") &
#                  !sample_sirius$run_sample_id %in% c("B00249303"))
# m_trn <- apply(Q[i_trn,], 2, mean, na.rm = T)
# sd_trn <- apply(Q[i_trn,], 2, sd, na.rm = T)
# Q_prime <- t(apply(Q, 1, function(x) (x-m_trn)/sd_trn))
# colnames(Q_prime) <- colnames(Q)

# ##  boxplot targeted region vs batch
# ##  and scater vs cg 6:12 scores of pos control regions

# i_crc <- which(sample_sirius$cancer_type == "CRC")
# sd_crc <- apply(counts_scaled_hyper_epi_bip[,i_crc], 1, mad, na.rm = T)
# m_crc <- apply(counts_scaled_hyper_epi_bip[,i_crc], 1, median, na.rm = T)
# i_use <- which(m_crc > -1)
# i_plot <- i_use[order(sd_crc[i_use], decreasing = T)][1:3]

# pdf(file = file.path(saveDir, paste0("REGIONS_SCORES_VS_POS_CTRL_STRATA_Q_SCORE.pdf")), height = 8, width = 12)
# for(i_target in i_plot) {
#   df <- data.frame(runid = rep(sample_sirius$runid[i_crc], ncol(Q_prime)),
#                    run_sample_id = rep(sample_sirius$run_sample_id[i_crc], ncol(Q_prime)),
#                    target_score = rep(counts_scaled_hyper_epi_bip[i_target,i_crc], ncol(Q_prime)),
#                    pos_ctrl_score =as.vector( Q_prime[i_crc,]),
#                    cg_strata = factor(rep(colnames(Q_prime), each = length(i_crc)), levels = colnames(Q_prime)))
  
#   g <- ggplot(df, aes(x = runid, y = target_score, colour = runid)) +
#     geom_boxplot(show.legend = F) +
#     ylab("target region BIP score") +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target], "CRC samples")) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 6, angle = 90),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df, aes(x = target_score, y = pos_ctrl_score, colour = runid)) +
#     geom_point(show.legend = FALSE) +
#     ylim(-1.5,1.5) +
#     xlab("target region BIP score") + ylab("Total pos control Q' score") +
#     facet_wrap(~ cg_strata, nrow = 2) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target], "CRC samples")) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 18, angle = 0),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   df <- data.frame(runid = rep(sample_sirius$runid[i_trn], ncol(Q_prime)),
#                    run_sample_id = rep(sample_sirius$run_sample_id[i_trn], ncol(Q_prime)),
#                    target_score = rep(counts_scaled_hyper_epi_bip[i_target,i_trn], ncol(Q_prime)),
#                    pos_ctrl_score =as.vector( Q_prime[i_trn,]),
#                    cg_strata = factor(rep(colnames(Q_prime), each = length(i_trn)), levels = colnames(Q_prime)))
  
#   g <- ggplot(df, aes(x = runid, y = target_score, colour = runid)) +
#     geom_boxplot(show.legend = F) +
#     ylab("target region BIP score") +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target], "NORMALS+TND samples")) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 6, angle = 90),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df, aes(x = target_score, y = pos_ctrl_score, colour = runid)) +
#     geom_point(show.legend = FALSE) +
#     ylim(-1.5,1.5) +
#     xlab("target region BIP score") + ylab("Total pos control Q' score") +
#     facet_wrap(~ cg_strata, nrow = 2) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target], "NORMALS+TND samples")) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 18, angle = 0),
#           axis.text.y=element_text(size = 18))
#   print(g)
# }
# dev.off()


# #############################################
# ##  normalize targeted regions to
# ##  cg 6:12 scores of pos control regions
# #############################################

# library(fastglm)

# #i_trn <- which((sample_sirius$cancer_type == "normal" | sample_sirius$source == "TND") &
# #                 !sample_sirius$run_sample_id %in% c("B00249303"))
# set.seed(12345)
# i_trn <- sample(which(sample_sirius$QC_call_merged == "PASS"), round(sum(sample_sirius$QC_call_merged == "PASS")/2))

# ##  regression lines
# x <-  Q_prime[i_trn,]
# target_region_coeffs <- mclapply(1:nrow(counts_scaled_hyper_epi_bip), function(i) {
#   y <- counts_scaled_hyper_epi_bip[i, i_trn];
#   l <- fastglm(x = x, y = y, family =  gaussian())
#   #l <- lm(y~x)
#   as.data.frame(summary(l)$coefficients)
# }, mc.cores = 12)
# intercept <- apply(counts_scaled_hyper_epi_bip[,i_trn], 1, mean, na.rm = T)

# region_bias <- t(sapply(1:nrow(counts_scaled_hyper_epi_bip), function(i) Q_prime %*% target_region_coeffs[[i]]$Estimate + intercept[i]))
# counts_scaled_hyper_epi_strata <- counts_scaled_hyper_epi_bip - region_bias
# save(counts_scaled_hyper_epi_strata, file = file.path(saveDir, "counts_scaled_hyper_epi_strata.230424.RData"))

# pdf(file = file.path(saveDir, paste0("REGIONS_SCORES_NORMALIZED.230424.pdf")), height = 8, width = 12)
# for(i_target in i_plot) {
#   df_norm <- data.frame(runid = rep(sample_sirius$runid, ncol(Q_prime)),
#                         run_sample_id = rep(sample_sirius$run_sample_id, ncol(Q_prime)),
#                         target_score = rep(counts_scaled_hyper_epi_strata[i_target,], ncol(Q_prime)),
#                         pos_ctrl_score =as.vector( Q_prime),
#                         cg_strata = factor(rep(colnames(Q_prime), each = nrow(sample_sirius)), levels = colnames(Q_prime)))
  
#   df_bip <- data.frame(runid = rep(sample_sirius$runid, ncol(Q_prime)),
#                        run_sample_id = rep(sample_sirius$run_sample_id, ncol(Q_prime)),
#                        target_score = rep(counts_scaled_hyper_epi_bip[i_target,], ncol(Q_prime)),
#                        pos_ctrl_score =as.vector( Q_prime),
#                        cg_strata = factor(rep(colnames(Q_prime), each = nrow(sample_sirius)), levels = colnames(Q_prime)))
  
#   df_batch_norm <- data.frame(runid = rep(sample_sirius$runid[i_samples], ncol(Q_prime)),
#                               run_sample_id = rep(sample_sirius$run_sample_id[i_samples], ncol(Q_prime)),
#                               target_score = rep(counts_scaled_hyper_epi_bip[i_target, i_samples], ncol(Q_prime)),
#                               pos_ctrl_score =as.vector( Q_prime[i_samples,]),
#                               cg_strata = factor(rep(colnames(Q_prime), each = length(i_samples)), levels = colnames(Q_prime)))
  
#   g <- ggplot(df_norm, aes(x = runid, y = target_score, colour = runid)) +
#     geom_boxplot(show.legend = F) +
#     ylab("target region STRATA normalized score") +
#     ylim(-5,10) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 6, angle = 90),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df_bip, aes(x = runid, y = target_score, colour = runid)) +
#     geom_boxplot(show.legend = F) +
#     ylim(-5,10) +
#     ylab("target region BIP score") +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 6, angle = 90),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df_batch_norm, aes(x = runid, y = target_score, colour = runid)) +
#     geom_boxplot(show.legend = F) +
#     ylab("target region batch norm score") +
#     ylim(-5,10) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 6, angle = 90),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df_norm, aes(x = target_score, y = pos_ctrl_score, colour = runid)) +
#     geom_point(show.legend = FALSE) +
#     ylim(-1.5,1.5) +
#     xlab("target region STRATA normalized score") + ylab("Total pos control Q' score") +
#     facet_wrap(~ cg_strata, nrow = 2) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 18, angle = 0),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df_bip, aes(x = target_score, y = pos_ctrl_score, colour = runid)) +
#     geom_point(show.legend = FALSE) +
#     ylim(-1.5,1.5) +
#     xlab("target region BIP score") + ylab("Total pos control Q' score") +
#     facet_wrap(~ cg_strata, nrow = 2) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 18, angle = 0),
#           axis.text.y=element_text(size = 18))
#   print(g)
  
#   g <- ggplot(df_batch_norm, aes(x = target_score, y = pos_ctrl_score, colour = runid)) +
#     geom_point(show.legend = FALSE) +
#     ylim(-1.5,1.5) +
#     xlab("target region BATCH normalized score") + ylab("Total pos control Q' score") +
#     facet_wrap(~ cg_strata, nrow = 2) +
#     ggtitle(paste("Target region", epi_panel$region_id[i_target])) +
#     theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
#           text = element_text(size = 14), # , family = "Tahoma"
#           axis.title = element_text(face="bold"),
#           axis.text.x=element_text(size = 18, angle = 0),
#           axis.text.y=element_text(size = 18))
#   print(g)
# }
# dev.off()
