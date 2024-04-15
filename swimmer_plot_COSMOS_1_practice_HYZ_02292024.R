library(readxl)
library(dplyr)
library(janitor)
library(lubridate)
library(ggplot2)
library(ggnewscale)
library(ggsci)
library(stringr)

setwd("/ghsfa/projects/pharma/projects/sirius_pharma/hazhang_projects/swimmer_plot_COSMOS_1_practice_HYZ_02292024")

input_file_path = "COS_01_Combined Report (Batch1,2&Enrollment&cfDNA retains)_MAFBAND UPDATE_TR working copy.xlsx"
sheet_name = "COS_01 DATA MERGE"
output_prefix = "COSMOS_Swimmer_plot_1"

target_patient_df =  readxl::read_excel("Patient Lists for SP_8.17.23.xlsx", sheet="Sheet1",
                                        na = c("", "#VALUE", "n/a", "NA", "/", "\\", "?")) %>% 
  janitor::clean_names()

target_patients = target_patient_df[["colon_relapse"]]
target_patients = target_patients[!is.na(target_patients)]
target_patients = gsub('COSMOS', '', target_patients)

dates_to_relative_times <- function(data_in, surg_column, relevant_columns) {
  # Takes lubridate "date" columns and creates matching "time" column, relative to surgery date
  
  for(date_column in relevant_columns) {
    col_name <- gsub("date", "time", date_column)
    # data_in[,col_name] <- difftime(as.Date(data_in[[date_column]]), as.Date(data_in[[surg_column]]), units="days")
    data_in[,col_name] <- difftime(parse_date_time(data_in[[date_column]], c("mdy", "ymd")), 
                                   parse_date_time(data_in[[surg_column]], c("mdy", "ymd")), 
                                   units="days")
  }
  
  # for some reason I can't just 'as.numeric()' the calculation above, so just sapply it after all of them have been made
  new_time_cols <- colnames(data_in)[grepl("time", colnames(data_in))]
  data_in[,new_time_cols] <- sapply(data_in[,new_time_cols], as.numeric)
  
  return(data_in)
}


original_df <- readxl::read_excel(input_file_path, sheet=sheet_name,
                                  na = c("", "#VALUE", "n/a", "NA", "/", "\\", "?")) %>% 
  janitor::clean_names() %>% 
  filter(is.na(patient_id) == F) %>% 
  filter(!visit_name %in% c("Enrollment", "7 days", "14 days")) %>%
  distinct() %>%
  rename("date_last_act_dose" = "last_act_dose") %>%
  mutate(patient_id = gsub("COSMOS", "", patient_id),
         relapse_site_site = str_replace(relapse_site_site, "Peritonal Only", "Peritoneal Only"))
  


test_df3 = read.table("test_3.csv", sep=",", stringsAsFactors = F, header = T) %>%
  janitor::clean_names() %>%
  mutate(patient_id = gsub("COSMOS", "", patient_id),
         date_of_recurrence = parse_date_time(date_of_recurrence, c("mdy", "ymd")))

original_df["date_of_recurrence"] = NULL
original_df["relapse_site_site"] = NULL
original_df = merge(original_df, test_df3, by.x = "patient_id", 
                    by.y = "patient_id", all.x = TRUE)


original_df = original_df %>%
  filter(patient_id %in% target_patients)

# order
patient_levels <- original_df %>% 
  select(c("patient_id")) %>%
  distinct() %>%
  pull(patient_id)

original_df <- original_df %>% mutate(patient_id = factor(patient_id, levels = rev(patient_levels)))


# 1. assay_df for ctdna_status
# assay_df = data.frame()
# for (timepoint in c("x28_days", "x3_months", "x6_months", "x9_months", 
#                     "x1_year", "x1_year_and_6_months", "x2_years", "x2_years_and_6_months",
#                     "recurrence_sample")) {
#   timepoint_df = original_df %>%
#     select(c("patient_id", "date_of_definitive_resection", paste0(timepoint, "_date"), paste0(timepoint, "_result"))) %>%
#     rename("bloodcoll_date" = paste0(timepoint, "_date") , 
#            "ctdna_status" = paste0(timepoint, "_result")) %>%
#     filter(is.na(bloodcoll_date) == F, is.na(ctdna_status) == F) %>%
#     dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("bloodcoll_date")) %>%
#     rename("days" = "bloodcoll_time") %>%
#     mutate(ctdna_detected = ifelse(ctdna_status == 0, "Not Detected", 
#                                    ifelse(ctdna_status == 1, "Detected", "QC Failed")))
#   assay_df = dplyr::bind_rows(assay_df, timepoint_df)
# }
assay_df =  original_df %>%
  select(c("patient_id", "date_of_definitive_resection", "blood_collection_date", "ct_dna_mafband_ruo")) %>%
  filter(is.na(blood_collection_date) == F, is.na(ct_dna_mafband_ruo) == F) %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("blood_collection_date")) %>%
  rename("days" = "blood_collection_time") %>% 
  mutate(ctdna_detected = ifelse(ct_dna_mafband_ruo == 0, "Not Detected",
                                     ifelse(ct_dna_mafband_ruo == 1, "Detected", "QC Failed"))) %>%
  distinct()
           


# 2.event_df for date_of_definitive_resection, date_of_last_fu

first_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("date_of_definitive_resection")) %>%
  select(c("patient_id", "time_of_definitive_resection")) %>%
  rename("days" = "time_of_definitive_resection")

last_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("date_of_last_fu")) %>%
  select(c("patient_id", "time_of_last_fu")) %>%
  rename("days" = "time_of_last_fu")

event_df = dplyr::bind_rows(first_df, last_df) %>%
  distinct()

# 3. treatment_df for 
treatment_start_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("start_date_act")) %>%
  select(c("patient_id", "start_time_act")) %>%
  rename("days" = "start_time_act")

treatment_end_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("date_last_act_dose")) %>%
  select(c("patient_id", "time_last_act_dose")) %>%
  rename("days" = "time_last_act_dose")

treatment_df = dplyr::bind_rows(treatment_start_df, treatment_end_df) %>%
  mutate(treatment="Adjuvant Chemotherapy") %>%
  distinct()

# 4. recurrence
recur_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("date_of_recurrence")) %>%
  select(c("patient_id", "time_of_recurrence", "relapse_site_site")) %>%
  rename("days" = "time_of_recurrence") %>%
  filter(is.na(days) == F) %>%
  mutate(Relapse = "Relapse") %>%
  distinct()

# 5. outcome
# 2.event_df for date_of_definitive_resection, date_of_last_fu
last_fu_df = original_df %>%
  dates_to_relative_times(surg_column = "date_of_definitive_resection", relevant_columns = c("date_of_last_fu")) %>%
  select(c("patient_id", "time_of_last_fu", "survival")) %>%
  rename("days" = "time_of_last_fu") %>%
  distinct()

outcome_df = last_fu_df%>%rename("outcome"="survival")

# 6. meta
meta_df = original_df %>%
  select(c("patient_id", "tumor_stage", "relapse_site_site")) %>%
  mutate(days_tumor_stage = -120, days_relapse_site = -60) %>%
  rename("Tumor Stage" = "tumor_stage", "Relapse Site" = "relapse_site_site")

###### 
# plot
dat <- list(
  assay = assay_df,
  background_line = event_df,
  tx = treatment_df,
  relapse = recur_df,
  outcome = outcome_df,
  meta = meta_df)

# Guardant health marketing color palette, with names for reference
gh_palette <- c("#225fac", "#e52e2c", "#ed2290", "#23a0db", "#a3cd3b", "#f89d4e")
names(gh_palette) <- c("dark_blue", "red", "pink", "light_blue", "green", "orange")

ggplot(data = NULL, aes(x = days, y = patient_id)) +
  
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.5) + # Vertical line at t=0
  
  # Scales
  scale_x_continuous("Years since surgery", breaks = c(-120, -60, seq(-365,1500,365)), labels = c("Tumor Stage", "Relapse Site", -1:4)) +
  labs(y = "Patient ID") +
  scale_alpha_manual(values = c(1,1)) +
  
  # Core data
  geom_line(data = dat$background_line, aes(x = days, y = patient_id, group = patient_id), alpha = 0.25) + # Background line
  
  # ctdna
  geom_point(data = dat$assay, aes(x = days, y = patient_id, fill = ctdna_detected), size = 3, shape = 21) + # ctDNA detected
  scale_fill_manual("ctDNA Result", values = c(gh_palette["red"][[1]], "white", "grey"),
                    guide=guide_legend(order=6)) +

  # treatment
  geom_line(data = dat$tx, aes(color = treatment), size = 5, alpha = 0.4) + # Treatment lines
  scale_color_manual("Treatment", values = c(gh_palette["light_blue"][[1]]),
                     guide=guide_legend(order=3)) +

  # outcome
  new_scale_fill() +
  geom_point(data = dat$outcome, aes(x = days, y = patient_id, fill = outcome), size = 3, shape = 24, color = "white") +
  scale_fill_manual("Outcome", values = c("Death" = gh_palette["pink"][[1]], 
                                          "Survival" = gh_palette["green"][[1]]),
                    guide=guide_legend(order=4)) +

  # recurrence
  new_scale_fill() +
  # Mapping recurrence
  geom_point(data = dat$relapse,
              fill = NA, size = 5, aes(alpha = Relapse), shape = "|", color = "red") +
  
  geom_tile(data=dat$meta, aes(x=days_tumor_stage, fill=`Tumor Stage`), width=50, color="white", lwd=1) +
  # scale_fill_lancet() +
  scale_fill_brewer(palette="Set3",
                    guide=guide_legend(order=1, ncol=2)) +
  new_scale_fill() +
  geom_tile(data=dat$meta %>% filter(is.na(`Relapse Site`) == F), aes(x=days_relapse_site, fill=`Relapse Site`), width=50, color="white", lwd=1) +
  scale_fill_d3(guide=guide_legend(order=2, ncol=2)) +
  
  # Use guides to make the legends work, and set legend order
  # guides(alpha = guide_legend(order = 3, title = element_blank(), override.aes = list(color = c("red","black"))),
  #        color = guide_legend(order = 2),
  #        fill = guide_legend(order = 1),
  #        shape = guide_legend(order = 4)) +
  
  guides(alpha = guide_legend(order = 5)) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = c(90, 90, 0, 0, 0, 0)))
# 
# hjust = c(1, 1, 0.5, 0.5, 0.5, 0.5),
# vjust = c(0.5, 0.5, 0, 0, 0, 0)
ggsave(paste0(output_prefix, ".png"), width = 9, height = 8)
