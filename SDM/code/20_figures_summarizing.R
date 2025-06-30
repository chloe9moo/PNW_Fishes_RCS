### PNW FISHES SDMs: MAKE FIGURES, TABLES, SUMMARIZE RESULTS ###

#set up
library(tidyverse)
options(readr.show_col_types = FALSE)

PATH <- "~/Documents/projects/PNW_fishes/SDM"

#VARIABLE IMPORTANCE ----
## read in all spp var imp results and combine into single dataframe ----
vi.files <- list.files(paste0(PATH, "/Output/Results_Tables/"), pattern = "*_VarImp.csv")

vi.res <- lapply(vi.files, function(x) {
  
  spp <- gsub("_", " ", gsub("_SpCont_VarImp.csv", "", x))
  df <- read_csv(paste0(PATH, "/Output/Results_Tables/", x))
  
  df$species <- spp
  
  return(df)
  
})

vi.res <- do.call(rbind, vi.res)

### add in species groups, correct names ----
vi.res <- vi.res %>% 
  filter(species != "Prosopium williamsoni") %>% ##salmonid snuck through
  mutate(species = if_else(grepl("Acrocheilus", species), "Gila alutacea", species)) %>%
  mutate(species_group = case_when(grepl("Gila|Mylocheilus|Ptychocheilus|Siphateles", species) ~ "Large minnows",
                                   grepl("Rhinichthys|Richardsonius|Oregonichthys", species) ~ "Small minnows",
                                   grepl("Cottus", species) ~ "Sculpin",
                                   grepl("Catostomus", species) ~ "Suckers",
                                   grepl("Percopsis", species) ~ "Sand roller",
                                   T ~ NA))

write_csv(vi.res, paste0(PATH, "/Output/variable_importance_allspp_nosumm.csv"))

## wide table  ----
vi.wide <- vi.res %>%
  # select(-percent.contribution) %>%
  pivot_wider(id_cols = c(species, species_group), names_from = variable, values_from = permutation.importance)

write_csv(vi.wide, paste0(PATH, "/Output/variable_importance_species_x_permimp.csv"))

## manuscript figure: var imp for each species group ----
var.names <- readxl::read_excel(paste0(PATH, "/env_preds_for_pnw_sdm.xlsx"))
var.names <- var.names %>% rename(variable = `Variable Name`, var_name = `Variable`) %>% select(variable, var_name)
vi.res.summ <- vi.res %>%
  group_by(species_group) %>% mutate(n_spp = n_distinct(species)) %>% ungroup() %>%
  group_by(variable) %>%
  mutate(mn_overall = mean(permutation.importance),
         species_gr_label = paste0(species_group, " (n = ", n_spp, ")")) %>%
  ungroup() %>%
  group_by(species_group, variable) %>%
  summarise(mn_pi = mean(permutation.importance),
            sd_pi = sd(permutation.importance),
            mn_overall = unique(mn_overall),
            species_gr_label = unique(species_gr_label)) %>%
  ungroup() %>%
  left_join(var.names) %>% 
  mutate(var_name = if_else(is.na(var_name), "Mean Sum Annual Flow", var_name),
         species_group = factor(species_group, levels = unique(species_group)),
         var_name = fct_reorder(var_name, mn_overall))

ggplot(data = vi.res.summ) +
  geom_col(aes(x = mn_pi, y = var_name)) +
  geom_errorbar(aes(y = var_name, xmin = pmax(mn_pi - sd_pi, 0), xmax = mn_pi + sd_pi), width = 0.4) +
  facet_wrap(~ species_gr_label) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 75)) +
  labs(x = "Mean permutation importance", y = "") +
  theme_bw()
ggsave(paste0(PATH, "/Output/variable_importance_species_grp_barchart.png"), bg = "white", width = 12, height = 6)

### alternative boxplot ----
vi.box <- vi.res %>%
  left_join(var.names) %>%
  group_by(species_group) %>% mutate(n_spp = n_distinct(species)) %>% ungroup() %>%
  group_by(variable) %>%
  mutate(mn_overall = mean(permutation.importance),
         species_gr_label = paste0(species_group, " (n = ", n_spp, ")"),
         var_name = if_else(is.na(var_name), "Mean Sum Annual Flow", var_name)) %>%
  ungroup() %>%
  mutate(var_name = fct_reorder(var_name, mn_overall))

ggplot(data = vi.box, aes(x = permutation.importance, y = var_name)) +
  geom_boxplot(fill = "lightgrey", outliers = F) +
  geom_point() +
  facet_wrap(~ species_gr_label) +
  labs(x = "Permutation importance", y = "") +
  theme_bw()
ggsave(paste0(PATH, "/Output/variable_importance_species_grp_boxplot.png"), bg = "white", width = 12, height = 6)



## tile plot of variable importance across species ----
ggplot() +
  geom_tile(data = vi.res, aes(x = variable, y = species, fill = permutation.importance)) +
  scale_fill_gradientn(colors = c("gray", viridis::viridis(100)), name = "permutation\nimportance") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(paste0(PATH, "/Output/variable_importance_species_x_permimp.png"), bg = "white", width = 7, height = 8)

## species summary table (species specific top contributing variables) ----
vi.spp <- vi.res %>%
  group_by(species) %>%
  summarise(top_var = variable[which.max(permutation.importance)],
            top_var_value = max(permutation.importance),
            n_vars_above_zero = sum(permutation.importance > 0),
            med_perm_imp = median(permutation.importance),
            sd_perm_imp = sd(permutation.importance))

write_csv(vi.spp, paste0(PATH, "/Output/variable_importance_species_summary.csv"))

## find which variables were most often the most important for a species ----
table(vi.spp$top_var)

## variable summary table (min, max, mean permutation importance) ----
vi.var <- vi.res %>%
  group_by(variable) %>%
  summarise(top_spp = species[which.max(permutation.importance)],
            max_perm_imp = max(permutation.importance),
            min_perm_imp = min(permutation.importance),
            mn_perm_imp = mean(permutation.importance),
            sd_perm_imp = sd(permutation.importance),
            n_zeros_per_spp = sum(permutation.importance == 0))
write_csv(vi.var, paste0(PATH, "/Output/variable_importance_variable_summary.csv"))



# table request for paper: eval results ----
##col1 = species, col2-X = permutation importance for spp x variable, last col = CBI
##read in cbi values
eval.res <- read_csv(paste0(PATH, "/DataRelease/PNWfishes_SDM_model_evaluation_results.csv")) ##from data release code

paper.table <- vi.wide %>%
  left_join(eval.res %>% select(species, CBI_val)) %>%
  arrange(-CBI_val) %>%
  relocate(species, CBI_val) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3)))

write_csv(paper.table, paste0(PATH, "/Output/sdm_result_varimp_summary_table.csv"))












