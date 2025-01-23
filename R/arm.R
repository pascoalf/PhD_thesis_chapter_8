# association rule mining 
set.seed(123)
# load env and ASVs data
env_cat <- readRDS("./data/environmental_data_category.rds")
env_data <- readRDS("./data/env_data.rds")

ASVs_df <- readRDS("./data/ASV_clean_full_df.rds")

# remove outlier sample (mislabeled, we don't know real metadata)
ASVs_df <- ASVs_df %>% filter(Sample != "M19_88")

# merge relevant data
ASVs_cat_df <- ASVs_df %>% 
  select(Sample, Classification, taxon) %>% 
  mutate(Sample = factor(Sample),
         taxon = factor(taxon)) %>% 
  left_join(env_cat) %>% ## env_cat includes the discretized numerical variables (made in "env_data_discretization.R")
  ungroup() %>% 
  select(-Sample, -Station) # discuss this (features that are redundant)

# make transactions
transactions_full <- ASVs_cat_df %>% transactions()

#
summary(transactions_full)

# get rhs items (set the consequent)
classifications_rhs <- grep("Classification=", 
                            itemLabels(transactions_full), 
                            value = TRUE)

# mine rules
full_rules <- apriori(transactions_full,
                      parameter = list(support = 0.001,minlen = 2, maxlen = 14),
                      appearance = list(rhs = classifications_rhs))

#plot(full_rules)
reds <- RColorBrewer::brewer.pal(9, "Reds")

#
full_rules_df <- DATAFRAME(full_rules)

#
full_rules_df %>% 
  ggplot(aes(support, confidence, fill = lift)) + 
  geom_jitter(shape = 21, col = "grey", size = 2) + 
  scale_fill_gradient(low = reds[1], high = reds[9]) + 
  # labs(title = paste(length(full_rules_df[,1]), "rules")) + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  labs(x = "Conviction",
       y = "Confidence",
       fill = "Lift: ")

# alternative
full_rules_df %>% 
  ggplot(aes(support, confidence, col = lift)) + 
  geom_jitter(size = 2) + 
  scale_color_gradient(low = reds[1], high = reds[9]) + 
  # labs(title = paste(length(full_rules_df[,1]), "rules")) + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.background = element_rect(fill = "grey80")) + 
  labs(x = "Conviction",
       y = "Confidence",
       fill = "Lift: ")

full_rules %>% summary()


# nutrients reference values
env_data %>% 
  pivot_longer(cols = c("PO4", "NO2", "NO3",
                        "NH4", "Si", "Chl", "Salinity"),
               names_to = "Metric",
               values_to = "Score") %>% 
  group_by(Metric) %>% 
  summarise(`very low` = paste0("]0,",
                                round(quantile(Score, na.rm = TRUE)[2],2),"]"),
            low = paste0("]", 
                         round(quantile(Score, na.rm = TRUE),2)[2],",",
                         round(median(Score, na.rm = TRUE),2), "]"),
            high = paste0("]", round(median(Score, na.rm = TRUE),2),
                          ",", round(quantile(Score, na.rm = TRUE),2)[4],
                          "]"),
            `very high` = paste0("]", 
                                 round(quantile(Score, na.rm = TRUE),2)[4], ",", 
                                 round(max(Score, na.rm = TRUE),2), "]")) %>% 
  knitr::kable()

# temperature reference values
env_data %>% 
  pivot_longer(cols = c("Temperature"),
               names_to = "Metric",
               values_to = "Score") %>% 
  group_by(Metric) %>% 
  summarise(`very low` = paste0("]", round(min(Score, na.rm = TRUE),2),",",
                                round(quantile(Score, na.rm = TRUE)[2],2),"]"),
            low = paste0("]", 
                         round(quantile(Score, na.rm = TRUE),2)[2],",",
                         round(median(Score, na.rm = TRUE),2), "]"),
            high = paste0("]", round(median(Score, na.rm = TRUE),2),
                          ",", round(quantile(Score, na.rm = TRUE),2)[4],
                          "]"),
            `very high` = paste0("]", 
                                 round(quantile(Score, na.rm = TRUE),2)[4], ",", 
                                 round(max(Score, na.rm = TRUE),2), "]")) %>% 
  knitr::kable()


# check top10 by lift
top100_rules_by_lift <- full_rules_df %>% 
  arrange(desc(lift)) %>% 
  head(100) 
# print in report
top100_rules_by_lift %>% 
  knitr::kable()
# save to file
top100_rules_by_lift %>% 
  write.csv("top100_rules_by_lift.csv")


# Apply independence test
full_rules_df$chi_p_adj <- 
  full_rules %>% 
  interestMeasure(measure = "chiSquared", 
                  significance = TRUE) %>% 
  p.adjust(method = "bonferroni") ## adjust p-values

# Calculate conviction
full_rules_df$Conviction <- full_rules %>% 
  interestMeasure(measure = "Conviction")

# update dataframe
full_rules_df_dependent <- 
  full_rules_df %>% 
  mutate(Dependence = ifelse(chi_p_adj <0.05, "Dependent", "Independent")) %>% # add test result 
  filter(Dependence == "Dependent")# select dependent rules

#
quality(full_rules)$conviction <- full_rules_df$Conviction
quality(full_rules)$chi_p_adj <- full_rules_df$chi_p_adj

#
full_rules_dependent <- full_rules[quality(full_rules)$chi_p_adj<0.05,]

#
plot(full_rules_dependent, 
     measure = c("conviction", "confidence"), 
     shading = "lift")


full_rules_dependent %>% 
  DATAFRAME() %>% 
  filter(!is.infinite(conviction)) %>% 
  ggplot(aes(conviction, confidence, col = lift)) + 
  geom_point(size = 4) + 
  scale_color_gradient(low = reds[1], high = reds[9]) + 
  labs(x = "Conviction",
       y = "Confidence",
       col = "Lift: ") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.background = element_rect(fill = "grey90")) + 
  lims(y = c(0.80, 1.0))

dependent_rules_high_lift <- full_rules_dependent %>% 
  DATAFRAME() %>% 
  filter(lift > 15) %>% 
  arrange(desc(conviction))
#
dependent_rules_high_lift %>% 
  knitr::kable()
#
dependent_rules_high_lift %>% 
  write.csv("dependent_rules_high_lift.csv")

## repeat, but for any value of Lift
dependent_rules_all <- full_rules_dependent %>% 
  DATAFRAME() %>% 
  arrange(desc(conviction))
#
dependent_rules_all %>% 
  knitr::kable()
#
dependent_rules_all %>% 
  write.csv("dependent_rules_all.csv")


quality(full_rules_dependent)$mutualInfo <- interestMeasure(full_rules_dependent,
                                                            measure = "mutualInformation")
quality(full_rules_dependent)$improvement <- interestMeasure(full_rules_dependent,
                                                             measure = "improvement")
#
plot(full_rules_dependent, 
     measure = c("mutualInfo", "confidence"), 
     shading = "improvement")
#
gridExtra::grid.arrange(
  # improvement
  quality(full_rules_dependent) %>% 
    ggplot(aes(conviction, confidence))+
    geom_point(aes(size = improvement, col = lift))+
    scale_color_gradient(low = reds[1], high = reds[9]) + 
    theme_classic() + 
    theme(#legend.position = "top",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)) + 
    labs(x = "Conviction",
         y = "Confidence",
         col = "Lift: ",
         size = "Improvement: ")# + 
  #guides(size  = FALSE)
  ,
  # mutual information
  quality(full_rules_dependent) %>% 
    ggplot(aes(conviction, confidence))+
    geom_point(aes(size = mutualInfo, col = lift))+
    scale_color_gradient(low = reds[1], high = reds[9]) + 
    theme_classic() + 
    theme(#legend.position = "top",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)) + 
    labs(x = "Conviction",
         y = "Confidence",
         col = "Lift: ",
         size = "Mutual\nInformation: ") + 
    guides(col = FALSE), 
  nrow=2
)


quality(full_rules_dependent) %>% 
  filter(!is.infinite(conviction)) %>% 
  pivot_longer(cols = c("mutualInfo", "improvement"),
               values_to = "Score",
               names_to = "Metric") %>% 
  mutate(Metric = ifelse(Metric == "mutualInfo", "Mutual Information", "Improvement")) %>% 
  ggplot(aes(conviction, 
             confidence, 
             col = lift,
             size = Score)) + 
  geom_point() + 
  scale_color_gradient(low = reds[1], high = reds[9]) + 
  theme_classic() + 
  facet_grid(~Metric) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "grey89")) + 
  labs(x = "Conviction",
       y = "Confidence",
       col = "Lift: ",
       size = "Score: ")
#

full_rules_dependent %>% 
  DATAFRAME() %>% 
  filter(lift > 15) %>% 
  arrange(desc(improvement)) %>% 
  knitr::kable()

# extract taxonomic side of antecedent
full_rules_dependent %>% 
  DATAFRAME() %>%
  filter(lift > 10) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = factor(str_remove(taxon, "\\{taxon="))) %>% 
  group_by(taxon) %>% 
  arrange(desc(improvement)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  knitr::kable()

full_rules_dependent %>% 
  DATAFRAME() %>%
  filter(str_detect(LHS, 
                    "SAR11 clade_Clade I_Clade Ia")) %>% 
  arrange(desc(improvement)) %>% 
  select(LHS, RHS, support,	confidence,	coverage,	lift,	count,	conviction, improvement) %>% 
  knitr::kable()

## dependent rules - nonredundant
# extract taxonomic side of antecedent
dependent_rules_non_redundant <- full_rules_dependent %>% 
  DATAFRAME() %>%
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = factor(str_remove(taxon, "\\{taxon="))) %>% 
  group_by(taxon) %>% 
  arrange(desc(improvement)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) 

dependent_rules_non_redundant %>%
  knitr::kable()

# save for supplementary table 
dependent_rules_non_redundant %>% 
  write.csv("./outputs/dependent_rules_non_redundant.csv")

# add mutual information and improvement
quality(full_rules)$mutualInfo <- interestMeasure(full_rules, measure = "mutualInformation")
quality(full_rules)$improvement <- interestMeasure(full_rules, measure = "improvement")

# n
full_rules %>% 
  DATAFRAME() %>% 
  filter(!is.infinite(conviction)) %>%
  dim()
# improvement
full_rules %>% 
  DATAFRAME() %>% 
  filter(!is.infinite(conviction)) %>%
  ggplot(aes(conviction, confidence, col = lift))+
  #geom_jitter(width = 0.1, size = 2) + 
  geom_point(size = 2) +
  scale_color_gradient(low = reds[1], high = reds[9]) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.background = element_rect(fill = "grey85")) + 
  labs(x = "Conviction",
       y = "Confidence",
       col = "Lift: ")


top100_all_high_quality_by_improvement <- full_rules %>% 
  DATAFRAME() %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>% 
  arrange(desc(improvement)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% 
  top_n(100) 
#
top100_all_high_quality_by_improvement %>% #only best 100
  knitr::kable()
#
top100_all_high_quality_by_improvement %>% 
  write.csv("top100_all_high_quality_by_improvement.csv")

top100_all_non_redundant_by_complexity <- 
  full_rules %>% 
  DATAFRAME() %>% 
  mutate(LHSsize = str_count(LHS, ",")) %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>%
  arrange(desc(LHSsize)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% #use lift again to order the rules by relevance
  top_n(100)

top100_all_non_redundant_by_complexity %>% #only best 100
  knitr::kable()
top100_all_non_redundant_by_complexity %>% 
  write.csv("top100_all_non_redundant_by_complexity.csv")

# View all rules
# hard to read
full_rules %>% 
  plot(method = "grouped", measure = "conviction")


full_rules_dependent %>% 
  head(by = "lift", n = 100) %>% 
  plot(method = "grouped", 
       measure = "conviction")

full_rules_dependent %>% 
  DATAFRAME() %>% 
  ggplot(aes(RHS, lift))+
  geom_boxplot()+ 
  theme_classic() +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_y_continuous(breaks = seq(1,25, by = 4)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) + 
  labs(y = "Lift",
       x = "Classification") + 
  coord_flip()

## make the specific table for paper discussion
full_rules %>% 
  DATAFRAME() %>% 
  mutate(LHSsize = str_count(LHS, ",")) %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>%
  arrange(desc(LHSsize)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% #use lift again to order the rules by relevance
  top_n(100) %>% #only best 100
  knitr::kable()

full_rules_dependent %>% 
  DATAFRAME() %>% 
  #filter(!is.infinite(conviction)) %>% 
  mutate(LHSsize = str_count(LHS, ",")) %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>%
  arrange(desc(LHSsize)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support, confidence, coverage,   lift,   count,  conviction, mutualInfo, improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% 
  knitr::kable()

# save
full_rules_dependent %>% 
  DATAFRAME() %>% 
  #filter(!is.infinite(conviction)) %>% 
  mutate(LHSsize = str_count(LHS, ",")) %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>%
  arrange(desc(LHSsize)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% #use lift again to order the rules by relevance
  #top_n(100) %>%
  write.table("fine-grained_taxonomy.txt")


full_rules %>% DATAFRAME() %>% 
  filter(str_detect(LHS, "Cyanobacteria")) %>% 
  mutate(LHSsize = str_count(LHS, ",")) %>% 
  filter(conviction > 5, confidence > 0.9, lift > 1) %>% 
  separate_wider_delim(LHS, ",", names = "taxon", too_many = "debug") %>% 
  select(-LHS_ok, -LHS_pieces) %>% 
  mutate(taxon = str_remove(taxon, "\\{taxon="),
         taxon = str_remove(taxon, "\\}"),
         taxon = str_remove(taxon, ","),
         taxon = factor(taxon)) %>%
  group_by(taxon) %>%
  arrange(desc(LHSsize)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(LHS, RHS, support,	confidence, coverage,	lift,	count,	conviction,	mutualInfo,	improvement) %>% 
  filter(str_detect(LHS, "taxon")) %>% # to focus only on LHS with taxonomy
  arrange(desc(lift)) %>% 
  knitr::kable()



