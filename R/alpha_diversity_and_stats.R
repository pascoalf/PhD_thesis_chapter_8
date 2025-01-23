# Alpha diversit and stats
# Group by pelagic layer and year
taxa_summary <- taxa_overview_df %>% 
  group_by(Phylum, PelagicLayer, year) %>%
  summarise(Abundance = sum(Abundance)) %>% 
  group_by(PelagicLayer, year) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance))

# for table with top phyla:
(phyla_relative_abundance_summary <- taxa_overview_df %>% 
    group_by(Phylum) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(Percentage = round(Abundance*100/sum(Abundance),10)) %>% 
    arrange(desc(Percentage)))

# save
phyla_relative_abundance_summary %>% write.csv("phya_relative_abundance_summary.csv")

# select top 5 phyla
top_5_phyla <- taxa_overview_df %>% 
  group_by(Phylum) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  arrange(desc(Abundance)) %>%
  head(5) %>% 
  pull(Phylum)

# overview
taxa_overview_df %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(taxonomy = ifelse(Phylum %in% top_5_phyla, Phylum, "Other")) %>%
  mutate(taxonomy = factor(taxonomy, levels = c(top_5_phyla, "Other"))) %>% 
  ggplot(aes(year, RelativeAbundance, col = taxonomy)) + 
  #geom_jitter(width = 0.15) + 
  stat_summary(position = position_dodge(width = 0)) + 
  stat_summary(fun.y = mean, geom = "line") + 
  facet_grid(~PelagicLayer) + 
  scale_color_manual(values = c(qualitative_colors[1:5],"grey82")) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, vjust =0.8),
        legend.position = "top") +
  labs(fill = "Top Phyla: ",
       y = "mean (\U002B sd) of relative abundance (%) \n in Log10 scale",
       x = "")+ 
  scale_y_log10()


# class level summary
class_summary <- taxa_overview_df %>% 
  group_by(Class, PelagicLayer, year) %>%
  summarise(Abundance = sum(Abundance)) %>% 
  group_by(PelagicLayer, year) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance))

# class report for specific values
(class_relative_abundance_summary <- taxa_overview_df %>% 
    group_by(Class) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(Percentage = round(Abundance*100/sum(Abundance),10)) %>% 
    arrange(desc(Percentage)))
#
class_relative_abundance_summary %>% write.csv("class_relative_abundance.csv")

# select top 5 phyla
top_5_class <- taxa_overview_df %>% 
  group_by(Class) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  arrange(desc(Abundance)) %>%
  head(5) %>% 
  pull(Class)

# class level taxonomy
taxa_overview_df %>%
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(taxonomy = ifelse(Class %in% c("", "2", "4", "54", "j", NA), "Unknown", Class)) %>%
  mutate(taxonomy = ifelse(taxonomy %in% c(top_5_class, "Unknown"), taxonomy, "Other")) %>%
  mutate(taxonomy = factor(taxonomy, levels = c(top_5_class, "Unknown","Other"))) %>% 
  ggplot(aes(year, RelativeAbundance, col = taxonomy)) + 
  stat_summary() +
  stat_summary(method = mean, geom = "line") +
  #geom_col(position = "dodge") + 
  facet_grid(~PelagicLayer) + 
  scale_color_manual(values = c(qualitative_colors[1:5],"grey42","grey82")) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, vjust =0.8),
        legend.position = "top") +
  labs(col = "Top Class: ",
       y = "mean (\U002B sd) of relative abundance (%) \n in Log10 scale",
       x = "") + 
  scale_y_log10()

# Quick check of samples with cyanobacteria in bathypelagic layer:
ASVs_df %>% 
  filter(str_detect(taxon, "Cyanobacteria"),
         year == 2019,
         PelagicLayer == "Bathypelagic")

## define rare, undetermined and abundant taxa
ASVs_df <- ASVs_df %>% define_rb() ## the clustering result is now better

# Calculate taxa richness
ASVs_diversity <- ASVs_df %>% 
  group_by(Sample, Classification) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data) # recover environmental data

# total different taxa
ASVs_df$taxon %>% length()

#
ASVs_diversity %>% 
  group_by(Sample, year, PelagicLayer) %>% 
  summarise(taxaRichness = sum(taxaRichness)) %>% 
  ungroup() %>% 
  group_by(year, PelagicLayer) %>% 
  summarise(totalRichness = sum(taxaRichness),
            meanRichness = mean(taxaRichness),
            sdRichness = sd(taxaRichness),
            medianRichness = median(taxaRichness),
            minRichness = min(taxaRichness),
            maxRichness = max(taxaRichness),
            n= n())

# Illustrate overall diversity

ASVs_diversity %>% 
  group_by(Sample, year, PelagicLayer) %>% 
  summarise(taxaRichness = sum(taxaRichness)) %>%
  ggplot(aes(x = as.factor(year), y = taxaRichness)) + 
  geom_boxplot(aes(fill = PelagicLayer)) + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) + 
  labs(y = "Number of taxa",
       x = "Year",
       fill = "Pelagic layer: ") + 
  scale_fill_manual(values = c("#C6DBEF", "#6BAED6", "#08519C")) 

# Inspect outlier sample
ASVs_diversity %>% 
  group_by(Sample, year, PelagicLayer) %>% 
  summarise(taxaRichness = sum(taxaRichness)) %>%
  filter(PelagicLayer == "Bathypelagic")

#ANOVA test for years, grouped by pelagic layers
years_by_pelagic_layer <- ASVs_diversity %>% 
  group_by(Sample, year, PelagicLayer) %>% 
  summarise(taxaRichness = sum(taxaRichness))

# identify outliers
years_by_pelagic_layer %>% 
  group_by(year, PelagicLayer) %>% 
  identify_outliers(taxaRichness)

lm_years_and_pelagic_vs_taxa <- lm(taxaRichness ~ year*PelagicLayer, data = years_by_pelagic_layer)
#
shapiro_test_years_and_layers <- shapiro_test(residuals(lm_years_and_pelagic_vs_taxa))
write.csv(shapiro_test_years_and_layers, "./outputs/shapiro_test_years_and_layers.csv")

# repeat, by groups
shapiro_test_groups <- years_by_pelagic_layer %>% 
  group_by(year, PelagicLayer) %>%filter(PelagicLayer != "Bathypelagic") %>% 
  shapiro_test(taxaRichness) %>% 
  mutate(Significance = ifelse(p < 0.05, "Significant", "ns")) %>% 
  mutate("Normal distribution" = ifelse(Significance == "ns", "Yes", "No"))
#
shapiro_test_groups %>% 
  write.csv("./outputs/shapiro_test_groups.csv")

# homogeneity of variance
years_by_pelagic_layer %>% ungroup() %>% 
  filter(PelagicLayer != "Bathypelagic") %>%
  mutate(year = as.factor(year)) %>% 
  levene_test(formula = taxaRichness ~ year*PelagicLayer)

# Anova 
years_by_pelagic_layer %>% 
  ungroup() %>% 
  mutate(year = as.factor(year)) %>% 
  anova_test(taxaRichness ~ year*PelagicLayer)

# pelagic layer
years_by_pelagic_layer %>% ungroup() %>% 
  pairwise_t_test(taxaRichness ~ PelagicLayer,
                  p.adjust.method = "bonferroni") %>% 
  write.csv("./outputs/pairwise_t_test_pelagic_layers.csv")
# years
years_by_pelagic_layer %>% ungroup() %>% 
  pairwise_t_test(taxaRichness ~ year,
                  p.adjust.method = "bonferroni") %>% 
  write.csv("./outputs/pairwise_t_test_years.csv")

# classifications overview
ASVs_diversity %>% 
  group_by(Sample) %>% 
  mutate(taxaRichness = taxaRichness*100/sum(taxaRichness)) %>%
  ungroup() %>% 
  group_by(Sample, Classification) %>% 
  summarise(taxaRichness = sum(taxaRichness)) %>% 
  ungroup() %>% 
  group_by(Classification) %>% 
  summarise(mean(taxaRichness),
            sd = sd(taxaRichness),
            n = n())

## Temperature ##
gridExtra::grid.arrange(
  ## all taxa
  ASVs_df %>% 
    group_by(Sample) %>% 
    summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
    left_join(env_data) %>% 
    ggplot(aes(Temperature, taxaRichness))+
    geom_point(size = 3) + 
    stat_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) +
    stat_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson")) + 
    labs(y = "Number of taxa",
         x = "Temperature (ºC)") + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)),
  # diversity as a function of temperature
  ASVs_diversity %>% 
    ggplot(aes(Temperature, taxaRichness, col = Classification))+
    geom_point(size = 3) + 
    stat_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) + 
    stat_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson")) + 
    labs(y = "Number of taxa",
         x = "Temperature (ºC)",
         col = "Classification: ") +
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), 
  ncol= 2 )


## Overall community
general_diversity <- ASVs_df %>% 
  group_by(Sample) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data)
#
m1 <- lm(taxaRichness ~ Temperature, data = general_diversity)
summary(m1)
cor.test(general_diversity$taxaRichness, general_diversity$Temperature)
cor.test(general_diversity$taxaRichness, general_diversity$Temperature, 
         method = "kendall")
# glm - poisson
m2 <- glm(taxaRichness ~ Temperature, data = general_diversity, family = "poisson")
summary(m2)
#
lm(taxaRichness ~ Temperature,
   data = ASVs_diversity[ASVs_diversity$Classification == "Rare",]) 
#summary(rb_v_t_model)

# poisson - rare vs temperature
glm(taxaRichness ~ Temperature,
    data = ASVs_diversity[ASVs_diversity$Classification == "Rare",], family = "poisson") %>% 
  summary()

# corr - rare vs temperature - pearson
# 
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Rare",]$taxaRichness, 
         ASVs_diversity[ASVs_diversity$Classification == "Rare",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Rare",]$taxaRichness, 
         ASVs_diversity[ASVs_diversity$Classification == "Rare",]$Temperature,
         method = "kendall")

# linear model
und_vs_t_model <- lm(taxaRichness ~ Temperature, 
                     data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]) 
summary(und_vs_t_model)
# poisson model
glm(taxaRichness ~ Temperature, 
    data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",],
    family = "poisson") %>% 
  summary()

# pearson
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$Temperature,
         method = "kendall")

# Abundant vs temperature
# linear least squares model
abun_vs_t_model <- lm(taxaRichness ~ Temperature, 
                      data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",]) 
summary(abun_vs_t_model)
# poisson
glm(taxaRichness ~ Temperature,
    data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",],
    family = "poisson") %>% 
  summary()

# pearson
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$Temperature,
         method = "kendall")

### richness as a function of temperature, by year
ASVs_df %>% 
  group_by(Sample) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data) %>% 
  ggplot(aes(Temperature, taxaRichness))+
  geom_point(aes(col = as.factor(year)),
             size = 2) + 
  geom_smooth(se = FALSE, method = "lm",
              lty = "dashed", lwd = 0.5) + 
  labs(y = "Number of taxa",
       x = "Temperature (ºC)",
       col = "Year: ") + #,
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  scale_color_manual(values = qualitative_colors)

# year must be interpreted as a categorical variable
general_diversity$catYear <- factor(general_diversity$year, levels = c(2016:2020))
summary(lm(taxaRichness ~ Temperature %in% catYear, general_diversity))

### Check water masses in same context
ASVs_df %>% 
  group_by(Sample) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data) %>% 
  ggplot(aes(Temperature, taxaRichness))+
  geom_point(aes(col = WaterMass), size = 2) + 
  stat_smooth(se = FALSE, method = "lm",
              lty = "dashed", lwd = 0.5) + 
  labs(y = "Number of taxa",
       x = "Temperature (ºC)",
       col = "Water mass: ") +
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  scale_color_manual(values = qualitative_colors)


##  Only water masses with more than 10 samples
selected_water_masses <- ASVs_diversity %>% 
  group_by(WaterMass) %>% 
  count() %>% 
  filter(n >= 10) %>% 
  pull(WaterMass)
selected_water_masses

ASVs_df %>% 
  group_by(Sample) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data) %>% 
  filter(WaterMass %in% selected_water_masses,
         WaterMass != "Unknown") %>%
  ggplot(aes(Temperature, taxaRichness))+
  geom_point(aes(col = WaterMass), size = 2) + 
  stat_smooth(aes(col = WaterMass), se = FALSE, method = "lm",
              lty = "dashed", lwd = 0.5) + 
  stat_smooth(aes(col = WaterMass), se = FALSE, 
              method = "glm", method.args = list(family = "poisson"),lwd = 0.5) + 
  labs(y = "Number of taxa",
       x = "Temperature (ºC)",
       col = "Water mass: ") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  scale_color_manual(values = qualitative_colors)

# AW 
AW_div <- general_diversity %>% 
  filter(WaterMass == "AW")

# lm
AW_div %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
AW_div %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()

# IW
IW_div <- general_diversity %>% 
  filter(WaterMass == "IW")

# lm
IW_div %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
IW_div %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()

# SW
SW_div <- general_diversity %>% 
  filter(WaterMass == "SW")

# lm
SW_div %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
SW_div %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()

# TAW
TAW_div <- general_diversity %>% 
  filter(WaterMass == "TAW")

# lm
TAW_div %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
TAW_div %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()

# WCW
WCW_div <- general_diversity %>% 
  filter(WaterMass == "WCW")

# lm
WCW_div %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
WCW_div %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()

#
test_water_mass <- function(x, water){x %>% filter(WaterMass == water)}

## Pearson tests
cor.test(test_water_mass(general_diversity, water = "AW")$taxaRichness,
         test_water_mass(general_diversity, water = "AW")$Temperature)
#
cor.test(test_water_mass(general_diversity, water = "IW")$taxaRichness,
         test_water_mass(general_diversity, water = "IW")$Temperature)
#
cor.test(test_water_mass(general_diversity, water = "LW")$taxaRichness,
         test_water_mass(general_diversity, water = "LW")$Temperature)
#
cor.test(test_water_mass(general_diversity, water = "SW")$taxaRichness,
         test_water_mass(general_diversity, water = "SW")$Temperature)
#
cor.test(test_water_mass(general_diversity, water = "TAW")$taxaRichness,
         test_water_mass(general_diversity, water = "TAW")$Temperature)
#
#
cor.test(test_water_mass(general_diversity, water = "WCW")$taxaRichness,
         test_water_mass(general_diversity, water = "WCW")$Temperature)

## Kendall tests
cor.test(test_water_mass(general_diversity, water = "AW")$taxaRichness,
         test_water_mass(general_diversity, water = "AW")$Temperature,
         method = "kendall")
#
cor.test(test_water_mass(general_diversity, water = "IW")$taxaRichness,
         test_water_mass(general_diversity, water = "IW")$Temperature,
         method = "kendall")
#
cor.test(test_water_mass(general_diversity, water = "SW")$taxaRichness,
         test_water_mass(general_diversity, water = "SW")$Temperature,
         method = "kendall")
#
cor.test(test_water_mass(general_diversity, water = "TAW")$taxaRichness,
         test_water_mass(general_diversity, water = "TAW")$Temperature,
         method = "kendall")
#
cor.test(test_water_mass(general_diversity, water = "WCW")$taxaRichness,
         test_water_mass(general_diversity, water = "WCW")$Temperature,
         method = "kendall")

### Coast vs open ocean
ASVs_df %>% 
  group_by(Sample) %>% 
  summarise(taxaRichness = vegan::specnumber(Abundance)) %>% 
  left_join(env_data) %>% 
  mutate(Region = case_when(Station %in% c("KB0", "KB3", "KB6",
                                           "KB7", "KB5", "KB2",
                                           "KB1", "KH") ~ "Fjord",
                            Station %in% c("HG-IV", "V12", "V6",
                                           "HG-I", "V10") ~ "Open ocean")) %>% 
  ggplot(aes(Temperature, taxaRichness))+
  geom_point(aes(col = Region),
             size = 2) + 
  stat_smooth(aes(col = Region), se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) + 
  stat_smooth(aes(col = Region), se = FALSE, method = "glm", method.args = list(family = "poisson"), lwd = 0.5) + 
  labs(y = "Number of taxa",
       x = "Temperature (ºC)",
       col = "Region: ") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  scale_color_manual(values = qualitative_colors)


general_diversity <- general_diversity %>%
  mutate(Region = case_when(Station %in% c("KB0", "KB3", "KB6",
                                           "KB7", "KB5", "KB2",
                                           "KB1", "KH") ~ "Fjord",
                            Station %in% c("HG-IV", "V12", "V6",
                                           "HG-I", "V10") ~ "Open ocean"))
# lm
general_diversity %>% 
  filter(Region == "Fjord") %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
#
general_diversity %>% 
  filter(Region == "Open ocean") %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# glm
general_diversity %>% 
  filter(Region == "Fjord") %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()
#
general_diversity %>% 
  filter(Region == "Open ocean") %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()
# Pearson corr
cor.test(general_diversity[general_diversity$Region == "Fjord",]$taxaRichness,
         general_diversity[general_diversity$Region == "Fjord",]$Temperature)
#
cor.test(general_diversity[general_diversity$Region != "Fjord",]$taxaRichness,
         general_diversity[general_diversity$Region != "Fjord",]$Temperature)
# Kendall corr
cor.test(general_diversity[general_diversity$Region == "Fjord",]$taxaRichness,
         general_diversity[general_diversity$Region == "Fjord",]$Temperature, 
         method = "kendall")
#
cor.test(general_diversity[general_diversity$Region != "Fjord",]$taxaRichness,
         general_diversity[general_diversity$Region != "Fjord",]$Temperature, 
         method = "kendall")


# temperature vs year 
# first, select water masses with more than 10 samples, and also remove the unknown
selected_water_masses <- ASVs_diversity %>% 
  group_by(WaterMass) %>% 
  count() %>% 
  filter(n >= 10) %>% 
  pull(WaterMass)
#
ASVs_diversity %>% 
  filter(WaterMass %in% selected_water_masses, WaterMass != "Unknown") %>% 
  ggplot(aes(year,Temperature, col = WaterMass))+
  geom_jitter(size = 2.5, width = 0.15, alpha = 0.49) +
  facet_grid(~PelagicLayer) + 
  geom_smooth(se = FALSE) + 
  #geom_smooth(method = "lm", se = FALSE, lty = "dashed", lwd = 0.5) +
  labs(x = "Year",
       y = "Temperature (ºC)",
       fill = "Pelagic layer: ") + #,
  # title = "Temperature as a function of sampling year, divided by pelagic layer",
  # subtitle = paste("range: 2016-2020")) + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 14, angle = 90),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 16))+
  scale_color_manual(values = qualitative_colors)


gridExtra::grid.arrange(
  general_diversity %>% 
    ggplot(aes(Temperature, taxaRichness, col = PelagicLayer))+
    geom_jitter(size = 2, width = 0.1) +
    geom_smooth(method = "lm", se = FALSE, lty = "dashed", lwd = 0.75) +
    geom_smooth(method = "glm", method.args = list(family = "poisson") , se = FALSE, lwd = 0.755) +
    labs(x = "Temperature (ºC)",
         y = "Number of taxa",
         col = "Pelagic layer: ") + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.background = element_rect(fill = "grey94")) + 
    scale_color_manual(values = c("#C6DBEF", "#6BAED6", "#08519C")),
  
  ## as a function of depth
  general_diversity %>% 
    ggplot(aes(Depth, taxaRichness, col = PelagicLayer))+
    geom_jitter(size = 2, width = 0.1) +
    geom_smooth(method = "lm", se = FALSE, lty = "dashed", lwd = 0.75) +
    geom_smooth(method = "glm", method.args = list(family = "poisson") , se = FALSE, lwd = 0.755) +
    labs(x = "Depth (m)(Log10 scale)",
         y = "Number of taxa") + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.background = element_rect(fill = "grey94")) + 
    scale_color_manual(values = c("#C6DBEF", "#6BAED6", "#08519C")) + 
    scale_x_log10(), ncol = 2) 

### Epipelagic
## Temperature
# lm
general_diversity %>% 
  filter(PelagicLayer == "Epipelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Epipelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()
#
epi_div <-general_diversity %>% filter(PelagicLayer == "Epipelagic") 
# poisson correlation
cor.test(epi_div$taxaRichness, epi_div$Temperature)
# Kendall correlation
cor.test(epi_div$taxaRichness, epi_div$Temperature,
         method = "kendall")

## Depth
# lm
general_diversity %>% 
  filter(PelagicLayer == "Epipelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Depth) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Epipelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Depth, family = "poisson") %>% 
  summary()
#
epi_div <-general_diversity %>% filter(PelagicLayer == "Epipelagic") 
# poisson correlation
cor.test(epi_div$taxaRichness, epi_div$Depth)
# Kendall correlation
cor.test(epi_div$taxaRichness, epi_div$Depth,
         method = "kendall")

### Mesopelagic
## Temperature
# lm 
general_diversity %>% 
  filter(PelagicLayer == "Mesopelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Mesopelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()
#
meso_div <- general_diversity %>% filter(PelagicLayer == "Mesopelagic")
# poisson correlation
cor.test(meso_div$taxaRichness, meso_div$Temperature)
# Kendall correlation
cor.test(meso_div$taxaRichness, meso_div$Temperature,
         method = "kendall")

## Depth
# lm 
general_diversity %>% 
  filter(PelagicLayer == "Mesopelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Depth) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Mesopelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Depth, family = "poisson") %>% 
  summary()
#
meso_div <- general_diversity %>% filter(PelagicLayer == "Mesopelagic")
# poisson correlation
cor.test(meso_div$taxaRichness, meso_div$Depth)
# Kendall correlation
cor.test(meso_div$taxaRichness, meso_div$Depth,
         method = "kendall", exact = FALSE)


### Bathypelagic
## Temperature
# lm 
general_diversity %>% 
  filter(PelagicLayer == "Bathypelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Temperature) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Bathypelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Temperature, family = "poisson") %>% 
  summary()
#
bathy_div <- general_diversity %>% filter(PelagicLayer == "Bathypelagic")
# poisson correlation
cor.test(bathy_div$taxaRichness, bathy_div$Temperature)
# Kendall correlation
cor.test(bathy_div$taxaRichness, bathy_div$Temperature,
         method = "kendall")

## Depth
# lm 
general_diversity %>% 
  filter(PelagicLayer == "Bathypelagic") %>% 
  lm(data = ., formula = taxaRichness ~ Depth) %>% 
  summary()
# poisson
general_diversity %>% 
  filter(PelagicLayer == "Bathypelagic") %>% 
  glm(data = ., formula = taxaRichness ~ Depth, family = "poisson") %>% 
  summary()
#
bathy_div <- general_diversity %>% filter(PelagicLayer == "Bathypelagic")
# poisson correlation
cor.test(bathy_div$taxaRichness, bathy_div$Depth)
# Kendall correlation
cor.test(bathy_div$taxaRichness, bathy_div$Depth,
         method = "kendall", exact = FALSE)

# Divide by pelagic layer AND classifications
# rare
summary(lm(taxaRichness ~ Temperature %in% PelagicLayer, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Rare",]))


# undetermined
summary(lm(taxaRichness ~ Temperature %in% PelagicLayer, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]))
# abundant
summary(lm(taxaRichness ~ Temperature %in% PelagicLayer, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",]))

# diversity as a function of years
ASVs_diversity %>% 
  ggplot(aes(year, taxaRichness, col = Classification))+
  geom_jitter(width = 0.15) + 
  stat_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 1.25, alpha = 0.5) + 
  stat_smooth(se = FALSE, method = "glm",  
              method.args = list(family = "poisson"),
              lwd = 0.5) + 
  labs(y = "Number of taxa",
       x = "Years",
       col = "Classification: ") + 
  theme_classic() + 
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title  = element_text(size = 12),
        axis.title = element_text(size = 12)) + 
  facet_grid(~PelagicLayer)+
  scale_color_manual(values = qualitative_colors[1:3])

##
general_diversity$corrected_year <- general_diversity$year - 2016

# statistical tests for richness vs year by pelagic layer
summary(lm(taxaRichness ~ corrected_year, data = general_diversity))
summary(glm(taxaRichness ~ corrected_year, data = general_diversity, family = "poisson"))
cor.test(general_diversity$taxaRichness, general_diversity$corrected_year)
cor.test(general_diversity$taxaRichness, general_diversity$corrected_year, method = "kendall")
#

# correct year
ASVs_diversity$year_corrected <- ASVs_diversity$year - 2016
# rare
## linear model
ASVs_diversity %>% filter(PelagicLayer == "Epipelagic", Classification == "Rare") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
#
ASVs_diversity %>% filter(PelagicLayer == "Mesopelagic", Classification == "Rare") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
#
ASVs_diversity %>% filter(PelagicLayer == "Bathypelagic", Classification == "Rare") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()

ASVs_diversity %>% filter(PelagicLayer == "Epipelagic", Classification == "Rare") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
#
ASVs_diversity %>% filter(PelagicLayer == "Mesopelagic", Classification == "Rare") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
#
ASVs_diversity %>% filter(PelagicLayer == "Bathypelagic", Classification == "Rare") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()

### Correlation tests
# rare
rares_epi <- ASVs_diversity %>% filter(PelagicLayer == "Epipelagic", Classification == "Rare")
cor.test(rares_epi$taxaRichness, rares_epi$year)
cor.test(rares_epi$taxaRichness, rares_epi$year, method = "kendall")
rares_meso <- ASVs_diversity %>% filter(PelagicLayer == "Mesopelagic", Classification == "Rare")
cor.test(rares_meso$taxaRichness, rares_meso$year)
cor.test(rares_meso$taxaRichness, rares_meso$year, method = "kendall")
rares_bathy <- ASVs_diversity %>% filter(PelagicLayer == "Bathypelagic", Classification == "Rare")
cor.test(rares_bathy$taxaRichness, rares_bathy$year)
cor.test(rares_bathy$taxaRichness, rares_bathy$year, method = "kendall")

# undetermined
# epipelagic layer
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Epipelagic") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
# meso
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Mesopelagic") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
# bathy
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Bathypelagic") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()

# Poisson model for undetermined
# epipelagic layer
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Epipelagic") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
# meso
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Mesopelagic") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
# bathy
ASVs_diversity %>% 
  filter(Classification == "Undetermined", PelagicLayer == "Bathypelagic") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()

## correlations
und_epi <- ASVs_diversity %>% filter(PelagicLayer == "Epipelagic", Classification == "Undetermined")
cor.test(und_epi$taxaRichness, und_epi$year)
cor.test(und_epi$taxaRichness, und_epi$year, method = "kendall")
und_meso <- ASVs_diversity %>% filter(PelagicLayer == "Mesopelagic", Classification == "Undetermined")
cor.test(und_meso$taxaRichness, und_meso$year)
cor.test(und_meso$taxaRichness, und_meso$year, method = "kendall", exact = FALSE)
und_bathy <- ASVs_diversity %>% filter(PelagicLayer == "Bathypelagic", Classification == "Undetermined")
cor.test(und_bathy$taxaRichness, und_bathy$year)
cor.test(und_bathy$taxaRichness, und_bathy$year, method = "kendall")

# abundant
# epipelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Epipelagic", Classification == "Abundant") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
# mesopelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Mesopelagic", Classification == "Abundant") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()
# bathypelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Bathypelagic", Classification == "Abundant") %>% 
  lm(data = ., formula = taxaRichness ~ year_corrected) %>% 
  summary()

# abundant
# epipelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Epipelagic", Classification == "Abundant") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
# mesopelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Mesopelagic", Classification == "Abundant") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()
# bathypelagic
ASVs_diversity %>% 
  filter(PelagicLayer == "Bathypelagic", Classification == "Abundant") %>% 
  glm(data = ., formula = taxaRichness ~ year_corrected, family = "poisson") %>% 
  summary()

## correlations
abun_epi <- ASVs_diversity %>% filter(PelagicLayer == "Epipelagic", Classification == "Abundant")
cor.test(abun_epi$taxaRichness, abun_epi$year)
cor.test(abun_epi$taxaRichness, abun_epi$year, method = "kendall")
abun_meso <- ASVs_diversity %>% filter(PelagicLayer == "Mesopelagic", Classification == "Abundant")
cor.test(abun_meso$taxaRichness, abun_meso$year)
cor.test(abun_meso$taxaRichness, abun_meso$year, method = "kendall")
abun_bathy <- ASVs_diversity %>% filter(PelagicLayer == "Bathypelagic", Classification == "Abundant")
cor.test(abun_bathy$taxaRichness, abun_bathy$year)
cor.test(abun_bathy$taxaRichness, abun_bathy$year, method = "kendall")

# diversity as a function of salinity
gridExtra::grid.arrange(
  general_diversity %>% 
    ggplot(aes(Salinity, taxaRichness))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Salinity (PSU)",
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20))
  ,
  ASVs_diversity %>% 
    ggplot(aes(Salinity, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Salinity (PSU)",
         col = "Classification:"
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), ncol = 2)


## rare
summary(lm(taxaRichness ~ Salinity, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Rare",]))
## undetermined
summary(lm(taxaRichness ~ Salinity, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]))
## abundant
summary(lm(taxaRichness ~ Salinity, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",]))


# diversity as a function of NO2
gridExtra::grid.arrange(
  general_diversity %>% 
    filter(NO2 > 0) %>% 
    ggplot(aes(NO2, taxaRichness))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Nitrite [\U00B5M]",
         #title = "NO2 (by classification)",
         #subtitle = "range: 2016-2020\nn=117", 
         #sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NO2")])$NO2))
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20)), 
  ASVs_diversity %>%
    filter(NO2 > 0) %>% 
    ggplot(aes(NO2, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Nitrite [\U00B5M]",
         col = "Classification:"
         # title = "NO2 (by classification)",
         #subtitle = "range: 2016-2020\nn=117", 
         #sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NO2")])$NO2))
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), 
  ncol = 2)


## rare
summary(lm(taxaRichness ~ NO2, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Rare", ]))
## undetermined
summary(lm(taxaRichness ~ NO2, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined", ]))
## abundant
summary(lm(taxaRichness ~ NO2, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Abundant", ]))

# diversity as a function of NO3
gridExtra::grid.arrange(
  general_diversity %>% 
    filter(NO3>0) %>% 
    ggplot(aes(NO3, taxaRichness))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Nitrate [\U00B5M]",
         #title = "NO3 by classification",
         #subtitle = paste0("range: 2016-2020\nn=", 
         #                 sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NO3")])$NO3)))
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20)),
  ASVs_diversity %>% 
    filter(NO3>0) %>% 
    ggplot(aes(NO3, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Nitrate [\U00B5M]",
         col = "Classification:"
         #title = "NO3 by classification",
         #subtitle = paste0("range: 2016-2020\nn=", 
         #                 sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NO3")])$NO3)))
    ) + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), ncol = 2)

## rare
summary(lm(taxaRichness ~ NO3, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Rare", ]))
## undetermined
summary(lm(taxaRichness ~ NO3, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined", ]))
## abundant
summary(lm(taxaRichness ~ NO3, 
           data = ASVs_diversity[ASVs_diversity$Classification == "Abundant", ]))

# diversity as a function of NH4
gridExtra::grid.arrange(
  general_diversity %>% 
    filter(NH4 > 0) %>% 
    ggplot(aes(NH4, taxaRichness))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Ammonium [\U00B5M]",
         #title = "NH4",
         #subtitle = paste0("range: 2016-2020\nn=", 
         #                 sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NH4")])$NH4)))
    ) + # count valid samples (after NaN removal) 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20)),  
  ASVs_diversity %>% 
    filter(NH4 > 0) %>% 
    ggplot(aes(NH4, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    labs(y = "Number of taxa",
         x = "Ammonium [\U00B5M]",
         col = "Classification:"
         #title = "NH4",
         #subtitle = paste0("range: 2016-2020\nn=", 
         #                 sum(!is.na(distinct(ASVs_diversity[, c("Sample", "NH4")])$NH4)))
    ) + # count valid samples (after NaN removal) 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), ncol = 2
)

# diversity as a function of Si
gridExtra::grid.arrange(
  general_diversity %>% 
    filter(Si > 0) %>% 
    ggplot(aes(Si, taxaRichness))+
    geom_point(size = 2) + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) +
    geom_smooth(se = FALSE, method = "glm", lwd = 0.5, method.args = list(family = "poisson"))+
    labs(y = "Number of taxa",
         x = "Silicates [\U00B5M]") + 
    theme_classic() + 
    theme(#legend.position = "top",
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)),
  ASVs_diversity %>% 
    filter(Si > 0) %>% 
    ggplot(aes(Si, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) +
    geom_smooth(se = FALSE, method = "glm", lwd = 0.5, method.args = "poisson") +
    labs(y = "Number of taxa",
         x = "Silicates [\U00B5M]",
         color = "Classification:") + 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]), 
  ncol = 2)


# overall 
# lm
summary(lm(taxaRichness ~ Si, data = general_diversity))

# poisson
summary(glm(taxaRichness ~ Si, data = general_diversity, family = "poisson"))
#

#plot(lm(taxaRichness ~ Si, data = general_diversity))
# pearson 
cor.test(general_diversity$taxaRichness, general_diversity$Si)
# Kendall
cor.test(general_diversity$taxaRichness, general_diversity$Si, 
         method = "kendall")

## rare
# lm
summary(lm(taxaRichness ~ Si, data = ASVs_diversity[ASVs_diversity$Classification == "Rare",]))
# poison
summary(glm(taxaRichness ~ Si, 
            data = ASVs_diversity[ASVs_diversity$Classification == "Rare",],
            family = "poisson"))

# pearson
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Rare",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Rare",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Rare",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Rare",]$Temperature,
         method = "kendall")

## undetermined
# lm
summary(lm(taxaRichness ~ Si, data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]))
# poisson
summary(glm(taxaRichness ~ Si, 
            data = ASVs_diversity[ASVs_diversity$Classification == "Undetermined",],
            family = "poisson"))

# pearson
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Undetermined",]$Temperature,
         method = "kendall")

## abundant
# pearson
summary(lm(taxaRichness ~ Si, data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",]))
# Kendall
summary(glm(taxaRichness ~ Si, 
            data = ASVs_diversity[ASVs_diversity$Classification == "Abundant",],
            family = "poisson"))

# pearson
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$Temperature)
# Kendall
cor.test(ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$taxaRichness,
         ASVs_diversity[ASVs_diversity$Classification == "Abundant",]$Temperature,
         method = "kendall")

#
general_diversity %>% 
  filter(Si > 0) %>% 
  ggplot(aes(Depth, Si)) +
  geom_point(size = 2) +
  theme_classic() + 
  labs(x = "Depth (m) in Log10 scale",
       y = "Silicates [\U00B5M]") + 
  scale_x_log10() +
  theme(legend.position = "top",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# diversity as a function of Chl
gridExtra::grid.arrange(
  general_diversity %>% 
    filter(Chl > 0) %>% 
    ggplot(aes(Chl, taxaRichness))+
    geom_point(size = 2) + 
    #  geom_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) +
    labs(y = "Number of taxa",
         x = "Chl a [mg/m]") + #,
    #       title = "Taxonomic richness as a function of Chl",
    #       subtitle = paste0("range: 2016-2020\nn=", 
    #                         sum(!is.na(distinct(ASVs_diversity[, c("Sample", "Chl")])$Chl)))) + # count valid samples (after NaN removal) 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)),
  ASVs_diversity %>% 
    filter(Chl > 0) %>% 
    ggplot(aes(Chl, taxaRichness, col = Classification))+
    geom_point(size = 2) + 
    #  geom_smooth(se = FALSE, method = "lm", lty = "dashed", lwd = 0.5) +
    labs(y = "Number of taxa",
         x = "Chl a [mg/m³]",
         col = "Classification: ") + 
    #      title = "Taxonomic richness as a function of Chl",
    #       subtitle = paste0("range: 2016-2020\nn=", 
    #                        sum(!is.na(distinct(ASVs_diversity[, c("Sample", "Chl")])$Chl)))) + # count valid samples (after NaN removal) 
    theme_classic() + 
    theme(legend.position = "top",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    scale_color_manual(values = qualitative_colors[1:3]),
  ncol = 2)
