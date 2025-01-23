# Data preparation

# Load ASVs table from ulrb paper
ASVs_env <- read.csv("./data/AO_ASVs_env_clean.csv", row.names = 1)

# in case you want to replace with GTDB taxonomy
#load("unique_taxonomy_gtdb")
#unique_taxonomy_gtdb_df <- unique_taxonomy_gtdb %>% as.data.frame() %>% mutate(Sequence = rownames(.))
#
#ASVs_gtdb <- ASVs_env %>% 
# select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>% # remove silva taxonomy 
#left_join(unique_taxonomy_gtdb_df, by = "Sequence") %>% 
#select(-Species) # remove species level


# Order Station levels by longitude
Stations_by_long <- ASVs_env %>% 
  select(Station, Longitude, Latitude) %>% 
  distinct() %>% 
  arrange(desc(Longitude)) %>% 
  pull(Station)

# correct format
ASVs_env_clean <- 
  ASVs_env %>% 
  as_tibble() %>% 
  mutate(PelagicLayer = case_when(Depth < 200 ~ "Epipelagic",
                                  Depth < 1000 ~ "Mesopelagic",
                                  Depth >=1000 ~ "Bathypelagic")) %>% 
  mutate(Kingdom = factor(Kingdom),
         Phylum = factor(Phylum),
         Class = factor(Class),
         Order = factor(Order),
         Family = factor(Family),
         Genus = factor(Genus),
         Station = factor(Station, levels = Stations_by_long),
         PelagicLayer = factor(PelagicLayer, 
                               levels = c("Epipelagic", "Mesopelagic", "Bathypelagic")),
         WaterMass = factor(WaterMAss))


# Group ASVs at genus level
taxa_table <- ASVs_env_clean %>% 
  filter(!is.na(Kingdom)) %>% 
  mutate(taxon = 
           str_remove_all(
             paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "_"), c("NA"))) %>%
  mutate(taxon = str_remove(taxon, "__")) %>% 
  group_by(Sample, taxon) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  filter(Abundance > 0) # remove zeros

# get just env data in a separate table
env_data <- ASVs_env_clean %>% 
  select(Sample, year, Station, Depth, 
         PelagicLayer, Longitude, 
         Latitude,PO4, NO2, NO3, 
         NH4, Si, Chl, Temperature, 
         Salinity, WaterMass) %>% 
  mutate(WaterMass = case_when(is.na(WaterMass) ~ "Unknown",
                               WaterMass == "[]" ~ "Unknown",
                               TRUE ~ WaterMass)) %>% 
  mutate(
    across(
      where(is.double), ~round(.x, digits = 3))) %>% 
  distinct()

# join collapsed taxa table with env data
full_table <- taxa_table %>% left_join(env_data)

# Remove outlier sample: M19_88
full_table <- full_table %>% 
  filter(Sample != "M19_88")


taxa_singles <- full_table %>% 
  group_by(taxon) %>% 
  count() %>% 
  summarise(category = ifelse(n == 1, "single sample", "more than 1 sample")) 

# by grouping, only 28% appear in a single sample
taxa_singles %>% 
  pull(category) %>% 
  table()

# verify relation with abundance
singles_df <- full_table %>% 
  left_join(taxa_singles) %>% 
  mutate(category = factor(category, levels = c("single sample", "more than 1 sample"))) %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>%
  mutate(isSingleton = ifelse(Abundance == 1, "yes", "no"))
#
singles_df %>% 
  mutate(isSingleton = factor(isSingleton, levels = c("yes", "no"))) %>% 
  ggplot(aes(Sample, RelativeAbundance, col = isSingleton)) + 
  geom_point() + 
  scale_y_log10() +
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "top") + 
  facet_grid(~category) +
  geom_hline(yintercept = 0.01, col = "red") +
  labs(y = "Relative abundance (%)",
       col = "Is singleton: ",
       x = "Taxa (ordered by sample)") +
  scale_color_manual(values = c("red","grey")) 
#
singles_df %>% 
  mutate(Status = ifelse(RelativeAbundance < 0.01, "remove", "keep")) %>% 
  ungroup() %>% 
  count(Status)
#
2425/13088


total_reads_df <- full_table %>% 
  group_by(Sample, year) %>% 
  summarise(Total_reads = sum(Abundance))

# sample with less than 10 000 reads (to be removed)
samples_to_remove <- total_reads_df %>% filter(Total_reads <10000) %>% pull(Sample)

ASVs_df <- full_table %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>% 
  filter(RelativeAbundance > 0.01,
         !Sample %in% samples_to_remove,
         !taxon %in% c("Bacteria___", "Archaea___"),
         !str_detect(taxon, "Chloroplast"),
         !str_detect(taxon, "Eukaryota___"))

taxa_overview_df <- ASVs_df %>% 
  separate(col = taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))




