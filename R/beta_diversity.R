## beta diversity overview
# load ASV table
ASVs_df <- readRDS("./data/ASV_clean_full_df.rds")
ASVs_df <- ASVs_df %>% filter(Sample != "M19_88") 
#
abundance_table_df <- ASVs_df %>% 
  select(Sample, taxon, Abundance)

#
abundance_table_wide <- abundance_table_df %>% 
  ungroup() %>% 
  pivot_wider(names_from = "taxon",
              values_from = "Abundance") %>% 
    mutate(placeholder = paste(Sample, Classification)) %>% 
  select(-Classification, -Sample)
#
matrix_sample_names <- abundance_table_wide$placeholder
#
abundance_matrix <- abundance_table_wide %>% select(-placeholder) %>% as.matrix()
#
rownames(abundance_matrix) <- matrix_sample_names

# rmeove NAs
abundance_matrix[is.na(abundance_matrix)] <- 0

# prepare env aesthetics
env_aes <- ASVs_df %>% select(Classification, Sample, year,
                              Station, Depth, PelagicLayer, 
                              PO4, NO2, NO3, NH4, Si, WaterMass,
                              Salinity, Salinity, Chl) %>% 
  mutate(year_col = case_when(year == 2016 ~ qualitative_colors[1],
                              year == 2017 ~ qualitative_colors[2],
                              year == 2018 ~ qualitative_colors[3],
                              year == 2019 ~ qualitative_colors[4],
                              year == 2020 ~ qualitative_colors[5]),
         water_col = case_when(WaterMass == "AW" ~ qualitative_colors[1],
                               WaterMass == "TAW" ~ qualitative_colors[2],
                               WaterMass == "SW" ~ qualitative_colors[3],
                               WaterMass == "Unknown" ~ "grey81",
                               WaterMass == "WCW" ~ qualitative_colors[4],
                               WaterMass == "IW" ~ qualitative_colors[5],
                               WaterMass == "LW" ~ qualitative_colors[6]),
         pelagic_color = case_when(PelagicLayer == "Epipelagic" ~ "#C6DBEF",
                                   PelagicLayer == "Mesopelagic" ~ "#6BAED6",
                                   PelagicLayer == "Bathypelagic" ~ "#08519C"),
         PelagicLayer_shape = case_when(PelagicLayer == "Epipelagic" ~ 21,
                                  PelagicLayer == "Mesopelagic" ~ 22,
                                  PelagicLayer == "Bathypelagic" ~ 24),
         Classification_color = case_when(Classification == "Rare" ~ qualitative_colors[1],
                                          Classification == "Undetermined" ~ qualitative_colors[2],
                                          Classification == "Abundant" ~ qualitative_colors[3]),
         placeholder = paste(Sample, Classification))

#
env_aes <- env_aes %>% distinct() %>% filter(placeholder %in% rownames(abundance_matrix))

#
dim(env_aes)
dim(abundance_matrix)

#nMDS
nMDS <- metaMDS(abundance_matrix)

##
## plain analysis of classifications ##
classification_color_labels <- data.frame(Classification = c("Rare", "Undetermined", "Abundant"),
                                          Color = qualitative_colors[1:3])
pelagic_layers_labels <- data.frame(PelagicLayer = c("Epipelagic", "Mesopelagic", "Bathypelagic"),
                                    Shape = c(21, 22, 24),
                                    Color = c("#C6DBEF", "#6BAED6", "#08519C"))

### beta diversity for overall community ###
abundance_table_wide_overall <- abundance_table_df %>%
  group_by(Sample, taxon) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "taxon",
              values_from = "Abundance")
#
matrix_sample_names_overall <- abundance_table_wide_overall$Sample
#
abundance_matrix_overall <- abundance_table_wide_overall %>% 
  select(-Sample) %>% 
  as.matrix()
#
rownames(abundance_matrix_overall) <- matrix_sample_names_overall

# rmeove NAs
abundance_matrix_overall[is.na(abundance_matrix_overall)] <- 0

# calculate nMDS
overall_nmds <- metaMDS(abundance_matrix_overall) 

# env info for overall plot
env_overall <- env_aes %>% 
  ungroup() %>% 
  select(-Classification, -placeholder, -Classification_color) %>% 
  distinct() %>% 
  filter(Sample %in% rownames(abundance_matrix_overall))

# legend for years
year_color_labels <- data.frame(years = 2016:2020,
                                Color = qualitative_colors[1:5])

par(mfrow = c(1,2))
# plot
plot(overall_nmds$points, 
     display = "sites", 
     type = "p",
     xlab = "nMDS", 
     ylab = "nMDS")
#
points(overall_nmds, 
       bg = env_overall$pelagic_color,
       pch = 21,
       #pch =env_aes$PelagicLayer_shape,
       col = "grey", 
       cex = 1.5)
#
with(env_overall,
     ordihull(overall_nmds, PelagicLayer, 
                 lty = "dashed", 
                 label = FALSE))
with(pelagic_layers_labels, 
     legend("topleft", legend = PelagicLayer, bty = "n",
            col = Color, pch = 19, cex = 1,
            title = "Pelagic Layer"))

## check years
plot(overall_nmds$points, 
     display = "sites", 
     type = "p",
     xlab = "nMDS", 
     ylab = "nMDS")
#
points(overall_nmds, 
       bg = env_overall$year_col,
       pch = 21,
       #pch =env_aes$PelagicLayer_shape,
       col = "grey", 
       cex = 1.5)
#
with(year_color_labels, 
     legend("topleft", legend = years, bty = "n",
            col = Color, pch = 19, cex = 1,
            title = "Year"))
##

