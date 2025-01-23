## Sequencing summary

# final reads and ASVs used

ASVs_df %>% 
  group_by(Sample, year) %>% 
  summarise(Reads = sum(Abundance)) %>% 
  group_by(year) %>% 
  summarise(meanReads = mean(Reads),
            sdReads = sd(Reads), n())

ASVs_df %>% 
  group_by(Sample, year) %>% 
  count() %>% 
  group_by(year) %>% 
  summarise(meanASVs = mean(n),
            sd = sd(n))
