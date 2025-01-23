## study of discretization of quantitative variables
library(dplyr)
library(stringr)


## temperature
env_data %>% summary()

## Categories - Temperature and salinity:
# very low TºC (min to 1º quartile);
# low TºC (1º to median)
# median TºC (median) ## removed
# high TºC (median to 3º)
# very high TºC (3º quartile to max)

## Nutrients
# Unknwon (if NA)
# Absent (if zero)
# very low TºC (min to 1º quartile);
# low TºC (1º to median)
# median TºC (median) ## removed
# high TºC (median to 3º)
# very high TºC (3º quartile to max)

## Depth by pelagic layers (also, there is another column with water mass, so this should be enough)

## make categorical env_data
#
nutrient_transformer <- function(x, y = "label"){
  case_when(is.na(x) ~ "Unknown",
            x == 0 ~ "Undetected",
            x <= quantile(x, na.rm = TRUE)[2] ~ "very low",
            x <= median(x, na.rm = TRUE) ~ "low",
            x <= quantile(x, na.rm = TRUE)[4] ~ "high",
            x > quantile(x, na.rm = TRUE)[4] ~ "very high")  
}

env_cat <- env_data %>% 
  mutate(NO2 = nutrient_transformer(NO2),
         NO3 = nutrient_transformer(NO3),
         NH4 = nutrient_transformer(NH4),
         Si = nutrient_transformer(Si),
         Chl = nutrient_transformer(Chl),
         Salinity = nutrient_transformer(Salinity),
         PO4 = nutrient_transformer(PO4)) %>% 
  mutate(Temperature = case_when(is.na(Temperature) ~ "Unknown", # temperature is transformed differently
                                 Temperature <= quantile(Temperature, na.rm = TRUE)[2] ~ "very low",
                                 Temperature <= median(Temperature, na.rm = TRUE) ~ "low",
                                 Temperature <= quantile(Temperature, na.rm = TRUE)[4] ~ "high",
                                 Temperature > quantile(Temperature, na.rm = TRUE)[4] ~ "very high")) %>% 
  select(Sample, Station, PelagicLayer, PO4, NO2, NH4, Si, Chl, Temperature, Salinity, WaterMass) %>% 
  mutate(across(!where(is.factor), as.factor))
##

saveRDS(env_cat, "environmental_data_category.rds")


