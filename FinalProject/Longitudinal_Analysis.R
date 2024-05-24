setwd("~/GISC 38400")

###### LOAD PACKAGES ###########

library(labelled)   # labeling data
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(GGally)     # advanced plot
library(car)        # useful for anova/wald test
library(Epi)        # easy getting CI for model coef/pred
library(lme4)       # linear mixed-effects models
library(lmerTest)   # test for linear mixed-effects models
library(emmeans)    # marginal means
library(multcomp)   # CI for linear combinations of model coef
library(geepack)    # generalized estimating equations
library(ggeffects)  # marginal effects, adjusted predictions
library(gt)         # nice tables

install.packages("areal")

library(tidyverse)  # for everything (data manipulation, visualization, coding, and more)
theme_set(theme_minimal() + theme(legend.position = "bottom")) # theme for ggplot

library(sf)
library(dplyr)
library(tmap)
library(units)
library(areal)

###### LOAD DATA ###########

fp_2024 <- st_read("~/GISC 38400/Footprints.csv")
blocks_2020 <- st_read("~/GISC 38400/Census_Blocks_2020/tl_2020_17_tabblock20.shp")
urp_blocks_1970 <- st_read("~/GISC 38400/Blocks_URP_1970_fixed.shp") %>%
  st_transform(crs = 3435)
urp_blocks_2020 <- st_read("~/GISC 38400/Blocks_URP_2020.shp") %>%
  st_transform(crs = 3435)
urp_areas <- st_read("~/GISC 38400/URP_Area.shp")
urp_boundary <- st_read("~/GISC 38400/URP_Boundary.shp")

###### WRANGLE DATA ###########

# calculate % of URP boundary that is URP project
st_area(urp_area_dissolve)/st_area(urp_boundary)

# add 1970 census data
census_data_1970 <- read_csv("~/GISC 38400/NHGIS Data/1970_Block_Revised.csv")

# add 2020 census data and make it workable
census_data_2020 <- read_csv("~/GISC 38400/NHGIS Data/2020_Block.csv")
census_data_2020 <- census_data_2020 %>%
  rename("GEOID20" = GEOCODE) %>%
  mutate(GEOID20 = as.character(GEOID20))

# join census data to 1970 blocks in URP
urp_blocks_1970_join <- left_join(urp_blocks_1970, census_data_1970, by = "GISJOIN")
urp_blocks_1970_join$AREA <- drop_units(st_area(urp_blocks_1970_join)) 
urp_blocks_1970_join <- urp_blocks_1970_join %>%
  rename("GISJOIN70" = GISJOIN)

# join census data to 2020 blocks in URP
urp_blocks_2020_join <- left_join(urp_blocks_2020, census_data_2020, by = "GEOID20")
urp_blocks_2020_join$AREA <- drop_units(st_area(urp_blocks_2020_join))
urp_blocks_2020_join <- urp_blocks_2020_join %>%
  rename("GISJOIN20" = GISJOIN)

###### INTERPOLATION: AREAL PACKAGE ######
aw_interpolate(urp_blocks_1970_join %>% dplyr::select(c("GISJOIN70", "CK0001")), GISJOIN70, urp_blocks_2020_join %>% dplyr::select("GISJOIN20", "U9V001") %>% rename("CK0001" = U9V001), GISJOIN20, weight = "sum", output = "sf", extensive = "CK0001")
ar_validate(urp_blocks_1970_join, urp_blocks_2020_join, varList = c("CK0001"), verbose = TRUE)
View(urp_blocks_2020_join)
###### INTERPOLATION ###########

# Ensure the GeoDataFrames have a common CRS (Coordinate Reference System)
st_crs(urp_blocks_2020_join) == st_crs(urp_blocks_1970_join)

# Rename GISJOIN columns for clarity
urp_blocks_1970_join <- urp_blocks_1970_join %>%
  rename(GISJOIN_1970 = GISJOIN)
urp_blocks_2020_join <- urp_blocks_2020_join %>%
  rename(GISJOIN_2020 = GISJOIN)

# Overlay the 2020 data to the 1970 data using intersection
interpolated_data <- st_intersection(urp_blocks_1970_join, urp_blocks_2020_join %>% dplyr::select(c(GISJOIN20, U7H001, U9V001, VAL001)))

interpolated_data$GEOID20 <- as.numeric(interpolated_data$GEOID20)

interpolated_data <- interpolated_data %>%
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"))

# Calculate weights for areal interpolation based on the intersection area
interpolated_data <- interpolated_data %>%
  mutate(area_intersection = drop_units(st_area(geometry))) %>%
  group_by(GISJOIN20) %>%
  mutate(weight = area_intersection / sum(area_intersection)) %>%
  ungroup()

View(interpolated_data)

# Perform areal interpolation for the target variable 'CA5001'
interpolated_data <- interpolated_data %>%
  mutate(pop20 = U7H001 * weight) %>%
  mutate(hou20 = U9V001 * weight) %>%
  mutate(vac20 = VAL001 * weight)

# Aggregate the interpolated values by the 1970 tract boundaries
pop20_70 <- interpolated_data %>%
  group_by(GISJOIN70) %>%
  summarise(POP20 = sum(pop20, na.rm = TRUE))

hou20_70 <- interpolated_data %>%
  group_by(GISJOIN70) %>%
  summarise(HOU20 = sum(hou20, na.rm = TRUE))

vac20_70 <- interpolated_data %>%
  group_by(GISJOIN70) %>%
  summarise(VAC20 = sum(vac20, na.rm = TRUE))

blk20_70 <- interpolated_data %>%
  group_by(GISJOIN70) %>%
  summarise(VAC20 = sum(vac20, na.rm = TRUE))

interp_70 <- pop20_70 %>% inner_join(st_drop_geometry(hou20_70), by = "GISJOIN70") %>% inner_join(st_drop_geometry(vac20_70), by = "GISJOIN70")

interp_70 <- urp_blocks_1970_join %>% left_join(st_drop_geometry(interp_70), by = 'GISJOIN70')

interp_70 <- interp_70 %>%
  mutate(across(11:47, ~ if_else(. < 0, NA_real_, .)))

# calculate statistics from data
urp_blocks_2020_join <- urp_blocks_2020_join %>%
  mutate(
    HOUS_DEN = U9V001 / drop_units(st_area(geometry)),
    VAC_PCT = VAL001 / (VAP001 + VAL001),
    OWN_M_PCT = U9Y002 / U9V001,
    OWN_F_PCT = U9Y003 / U9V001,
    RENTER_PCT = U9Y004 / U9V001,
    BLK_PCT = U7L004 / U7L001
  )

tm_shape(urp_blocks_2020_join)+tm_polygons('BLK_PCT')

interp_70$POP70 <- interp_70$CM5001+interp_70$CM5002
interp_70$HOU70 <- interp_70$CK0001
interp_70$VAC70 <- interp_70$CK1002+interp_70$CK1003+interp_70$CK1004

# Calculate the percentage change in housing units
interp_70 <- interp_70 %>%
  mutate(PCT_POP_CHG = (POP20 - POP70) / POP70)%>%
  mutate(PCT_HOU_CHG = (HOU20 - HOU70) / HOU70) %>%
  mutate(PCT_VAC_CHG = (VAC20 - VAC70) / VAC70)

###### LONGITUDINAL ANALYSIS ###########

# dissolve URP direct project area, get portions of each block in URP area
urp_area_dissolve <- st_union(urp_areas)
urp_clip_70 <- st_intersection(urp_blocks_1970_join %>% dplyr::select(GISJOIN70), urp_area_dissolve)

# calculate area of direct URP effect portions
urp_clip_70$URP_AREA <- st_drop_geometry(st_area(urp_clip_70))

# join URP portions to blocks, calculate %URP
interp_70 <- left_join(interp_70, st_drop_geometry(urp_clip_70), by = "GISJOIN70")
interp_70$PCT_URP <- drop_units(interp_70$URP_AREA/st_drop_geometry(st_area(interp_70)))
interp_70$PCT_URP[is.na(interp_70$PCT_URP)] <- 0

# map it out
tm_shape(interp_70) + tm_polygons(col = "VAC20")
tm_shape(interp_70) + tm_polygons("PCT_POP_CHG", style = "jenks")
tm_shape(interp_70) + tm_polygons("PCT_URP", style = "jenks")
plot(interp_70$VAC20/interp_70$HOU20,(interp_70$CK1002+interp_70$CK1003+interp_70$CK1004)/interp_70$CK0001)

summary(lm(PCT_URP ~ PCT_VAC_CHG, data = interp_70[,c('PCT_URP','PCT_VAC_CHG')] %>%
             st_drop_geometry() %>%
             na.omit() %>%
             filter(across(everything(), ~ is.finite(.)))))

# write files
st_write(interp_70 %>% st_transform("WGS84"), "~/GISC 38400/interp_70_wgs.geojson")
st_write(urp_areas %>% st_transform("WGS84"), "~/GISC 38400/urp_areas.geojson")
st_write(urp_boundary %>% st_transform("WGS84"), "~/GISC 38400/urp_boundary.geojson")
tm_shape(urp_boundary)+tm_borders("black")
