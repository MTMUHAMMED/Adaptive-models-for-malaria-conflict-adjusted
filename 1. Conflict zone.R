#################################################################
######## ACLED DATA A YEAR BEFORE NMIS 2021 (12TH OCT 2021)
library(tidyverse)
library(readr)
library(lubridate)
library(sf)
library(spdep)
ACLED_Data <- read_delim("ACLED_2000_2024.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
#View(ACLED_Data)
dim(ACLED_Data)
unique(ACLED_Data$admin1)

acled_nmis <- ACLED_Data |>  filter(country=='Nigeria' & year==c(2020,2021)) |> 
  mutate(admin1 = recode(admin1, 
                         "Federal Capital Territory" = "FCT Abuja",
                         "Nassarawa" = "Nasarawa"))

unique(acled_nmis$admin1)
#View(acled_nmis)
# NMIS 2021 survey date
survey_date <- as.Date("2021-10-12")

# Start date: exactly one year before the survey
start_date <- survey_date %m-% years(1)
# Select conflict events occurring within the one-year period
ACLED_2021 <- acled_nmis |> mutate(
  event_date = as.Date(event_date)
) |> 
  filter(
    event_date >= start_date,
    event_date <= survey_date
  ) |> filter(event_type %in% c("Violence against civilians","Battles","Explosions/Remote violence"))
#View(ACLED_2021)
# Overall number of conflict events and total fatalities
ACLED_2021_Summary <- ACLED_2021 |> group_by(admin1) |> 
  summarise(
    conflict_count = n(),
    total_fatalities = sum(fatalities, na.rm = TRUE),
    log_fatality = log(total_fatalities)
  )

#View(ACLED_2021_Summary)

# Dividing into groups
q <- quantile(
  ACLED_2021_Summary$total_fatalities,
  #ACLED_2021_Summary$log_fatality,
  probs = c(1/3, 2/3),
  na.rm = TRUE
)

# Classify fatalities into three groups for states
conflict_2021 <- ACLED_2021_Summary |>  mutate(
  fatality_group = case_when(
    total_fatalities <= q[1] ~ "Low",
    total_fatalities <= q[2] ~ "Medium",
    total_fatalities >  q[2] ~ "High"
  ),
  fatality_group = factor(
    fatality_group,
    levels = c("Low", "Medium", "High"),
    ordered = TRUE
  )
) |> rename(zone = fatality_group, REGNAME = admin1)

table(conflict_2021$zone)

# Registration code for states in Nigeria
reg_data <- data.frame(
  REGCODE = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370),
  REGNAME = c('Sokoto', 'Zamfara', 'Katsina', 'Jigawa', 'Yobe', 'Borno', 'Adamawa', 'Gombe', 'Bauchi', 'Kano', 'Kaduna', 'Kebbi', 'Niger', 'FCT Abuja', 'Nasarawa', 'Plateau', 'Taraba', 'Benue', 'Kogi', 'Kwara', 'Oyo', 'Osun', 'Ekiti', 'Ondo', 'Edo', 'Anambra', 'Enugu', 'Ebonyi', 'Cross River', 'Akwa Ibom', 'Abia', 'Imo', 'Rivers', 'Bayelsa', 'Delta', 'Lagos', 'Ogun')
)
ACLED_ZONE_2021 <- merge(conflict_2021, reg_data,by = "REGNAME", all.x = TRUE)
ACLED_ZONE_2021 <- ACLED_ZONE_2021 |> arrange(REGCODE) |> select(
  REGNAME,REGCODE, zone
)
#View(ACLED_ZONE_2021)

# Load Nigeria shapefile
map <- st_read("Nigeria 2018 DHS sdr shapefiles/sdr_subnational_boundaries2.shp")
#View(map)

Zone_map <- left_join(map, ACLED_ZONE_2021, by="REGCODE")
#View(Zone_map)

#map
######################################################
# ADJACENCY GRAPH PLOT
# Centroids or representative points of each polygon
map_aux <- st_set_crs(Zone_map, 4326)

coords <- st_coordinates(st_point_on_surface(map_aux))

# Convert neighbor object to line segments
nb <- poly2nb(map)
adjacency_lines <- spdep::nb2lines(
  nb,
  coords = coords,
  as_sf = TRUE
)

# IMPORTANT: assign CRS to the adjacency lines
adjacency_lines <- st_as_sf(adjacency_lines)
st_crs(adjacency_lines) <- st_crs(map_aux)


# Plot map + conflict zones + adjacencies
adjacency_map <- ggplot() +
  
  # Shade conflict zones
  geom_sf(
    data = map_aux,
    aes(fill = zone),
    color = "black",
    linewidth = 0.4
  ) +
  
  # Add adjacency network
  geom_sf(
    data = adjacency_lines,
    color = "black",
    linewidth = 0.5,
    alpha = 0.6
  ) +
  
  # Conflict zone colours
  scale_fill_manual(
    values = c(
      "Low" = "#01ba38",
      "Medium" = "#619dff",
      "High" = "#f9766e"
    ),
    name = "Conflict Zone"
  ) +
  
  labs(
    title = "Conflict Zones and Neighbourhood Structure"
  ) 

adjacency_map


################################################################
# Plot

Zones <- ggplot(conflict_2021, 
                aes(x = log(total_fatalities + 1),
                    y = REGNAME,
                    color = as.factor(zone))) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c(
      "#00BA38",   # Green
      "#619CFF",  # Blue
      "#F8766D"  # Salmon red
      
    )
  ) +
  labs(
    title = "Tertile classification of armed conflict fatality count",
    x = "Log(Total fatality + 1)",
    y = "State",
    color = "Conflict zone"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

Zones

