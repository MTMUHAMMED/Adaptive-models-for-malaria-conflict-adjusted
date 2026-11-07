# ============================================================
# ACLED CONFLICT EXPOSURE FOR THE 2021 NIGERIA MIS
# Exposure period: 12 October 2020 to 11 October 2021
# Spatial exposure: conflict events within 10 km of each cluster
# ============================================================
library(sf)
library(dplyr)
library(readr)
library(data.table)
library(lubridate)
library(janitor)
library(purrr)
library(tidyr)
library(units)

# Use spherical geometry for longitude-latitude distance calculations
sf_use_s2(TRUE)
# ------------------------------------------------------------
# Define the survey and conflict exposure dates
# ------------------------------------------------------------
survey_start <- ymd("2021-10-12")
survey_end   <- ymd("2021-12-04")
# One complete year before NMIS fieldwork began
conflict_start <- survey_start %m-% years(1)
conflict_end   <- survey_start - days(1)

conflict_start
conflict_end
# ------------------------------------------------------------
# Read the downloaded ACLED Nigeria dataset
# ------------------------------------------------------------
acled_raw <- read_delim("ACLED_2000_2024.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(acled_raw)
acled <- acled_raw |>  clean_names()
# Examine the available ACLED variables
names(acled)
# ------------------------------------------------------------
# Prepare the ACLED event date and variables
acled <- acled |>  mutate(
  event_date = as.Date(
    parse_date_time(
      event_date,orders = c("ymd", "dmy", "mdy")
    )
  ),
  latitude  = as.numeric(latitude),
  longitude = as.numeric(longitude),
  fatalities = replace_na(
    as.numeric(fatalities),
    0
  )
)

# ------------------------------------------------------------
# Extract violent conflict events occurring in Nigeria
violent_event_types <- c(
  "Battles",
  "Explosions/Remote violence",
  "Violence against civilians"
)

acled_one_year <- acled %>%
  filter(
    country == "Nigeria",
    event_date >= conflict_start,
    event_date <= conflict_end,
    event_type %in% violent_event_types,
    !is.na(latitude),
    !is.na(longitude)
  )

# Check the extracted date range
range(acled_one_year$event_date, na.rm = TRUE)

# Check event frequencies
acled_one_year |> 
  count(event_type, sort = TRUE)

# Number of extracted conflict events
nrow(acled_one_year)

# Total reported fatalities
sum(acled_one_year$fatalities, na.rm = TRUE)


# ------------------------------------------------------------
# Read the 2021 NMIS GPS cluster shapefile
# ------------------------------------------------------------

nmis_gps_raw <- st_read("NGGE81FL_NMIS_2021/NGGE81FL.shp",quiet = TRUE) |> 
  clean_names()

# Examine the variable names
names(nmis_gps_raw)

nmis_clusters <- nmis_gps_raw |> 
  transmute(
    cluster_id = as.integer(dhsclust),
    
    # Keep urban-rural classification when available
    residence = if ("urban_rura" %in% names(nmis_gps_raw)) {
      as.character(urban_rura)
    } else {
      NA_character_
    },
    
    geometry = geometry
  ) %>%
  filter(!st_is_empty(geometry))

# Confirm the number of clusters
nrow(nmis_clusters)

# Number of distinct clusters
n_distinct(nmis_clusters$cluster_id)


# ------------------------------------------------------------
# Convert the ACLED events to spatial points
acled_sf <- st_as_sf(
  acled_one_year,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# Ensure both datasets use the same coordinate reference system
nmis_clusters <- st_transform(nmis_clusters, 4326)
acled_sf      <- st_transform(acled_sf, 4326)


# ------------------------------------------------------------
# Function for aggregating conflict within a chosen buffer
extract_cluster_conflict <- function(
    clusters,
    conflict_events,
    buffer_km = 10
) {
  
  # Identify all conflict events within the specified buffer
  nearby_events <- st_is_within_distance(
    clusters,
    conflict_events,
    dist = units::set_units(buffer_km, "km")
  )
  
  conflict_summary <- map_dfr(
    seq_len(nrow(clusters)),
    function(i) {
      
      cluster_number <- clusters$cluster_id[i]
      event_indices  <- nearby_events[[i]]
      
      # Return zeros where no conflict event was recorded
      if (length(event_indices) == 0) {
        
        return(
          tibble(
            cluster_id = cluster_number,
            buffer_km = buffer_km,
            conflict_events = 0L,
            total_fatalities = 0,
            battles = 0L,
            remote_violence = 0L,
            violence_against_civilians = 0L,
            civilian_targeting_events = 0L,
            days_since_last_event = NA_real_
          )
        )
      }
      
      cluster_events <- conflict_events[
        event_indices,
        ,
        drop = FALSE
      ] %>%
        st_drop_geometry()
      
      # Check whether civilian_targeting exists
      civilian_events <- if (
        "civilian_targeting" %in% names(cluster_events)
      ) {
        sum(
          !is.na(cluster_events$civilian_targeting) &
            cluster_events$civilian_targeting != "",
          na.rm = TRUE
        )
      } else {
        sum(
          cluster_events$event_type ==
            "Violence against civilians",
          na.rm = TRUE
        )
      }
      
      tibble(
        cluster_id = cluster_number,
        buffer_km = buffer_km,
        
        # Total number of distinct conflict events
        conflict_events = n_distinct(
          cluster_events$event_id_cnty
        ),
        
        # Total reported fatalities
        total_fatalities = sum(
          cluster_events$fatalities,
          na.rm = TRUE
        ),
        
        # Conflict events by event type
        battles = sum(
          cluster_events$event_type == "Battles",
          na.rm = TRUE
        ),
        
        remote_violence = sum(
          cluster_events$event_type ==
            "Explosions/Remote violence",
          na.rm = TRUE
        ),
        
        violence_against_civilians = sum(
          cluster_events$event_type ==
            "Violence against civilians",
          na.rm = TRUE
        ),
        
        civilian_targeting_events = civilian_events,
        
        # Recency of the latest event before fieldwork
        days_since_last_event = as.numeric(
          conflict_end -
            max(cluster_events$event_date, na.rm = TRUE)
        )
      )
    }
  )
  
  conflict_summary %>%
    mutate(
      # Binary conflict shock
      any_conflict = as.integer(conflict_events > 0),
      
      # Log transformations for highly skewed variables
      log_conflict_events = log1p(conflict_events),
      log_fatalities = log1p(total_fatalities),
      
      # Ordered conflict intensity category
      conflict_intensity = case_when(
        conflict_events == 0 ~ "No conflict",
        conflict_events == 1 ~ "One event",
        conflict_events <= 3 ~ "Two to three events",
        conflict_events >= 4 ~ "Four or more events"
      ),
      
      conflict_intensity = factor(
        conflict_intensity,
        levels = c(
          "No conflict",
          "One event",
          "Two to three events",
          "Four or more events"
        ),
        ordered = TRUE
      )
    )
}


# ------------------------------------------------------------
# Extract conflict exposure within 10 km
# ------------------------------------------------------------

cluster_conflict_10km <- extract_cluster_conflict(
  clusters = nmis_clusters,
  conflict_events = acled_sf,
  buffer_km = 10
)

# Inspect the results
head(cluster_conflict_10km)
#View(cluster_conflict_10km)
summary(cluster_conflict_10km$conflict_events)
summary(cluster_conflict_10km$total_fatalities)

# Number of clusters exposed to at least one conflict event
table(cluster_conflict_10km$any_conflict)

# Percentage of exposed clusters
mean(cluster_conflict_10km$any_conflict) * 100

#################################################################
# Read in malaria data
load("Malaria_2010_2021_covariates.RData")
nmis_unit <- cleaned_data |> filter(year==2021)
View(nmis_unit)
length(unique(nmis_unit$cluster))
# Ensure the cluster identification variables have the same type
nmis_unit <- nmis_unit |> 
  mutate(
    cluster_id = as.integer(cluster)
  )

cluster_conflict_10km <- cluster_conflict_10km |> 
  mutate(
    cluster_id = as.integer(cluster_id)
  )

# Merge ACLED conflict exposure with each sampled child
nmis_acled_model_data <- nmis_unit |> 
  left_join(
    cluster_conflict_10km,
    by = "cluster_id"
  )
#View(nmis_acled_model_data)

#------------------------------------------------------------#
# Categorizing total fatalities into groups based on quantile
positive_fatality <- nmis_acled_model_data |> 
  filter(total_fatalities > 0)


q <- quantile(
  positive_fatality$total_fatalities,
  probs=c(0.33,0.66),
  na.rm=TRUE
)

# Grouping conflict based on itensity
conflict_data <- nmis_acled_model_data |> 
  mutate(
    conflict_intensity = case_when(
      
      total_fatalities == 0 ~ "No conflict",
      
      total_fatalities <= q[1] ~ "Low conflict",
      
      total_fatalities <= q[2] ~ "Medium conflict",
      
      total_fatalities > q[2] ~ "High conflict"
      
    )
  )
#View(conflict_data)
table(conflict_data$conflict_intensity)

NMIS2021_Conflict_data <- conflict_data |> select(
  REGNAME, REGCODE, cluster_id,year, malaria,NoOfHouseHoldMember,
  NoOfChildren,residence,noOfRooms,ageOfHouseholdHead,MosquitoNet,
  NoOfChildSplepUnderMN,WealthIndex,Gender,Age,MotherEduLevel,any_conflict,
  total_fatalities, conflict_intensity
)
save(NMIS2021_Conflict_data, file = "NMIS2021_Conflict_data.RData")

