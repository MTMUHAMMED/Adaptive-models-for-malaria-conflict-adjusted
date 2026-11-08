#------------------------------------------------------------#
#----------NMIS 2021 and Conflict intensity data models--------#
library(tidyverse)
library(SUMMER)
library(sf)
library(spdep)
library(INLA)
library(pacman)
library(patchwork)

# Data
load("NMIS2021_Conflict_data.RData")
#View(NMIS2021_Conflict_data)

#--------------------------------------------------------------#
# CONVERTING THE VARIABLES TO FACTORS
data <- NMIS2021_Conflict_data |> 
  mutate(
    residence = factor(residence, levels =c(1,2),labels=c("Urban","Rural")),
    MosquitoNet = factor(MosquitoNet, levels=c(0,1), labels=c("No","Yes")),
    WealthIndex = factor(WealthIndex, levels =1:5, labels =c("Poorest","Poorer","Middle","Richer","Richest")),
    MotherEduLevel=factor(MotherEduLevel,levels = 0:3, labels = c("No Education","Primary", "Secondary", "Higher")),
    Gender = factor(Gender, levels = c(1,2),labels = c("Male","Female")),
    any_conflict = factor(any_conflict, levels=c(0,1), labels=c("No","Yes")),
    conflict_intensity = factor(conflict_intensity,
                                levels=c("No conflict","Low conflict","Medium conflict","High conflict"),
                                ordered=TRUE)
  )

#View(data)

#--------------------------------------------------------------#
################# DESCRIPTIVE STATISTICS ##################
#   SUMMARY OF OBSERVATION DATA  ######################################
library(gtsummary)
data_descrip <- data |> select(-c(REGNAME, REGCODE,cluster_id, year))

data_pcent <- data_descrip |> 
  tbl_summary(
    by = conflict_intensity,                                    
    statistic = list(all_continuous() ~ "{mean}±{sd}",   # mean ± SD for numeric
                     all_categorical() ~ "{n}({p}%)"), # counts and percentages
    digits = all_continuous() ~ 2
  ) |> 
  add_n() |>              # adds total N
  add_overall() |>       # adds overall column
  modify_header(label = "Characteristics") |> 
  bold_labels()

as_gt(data_pcent) |>  gt::gtsave("data_description.docx")

###################################################################
# Neighborhood matrix (from Nigeria shapefile)
map <- st_read("Nigeria 2018 DHS sdr shapefiles/sdr_subnational_boundaries2.shp")

# From SUMMER package
amat <- getAmat(map, map$REGCODE) # amat == adjacency matrix

# From SPDEP
nb <- poly2nb(map)
#head(nb)

nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

# ADJACENCY GRAPH PLOT
# Centroids or representative points of each polygon
map_aux <- st_set_crs(map, 4326)
coords <- st_coordinates(st_point_on_surface(map_aux))

# Convert neighbor object to line segments
adjacency_lines <- spdep::nb2lines(
  nb,
  coords = coords,
  as_sf = TRUE
)

# IMPORTANT: assign CRS to the adjacency lines
adjacency_lines <- st_as_sf(adjacency_lines)
st_crs(adjacency_lines) <- st_crs(map_aux)

# Plot map + adjacencies
adjacency_map <- ggplot() +
  geom_sf(data = map_aux, fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_sf(data = adjacency_lines, color = "red", linewidth = 0.3, alpha = 0.7) +
  labs(
    title = " ",
  ) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  ) 

adjacency_map

###################################################################
# prior
prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)),
  phi = list(
    prior = "pc",
    param = c(0.5, 2 / 3))
)

precprior <- list(prec = list(prior = "pc.prec", param = c(0.5/0.31, 0.01)))
# Models' formula
data$idx_state <- match(data$REGCODE, map$REGCODE)
data$idx_conflict <-match(data$conflict_intensity, unique(data$conflict_intensity))

#-----------------------------------------------------------#
### SETUP FOR ADAPTIVE BYM 2 (General AICAR Model) 
# General AICAR Model
# Using 2021 ACLED only
load("zonal_admin_table_2021.RData")
zonal_admin_table
num_state <- nrow(zonal_admin_table)
zone <- unique(zonal_admin_table$zone)
zone_table <- zonal_admin_table
num_zone <- length(zone)
W_list <- vector(mode = "list", length = num_zone + 1)
for(m in 1:(num_zone + 1)){
  W_list[[m]] <- amat
}
for(i in 1:num_state){
  for(j in i:num_state){
    if(amat[i, j] == 1){
      
      if(zone_table[i, "zone"] == 
         zone_table[j, "zone"]){
        W_list[[num_zone + 1]][i, j] = 0
        W_list[[num_zone + 1]][j, i] = 0
        for(m in 1:num_zone){
          if(zone_table[i, "zone"] != zone[m]){
            W_list[[m]][i, j] = 0
            W_list[[m]][j, i] = 0
          }
        }
      }
      else{
        for(m in 1:num_zone){
          W_list[[m]][i, j] = 0
          W_list[[m]][j, i] = 0
        }
      }
    }
  }
}

D_list <- vector(mode = "list", length = num_zone + 1)
R_list <- vector(mode = "list", length = num_zone + 1)
for(m in 1:(num_zone + 1)){
  D_list[[m]] <- diag(rowSums(W_list[[m]]))
  R_list[[m]] <- D_list[[m]] - W_list[[m]]
}

R_sum <- R_list[[num_zone + 1]]
for(m in 1:num_zone){
  R_sum <- R_sum + R_list[[m]]
}
scaled_Q <- INLA:::inla.scale.model.bym.internal(R_sum, 
                                                 adjust.for.con.comp = TRUE)$Q
gv <- scaled_Q[1, 1] / (R_sum)[1, 1]

R_star_list <- vector(mode = "list", length = num_zone + 1)
for(m in 1:(num_zone + 1)){
  R_star_list[[m]] <- gv * (D_list[[m]] - W_list[[m]])
}

vals <- (1:num_state)[-num_state]
R_1_star <- gv * (R_sum - R_list[[num_zone + 1]])
R_2_star <- R_star_list[[num_zone + 1]]
R_1_star_hat <- R_1_star[vals, vals]
R_2_star_hat <- R_2_star[vals, vals]
eps_multi <- eigen(solve(R_1_star_hat + R_2_star_hat) %*% R_2_star_hat)$values
gamma_tilde_multi <- c(1 / eigen(R_1_star + R_2_star)$values[1:(num_state - 1)], 0)
save(R_star_list, eps_multi, gamma_tilde_multi,
     file = "Gen_AICAR_2021_structure_matrices.RData")

load("Gen_AICAR_2021_structure_matrices.RData")

# rgeneric code for General AICAR model
source("multizone_rgeneric.R")
num_state <- 37
num_region <- num_state
num_zone <- 3
pc.u <- 0.5 /0.31
pc.alpha <- 0.01

pc.u.phi <- 0.5
pc.alpha.phi <- 2 / 3

pc.u.theta <- 0.75
pc.alpha.theta <- 0.75

# Defining the General AICAR model
General_AICAR_model <- INLA::inla.rgeneric.define(
  model = inla.rgeneric.multizone.model, 
  n = num_state,
  M = num_zone,
  R_star_list = R_star_list,
  gamma_tilde = gamma_tilde_multi, 
  eps = eps_multi, U_prec = pc.u, alpha_prec = pc.alpha, 
  U_phi = pc.u.phi, alpha_phi = pc.alpha.phi,
  U_theta = pc.u.theta, 
  alpha_theta = pc.alpha.theta,
  psi_prior_setting = 1)



constr <- list(A = matrix(c(rep(0, num_region),rep(1, num_region)),nrow = 1,
                          ncol = 2 * num_region),e = 0)

options <- list(dic = TRUE,mlik = TRUE,cpo = TRUE,openmp.strategy = "default")

control.inla <- list(strategy = "adaptive",int.strategy = "auto")


################### MODELS ######################################
#################################################################
#### BYM2
formula1 <- malaria ~ residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel + 
  f(idx_state, model = "bym2", graph = g, hyper = prior)

formula2 <- malaria ~ any_conflict + residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel + f(idx_state, model = "bym2", graph = g, hyper = prior)

formula3 <- malaria ~ conflict_intensity + residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel + 
  f(idx_state, model = "bym2", graph = g, hyper = prior)

formula4 <- malaria ~ total_fatalities + residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel + f(idx_state, model = "bym2", graph = g, hyper = prior)

formula5 <- malaria ~ residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel + 
  f(idx_state, model = "bym2", graph = g, hyper = prior) +
  f(idx_conflict,model = "iid", hyper = precprior)

formula6 <- malaria ~ residence + MosquitoNet + WealthIndex + 
  NoOfHouseHoldMember + NoOfChildren + noOfRooms + ageOfHouseholdHead +
  NoOfChildSplepUnderMN + Gender + Age + MotherEduLevel +
  f(idx_state,model = General_AICAR_model,diagonal = 1e-06,extraconstr = constr,n = 2 * num_region)

###################################################################
# Models fitting
ccomp <- list(dic = TRUE, waic = TRUE, cpo = TRUE)

# Model 1: No conflict variable + spatial effect
model1 <- inla(formula1, family = "binomial", data = data,
               control.compute = ccomp, verbose = FALSE)

# Model 2: conflict:yes/no + covariates + spatial effect
model2 <- inla(formula2, family = "binomial", data = data,
               control.compute = ccomp, verbose = FALSE)

# Model 3: conflict intensity + covariates + spatial effect
model3 <- inla(formula3, family = "binomial", data = data,
               control.compute = ccomp, verbose = FALSE)

# Model 4: total_fatalities + covariates + spatial effect
model4 <- inla(formula4, family = "binomial", data = data,
               control.compute = ccomp, verbose = FALSE)

# Model 5: covariate + spatial effect + conflict effect
model5 <- inla(formula5, family = "binomial", data = data,
               control.compute = ccomp, verbose = FALSE)

model6 <- inla(formula6, family = "binomial", data = data, 
               control.compute = ccomp,verbose = FALSE)


#------------------------------------------------------------------------#
# Comparing models
#############################################
# Logarithmic Score
logscore <- data.frame(
  Model = c("Model 1","Model 2","Model 3","Model 4","Model 5",
            "Model 6"),
  LogScore = c(
    -mean(log(model1$cpo$cpo), na.rm=TRUE),
    -mean(log(model2$cpo$cpo), na.rm=TRUE),
    -mean(log(model3$cpo$cpo), na.rm=TRUE),
    -mean(log(model4$cpo$cpo), na.rm=TRUE),
    -mean(log(model5$cpo$cpo), na.rm=TRUE),
    -mean(log(model6$cpo$cpo), na.rm=TRUE)
  )
)

logscore

#DIC
DIC_values <- data.frame(
  Model = c("Model 1","Model 2","Model 3","Model 4","Model 5",
            "Model 6"),
  DIC = c(
    model1$dic$dic,
    model2$dic$dic,
    model3$dic$dic,
    model4$dic$dic,
    model5$dic$dic,
    model6$dic$dic
    
  )
)
DIC_values

# WAIC
WAIC_values <- data.frame(
  Model = c("Model 1","Model 2","Model 3","Model 4","Model 5",
            "Model 6"),
  WAIC = c(
    model1$waic$waic,
    model2$waic$waic,
    model3$waic$waic,
    model4$waic$waic,
    model5$waic$waic,
    model6$waic$waic
    
  )
)
WAIC_values
#############################################################
#-------------------------------------------------------------#
# GENERAL AICAR MODEL PARAMETERS
# Intercept (mean)
fit <- inla.tmarginal(function(x) expit(x),
                      model6$marginals.fixed[[1]])

inla.zmarginal(fit)

# Adaptive
tau_post <- inla.tmarginal(fun = exp, 
                           marginal = model6$marginals.hyperpar$`Theta1 for idx_state`)
tau_sum <- inla.zmarginal(tau_post)

phi_post <- inla.tmarginal(fun = SUMMER::expit, 
                           marginal = model6$marginals.hyperpar$`Theta2 for idx_state`)
phi_sum <- inla.zmarginal(phi_post)

theta_post <- inla.tmarginal(fun = SUMMER::expit, 
                             marginal = model6$marginals.hyperpar$`Theta3 for idx_state`)
theta_sum <- inla.zmarginal(theta_post)

theta_4 <- inla.tmarginal(fun = SUMMER::expit, 
                          marginal = model6$marginals.hyperpar$`Theta4 for idx_state`)
theta_sum <- inla.zmarginal(theta_4)

theta_5 <- inla.tmarginal(fun = SUMMER::expit, 
                          marginal = model6$marginals.hyperpar$`Theta5 for idx_state`)
theta_sum <- inla.zmarginal(theta_5)

#--------------------------------------------------------------------#
#########################################################
# for model 4 (Best model)
fixed_effects <- model4$summary.fixed

OR_table <- data.frame(
  Variable = rownames(fixed_effects),
  OR = exp(fixed_effects$mean),
  Lower_CI = exp(fixed_effects$`0.025quant`),
  Upper_CI = exp(fixed_effects$`0.975quant`)
)

OR_table


######################################
# OR graph
OR_plot <- OR_table |> 
  filter(Variable != "(Intercept)") |> 
  mutate(
    Variable = recode(
      Variable,
      "total_fatalities" = "Conflict fatalities",
      "residenceRural" = "Rural residence",
      "MosquitoNetYes" = "Mosquito net use",
      "WealthIndexPoorer" = "Poorer wealth",
      "WealthIndexMiddle" = "Middle wealth",
      "WealthIndexRicher" = "Richer wealth",
      "WealthIndexRichest" = "Richest wealth",
      "NoOfHouseHoldMember" = "Household members",
      "NoOfChildren" = "Number of children",
      "noOfRooms" = "Number of rooms",
      "ageOfHouseholdHead" = "Age of household head",
      "NoOfChildSplepUnderMN" = "Children sleeping under net",
      "GenderFemale" = "Female child",
      "Age" = "Child age",
      "MotherEduLevelPrimary" = "Primary education",
      "MotherEduLevelSecondary" = "Secondary education",
      "MotherEduLevelHigher" = "Higher education"
    )
  )

# plot
ggplot(OR_plot, 
       aes(x = OR, 
           y = reorder(Variable, OR))) +
  
  geom_point(size = 3) +
  
  geom_errorbarh(
    aes(xmin = Lower_CI,
        xmax = Upper_CI),
    height = 0.25
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  geom_text(
    aes(label = paste0(
      round(OR,2),
      " (",
      round(Lower_CI,2),
      "-",
      round(Upper_CI,2),
      ")"
    )),
    hjust = -0.1,
    size = 3
  ) +
  
  scale_x_log10() +
  
  labs(
    x = "Odds Ratio (95% CrI)",
    y = NULL
  ) +
  
  theme_classic()

#########################################
forest_data <- OR_table |> 
  filter(Variable != "(Intercept)") |> 
  mutate(
    Group = case_when(
      Variable == "total_fatalities" ~ "Conflict exposure",
      Variable %in% c("WealthIndexPoorer",
                      "WealthIndexMiddle",
                      "WealthIndexRicher",
                      "WealthIndexRichest") ~ "Wealth",
      Variable %in% c("MotherEduLevelPrimary",
                      "MotherEduLevelSecondary",
                      "MotherEduLevelHigher") ~ "Education",
      Variable %in% c("residenceRural") ~ "Location",
      Variable %in% c("GenderFemale",
                      "Age") ~ "Child characteristics",
      TRUE ~ "Household characteristics"
    )
  )

# Renaming
forest_data <- forest_data |> 
  mutate(
    Variable = recode(
      Variable,
      "total_fatalities"="Conflict fatalities",
      "residenceRural"="Rural residence",
      "MosquitoNetYes"="Mosquito net use",
      "WealthIndexPoorer"="Poorer",
      "WealthIndexMiddle"="Middle",
      "WealthIndexRicher"="Richer",
      "WealthIndexRichest"="Richest",
      "NoOfHouseHoldMember"="Household members",
      "NoOfChildren"="Number of children",
      "noOfRooms"="Number of rooms",
      "ageOfHouseholdHead"="Age of household head",
      "NoOfChildSplepUnderMN"="Children sleeping under net",
      "GenderFemale"="Female",
      "Age"="Child age",
      "MotherEduLevelPrimary"="Primary",
      "MotherEduLevelSecondary"="Secondary",
      "MotherEduLevelHigher"="Higher"
    ))

# create forest plot
ggplot(forest_data,
       aes(
         x = OR,
         y = reorder(Variable, OR),
         xmin = Lower_CI,
         xmax = Upper_CI,
         colour = Group
       )) +
  
  geom_errorbarh(
    height = 0.25,
    linewidth = 0.8
  )+
  
  geom_point(
    size = 3
  )+
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey40"
  )+
  
  scale_x_log10(
    breaks=c(0.1,0.25,0.5,1,2)
  )+
  
  facet_grid(
    Group ~ .,
    scales="free_y",
    space="free_y"
  )+
  
  labs(
    x="Odds Ratio (95% Credible Interval)",
    y=NULL,
    colour=NULL
  )+
  
  theme_classic() +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face="bold",
      size=11
    ),
    axis.text.y = element_text(
      size=10
    ),
    legend.position="none"
  )

#-------------------------------------------------------------------#
######################################################
## Mapping
dim(data)
data <- data |> mutate(
  prevalence1 = model1$summary.fitted.values[,"mean"],
  ci_width1  = model1$summary.fitted.values[, "0.975quant"] - 
    model1$summary.fitted.values[, "0.025quant"],
  
  prevalence2 = model2$summary.fitted.values[,"mean"],
  ci_width2  = model2$summary.fitted.values[, "0.975quant"] - 
    model2$summary.fitted.values[, "0.025quant"],
  
  prevalence3 = model3$summary.fitted.values[,"mean"],
  ci_width3  = model3$summary.fitted.values[, "0.975quant"] - 
    model3$summary.fitted.values[, "0.025quant"],
  
  prevalence4 = model4$summary.fitted.values[,"mean"],
  ci_width4  = model4$summary.fitted.values[, "0.975quant"] - 
    model4$summary.fitted.values[, "0.025quant"],
  
  prevalence5 = model5$summary.fitted.values[,"mean"],
  ci_width5  = model5$summary.fitted.values[, "0.975quant"] - 
    model5$summary.fitted.values[, "0.025quant"],
  
  prevalence6 = model6$summary.fitted.values[,"mean"],
  ci_width6  = model6$summary.fitted.values[, "0.975quant"] - 
    model6$summary.fitted.values[, "0.025quant"],
)

#View(data)

prev_data <- data |> group_by(REGCODE) |> 
  summarise(
    prev1 = mean(prevalence1),
    ci_width_1 = mean(ci_width1),
    
    prev2 = mean(prevalence2),
    ci_width_2 = mean(ci_width2),
    
    prev3 = mean(prevalence3),
    ci_width_3 = mean(ci_width3),
    
    prev4 = mean(prevalence4),
    ci_width_4 = mean(ci_width4),
    
    prev5 = mean(prevalence5),
    ci_width_5 = mean(ci_width5),
    
    prev6 = mean(prevalence6),
    ci_width_6 = mean(ci_width6)
  )
#View(prev_data)
#View(map)
Prev_map <- map |> left_join(prev_data, by="REGCODE")
#View(Prev_map)
map1 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev1)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 1'
  )
#map1

map2 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev2)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 2'
  )
#map2

map3 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev3)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 3'
  )
#map3

map4 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev4)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 4'
  )
#map4

map5 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev5)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 5'
  )
#map5

map6 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = prev6)) + 
  scale_fill_distiller(palette = "RdYlBu", name = 'Prevalence') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 6'
  )
#map6


combined_maps <- (map1 | map2 | map3) /(map4 | map5 | map6)

prevalence_maps <- combined_maps +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

prevalence_maps
####################################################################
####################################################################
# 95% Uncertainty maps
width1 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_1)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 1'
  )
width1

width2 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_2)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 2'
  )
width2

width3 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_3)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 3'
  )
width3

width4 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_4)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 4'
  )
width4

width5 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_5)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 5'
  )
width5

width6 <-  ggplot(data = Prev_map) +
  geom_sf(aes(fill = ci_width_6)) + 
  scale_fill_distiller(palette = "RdYlBu", name = '95% CI width') +
  labs(title = NULL) + 
  theme_minimal() + 
  theme(legend.position = "right")+
  labs(
    title = 'Model 6'
  )
width6


uncertainty_maps <- (width1 | width2 | width3) /(width4 | width5 | width6)

Width_maps <- uncertainty_maps +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

Width_maps

#######################################################
# Random effects
Prev_map$spatial4 <-model4$summary.random$idx_state$mean[1:37]
Prev_map$spatial6 <-model6$summary.random$idx_state$mean[1:37]

# Find common range across both models
common_range <- range(c(Prev_map$spatial4, Prev_map$spatial6), na.rm = TRUE)
common_range

# Maps
spatial_map4 <- ggplot(data = Prev_map) +
  geom_sf(aes(fill = spatial4)) + 
  scale_fill_distiller(
    palette = "RdYlBu",
    name = "Random Effect",
    limits = common_range
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Spatial Random Effect for model 4")


spatial_map6 <- ggplot(data = Prev_map) +
  geom_sf(aes(fill = spatial6)) + 
  scale_fill_distiller(
    palette = "RdYlBu",
    name = "Random Effect",
    limits = common_range
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Spatial Random Effect for model 6")


spatial_map4
spatial_map6

spatial_map4 | spatial_map6

