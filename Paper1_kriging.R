######################### load packages #####################################

library(gstat)
library(sp)
library(rgdal)
library(lattice)
library(ggplot2)
library(mapview)
library(sf)

################################################################################################
####################################### Hf #####################################################
################################################################################################

########################## load data ######################################
mydata <- read.csv("Hf_data")
df <- data.frame(mydata)
t <- df[] #U-Pb age column
EHFI <- df[] #eHf values column

########################### constant value ###########################################
m = -0.015 # average continental crust 176Lu/177Hf
tc <- seq(2600, 3000, by = 10) #each time-slice for every 10 Ma interval
########################### Temporal adjustment #####################################################
# Loop through the indices and calculate Hf for each
for (i in 1:41) {
  t_diff <- t - tc[i]
  a <- ifelse(t > tc[i], t_diff, numeric(0))
  Hf <- (m * a) + EHFI
  df[[paste0("Hf", i)]] <- Hf
}

############################## creating grid ###################################################

coordinates(df) <- c("LONG","LAT")
proj4string(df) <- CRS("+init=epsg:4326")
sp.df <- spTransform(df, CRSobj = "+proj=utm +zone=51 +south +datum=WGS84 +units=m
                     +no_defs") # use projection for coordinate system
sp_sp <- SpatialPoints(sp.df@coords, CRS("+proj=utm +zone=51 +south +datum=WGS84 +units=m
                                         +no_defs")) # Extract the coordinates from the transformed spatial data frame
grid <- SpatialPixels(SpatialPoints(makegrid(sp_sp, n=10000)),
                      proj4string = proj4string(sp_sp)) # Create a grid of spatial pixels based on the transformed spatial points


############################## extracting objects, i and locations, z for every time-slice ##################################################################

# Create an empty list to store the variables
Hf_list <- list() # creating the 'object' arguments, i for variogram()

# Loop through Hf variables
for (i in 1:41) {
  # Extract the ith Hf variable and remove NAs
  Hf_list[[i]] <- na.omit(df@data[[paste0("Hf", i)]])
}

# Assign variables to individual variables i1, i2, ..., i41 is equal to its respective time-slice 2600,2610 ...., 3000 Ma (object for variogram calculations)
for (i in 1:41) {
  assign(paste0("i", i), Hf_list[[i]])
}

##########################################################
# Create an empty list to store the variables
z_list <- list() # creating the 'locations' arguments, z for variogram()

# Define the list of indices based on the specific pattern
indices_list <- list(
  1:length(i1), 2:233, 5:233, 10:233, 25:233, 39:233, 58:233, 79:233,
  98:233, 118:233, 129:233, 140:233, 146:233, 153:233, 161:233,
  167:233, 171:233, 177:233, 179:233, 185:233, 188:233, 196:233,
  203:233, 204:233, 207:233, 207:233, 207:233, 207:233, 207:233,
  207:233, 207:233, 208:233, 210:233, 214:233, 216:233, 219:233,
  224:233, 227:233, 230:233, 232:233, 232:233
)

# Create a loop to generate z variables
for (i in 1:length(indices_list)) {
  z_list[[i]] <- sp_sp[indices_list[[i]], ]
}

# Assign variables to individual variables z1, z2, ..., z41 is equal to its respective time-slice 2600,2610 ...., 3000 Ma 
for (i in 1:length(indices_list)) {
  assign(paste0("z", i), z_list[[i]])
}

################### computing variogram #########################

vrgm <- function(i, z){
  var_plot <- variogram(i~x+y, z)
  return(var_plot)
}
plot_var <- vrgm(ix,zy) #replace x and y with the chosen timeslice (choose object i1,i2...i41 its respective locations z1,z2...z41)
plot(plot_var)

#variogram parameters
m_eHf = vgm(a, "Model", b, c) #replace a with sill, "Model" with chosen model to fit the variogram, b with range, and c with nugget value.

fit.v <- fit.variogram(plot_var, m_eHf)#fit the variogram using the sill, range and nugget parameters
plot(plot_var, fit.v)

################# kriging #################################

krig_pred_Hf = krige(ix~x+y, zy, grid, m_eHf)
mapview(krig_pred_Hf["var1.pred"]) #isotope map
mapview(krig_pred_Hf["var1.var"]) #variance map

################################################################################################
####################################### Nd #####################################################
################################################################################################

########################## load data ######################################
mydata_Nd <- read.csv("Nd_data")
df_Nd <- data.frame(mydata_Nd)
t_Nd <- df_Nd[] #U-Pb age column
ENDI <- df_Nd[] #eNd values column

########################### constant value ###########################################
m = -0.11 # average continental crust 147Sm/144Nd
tc <- seq(2600, 3000, by = 10)  #each time-slice for every 10 Ma interval

########################### Temporal adjustment #####################################################

# Loop to calculate values
for (i in seq_along(tc)) {
  t_diff_Nd <- t_Nd - tc[i]
  a_Nd <- ifelse(t_Nd > tc[i], t_diff_Nd, numeric(0))
  df_Nd[[paste0("Nd", i)]] <- (m * a_Nd) + ENDI
}

############################## creating grid ###################################################

coordinates(df_Nd) <- c("LONG","LAT")
proj4string(df_Nd) <- CRS("+proj=utm +zone=51 +south +datum=WGS84 +units=m +no_defs")

sp.df_Nd <- spTransform(df_Nd, CRSobj = "+proj=utm +zone=51 +south +datum=WGS84 +units=m
+no_defs")
sp_sp_Nd <- SpatialPoints(sp.df_Nd@coords, CRS("+proj=utm +zone=51 +south +datum=WGS84 +units=m
+no_defs"))
grid_Nd <- SpatialPixels(SpatialPoints(makegrid(sp_sp_Nd, n=10000)),
                         proj4string = proj4string(sp_sp_Nd))

############################## extracting objects, i_Nd and locations, z_Nd for every time-slice ##################################################################

# Create an empty list to store the variables
Nd_list <- list() # creating the 'object' arguments, i_Nd for variogram()

# Loop through Nd variables
for (i in 1:41) {
  # Extract the ith Nd variable and remove NAs
  Nd_list[[i]] <- na.omit(df_Nd@data[[paste0("Nd", i)]])
}

# Assign variables to individual variables i1_Nd, i2_Nd, ..., i41_Nd is equal to its respective time-slice 2600,2610 ...., 3000 Ma (object for variogram calculations)
for (i in 1:41) {
  assign(paste0("i", i, "_Nd"), Nd_list[[i]])
}

##########################################################
# Create an empty list to store the variables
z_Nd_list <- list() # creating the 'location' arguments for variogram()

# Define the list of indices based on the specific pattern
indices_list_Nd <- list(
  1:length(i1_Nd), 7:551, 23:551, 43:551, 84:551, 187:551, 232:551, 293:551,
  338:551, 373:551, 395:551, 405:551, 411:551, 421:551, 430:551,
  441:551, 457:551, 462:551, 463:551, 467:551, 470:551, 477:551,
  484:551, 484:551, 485:551, 485:551, 486:551, 486:551, 486:551,
  486:551, 486:551, 486:551, 490:551, 490:551, 495:551, 513:551,
  521:551, 524:551, 530:551
) 

# Create a loop to generate z variables
for (i in 1:length(indices_list_Nd)) {
  z_Nd_list[[i]] <- sp_sp_Nd[indices_list_Nd[[i]], ]
}

# Assign variables to individual variables z1_Nd, z2_Nd, ..., z41_Nd is equal to its respective time-slice 2600,2610 ...., 3000 Ma 
for (i in 1:length(indices_list_Nd)) {
  assign(paste0("z", i, "_Nd"), z_Nd_list[[i]])
}

################### computing variogram #########################

vrgm <- function(i, z){
  var_plot <- variogram(i~x+y, z)
  return(var_plot)
}
plot_var <- vrgm(ix_Nd,zy_Nd) #replace x and y with the chosen timeslice (choose object i1,i2...i41 its respective locations z1,z2...z41)
plot(plot_var, pl = T)

#variogram parameters
m_eNd = vgm(a, "Model", b, c) #replace a with sill, "Model" with chosen model to fit the variogram, b with range, and c with nugget value.

fit.v <- fit.variogram(plot_var, m_eNd)#fit the variogram using the sill, range and nugget parameters
plot(plot_var, fit.v)

################# kriging #################################

krig_pred_Nd = krige(i1_Nd~x+y, z1_Nd, grid, m_eNd)
mapview(krig_pred_Nd["var1.pred"]) #isotope map
mapview(krig_pred_Nd["var1.var"]) #variance map

################ cokriging ##################################

#start of co-kriging: cross variogram

g <- gstat(NULL, id = "eHf", form = ix ~ x+y, data=zy) #replace x and y with the chosen timeslice (choose object i1,i2...i41 its respective locations z1,z2...z41)
g <- gstat(g, id = "eNd", form = ix_Nd ~ x+y, data=zy_Nd) #replace x and y with the chosen timeslice (choose object i1,i2...i41 its respective locations z1,z2...z41)
v.cross <- variogram(g)
plot(v.cross)

#variograms parameters
m_fit <- vgm(a, "Model", b, c) #replace a with sill, "Model" with chosen model to fit the variogram, b with range, and c with nugget value.
m_eHf = vgm(a, "Model", b, c) #replace a with sill, "Model" with chosen model to fit the variogram, b with range, and c with nugget value.
m_eNd = vgm(a, "Model", b, c) #replace a with sill, "Model" with chosen model to fit the variogram, b with range, and c with nugget value.

g <- gstat(g, id = "eHf", model = m_eHf)
g <- gstat(g, id = "eNd", model = m_eNd)
g <- gstat(g, id = c("eHf", "eNd"), model = m_fit)
g <- fit.lmc(v.cross, g)

vcross.fit <- variogram(g)
plot(vcross.fit, model=g$model)

#co-Kriging
k.c <- predict(g, grid)
mapview(k.c["eHf.pred"]) 

################# eHf-time plot ##########################

library(dplyr)
library(tidyr)

eHf_time <- read.csv("C:/Users/23030542/Desktop/sampled_coHf_1000km.csv")
# Gather columns rvalue_1 to rvalue_39 into two columns: time and Hf_value
df_long <- eHf_time %>%
  gather(key = "time", value = "Hf_value", -distance) %>%
  mutate(time = as.numeric(sub("rvalue_", "", time)),
         time = 2600 + (time - 1) * 10) %>%
  arrange(time, distance)

# Create a new column for the repeated distances
df_long <- df_long %>%
  mutate(repeated_distance = rep(1:11, each = nrow(df_long)/11))

# Select and arrange the columns in the desired order
df_long <- df_long %>%
  select(time, distance, Hf_value)

# Create the plot with multiple distances
ggplot(df_long, aes(x = time, y = Hf_value, color = as.factor(distance))) +
  geom_line(linewidth = 1.5) +
  labs(x = "Time Slice (Ma)", y = "εHf Value", color = "Sampled") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  theme(
    axis.text = element_text(size = 15, color = "black", face = "bold"),
    axis.title = element_text(size = 15, face = "bold", margin = margin(t = 10, r = 10, b = 30, l = 10)),
    panel.background = element_rect(fill = F), 
    legend.text = element_text(size = 15, color = "black", face = "bold"),
  ) +
  scale_x_continuous(limits = c(2600, 3000), breaks = seq(2600, 3000, by = 100)) +
  scale_y_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_vline(xintercept = seq(2600, 3000, by = 100), linetype = "solid", color = "grey", size = 0.7)

