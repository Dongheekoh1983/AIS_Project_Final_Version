#################################################
### AIS Data Visualization Project Final Code ###
#################################################

# This project began in March 06th of 2023 at the request of Korean Navy
# The goal of this project was developing a Graphical User Interface
# that can help users to effectively visualize AIS (Automatic Identification System) data.
# The following work shows all the r source code that were used to develop the visualization application.
# FYI. After prototype was demonstrated using r, all the corresponding codes were translated into qt cpp.
#I am testing commit to github
######################################
### Setting up working environment ###
######################################
library(rgdal)
library(dplyr)
library(RSQLite)
library(sf)
library(ggplot2)
library(ggblend)
library(viridis)
library(ggmap)
library(gganimate)
library(osmdata)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(gganimate)
library(tmap)
library(units)
library(gifski)
library(DiagrammeR)
library(tidyverse)
library(rnaturalearth)
library(RMariaDB)
library(plotly)
library(move)
library(moveVis)
library(mapview)



#Set up working directory
setwd("/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization")


#bringing Korean peninsula shapefiles for visualization
south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                 quiet=TRUE, layer = "SouthKorea")
north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                 quiet=TRUE, layer = "NorthKorea")


###############################################
### Bringing AIS datasets from MySQL server ###
###############################################


storiesDb <- dbConnect(RMariaDB::MariaDB(), user = 'root',
                       password='gpryddl87',dbname = 'AIS_DATA_TABLES',
                       host='localhost')

dbListTables(storiesDb)


Dec_01 <- dbGetQuery(storiesDb, "SELECT *FROM AIS_DATA_ORIGINAL
                                 WHERE timestamp LIKE '2022-12-01%';")

BADA_NURI <- dbGetQuery(storiesDb, "SELECT *FROM AIS_DATA_ORIGINAL
                                 WHERE name = 'BADA NURI HO';")



Dec_01 <- read.csv("/Users/dongheekoh/Desktop/AIS_Dec01_Dec07/AIS_Korea_Dec_01.csv")
names(Dec_01)

########################################
####  AIS Trajectory Mining Workflow ###
########################################

# Point Pattern Analysis Visualization plan

DiagrammeR::grViz("               # All instructions are within a large character string
digraph surveillance_diagram {    # 'digraph' means 'directional graph', then the graph name 

  # graph statement
  #################

  graph [layout = dot,
         rankdir = LR,            # layout top-to-bottom
         fontsize = 10]

  
  # nodes (circles)
  #################

  node [shape = circle,           # shape = circle
       fixedsize = true
       width = 1.3]                      

  # Main tree
  Original  [label = 'Original\nPoint Data'] 
  Hexagon [label = 'Create hexagon\n(e.g.,20km)'] 
  Count  [label = 'CountPoints\nin polygon'] 
  Spatial_concentration [label = 'Concentration\nVisualization']
  Raster_surface [label = 'Kernel Density\nEstimation']
  Heat_map [label = 'Heat map', shape =  square, fontcolor=blue, color=blue]
  Isopleth [label = 'Isopleth map', fontcolor = darkgreen, shape = square,
  color=darkgreen]
  
  #Branch1 : Creating 3D and animation
  SplitMMSI [label = 'Split by\nMMSI']
  OrderTime [label = 'Order by\ntimestamp']
  MovingObject [label = 'Create\nmoving object']
  ThreeD [label = '3D time-space\ntrajectory\nvisualization', fontcolor = darkgreen,
  shape = square, color=darkgreen]
  Animation [label = 'Animation\n(e.g.,gif, mp4)', fontcolor = darkgreen,
  shape=square, color=darkgreen]

  #Branch2: Spatial autocorrelation
  Choropleth [label = 'Choropleth\nmap']
  Auto_cor  [label = 'Anselin\nLocal Morans I', 
  shape=square, color=orange, fontcolor=orange] 

  # edges
  ####### Main tree
  Original -> {Spatial_concentration Hexagon SplitMMSI}
  Spatial_concentration -> Raster_surface
  Raster_surface -> {Isopleth Heat_map}
  Hexagon -> Count                      
  Count -> Choropleth 
  
  #### Branch1: 3D and animation                        
  SplitMMSI -> OrderTime\
  OrderTime -> MovingObject[label = ' linear\n  interpolation', fontcolor=red]
  MovingObject -> ThreeD [style = dashed, color = darkgreen]
  MovingObject -> Animation [style = dashed, color = darkgreen]

  #### Branch2: Spatial Autocorrenation
  Choropleth -> Auto_cor[label = 'Statistical Validation\nof Spatial Clustering']
 }
")


# Line density analysis plan

DiagrammeR::grViz("               # All instructions are within a large character string
digraph surveillance_diagram {    # 'digraph' means 'directional graph', then the graph name 

  # graph statement
  #################
  graph [layout = dot,
         rankdir = LR,            # layout top-to-bottom
         fontsize = 10]

  # nodes (circles)
  #################
  node [shape = circle,           # shape = circle
       fixedsize = true
       width = 1.3]                      

  # Main tree
  Original  [label = 'Original\nPoint Data'] 
  SplitMMSI [label = 'Split by\nMMSI']
  OrderTime [label = 'Order by\ntimestamp']
  Threshold [label = 'Ti+1 - Ti >\nThreshold\n(e.g.,5hours)',
  shape=diamond, height=1.6, width=1.6, color=blue, fontcolor=blue]
  Separate_lines [label = 'CreateSeparate\nlines', color=blue]
  One_line [label = 'Create\ncontinuous line', color=blue, fontcolor=blue]
  GIS [label = 'Merge\ninto\nGIS\nlinestring']
  Line_Density [label =  'LineDensity\nVisualization',
  shape=square, height = 1.6, width = 1.6, color = orange]
  Detection_range [label = '레이더\n탐지범위\n가시화', 
  fontcolor = darkgreen, color=darkgreen, shape=square ]
  Join_attribute [label = 'JoinBy\nattribute', 
  fontcolor = black, color=black]
  Feature_blending [label = 'Feature\nblending\n효과', 
  fontcolor = darkgreen, color=darkgreen, shape = square]
  Traj_Clus [label =  'Trajectory\nClustering\n\n-DBSCAN\n-HDBSCAN\n-Kmedoids',
  shape=square, height = 1.6, width = 1.6, color = orange]
  Vis_Type [label = 'Visualization\nbyType', 
  shape=square, color=orange, fontcolor=orange]
  Filtering [label = 'Filtering\nby users']
  Overlay_point [label = 'Overlay\npointLayer']
  Speed [label='Visualization\n by speed', 
  shape=square, color=orange, fontcolor=orange]
  Heading[label='Visualization\n by heading',
  shape=square, color=orange, fontcolor=orange]

  # edges
  ####### Main tree
  Original ->  SplitMMSI 
  SplitMMSI -> OrderTime                      
  OrderTime -> Threshold
  Threshold -> Separate_lines[label = 'yes', fontcolor=red, style=dashed, color=blue]
  Threshold -> One_line [label = 'no', fontcolor=red, style=dashed, color=blue]
  {Separate_lines One_line} -> GIS[style=dashed, color=blue]  
  GIS -> {Line_Density Traj_Clus}
  Line_Density -> {Detection_range Join_attribute Feature_blending} [style=dashed,
  color=darkgreen]
  Join_attribute -> {Filtering Vis_Type}
  Filtering -> Overlay_point
  Overlay_point -> {Speed Heading} [style=dashed, color=orange]
  Vis_Type -> {Speed Heading} [style=dashed, color=orange]
  }
")


#############################################
### Developing AIS data cleaning function ###
#############################################

### Data cleaning process
# 1) delete duplicate timestamp within mmsi
# 2) delete duplicate xy coordinate within mmsi
# 3) delete NULL values if necessary
# 4) delete mmsi(s) whose total number of points are less than certain threshold (e.g., 10, 50)
# 5) also think about what to do for data whose speed, course, and heading values do not make sense
# 6) find out "outliers" and delete them if necessary

ais_data_cleaning <- function(df, threshold= 10,cutoff_mins=90  ,layer_name) {
  
  #extracting fields that are needed for visualization
  large_data <- df %>% dplyr::select(ship_and_cargo_type, name, mmsi, timestamp, course, speed, longitude, latitude)
  
  #creating "shiptype" column
  large_data <- large_data %>%
    mutate(SHIPTYPE = case_when(ship_and_cargo_type == "70" | ship_and_cargo_type == "71" |
                                ship_and_cargo_type == "79" ~ "Cargo",
                                ship_and_cargo_type == "80" ~ "Tanker",
                                ship_and_cargo_type == "52" ~ "Tug",
                                ship_and_cargo_type == "30" ~ "Fishing",
                                ship_and_cargo_type == "60" ~ "Passenger",
                                ship_and_cargo_type == "50" ~ "Pilot",
                                .default = "Other"))
  
  #creating a xy_coords field to delete duplicate coordinates
  large_data <- large_data %>%
    mutate(lon = round(longitude, 3),
           lat = round(latitude, 3),
           xy_combined = paste(as.character(lon), ", ",
                               as.character(lat))) %>% dplyr::select(-c(lon, lat))
  
  
  #split a dataset into a large list by mmsi
  split <- split(large_data, large_data$mmsi)
  
  #a function that removes duplicate timestamps
  f_time <- function(x) x[!duplicated(x[,c("timestamp")]),]
  
  #a function that removes duplicate cy_coords
  f_xy <- function(x) x[!duplicated(x[,c("xy_combined")]),]
  
  #applying timestamp removing function
  split <- lapply(split, f_time)
  
  #applying xy_coords removing function
  split <- lapply(split, f_xy)
  
  #converting a large list back into a dataframe again
  large_data <- do.call(what="rbind", split) %>% dplyr::select(-xy_combined)
  
  #selecting necessary fields in a right order
  large_data <- large_data %>% dplyr::select(SHIPTYPE, name, mmsi, timestamp, course, speed, longitude, latitude)
  
  #Rename columns
  colnames(large_data) <- c("SHIPTYPE","SHIPNAME", "MMSI", "TIMESTAMP", "COURSE", "SPEED", "LONGITUDE", "LATITUDE")
  
  #creating "HIGHER_TYPE" & "RADAR" & "SONAR" fields
  #the following fields are created as they are needed in the model developed for the NAVY
  large_data <- large_data %>% mutate(HIGHER_TYPES = case_when(SHIPTYPE == "Cargo" ~ "FIRST",
                                                               SHIPTYPE == "Tanker" ~ "SECOND",
                                                               SHIPTYPE == "Tug" ~ "THIRD",
                                                               SHIPTYPE == "Fishing" ~ "FOURTH",
                                                               SHIPTYPE == "Passenger" ~ "FIFTH",
                                                               SHIPTYPE == "Pilot" ~ "SIXTH",
                                                               .default = "OTHER")) %>%
    mutate(RADUIS = case_when(SHIPTYPE == "Cargo" ~ 15000,
                              SHIPTYPE == "Tanker" ~ 12000,
                              SHIPTYPE == "Tug" ~ 10000,
                              SHIPTYPE == "Fishing" ~ 8000,
                              SHIPTYPE == "Passenger" ~ 7000,
                              SHIPTYPE == "Pilot" ~ 6000,
                              .default = 5000)) %>%
    
    mutate(SONAR = case_when(SHIPTYPE == "Cargo" ~ 12000,
                             SHIPTYPE == "Tanker" ~ 9000,
                             SHIPTYPE == "Tug" ~ 7000,
                             SHIPTYPE == "Fishing" ~ 5000,
                             SHIPTYPE == "Passenger" ~ 4000,
                             SHIPTYPE == "Pilot" ~ 3000,
                             .default = 2000))
  
  #removing ', ' from the shipname field
  large_data$SHIPNAME <- gsub(","," ", large_data$SHIPNAME)
  
  # ChatGPT's suggestion to make previous code run faster (as of Feb/17/2024)
  library(data.table)
  
  # Convert to data.table
  setDT(large_data)

  # Convert TIMESTAMP to POSIXct
  large_data[, TIMESTAMP := as.POSIXct(TIMESTAMP)]
  
  # Sort by MMSI and TIMESTAMP
  setorder(large_data, MMSI, TIMESTAMP)
  
  # Compute time_lag
  large_data[, time_lag := TIMESTAMP - shift(TIMESTAMP, fill = first(TIMESTAMP)), by = MMSI]
  
  # Compute time_lag_exceeds_threshold
  cutoff_mins <- minutes(cutoff_mins)  # Adjust this threshold as needed
  large_data[, time_lag_exceeds_threshold := time_lag > cutoff_mins]
  
  # Compute group_id
  large_data[, group_id := cumsum(time_lag_exceeds_threshold), by = MMSI]
  
  # Compute MMSI_NEW
  large_data[, MMSI_NEW := paste0(MMSI, "_", group_id)]
  
  # Sort by MMSI_NEW and TIMESTAMP
  setorder(large_data, MMSI_NEW, TIMESTAMP)
  
  #counting number of points per MMSI
  count <- large_data %>% group_by(MMSI_NEW) %>% count()    
  
  #Inner_join count data to large data, then filtering mmsi whose total count is greater than or equal to 10 
  large_data <- large_data %>% inner_join(count, by=join_by(MMSI_NEW)) %>%
    mutate(MMSI = MMSI_NEW) %>% filter(n > threshold) %>% 
    dplyr::select(-c(n, time_lag, time_lag_exceeds_threshold, group_id, MMSI_NEW)) %>% 
    as.data.frame()

  #writing a cleaned ais data unto Global Environment with a new name
  assign(layer_name, large_data, envir = .GlobalEnv)
  
  #Compare before and after
  comparison <- function(data1, data2) {
    
    before <- dim(data1)[1]
    after <- dim(data2)[1]
    
    diff <- before - after
    
    cat("the total row number of the input dataset is", dim(data1)[1], '\n')
    cat("the total row number of the output dataset is", dim(data2)[1], '\n')
    cat("the difference between the two dataset is", diff)
  }
  
  comparison(df, large_data)
  
}


ais_data_cleaning(Dec_01, threshold = 10, cutoff_mins=180, layer_name = "Dec_01_NEW_DATA")
ais_data_cleaning(Dec_03, threshold = 10, cutoff_mins=180, layer_name = "Dec_03_NEW_DATA")
ais_data_cleaning(BADA_NURI, threshold = 2, cutoff_mins=180, layer_name = "BADA_NURI")

#measuring system processing time 
start <- Sys.time()
ais_data_cleaning(Dec_01, threshold = 10, cutoff_mins=180, layer_name = "Dec_01_NEW_DATA")
end <- Sys.time()
diff <- end - start
diff

Dec_01_NEW_DATA <- read.csv("Dec_01_Cleaned.csv")

#############################################
### How to write final output as csv file ###
#############################################

write.csv(Dec_01_NEW_DATA, 'Dec_01_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_03_NEW_DATA, 'Dec_03_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_04_NEW_DATA, 'Dec_04_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_05_NEW_DATA, 'Dec_05_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_07_NEW_DATA, 'Dec_07_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_08_NEW_DATA, 'Dec_08_Cleaned.csv', quote = FALSE, row.names = FALSE)
write.csv(Dec_10_NEW_DATA, 'Dec_10_Cleaned.csv', quote = FALSE, row.names = FALSE)


######################
### Patchwork Tool ###
######################

# the patchwork function combines and display multiple plots side by side for comparison purpose.
patchwork <- function(df1, df2=NULL, df3=NULL, df4=NULL, df5=NULL, df6=NULL, df7=NULL, ncol=2, title=NULL){
  library(patchwork)
  df1 + df2 + df3 + df4 + df5 + df6 + df7 +
    plot_layout(ncol=ncol) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(hjust=0.5)))
}


##################################
### Point map drawing function ###
##################################
point_map <- function(df, xmin=121, xmax=136, ymin=31, ymax=43) {
  
  library(sf)
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  ggplot(df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_point(mapping = aes(x=LONGITUDE, y=LATITUDE, colour=SHIPTYPE), size=0.5) +
    xlim(min(df$LONGITUDE - 0.5), max(df$LONGITUDE + 0.5)) +
    ylim(min(df$LATITUDE - 0.5), max(df$LATITUDE + 0.5)) + xlab(NULL) + ylab(NULL) + 
    xlim(xmin, xmax) + ylim(ymin, ymax) #+ theme(legend.position = "none")
  
}


fig1 <- point_map(Dec_01_NEW_DATA) # Visualizing Dec 3rd data
point_map(Dec_05_NEW_DATA) # visualizing Dec 5th data

jeju_point <- point_map(Dec_01_NEW_DATA, xmin=125, xmax=128, ymin=32.5, ymax=34 )

ggsave("./img/testing_1234.png", fig1, height = 5, width = 7, dpi = 320) # saving high resolution png image. 
ggsave("./img/jeju_point.png", jeju_point, height = 5, width = 7, dpi = 320) # saving high resolution png image. 

################################
### Heatmap Drawing Function ###
################################

# color presets that I can choose for drawing heat map
# "magma" (or "A")
# "inferno" (or "B")
# "plasma" (or "C")
# "viridis" (or "D")
# "cividis" (or "E")
# "rocket" (or "F")
# "mako" (or "G")
# "turbo" (or "H")

heatmap_custom <- function(df, bins = 150, option = "G") {
  
  library(viridis)
  library(sf)
  library(ggplot2)
  
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  ggplot() + geom_hex(data=df, aes(x=LONGITUDE, y=LATITUDE), bins = bins) +
    geom_sf(data=south) + geom_sf(data=north) +
    xlab("") + ylab("") +
    scale_fill_viridis(option=option, direction = -1, trans = "log") + #log scale for bin count
    theme(panel.background = element_rect("white"), #dark background
          axis.ticks = element_blank(),
          panel.grid = element_blank(), #remove panel grid
          axis.text = element_blank()) #remove x-axis value
  
}

heatmap_01 <- heatmap_custom(Dec_01_NEW_DATA, 1000)
heatmap_03 <- heatmap_custom(Dec_03_NEW_DATA, 1000)
heatmap_04 <- heatmap_custom(Dec_04_NEW_DATA, 1000)
heatmap_05 <- heatmap_custom(Dec_05_NEW_DATA, 1000)

patchwork(heatmap_01,heatmap_03,heatmap_04,heatmap_05, ncol=2)


#If we use ggplotly library, we make our map interactive.
ggplotly(heatmap_01)


#####################################
### Density map with contour line ###
#####################################
# biostas.w.uib.no/creating-a-2d-density-plot/.  <- refer to this source when writing portfolio
contour_map <- function(df, xmin=121, xmax=136, ymin=31, ymax=43) {
  
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  df <- df %>% filter(LONGITUDE > xmin & LONGITUDE < xmax & LATITUDE > ymin & LATITUDE < ymax)
  
  ggplot(data = df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_density_2d(mapping = aes(x=LONGITUDE, y=LATITUDE, alpha=0.7), lwd=0.1) +
    geom_density_2d_filled(mapping = aes(x=LONGITUDE, y=LATITUDE, alpha=0.7)) +
    xlim(xmin, xmax) + ylim(ymin, ymax) +
    theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
}

contour_01 <- contour_map(Dec_01_NEW_DATA)
contour_03 <- contour_map(Dec_03_NEW_DATA)

patchwork(contour_01, contour_03, ncol=2)

jeju <- contour_map(Dec_01_NEW_DATA, xmin=125, xmax=128, ymin=32.5, ymax=34) #around jeju subset
incheon <- contour_map(Dec_01_NEW_DATA, xmin=125, xmax=127.5, ymin=36.5, ymax=37.7) #around incheon subset

kde_contour <- patchwork(jeju, incheon, ncol=2)

ggsave("./img/jeju_kde_contour.png", jeju, height = 5, width = 7, dpi = 320) # saving high resolution png image. 


contour_map(Dec_01_NEW_DATA, xmin=127, xmax=130.5, ymin=33.5, ymax=36) #around Busan

#######################
### Points to Paths ###
#######################

# 1) points to paths using ggplot2.

point_line <- function(df, title = "AIS Vis", linewidth = 0.07, color="purple") {
  
  library(sf)
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  ggplot(data=df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_path(linewidth=linewidth, color=color,
              mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI)) +
    labs(title = title) + xlab(NULL) + ylab(NULL)
}

point_line_01 <- point_line(Dec_01_NEW_DATA, title = "Dec_01")
point_line_03 <- point_line(Dec_03_NEW_DATA, title = "Dec_03", color="orange")
point_line_04 <- point_line(Dec_04_NEW_DATA, title = "Dec_04", color="pink")

patchwork(point_line_01, point_line_03, point_line_04, ncol=3)

# 2) points to path : different approach
# The following approach is from
# "https://geobgu.xyz/r/processing-spatio-temporal-data.html" 11.4.4 which explains
# "To transform a point layer of object locations over time to a line layer of trajectories,
# we go through the following steps"
# 1) "split" the point layer to subset of points for each mmsi
# 2) "sort" each group of points chronologically, from the earliest to the latest
# 3) "transform" each of the point sequence to a line
# 4) "combine" separate lines back into a single polyline layer
# 5) Export polyline layer as a stanalone shapefile so as to be used for visualization
# 5) Make a function called points to path which include above algorithm
# 6) Also find out how I can apply "feature blending" effect to AIS paths

#1) developing a function that turns points to lines

pnt_data_frame <- Dec_01_NEW_DATA

points_to_paths <- function(pnt_data_frame, threshold, layer_name){
  
  #filtering an input data so each vessel should have at least no. of points greater than threhold
  mmsi_count <- pnt_data_frame %>% group_by(MMSI) %>% count()
  
  ais_df <- inner_join(pnt_data_frame, mmsi_count, by = "MMSI")
  
  ais_df <- ais_df %>% filter(n > threshold) #here is threshold parameter
  
  #creating sf point object
  pnt <- st_as_sf(ais_df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  
  #splitting the above 'pnt' layers by mmsi (each ship)
  lines <- split(pnt, pnt$MMSI)
  
  #timestamp ordering function
  f <- function(x) x[order(x$TIMESTAMP),]
  
  #apply time ordering function to all the list in line object
  lines <- lapply(lines, f)
  
  #next step is to combine all point geometries to a single multipoint geometry.
  #using "st_combine" to do that
  
  lines <- lapply(lines, st_combine)
  
  #casting multipoint geometry into linestring object
  lines <- lapply(lines, st_cast, to = "LINESTRING")
  
  # At this stage we have a list of 16,130 individual "linestring" geometries, one for each ship. The list can be combined back to an sfc geometry column using do.call
  geom <- do.call(c, lines)
  
  #transforming geom object into sf object
  layer_lines <- st_as_sf(geom, data.frame(id = names(lines)))
  
  # Assigning the result to a variable in the global environment
  layer_lines_global <- layer_lines
  assign(layer_name, layer_lines_global, envir = .GlobalEnv)
  
  #drawing the output on a tm_map
  #library(mapview)
  #mapview(layer_lines, lwd = 0.1, legend=FALSE)
}

points_to_paths(Dec_01_NEW_DATA, 50, "linestring_amazing")

library(mapview)
mapview(linestring_amazing, lwd=0.1, legend=FALSE)


#######################################
### Adding feature blending effects ###
#######################################

# Here, I am adding feature blending effects to AIS trajectories
# An obvious advantage of employing this approach is that it highlights areas where trajectory densities are higher,
# therefore giving more intuitive interpretation of the overall spatial distribution patterns of AIS data.????
# Currently available feature blending mode includes
#"overlay", "add", "saturate", "multiply", "screen", "overlay", "darken", "lighten", "color.dodge",??
#"color.burn", "hard.light", "soft.light", "difference", and "exclusion"
# source: https://mjskay.github.io/ggblend/


point_line_blending_effect <- function(df, title = "AIS Vis", linewidth = 0.05, color="purple" ,blend = "add") {
  
  library(sf)
  library(ggblend)
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  ggplot(data=df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_path(linewidth=linewidth, color=color,
              mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI)) * blend(blend=blend) +
    labs(title = title) + xlab(NULL) + ylab(NULL)
}


#comparing trajectory plots with or without feature blending effects (using ggplot)
plain_paths <- point_line(Dec_01_NEW_DATA, title="Plain", color="darkorange")
blending_path_multiply <- point_line_blending_effect(Dec_01_NEW_DATA, blend="multiply", title="Multiply", color="darkorange")
blending_path_add <- point_line_blending_effect(Dec_01_NEW_DATA, blend="add", title = "Add", color="darkorange")
feature_blending <- patchwork(plain_paths, blending_path_multiply, blending_path_add, ncol=3)

ggsave("./img/feature_blending.png", feature_blending, height = 2.5, width = 5.3, dpi = 320) # saving high resolution png image.

#applying a feature blending effect to a shapefile produced by "points_to_paths" tool developed above
ggplot() + geom_sf(data=south) + geom_sf(data=north) +
  geom_sf(data=linestring_amazing, lwd=0.07, color = "purple") * blend(blend="multiply")


##############################
### Time filtering example ###
##############################
# 0)First, defining TIMESTAMP field using as.POSIXct

Dec_01_NEW_DATA_m <- Dec_01_NEW_DATA

# 1) Second, filtering timestamp field using "ymd_hms"
time_filtering1 <- Dec_01_NEW_DATA_m %>%
  filter(TIMESTAMP > ymd_hms("2022-12-01 18:00:00"), TIMESTAMP < ymd_hms("2022-12-01 20:00:00")) %>%
  arrange(TIMESTAMP)

# 2) Second, filtering timestamp after turning timestamp field into "POSIXct"
time_filtering2 <- Dec_01_NEW_DATA_m %>%
  filter(TIMESTAMP > as.POSIXct("2022-12-01 13:00:00", tz="UTC"),
         TIMESTAMP < as.POSIXct("2022-12-01 18:00:00", tz="UTC")) %>%
  arrange(MMSI,TIMESTAMP)

# 3) creating time filtering function
time_filtering <- function(df, min, max) {
  df %>% filter(TIMESTAMP > as.POSIXct(min, tz="UTC"),
                TIMESTAMP < as.POSIXct(max, tz="UTC")) %>%
    arrange(MMSI, TIMESTAMP)
}

test <- time_filtering(Dec_01_NEW_DATA_m, min="2022-12-01 00:00:00", max = "2022-12-01 02:00:00")

####################################################################################################
### BADA NURI HO subset - to demonstrate the need for improving current points_to_path algorithm ###
#####################################################################################################

#Extracting "BADA NURI HO" for demonstration
validation <- Dec_01_NEW_DATA_m %>% filter(SHIPNAME=="BADA NURI HO") %>% arrange(TIMESTAMP)
#validation$TIMESTAMP <- as.POSIXct(validation$TIMESTAMP)
#extract xy extent values from the validation dataset

x_min <- min(validation$LONGITUDE) - 0.1
x_max <- max(validation$LONGITUDE) + 0.1
y_min <- min(validation$LATITUDE) - 0.1
y_max <- max(validation$LATITUDE) + 0.1

#point_line function modification
point_line_mod <- function(df, title) {
  
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  
  ggplot(data=df) + geom_sf(data=south) + geom_point(mapping = aes(x=LONGITUDE, y=LATITUDE)) +
    geom_path(mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI)) +
    theme(axis.text.x = element_blank(), axis.text.y=element_blank()) +
    xlim(x_min, x_max) + labs(title = title) +
    ylim(y_min, y_max)
  
}

#use the time_filtering function.
validation_01 <- time_filtering(validation, min = "2022-12-01 04:00:00", max = "2022-12-01 08:00:00")
validation_02 <- time_filtering(validation, min = "2022-12-01 21:00:00", max = "2022-12-01 24:00:00")

#Drawing maps by time period
first <- point_line_mod(validation, "Entire")
second <- point_line_mod(validation_01, "Between 4am and 8am")
third <- point_line_mod(validation_02, "Between 9pm and 12am")


#displaying plots together
patchwork(first, second, third, ncol=2)

#(Summary) As can be clearly seen from the plot above, there is actually a discontinuation in the path of "BADA NURI HO".
# I believe we can easily observe this problem in any other trajectories in our dataset.
# In this project so far, I have not really handled this issue.
# But knowing this now, I need to improve upon this limitation in the future development.

########################################################
### Filtering by type and visualization using paths  ###
########################################################
# In this project, I have created following SHIPTYPE categories
# "Other", "Cargo", "Fishing", "Tanker", "Passenger", "Tug", "Pilot"
# In this section, I am developing a function that can visualize AIS data by ships' major types as defined above.

filtering_vis_type <- function(df, type, color = "purple") {
  
  df <- df
  df <- df %>% filter(SHIPTYPE == type)
  df <- df %>% arrange(MMSI, TIMESTAMP)
  
  library(sf)
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  
  ggplot(data=df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_path(linewidth=0.09, color=color, mapping = aes(x=LONGITUDE, y=LATITUDE, group=MMSI))
  
}


Cargo <- filtering_vis_type(Dec_01_NEW_DATA, "Cargo", color="#0000a2")
Fishing <- filtering_vis_type(Dec_01_NEW_DATA, "Fishing", color="#e9c716")
Tug <- filtering_vis_type(Dec_01_NEW_DATA, "Tug", color="#bc272d")
Passenger <- filtering_vis_type(Dec_01_NEW_DATA, "Passenger", color="#50ad9f")


patchwork(Cargo, Fishing, Tug, Passenger)

####################################
### Visualization with Animation ###
####################################

# what will be added later in my final product.
#(comments regarding the advantages of animation plots when visualizing AIS data)

ais_animation <- function(df, xmin=121, xmax=136, ymin=31, ymax=43) {
  
  #Calling necessary libraries
  library(ggplot2)
  library(gganimate)
  
  #Defining a data frame and selecting necessary fields
  ais_df <- df
  tracks <- ais_df %>% dplyr::select(LONGITUDE, LATITUDE, MMSI, TIMESTAMP, SHIPNAME, SPEED, SHIPTYPE)
  
  #Casting a timestamp field from character to POSIXct
  tracks$TIMESTAMP <- parse_date_time(tracks$TIMESTAMP, "Ymd HMS")
  
  #Re-arranging tracks data by timestamp for each mmsi
  tracks_ordered <- tracks %>% arrange(MMSI, TIMESTAMP)
  
  #Deleting duplicate timestamp records
  tracks <- tracks_ordered[!duplicated(tracks_ordered[,c("TIMESTAMP")]),]
  
  #Casting MMSI field from integer to factor
  tracks$MMSI <- as.factor(tracks$MMSI)
  
  #Filtering data - ships that has greater than 10 location points
  count <- tracks %>% group_by(MMSI) %>% count() 
  tracks <- inner_join(tracks, count, by = "MMSI") %>% filter(n > 10)
  
  #Creating a Large MoveStack data
  mdata <- move(x=tracks$LONGITUDE, y=tracks$LATITUDE, time=tracks$TIMESTAMP, 
                data=tracks, proj = CRS("+proj=longlat +ellps=WGS84"),
                animal=tracks$MMSI)
  
  #Aligning movement data to a uniform time scale with a uniform temporal resolution throughout the data
  m <- align_move(mdata, res=5, digit=0, unit="mins", spaceMethod = "greatcircle")
  
  #Casting Large MoveStack data into a data.frame again
  tracks <- as.data.frame(m)
  
  #Producing an animation plot using ggplot and gganimate (original one)
  ggplot() + geom_sf(data=south) + geom_sf(data=north) +
    geom_point(data = tracks, aes(x=x,y=y, group=trackId, colour=trackId), alpha = 0.7, shape=20) +
    theme(legend.position = 'none') + transition_time(time) + labs(title = "Time: {frame_time}") +
    shadow_mark(alpha = 0.3, size = 0.5) + xlim(xmin, xmax) + ylim(ymin, ymax)
  
}


#Running and testing ais_animation function and cropping an animation according to your needs
ais_animation(Dec_01_NEW_DATA) #showing entire study area
ais_animation(Dec_01_NEW_DATA, xmin=126, xmax=132, ymin=32, ymax=36) #Showing only the south east region of South Korea
ais_animation(Dec_01_NEW_DATA, xmin=122, xmax=128, ymin=35, ymax=38) #Showing only the west region of South Korea


x_min <- min(BADA_NURI$LONGITUDE) - 0.1
x_max <- max(BADA_NURI$LONGITUDE) + 0.1
y_min <- min(BADA_NURI$LATITUDE) - 0.1
y_max <- max(BADA_NURI$LATITUDE) + 0.1
ais_animation(BADA_NURI, xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max)

#if you want to save your animation on your local drive!
anim_save("south_shadow_mark.gif")
anim_save("west_shadow_mark.gif")

#Issue I encountered as of Feb/19th/2024
# After refining the data_cleaning function, animation function started giving me an error message that says 
# length of 'dimnames' [1] not equal to array extent
# with regard to this, I found a thread which seemed to suggest potential solution which says
#> Temporal resolution of 5 [mins] is used to align trajectories.
#> Warning: The full temporal coverage of at least one trajectory is shorter than the specified resolution. 
#> You may want to choose a finer resolution. or you may want to remove trajectories that are two short
#> In my case, in the data_cleaning function, I changed the threhold = 100 from threshold = 10



###########################
#Spire like animation plot#
###########################
#https://spire.com/blog/maritime/the-red-sea-crisis-tracking-the-volatile-security-situation/?utm_campaign=maritime_monthly_newsletter&utm_medium=email&_hsmi=291862574&_hsenc=p2ANqtz-9kHA0Nr3rV--HLPySGU540cOQUZcjW72_JBnGgRQniAehN9Kq1_WRHz8jYHaYz9q4Ca--s_n8IUYTaFBhc3fO4EeRtkQ&utm_content=291841014&utm_source=hs_email

ais_df <- Dec_01_NEW_DATA
tracks <- ais_df %>% dplyr::select(LONGITUDE, LATITUDE, MMSI, TIMESTAMP, SHIPNAME, SPEED, SHIPTYPE)

#Casting a timestamp field from character to POSIXct
tracks$TIMESTAMP <- parse_date_time(tracks$TIMESTAMP, "Ymd HMS")

#Re-arranging tracks data by timestamp for each mmsi
tracks_ordered <- tracks %>% arrange(MMSI, TIMESTAMP)

#Deleting duplicate timestamp records
tracks <- tracks_ordered[!duplicated(tracks_ordered[,c("TIMESTAMP")]),]

#Casting MMSI field from integer to factor
tracks$MMSI <- as.factor(tracks$MMSI)

#Filtering data - ships that has greater than 10 location points
count <- tracks %>% group_by(MMSI) %>% count() %>% arrange(desc(n))
tracks <- inner_join(tracks, count, by = "MMSI") %>% filter(n > 50)

#Creating a Large MoveStack data
mdata <- move(x=tracks$LONGITUDE, y=tracks$LATITUDE, time=tracks$TIMESTAMP, data=tracks, proj = CRS("+proj=longlat +ellps=WGS84"),
              animal=tracks$MMSI)

#Aligning movement data to a uniform time scale with a uniform temporal resolution throughout the data
m <- align_move(mdata, res=5, digit=0, unit="mins", spaceMethod = "greatcircle")

#Casting Large MoveStack data into a data.frame again
tracks <- as.data.frame(m)

#Creating an animation 
ggplot() + geom_sf(data=south, fill="black") +
  geom_path(Dec_01_NEW_DATA, linewidth=0.1, color="#4A5D75", alpha=0.7 , mapping = aes(x=LONGITUDE, y=LATITUDE, group=MMSI)) +
  xlim(126,132) + ylim(32,36) +
  theme(panel.background=element_rect(fill="#282D34"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(data = tracks, aes(x=x,y=y, group=trackId, colour=trackId),alpha = 0.7, shape=20) +
  theme(legend.position = 'none') + transition_time(time) + labs(title = "Time: {frame_time}") +
  shadow_wake(wake_length = 0.2, size = 0.5)

anim_save("south_shadow_wake_colorful.gif")

##########################################################
### Time Space Trajectory 3D Plot for AIS trajectories ###
##########################################################
#What will be added later
#comments regarding the advantages of using time-space trajectory plots
# 1) Extract 4 paths
# 2) Draw 3D time space trajectory plot

#1) Extract 4 paths to be visualized
count <- Dec_01_NEW_DATA %>% group_by(MMSI) %>% count() %>% arrange(desc(n))
ais_3d <- Dec_01_NEW_DATA
ais_3d$MMSI <- as.character(ais_3d$MMSI)
tracks_3d <- ais_3d %>% filter(MMSI == "636020542_0"|MMSI == "352002128_0"|MMSI == "636012808_0"|MMSI == "432596000_0")

#2) Defininfg timestamp field as it will be used as a z-axis in the 3d visualization
tracks_3d$TIMESTAMP <- as.POSIXct(tracks_3d$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

#3) Arranging the dataset by MMSI and TIMESTAMP
tracks_3d <- tracks_3d %>% arrange(MMSI, TIMESTAMP)

#4) Finally 3D visualization
fig_4paths <- plot_ly(tracks_3d, x = ~LONGITUDE, y = ~LATITUDE, z = ~TIMESTAMP, type = 'scatter3d', color = ~MMSI,
                      mode = 'lines') %>% add_markers(color=~MMSI, size = I(15)) %>%
  layout(scene = list(xaxis = list(title = 'Longitude'),
                      yaxis = list(title = 'Latitude'),
                      zaxis = list(title = "Time")))

fig_4paths

ais_animation(tracks_3d, xmin=126, xmax=132, ymin=32, ymax=36) #comparing with animation results

####################################################
### Making hexagaon grid and draw choropleth map ###
####################################################

choropleth_map <- function(df, cell_size=50, square=FALSE, interactive = FALSE) {
  
  ais_df_test <- df %>% filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326 , remove = FALSE)
  
  #converting AIS dataframe into st_geometry() and transform its projection
  ais_df_test <- ais_df_test %>% st_set_crs(4326) %>% st_transform(3857)
  
  #Creating hex grid; for now I am using 70*1000 as a cell size; EPSG3857 is in meter, therefore 70*1000 means 70 km
  grd <- st_make_grid(ais_df_test, cell_size*1000, square=square, flat_topped = TRUE)
  
  # To sf and add grid ID
  fishnet_grid_sf = st_sf(grd) %>%
    # add grid ID
    mutate(grid_id = 1:length(lengths(grd)))
  
  # intersect fishnet grid and AIS point data, then counting how many points there are per cell
  fishnet_grid_sf$n_colli = lengths(st_intersects(fishnet_grid_sf, ais_df_test))
  
  # remove grids whose values are 0 (i.e. no points in side that grid)
  fishnet_count = filter(fishnet_grid_sf, n_colli > 0)
  
  if(interactive == FALSE){
    
    ggplot() +
      geom_sf(data=fishnet_count, aes(fill=n_colli), color="gray", lwd=0.1) +
      theme_bw() + scale_fill_viridis_c(option="D", alpha=0.5, direction=-1) +
      geom_sf(data=south) + geom_sf(data=north) + theme(legend.position='none')
    
  } else if (interactive == TRUE) {
    
    tmap_mode("view")
    
    map_fishnet = tm_shape(fishnet_count) +
      tm_fill(
        col = "n_colli",
        palette = "YlOrRd",
        style = "cont",
        title = "Number of AIS points",
        id = "grid_id",
        showNA = FALSE,
        alpha = 0.5,
        popup.vars = c("Number of AIS: " = "n_colli"),
        popup.format = list(n_colli = list(format = "f", digits = 0))) +
      tm_borders(col = "grey40", lwd = 0.1)
    
    map_fishnet
  }
}

choropleth_map(Dec_01_NEW_DATA, cell_size=30, square=FALSE, interactive=TRUE) # interactive=TRUE
choropleth_map(Dec_01_NEW_DATA, cell_size=30, square=FALSE, interactive=FALSE) # interactive=FALSE

size_30_km <- choropleth_map(Dec_01_NEW_DATA, cell_size=30, square=FALSE, interactive=FALSE) # interactive=FALSE
size_50_km <- choropleth_map(Dec_01_NEW_DATA, cell_size=50, square=FALSE, interactive=FALSE) # interactive=FALSE

choropleth_two_combined <- patchwork(size_30_km, size_50_km, ncol=2)
ggsave("./img/choropleth_three_combined.png", choropleth_two_combined, height = 4, width = 8, dpi = 320) # saving high resolution png image. 

#########################################
### Buffer (Detection Range) Analysis ###
#########################################

# Korean Navy wanted to 1) visualize how their vessels' radar/somar detection range accumulates as time progresses
# 2) to calculate total coverage areas per fixed time period
# So I am here applying buffer analysis for that end

buffer_union_analysis <- function(df, dist_or_type=NULL, buffer_dist=20){
  
  if(!is.numeric(dist_or_type) & !is.character(dist_or_type)) {
    
    print("Please input 'dist_or_type' parameter")
    
  } else {
    
    #1) Casting MMSI into character
    df$MMSI <- as.character(df$MMSI)
    
    #2) Extracting unique MMSIs from an original data
    inter <- df %>% distinct(MMSI, .keep_all = TRUE)
    
    #3) Points to paths
    points_to_paths(df, threshold = 50, layer_name = "linestring")
    
    #4) Joining linestgin polygon shapefile with attribute data
    join <- linestring %>% left_join(inter, by = c("id" = "MMSI"))
    
    #5) Re-projecting layers so as to have a meter as its unit
    join <- join %>% st_set_crs(4326) %>% st_transform(3857)
    
    #6) Calculating the lengths of trajectories using st_length
    join <- join %>% mutate(distance = as.numeric(st_length(.))) %>% arrange(desc(distance))
    
    
    if(is.numeric(dist_or_type)) {
      
      trajs <- join %>% filter(distance > dist_or_type)
      trajs_buffer <- st_buffer(trajs, buffer_dist*1000)
      
      #union and calculate area
      trajs_buffer_union <- trajs_buffer %>% st_union() %>% st_as_sf() %>%mutate(areas = as.numeric(st_area(.))/1000000)
      
      #draw both buffer and union map side by side
      
      buffer_plot <- ggplot() + geom_sf(data=south) + geom_sf(data=north) + geom_sf(data=trajs_buffer, fill="pink", alpha=0.5)
      union_plot <- ggplot() + geom_sf(data=south) + geom_sf(data=north) + geom_sf(data=trajs_buffer_union, fill="pink", alpha=0.5)
      
      patchwork(buffer_plot, union_plot,
                title=paste("the total area covered by chosen vessels is", round(trajs_buffer_union$areas[1],2), 'm\u00B2'))
      
      
    } else if(is.character(dist_or_type)) {
      
      ship_type <- join %>% filter(SHIPTYPE == dist_or_type)
      ship_type_buffer <- st_buffer(ship_type, buffer_dist*1000)
      
      #union and calculate area
      ship_type_buffer_union <- ship_type_buffer %>% st_union() %>% st_as_sf() %>%mutate(areas = as.numeric(st_area(.))/1000000)
      
      #draw both buffer and union map side by side
      
      buffer_plot <- ggplot() + geom_sf(data=south) + geom_sf(data=north) + geom_sf(data=ship_type_buffer, fill="pink", alpha=0.5)
      union_plot <- ggplot() + geom_sf(data=south) + geom_sf(data=north) + geom_sf(data=ship_type_buffer_union, fill="pink", alpha=0.5)
      
      patchwork(buffer_plot, union_plot,
                title=paste("the total area covered by chosen vessels is", round(ship_type_buffer_union$areas[1],2), 'km\u00B2'))
    }
    
  }
  
}

Tug <- buffer_union_analysis(Dec_01_NEW_DATA, dist_or_type = "Tug")
ggsave("./img/Tug_detection_range.png", Tug, height = 3.3, width = 5.3, dpi = 320) # saving high resolution png image.
buffer_union_analysis(Dec_01_NEW_DATA)

################################################
### Spatial Auto-correlation using Moran's I ###
################################################

# source: https://www.paulamoraga.com/book-spatial/spatial-autocorrelation.html
# https://carto.com/blog/spatial-hotspot-tools.   <- when writing portfolio refer to this site
# https://mgimond.github.io/Spatial/spatial-autocorrelation.html <- this one offers a through explanation about Morans' I anlaysis. 
# so refer to this source when writing a portfolio 
# 1) Turning data frame into sf object with also assigning crs=4326
# In order to see the spatial log effects more clearly, I tried log-transformation (as of Feb 27)
ais_df_moran <- Dec_01_NEW_DATA %>% filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326 , remove = FALSE)

# 2) re-projecting layer to EPSG 3857
ais_df_moran <- ais_df_moran %>% st_set_crs(4326) %>% st_transform(3857)

#Creating hex grid; for now I am using 70*1000 as a cell size; EPSG3857 is in meter, therefore 70*1000 means 70 km
grd <- st_make_grid(ais_df_moran, 50*1000, square=FALSE, flat_topped = TRUE)

# To sf and add grid ID
fishnet_grid_sf = st_sf(grd) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(grd))) 


# intersect fishnet grid and AIS point data, then counting how many points there are per cell
fishnet_grid_sf$n_colli = lengths(st_intersects(fishnet_grid_sf, ais_df_moran)) 
fishnet_grid_sf <- fishnet_grid_sf %>% filter(n_colli > 0)  #lengths(st_intersects(fishnet_grid_sf, ais_df_moran)) 
fishnet_grid_sf$n_colli_log = log(fishnet_grid_sf$n_colli) #log-transform

hist(fishnet_grid_sf$n_colli_log, breaks=100)

#fishnet_grid_sf <- fishnet_grid_sf %>% filter(n_colli > 0)
#visualize with ggplot
ggplot() + geom_sf(data=south) + geom_sf(data=north) +
  geom_sf(data=fishnet_grid_sf, aes(fill=n_colli_log), color="gray", lwd=0.1) +
  theme_bw() + scale_fill_viridis_c(option="D", alpha=0.5, direction=-1)

test_poly <- fishnet_grid_sf

library(spdep)

nb <- poly2nb(test_poly, queen=TRUE)
nbw <- nb2listw(nb, style="W", zero.policy = TRUE)

#drawing neighborhood connection map: spatial: neightborhood based on contiguity
plot(st_geometry(test_poly), border = "lightgray")
plot.nb(nb, st_geometry(test_poly) ,add=TRUE, lwd=0.4)

gmoran <- moran.test(test_poly$n_colli_log, nbw,
                     alternative = "greater", zero.policy = TRUE)

# If you want to do away with scientific notation type the following --> options(scipen=999)

gmoran

# the result show that Moran test statistics(=z-score) is 91.581 and associated p-value is near zero
# which indicates that the likelihood of observing this spatial patterns randomly is extremely small

# Monte Carlo approach creates random patterns by re-assigning the values among the fixed areas and calculates the Moran's I for each of these pattern
# As we can see both from the plot and p-value, the probability of observing tgiven spatial pattern is very small, thus our evidence suggests that
# there is a spatial clustering.
gmoranMC <- moran.mc(test_poly$n_colli_log, nbw, nsim = 4999, zero.policy = TRUE)
gmoranMC

gmoranMC$res

hist(gmoranMC$res, breaks=300, border="grey")
abline(v = gmoranMC$statistic, col = "red")

# we also observe a positive linear relationship between observations and their spatially lagged values.
moran.plot(test_poly$n_colli_log, nbw)

lmoran <- localmoran(test_poly$n_colli_log, nbw, alternative = "two.sided", zero.policy = TRUE)

library(tmap)
tmap_mode("view")


test_poly$lmI <- lmoran[,"Ii"]
test_poly$lmZ <- lmoran[,"Z.Ii"] #z-score
test_poly$lmp <- lmoran[,"Pr(z != E(Ii))"]

p1 <- tm_shape(test_poly) + tm_polygons(col = "n_colli_log", title="vble", style="quantile", alpha=0.5, lwd=0.05) + tm_layout(legend.outside = TRUE)
p2 <- tm_shape(test_poly) + tm_polygons(col = "lmI", title="Local Moran's I", style="quantile", alpha=0.5, lwd=0.05) + tm_layout(legend.outside = TRUE)
p3 <- tm_shape(test_poly) + tm_polygons(col = "lmZ", title="Z-score", breaks = c(-Inf, 1.65, Inf), alpha=0.5, lwd=0.05) + tm_layout(legend.outside = TRUE)
p4 <- tm_shape(test_poly) + tm_polygons(col = "lmp", title="p-value", breaks = c(-Inf, 0.05, Inf), alpha=0.5, lwd=0.05) + tm_layout(legend.outside = TRUE)

tmap_arrange(p1, p2, p3, p4)


tm_shape(test_poly) +
  tm_polygons(col = "lmZ",
              title = "Local Moran's I", style = "fixed",
              breaks = c(-Inf, -1.96, 1.96, Inf),
              labels = c("Negative SAC", "No SAC", " Positive SAC"),
              palette = c("blue", "white", "red"), alpha=0.5, lwd=0.05) +
  tm_layout(legend.outside = TRUE)

#######################
### Local Moran's I ###
#######################

# Local Moran's I alloows us to identify clusters of the following type.
# 1) High-High: areas of high values with neighbors of high values,
# 2) High-Low: areas of high values with neighbors of low values,
# 3) Low-High: areas of low values with neighbors of high values,
# 4) Low-Low: areas of low values with neighbors of low values.
# specifically, we identify the cluster types by using the quadrant of the scaled values(mp$x)
# their spatially lagged values (mp$wx), and the p-values obtained with the local Moran's I for each of the areas(mp$lmp)

lmoran <- 
  localmoran(test_poly$n_colli_log, nbw, alternative = "two.sided", zero.policy = TRUE)

test_poly$lmp <- lmoran[,5]

mp <- moran.plot(as.vector(scale(test_poly$n_colli_log)), nbw, zero.policy = TRUE)

test_poly$quadrant <- NA

#high-high
test_poly[(mp$x >= 0 & mp$wx >= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 1

#low-low
test_poly[(mp$x <= 0 & mp$wx <= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 2

#high-low
test_poly[(mp$x >= 0 & mp$wx <= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 3

#low-high
test_poly[(mp$x <= 0 & mp$wx >= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 4

#non-significant
test_poly[(test_poly$lmp > 0.05), "quadrant"] <- 5


#Drawing interactive map with tm package
tm_shape(test_poly) + tm_fill(col = "quadrant", title = "", breaks = c(1,2,3,4,5,6),
                              palette = c("red", "blue", "lightpink", "skyblue2", "white"),
                              labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"), alpha=0.4) +
  tm_legend(text.size=1) + tm_borders(alpha=0.5, lwd=0.05) + tm_layout(frame=FALSE, title = "Clucters") + tm_layout(legend.outside = TRUE)

#Drawing static map with ggplot
ggplot() + geom_sf(data=south) + geom_sf(data=north) +
  geom_sf(data = test_poly, alpha=0.5,aes( fill = as.factor(quadrant)), lwd=0.03) +
  scale_fill_manual("Local Morans'I",
                    values = c("1" = "red",
                               "2" = "blue",
                               "3" = "lightpink",
                               "4" = "skyblue",
                               "5" = "white"),
                    labels = c("High-high", "Low-low","Low-high" ,"non-significant"))


#################################
### Moran's I into a function ###
#################################

morans_I <- function(df, nsim=999, size=20, interactive = FALSE){
  
  ais_df_moran <- df %>% filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
    
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326 , remove = FALSE)
  
  ais_df_moran <- ais_df_moran %>% st_set_crs(4326) %>% st_transform(3857)
  
  grd <- st_make_grid(ais_df_moran, size*1000, square=FALSE, flat_topped = TRUE)
  
  fishnet_grid_sf = st_sf(grd) %>%
    mutate(grid_id = 1:length(lengths(grd)))
  
  fishnet_grid_sf$n_colli = lengths(st_intersects(fishnet_grid_sf, ais_df_moran))
  
  test_poly <- fishnet_grid_sf
  
  nb <- poly2nb(test_poly, queen=TRUE)
  nbw <- nb2listw(nb, style="W")
  
  gmoran <- moran.test(test_poly$n_colli, nbw,
                       alternative = "greater")
  print(gmoran)
  
  gmoranMC <- moran.mc(test_poly$n_colli, nbw, nsim = nsim)
  
  print(gmoranMC)
  
  lmoran <- localmoran(test_poly$n_colli, nbw, alternative = "two.sided")
  
  test_poly$lmp <- lmoran[,5]
  
  mp <- moran.plot(as.vector(scale(test_poly$n_colli)), nbw)
  
  test_poly$quadrant <- NA
  
  #high-high
  test_poly[(mp$x >= 0 & mp$wx >= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 1
  
  #low-low
  test_poly[(mp$x <= 0 & mp$wx <= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 2
  
  #high-low
  test_poly[(mp$x >= 0 & mp$wx <= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 3
  
  #low-high
  test_poly[(mp$x <= 0 & mp$wx >= 0) & (test_poly$lmp <= 0.05), "quadrant"] <- 4
  
  #non-significant
  test_poly[(test_poly$lmp > 0.05), "quadrant"] <- 5
  
  if(interactive==FALSE) {
    
    ggplot() + geom_sf(data=south) + geom_sf(data=north) +
      geom_sf(data = test_poly, alpha=0.5, aes(fill = as.factor(quadrant)), lwd=0.03) +
      scale_fill_manual("Local Morans'I",
                        values = c("1" = "red",
                                   "2" = "blue",
                                   "3" = "lightpink",
                                   "4" = "skyblue",
                                   "5" = "white"),
                        labels = c("High-high", "Low-high", "non-significant"))
    
  } else {
    
    tmap_mode("view")
    tm_shape(test_poly) + tm_fill(col = "quadrant", title = "", breaks = c(1,2,3,4,5,6),
                                  palette = c("red", "blue", "lightpink", "skyblue2", "white"),
                                  labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-significant"), alpha=0.4) +
      tm_legend(text.size=1) + tm_borders(alpha=0.5, lwd=0.05) + tm_layout(frame=FALSE, title = "Clusters") + tm_layout(legend.outside = TRUE)
    
  }
  
}

morans_I(Dec_01_NEW_DATA, size=30, interactive = FALSE)
morans_I(Dec_01_NEW_DATA, size=10, interactive = FALSE)


#subset (around jeju) Morans'I testing
subset <- Dec_01_NEW_DATA %>% filter(LONGITUDE > 125 & LONGITUDE < 128 & LATITUDE > 32 & LATITUDE < 34)
subset_fishing <- Dec_01_NEW_DATA %>% filter(LONGITUDE > 125 & LONGITUDE < 128 & LATITUDE > 32 & LATITUDE < 34, SHIPTYPE == "Fishing")
subset_cargo <- Dec_01_NEW_DATA %>% filter(LONGITUDE > 125 & LONGITUDE < 128 & LATITUDE > 32 & LATITUDE < 34, SHIPTYPE == "Cargo")
subset_passenger <- Dec_01_NEW_DATA %>% filter(LONGITUDE > 125 & LONGITUDE < 128 & LATITUDE > 32 & LATITUDE < 34, SHIPTYPE == "Passenger")
subset_tug <- Dec_01_NEW_DATA %>% filter(LONGITUDE > 125 & LONGITUDE < 128 & LATITUDE > 32 & LATITUDE < 34, SHIPTYPE == "Tug")


fishing <- morans_I(subset_fishing, size=3, interactive = TRUE)
cargo <- morans_I(subset_cargo, size=7, interactive = TRUE)
passenger <- morans_I(subset_passenger, size=3, interactive = TRUE)
tug <- morans_I(subset_tug, size=3, interactive = TRUE)

tmap_arrange(fishing, cargo, passenger, tug)

# why don't I compare plot Dec01 ~ Dec 14
# by time as well (extract data from DB by timestamp and draw the map)


###################################################
### Trajectory clustering Analysis using python ###
###################################################

  
######################################
##Bringing DBSCAN result from python##
######################################

dbscan <- read.csv("/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization/dbscan_result.csv")

dbscan <- dbscan %>% mutate(MMSI = paste0(X.1,"_", Cluster_ID)) %>% 
  rename(LONGITUDE=X, LATITUDE=Y, TIMESTAMP=Index) %>% mutate(Cluster_ID = as.factor(Cluster_ID)) 

names(dbscan)

point_line <- function(df, title = "AIS Vis", linewidth = 0.07, col = "purple") {
  
  library(sf)
  south <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "SouthKorea")
  north <- st_read('/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization',
                   quiet=TRUE, layer = "NorthKorea")
  

  ggplot(data=df) + geom_sf(data=south) + geom_sf(data=north) +
    geom_path(linewidth=linewidth, color = col,  
              mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI, colour=Cluster_ID)) +
    theme(legend.position = "none", text=element_text(size=7))+
    labs(title = title) + xlab(NULL) + ylab(NULL) + xlim(125,126.8) + ylim(36.2, 37.6)
}

point_line(dbscan, title="DBSCAN", linewidth=0.16)

cluster_0 <- dbscan %>% filter(Cluster_ID==0)
db_first <- point_line(cluster_0, title = "Cluster_0", linewidth = 0.3, col="#377eb8")

cluster_1 <- dbscan %>% filter(Cluster_ID==1)
db_second <- point_line(cluster_1, title = "Cluster_1", linewidth = 0.3, col= "#4daf4a")

cluster_4 <- dbscan %>% filter(Cluster_ID==4)
db_third <- point_line(cluster_4, title = "Cluster_4", linewidth = 0.3, col="#984ea3")

cluster_10 <- dbscan %>% filter(Cluster_ID==-1)
db_fourth <- point_line(cluster_10, title = "noise",linewidth = 0.3, col="#e41a1c")

dbscan_results_sep <- patchwork(db_first,db_second,db_third,db_fourth, ncol=2)

ggsave("./img/dbscan_results_sep.png", dbscan_results_sep, height = 5, width = 5, dpi = 320) # saving high resolution png image. 

#Creating a combined plot that shows four clusters together. 
dbscan_subset <- dbscan %>% filter(Cluster_ID == 0 | Cluster_ID == 1 | Cluster_ID == 4 | Cluster_ID == -1)

dbscan_results_all <-
  ggplot(data=dbscan_subset) + geom_sf(data=south) + geom_sf(data=north) +
  geom_path(linewidth=0.2,
            mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI, color=Cluster_ID)) +
  scale_color_manual(labels = c("noise", "cluster_0", "cluster_1", "cluster_4"), 
                    values = c("#e41a1c","#377eb8","#4daf4a","#984ea3")) +
  theme(text=element_text(size=7), legend.position = "bottom") +
  labs(title = "DBSCAN result") + xlab(NULL) + ylab(NULL) + xlim(125.3,126.8) + ylim(36.3, 37.55)

ggsave("./img/dbscan_results_all.png", dbscan_results_all, height = 4, width = 5, dpi = 320) 

#######################################
##Bringing Kmedoid result from python##
#######################################

kmedoid <- read.csv("/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization/k_medoid_result.csv")

kmedoid <-kmedoid %>% mutate(MMSI = paste0(X.1,"_", Cluster_ID)) %>% 
  rename(LONGITUDE=X, LATITUDE=Y, TIMESTAMP=Index) %>% mutate(Cluster_ID = as.factor(Cluster_ID)) 

names(dbscan)

count <- kmedoid %>% group_by(Cluster_ID) %>% count()


point_line(kmedoid, title="K_medoid", linewidth=0.16)

cluster_13 <- kmedoid %>% filter(Cluster_ID==13)
k_first <- point_line(cluster_13, title = "Cluster_13",linewidth = 0.3, col="#e41a1c")

cluster_7 <- kmedoid %>% filter(Cluster_ID==7)
k_second <- point_line(cluster_7, title = "Cluster_7",linewidth = 0.3, col= "#377eb8")

cluster_3 <- kmedoid %>% filter(Cluster_ID==3)
k_third <- point_line(cluster_3, title = "Cluster_3",linewidth = 0.3, col="#4daf4a")

cluster_6 <- kmedoid %>% filter(Cluster_ID==6)
k_fourth <- point_line(cluster_6, title = "Cluster_6",linewidth = 0.3, col="#984ea3")

kmedoids_resuls_sep <-patchwork(k_first, k_second, k_third, k_fourth, ncol=2)
ggsave("./img/kmedoids_results_sep.png", kmedoids_resuls_sep, height = 5, width = 5, dpi = 320)

#Creating a combined plot that shows four clusters together. 
kmedoids_subset <- kmedoid %>% filter(Cluster_ID == 13 | Cluster_ID == 7 | Cluster_ID == 3 | Cluster_ID == 6)

kmedoids_results_all <-
  ggplot(data=kmedoids_subset) + geom_sf(data=south) + geom_sf(data=north) +
  geom_path(linewidth=0.2,
            mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI, color=Cluster_ID)) +
  scale_color_manual(labels = c("cluster_3", "cluster_6", "cluster_7", "cluster_13"), 
                     values = c("#4daf4a","#984ea3","#377eb8","#e41a1c")) +
  theme(text=element_text(size=7), legend.position = "bottom") +
  labs(title = "Kmedoids result") + xlab(NULL) + ylab(NULL) + xlim(125.3,126.8) + ylim(36.3, 37.55)

ggsave("./img/kmedoids_results_all.png", kmedoids_results_all, height = 4, width = 5, dpi = 320) 


#######################################
##Bringing HDBSCAN result from python##
#######################################

hdbscan <- read.csv("/Users/dongheekoh/Documents/Data Science Training/portfolio/projects/AIS_visualization/hdbscan_result.csv")

hdbscan <-hdbscan %>% mutate(MMSI = paste0(X.1,"_", Cluster_ID)) %>% 
  rename(LONGITUDE=X, LATITUDE=Y, TIMESTAMP=Index) %>% mutate(Cluster_ID = as.factor(Cluster_ID)) 

count <- hdbscan %>% group_by(Cluster_ID) %>% count()

point_line(hdbscan, title="HDBSCAN", linewidth=0.16)

cluster_9 <- hdbscan %>% filter(Cluster_ID==9)
hdb_first <- point_line(cluster_9, title = "Cluster_9", linewidth = 0.3,col="#e41a1c")

cluster_1 <- hdbscan %>% filter(Cluster_ID==-1)
hdb_second <- point_line(cluster_1, title = "Noise", linewidth = 0.3, col= "#377eb8")

cluster_6 <- hdbscan %>% filter(Cluster_ID==6)
hdb_third <- point_line(cluster_6, title = "Cluster_6", linewidth = 0.3,col="#4daf4a")

cluster_8 <- hdbscan %>% filter(Cluster_ID==8)
hdb_fourth <- point_line(cluster_8, title = "Cluster_8", linewidth = 0.3,col="#984ea3")

hdbscan_resuls_sep <- patchwork(hdb_first, hdb_second, hdb_third, hdb_fourth, ncol=2)

ggsave("./img/hdbscan_resuls_sep.png", hdbscan_resuls_sep, height = 5, width = 5, dpi = 320)

#Creating a combined plot that shows four clusters together. 
hdbscan_subset <- hdbscan %>% filter(Cluster_ID == 9 | Cluster_ID == -1 | Cluster_ID == 6 | Cluster_ID == 8)

hdbscan_results_all <-
  ggplot(data=hdbscan_subset) + geom_sf(data=south) + geom_sf(data=north) +
  geom_path(linewidth=0.2,
            mapping=aes(x=LONGITUDE, y=LATITUDE, group=MMSI, color=Cluster_ID)) +
  scale_color_manual(labels = c("noise", "cluster_6", "cluster_8", "cluster_9"), 
                     values = c("#377eb8","#4daf4a","#984ea3","#e41a1c")) +
  theme(text=element_text(size=7), legend.position = "bottom") +
  labs(title = "HDBSCAN result") + xlab(NULL) + ylab(NULL) + xlim(125.3,126.8) + ylim(36.3, 37.55)

ggsave("./img/hdbscan_results_all.png", hdbscan_results_all, height = 4, width = 5, dpi = 320) 


###############
### The End ### The project analysis has been completed on Feb/20/2024! 
###############

# This has been a long journey, and it definitely wasn't the easy one. Coding was the fun yet frustrating at times. 
# A lot of trials and errors has not only enhanced my programming skills, but it also enhanced the quality of this project.
# Most of all, the Lord has given me wisdom. When I was in a quandary not knowing where to take my next steps, I felt frustrated.
# At one point, it felt like this whole project was completely failing.
# But whenever I was in such a situation, I prayed. I prayed that the Lord would give me wisdom and perseverance. 
# And the Lord every time I asked of Him -  this is neither hyperbole nor my imagination - faithfully answered all of my prayers.
# I knocked and it was given. I sought His Kingdom and righteousness first, and all these other things have been granted to me. 
# Therefore, I give all the glory to the Lord. It is not I who has completed this project but the Lord!!
# Obviously, compared to other big projects, this project may not seem much, but it does not matter. 
# What matters before the Lord is my motivation, devotion, sincerity, and hard works. 
# Lord, please continue to give me wisdom and opportunities so that I would glorify Your name with my works all the days of my life.
# In Jesus Name! Amen!





# as ofFeb/5th/2024
# creating animations like "spire" (done)
# Morans' I spatial autocorrelation (done)
# flowchart (done)
# Trajectory clustering (Done)




### Reference
#In the following link, https://bootswatch.com/journal/, I can find style source code
#For color choice, https://mycolor.space/?hex=%23072B6F&sub=1
#For Rmarkdown math equations https://rpruim.github.io/s341/S19/from-class/MathinRmd.html
#Tobler, Waldo R. 1970. “A Computer Movie Simulating Urban Growth in the Detroit Region.” Economic Geography 46 (2): 234???40. http://www.geog.ucsb.edu/~tobler/publications/pdf_docs/A-Computer-Movie.pdf.




#2024 salary. Be thankful to the Lord!! He provides!!
#67961544 + (67961544/12)


#(as a conclusion of this section)
#What is the purpose of trajectory clustering analysis? What are we trying to accomplish by this?
# - capturing motion characteristics of moving objects is an interest of many areas including navigation, surveillance, tracking and anomaly detection  
#- In maritime traffic applications, the clustering results may show the usual route and traffic volume distribution and regularity of marine environmental change
#- As a commonly used data mining method, ship trajectory clustering integrates the trajectory data of different ships into different categories or clusters.; 
#it is beneficial for the maritime traffic management stakeholders, such as Maritime Safety Administration (MSA), to obtain insights on the operation status 
#and characteristics of the regional traffic. In the meantime, ship trajectory clustering is one of the fundamental methods for trajectory prediction,
#anomaly detetion, and avoiding ship sollision, which draws much attention from academia






















