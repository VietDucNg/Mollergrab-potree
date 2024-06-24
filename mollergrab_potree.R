# Author: Viet Nguyen, faculty of Forest and Environment, HNEE
# Email: Duc.Nguyen@hnee.de

# setup workspace
getwd()

# install packages
install.packages('assertthat')
install.packages('dbscan')
install.packages('E:\\OneDrive_HNEE\\01_study\\03_Master_FIT\\03_semester_3\\LiDAR\\Material_LidarCourse_2022\\Software\\crownsegmentr-0.4.4.tar.gz')

# load packages
library(lidR)
library(sf)
library(crownsegmentr)

# Parallel computation in lidR
get_lidr_threads()
set_lidr_threads(0.8) # recommended to use 80% of total threads
get_lidr_threads()

# import las
las = readLAS('data/pc_Mollergrab.laz')

# load aoi shapefile
aoi.sf = st_read("data/aoi_Mollergrab.shp")


##############################
#### clip las file to aoi ####
##############################

# add buffer of 5m to aoi.sf to make sure entire trees included
aoi_buff = st_buffer(aoi.sf, dist = 5)
plot(aoi.sf[0])
plot(aoi_buff[0], add=T)

# transform aoi CRS to las CRS
st_crs(aoi_buff)
lidR::projection(las)
lidR::projection(las) = 25833

# clip point cloud by polygon
las = clip_roi(las,aoi_buff)


##################
#### thinning ####
##################

# check point density
density = grid_density(las, res = 1) #density in 1 square meter
plot(density)
density

# thin point cloud to 100points/m2
las = decimate_points(las, homogenize(500,res = 1))

# check point density after thinning
density = grid_density(las, res = 1) #density in 1 square meter
plot(density)
density


###############################
#### ground classification ####
###############################

# Classify ground points with the cloth simulation filter (CSF) algorithm
# Drop a simulated cloth on the inverted point cloud,
# Ground points classified by analysing
# the interactions between cloth nodes and the inverted surface.
las = classify_ground(las, algorithm=csf())
table(las$Classification)


#######################
#### normalization ####
#######################

# point cloud-based approach
# Inverse Distance Weighting (IDW) interpolation
las = normalize_height(las, knnidw())

hist(filter_ground(las)$Z, breaks = seq(-10, 10, 0.05), 
     main = "point cloud-based approach", xlab = "Z", xlim=c(-1, 1))


######################
#### segmentation ####
######################

# remove ground points may help to reduce significantly processing time
las = filter_poi(las, Classification != 2)

# AMS3D
las.ams3d = segment_tree_crowns(las, 
                                crown_diameter_to_tree_height = 0.3,
                                crown_height_to_tree_height = 0.4,
                                segment_crowns_only_above = 1)

# count tree detected
length(unique(las.ams3d@data[["crown_id"]]))

# Remove points which are not classified as tree
las.ams3d = filter_poi(las.ams3d, !is.na(crown_id))

# count tree detected
length(unique(las.ams3d@data[["crown_id"]]))

# visualization
plot(las.ams3d, size=2, color="crown_id")

# add treeid value to pointsource id
las.ams3d@data$PointSourceID = las.ams3d@data$crown_id

#### save preprocessed point cloud
writeLAS(las.ams3d,file = "data/pc_Mollergrab_segment.laz")
