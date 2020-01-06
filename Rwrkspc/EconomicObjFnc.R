library(rgdal)
#source("c2vsim_io.R")
landPriceTable <- read.csv(file = "ie_c2vsim_landuse_saleprice.csv")
#LU <- c2vsim.read.LandUse(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVlanduse.dat", 
#                          colNames = c("IE", "AG", "UR", "NV", "RV"))
m <- rgdal::readOGR("../gis_data/LandUseNASS_2018.geojson")
