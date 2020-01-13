library(pracma)
library("rgdal")
library("readxl")
source("Utilities.R")
landPriceTable <- read.csv(file = "ie_c2vsim_landuse_saleprice.csv")
ca_cdl_xwalk <- read_excel(path = "ca_cdl_xwalk.xlsx", sheet = 1, range = "A1:C102",col_names = T)
# get unique descriptions of land uses
lu_descr <- levels(factor(landPriceTable$davis_class_any))
LUNASS_codes <- read.csv(file = "LUNASS_categories.csv",header = F,col.names = c("id", "Name", "Code"))



#LU <- c2vsim.read.LandUse(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVlanduse.dat", 
#                          colNames = c("IE", "AG", "UR", "NV", "RV"))
m <- rgdal::readOGR("../gis_data/LandUseNASS_2017.geojson")
cv_mesh <- readOGR("../gis_data/C2Vsim_mesh.shp")
cv_mesh_poly <- slot(cv_mesh,"polygons")


LU_elem_crop_percent <- vector(mode = "list", length = length(m))
for (i in 1:length(m)) {
  m_index <- which(m$IE == i)
  LU_elem_crop_percent[[i]] <- readHistogramLine(m$histogram[m_index])
}

# Check if any element has the code 215
for (i in 1:length(m)){
  if(isempty(which(LU_elem_crop_percent[[i]][,1] == 215))){
    next
  }
  else{
    print(i)
  }
}


# make a unique list of codes
lucodelist <- c()
for (i in 1:length(m)){
  lucodelist <- c(lucodelist,LU_elem_crop_percent[[i]][,1])
}
lucodelist <- unique(lucodelist)
# add a bool field indicating whether the LU code exist in CV
cvexist <- vector(mode = "numeric", length = dim(LUNASS_codes)[1])
for (i in 1:length(cvexist)) {
  if (length(which(lucodelist == LUNASS_codes$Code[i]) == 0) != 0){
    cvexist[i] <- 1
  }
}
LUNASS_codes$CVexist  <-  cvexist

# extract the codes that exist in CV
CV_LUNASS_codes <- LUNASS_codes[which(LUNASS_codes$CVexist == 1),]
write.csv(CV_LUNASS_codes,file = "CV_LUNASS_codes.csv",)



# Assign a depth value for each element--------------
# generate CVstrat, headAll, MSH and XY from ExtractDataFromC2Vsim.R script
Depth_per_elem <- vector(mode = "numeric", length = dim(MSH)[1])

# put the head values of the last decade into a matrix
itm <- seq(dim(headAll[[2]])[1]-119,dim(headAll[[2]])[1])
head_tmp <- matrix(data = NA, nrow = dim(XY)[1], ncol = 120)
for (i in 1:length(itm)) {
  head_tmp[,i] <- headAll[[3]][[itm[i]]][,1]
}

head_tmp[which(head_tmp > 15000)] <- NA

# calculate element barycenter
elem_bc <- matrix(data = NA, nrow = dim(MSH)[1], ncol = 2)
for (i in 1:dim(MSH)[1]) {
  # find the node ids for the element i
  ids <- MSH[i,2:5]
  ids <- as.numeric(ids[which(ids != 0)])
  Depth_per_elem[i] <- mean(CVstrat$ELV[ids] - apply(head_tmp[ ids,], 1, mean, na.rm = T), na.rm = T)
  
  elem_bc[i,] <- c(mean(XY$X[ids]),mean(XY$Y[ids]))
}

# Calculate cost for each element ----------
elem_price <- vector(mode = "numeric", length = length(m))
for (iel in 1:length(m)) {
  print(paste("element:", iel))
  crop_el <- LU_elem_crop_percent[[iel]]
  ncrops <- dim(crop_el)[1]
  crop_el <- cbind(crop_el, vector(mode = "numeric", length = ncrops))
  for(icrop in 1:ncrops){
    davis_class <- ca_cdl_xwalk$davis_class[ which(ca_cdl_xwalk$class == crop_el[icrop])]
    if (isempty(davis_class)){
      print(paste("Not found", crop_el[icrop]))
      next
    }
    elem_with_class <- which(landPriceTable$davis_class_any == davis_class)
    if (isempty(elem_with_class)){
      print(paste(davis_class, "not found"))
      next
    }
    elem_with_class <- unique(landPriceTable$ie[elem_with_class])
    # Find the distance between the elements that contain price for that class with the current element
    dst_ml <- sqrt((elem_bc[iel,1] - elem_bc[elem_with_class,1])^2 + (elem_bc[iel,2] - elem_bc[elem_with_class,2])^2)*(1/1609.34)
    
    if (length(which(dst_ml<20)) > 1){
      # If there are more than one close elements then choose the one with the closet depth
      closest_elem <- sort(dst_ml,index.return=T)
      if (closest_elem$x[1] == 0){
        selected_elem <- elem_with_class[closest_elem$ix[1]]
      }
      else{
        close_picked_elem <- elem_with_class[closest_elem$ix[which(closest_elem$x<20)]]
        depth_elem_list <- Depth_per_elem[close_picked_elem]
        depth_current_el <- Depth_per_elem[iel]
        depth_sort <- sort(sqrt((depth_elem_list-depth_current_el)^2), index.return=T)
        selected_elem <- close_picked_elem[depth_sort$ix[1]]
      }
    }
    else{
      # if there is only one or none close elements choose the closest
      closest_elem <- sort(dst_ml,index.return=T)
      selected_elem <- elem_with_class[closest_elem$ix[1]]
      
    }
    iprice <- which(landPriceTable$davis_class_any == davis_class & landPriceTable$ie == selected_elem)
    crop_el[icrop,3] <- landPriceTable$mean_ie_lu_price[iprice]
  }
  
  elem_area <- slot(cv_mesh_poly[[iel]], "area")*(1/4046.86)
  crop_acres <- elem_area*crop_el[,2]/sum(crop_el[,2])
  elem_price[iel] <- sum(crop_acres*crop_el[,3])
}
