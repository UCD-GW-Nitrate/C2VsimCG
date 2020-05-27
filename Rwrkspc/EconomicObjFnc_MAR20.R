library("readxl")
library("rgdal")

# The file with the prices
IE_price <- read_xlsx('cv_land_price_022620.xlsx')

c2vsim_mesh <- readOGR(dsn = "../gis_data", layer = "C2Vsim_mesh")

c2vsim_mesh@data$price <- NA
c2vsim_mesh$price[IE_price$ie] <- IE_price$pred_price_county_fe

# Write prices to js file
c2vsim_mesh_4326 <-spTransform(c2vsim_mesh, CRS("+init=epsg:4326"))

filename = "../js_scripts/C2VsimElemPrice_MARCH20_mod1.js"
cat("var c2vsimPrices = [\n", file = filename)

min_price <- min(IE_price$pred_price_county_fe, na.rm = T)
max_price <- max(IE_price$pred_price_county_fe, na.rm = T)
{# min max based on modified prices 
  min_price <- min(c2vsim_mesh$price, na.rm = T)
  max_price <- max(c2vsim_mesh$price, na.rm = T)
}

for (i in 1:length(c2vsim_mesh_4326)) {
  clr_elem <- color_elem_price(c2vsim_mesh_4326$price[i], min_price, max_price)
  
  coords <- c2vsim_mesh_4326@polygons[[i]]@Polygons[[1]]@coords
  coords <- coords[-dim(coords)[1],]
  
  # price <- c2vsim_mesh_4326$price[i]
  # use the modified instead 
  price <- c2vsim_mesh$price[i]
  if (is.na(price)){
    price = 0
  }
  
  tmp <- paste0("\t{ id:", c2vsim_mesh_4326$IE[i], ",color:'", clr_elem,  "', cost:", price,  ", Polygon: [")
  for (j in 1:dim(coords)[1]) {
    tmp <- paste0(tmp, "[", coords[j,2], ", ", coords[j,1], "]")
    if (j < dim(coords)[1])
      tmp <- paste0(tmp, ",")
  }
  tmp <- paste(tmp, "]},\n")
  cat(tmp, file = filename, append = TRUE)
  
}
cat("];", file = filename, append = TRUE)


###### Interpolate missing price values #####
# Isolate the barycenters and the area of the elements that have price
cc_price <- matrix(data = NA, nrow = 0, ncol = 5) # COLUMNS: IE Barycenters X, Y, area, Price
for (i in 1:length(c2vsim_mesh)) {
  if (!is.na(c2vsim_mesh$price[i])){
    n <- dim(c2vsim_mesh@polygons[[i]]@Polygons[[1]]@coords)[1]-1
    cc_price <- rbind(cc_price,
    c(i,apply(c2vsim_mesh@polygons[[i]]@Polygons[[1]]@coords[1:n,],2,mean), 
      c2vsim_mesh@polygons[[i]]@Polygons[[1]]@area, c2vsim_mesh$price[i]))
  }
}

for (i in 1:length(c2vsim_mesh)) {
  if (!is.na(c2vsim_mesh$price[i])){
    next
  }
  n <- dim(c2vsim_mesh@polygons[[i]]@Polygons[[1]]@coords)[1]-1
  cc <- apply(c2vsim_mesh@polygons[[i]]@Polygons[[1]]@coords[1:n,],2,mean)
  dst <- sqrt((cc_price[,2] - cc[1])^2 + (cc_price[,3] - cc[2])^2 )
  id_sort <- sort.int(dst, index.return = T)
  
  # Select the 10 closest elements to contribute to interpolation
  d <- id_sort$x[1:10]
  w <- 1/(d^4)
  c2vsim_mesh$price[i] <- sum(w*cc_price[id_sort$ix[1:10],5])/sum(w)
}

{# Calculate element area
  area <- vector(mode = "numeric", length = length(c2vsim_mesh))
  for (i in 1:length(c2vsim_mesh)) {
    area[i] <- c2vsim_mesh@polygons[[i]]@area
  }
}

{# Print economic Objective function
  elemInfoFile <- "../OptimResults/inputFiles/ElemCost_MAY20_mod.dat"
  cat(length(c2vsim_mesh), "\n", sep = "",  file = elemInfoFile)
  write(x = t(cbind(c2vsim_mesh$IE, c2vsim_mesh$price/1000, area/(1000*1000))), file = elemInfoFile, ncolumns = 3, sep = " ", append = TRUE)
  
}
