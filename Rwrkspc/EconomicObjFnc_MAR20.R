library("readxl")
library("rgdal")

# The file with the prices
IE_price <- read_xlsx('cv_land_price_022620.xlsx')

c2vsim_mesh <- readOGR(dsn = "../gis_data", layer = "C2Vsim_mesh")

c2vsim_mesh@data$price <- NA
c2vsim_mesh$price[IE_price$ie] <- IE_price$pred_price_county_fe

# Write prices to js file
c2vsim_mesh_4326 <-spTransform(c2vsim_mesh, CRS("+init=epsg:4326"))

filename = "../js_scripts/C2VsimElemPrice_MARCH20.js"
cat("var c2vsimPrices = [\n", file = filename)

min_price <- min(IE_price$pred_price_county_fe, na.rm = T)
max_price <- max(IE_price$pred_price_county_fe, na.rm = T)

for (i in 1:length(c2vsim_mesh_4326)) {
  clr_elem <- color_elem_price(c2vsim_mesh_4326$price[i], min_price, max_price)
  
  coords <- c2vsim_mesh_4326@polygons[[i]]@Polygons[[1]]@coords
  coords <- coords[-dim(coords)[1],]
  
  price <- c2vsim_mesh_4326$price[i]
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

