library("rgdal")
source("c2vsim_io.R")
# Read the diversion specification data
divSpec <- c2vsim.readDivSpec(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdivspec.dat")

# Read the river nodes
rivernodes <- readOGR(dsn = "../gis_data/C2Vsim_riverNodes.shp")
rivernodes_4326 <-spTransform(rivernodes, CRS("+init=epsg:4326"))  



canalNames <- c()
canalNodes <- matrix(data = NA, nrow = 1, ncol = 2)

# Cross-Valley Canal
# https://www.anyplaceamerica.com/directory/ca/kern-county-06029/canals/cross-valley-canal-2680156/
canalNames <- c(canalNames, "Cross-Valley Canal")
canalNodes[1,] <- c(35.3466667, -119.1872222)

# Friant-Kern Canal
canalNames <- c(canalNames, "Friant-Kern Canal")
canalNodes <- rbind(canalNodes, c(36.998, -119.703))

# Madera Canal
canalNames <- c(canalNames, "Madera Canal")
canalNodes <- rbind(canalNodes, c(37.003, -119.707))

# San Luis Canal
canalNames <- c(canalNames, "San Luis Canal")
canalNodes <- rbind(canalNodes, c(36.386, -119.878))

# O'Neill Forebay
canalNames <- c(canalNames, "O'Neill Forebay")
canalNodes <- rbind(canalNodes, c(37.082, -121.044))

# Mendota Pool
canalNames <- c(canalNames, "Mendota Pool")
canalNodes <- rbind(canalNodes, c(36.784, -120.371))

# Delta-Mendota Canal
canalNames <- c(canalNames, "Delta-Mendota Canal")
canalNodes <- rbind(canalNodes, c(37.817, -121.556))

# Merced R
canalNames <- c(canalNames, "Merced R")
canalNodes <- rbind(canalNodes, c(37.525, -120.321))

# Turlock Canal
canalNames <- c(canalNames, "Turlock Canal")
canalNodes <- rbind(canalNodes, c(37.609, -120.766))

# Modesto Canal
canalNames <- c(canalNames, "Modesto Canal")
canalNodes <- rbind(canalNodes, c(37.656, -120.675))

# Stanislaus R
canalNames <- c(canalNames, "Stanislaus R")
canalNodes <- rbind(canalNodes, c(37.864, -120.629))

# Mokelumne R
canalNames <- c(canalNames, "Mokelumne R")
canalNodes <- rbind(canalNodes, c(38.225, -121.021))

# Folsom Lake
canalNames <- c(canalNames, "Folsom Lake")
canalNodes <- rbind(canalNodes, c(38.714, -121.143))

# Cache Creek
canalNames <- c(canalNames, "Cache Creek")
canalNodes <- rbind(canalNodes, c(38.712, -122.087))

# Cross Canal
canalNames <- c(canalNames, "Cross Canal")
canalNodes <- rbind(canalNodes, c(38.821, -121.545))

# Combie (Gold Hill) Canal
canalNames <- c(canalNames, "Combie (Gold Hill) Canal")
canalNodes <- rbind(canalNodes, c(38.987, -121.099))

# Boardman Canal
canalNames <- c(canalNames, "Boardman Canal")
canalNodes <- rbind(canalNodes, c(38.803, -121.164))

# Bear River Canal
canalNames <- c(canalNames, "Bear River Canal")
canalNodes <- rbind(canalNodes, c(39.036, -121.34))

# Feather River
canalNames <- c(canalNames, "Feather River")
canalNodes <- rbind(canalNodes, c(39.53, -121.545))

# Thermalito Afterbay
canalNames <- c(canalNames, "Thermalito Afterbay")
canalNodes <- rbind(canalNodes, c(39.459, -121.679))

# Bangor Canal
canalNames <- c(canalNames, "Bangor Canal")
canalNodes <- rbind(canalNodes, c(39.377, -121.419))

# Little Dry Creek
canalNames <- c(canalNames, "Little Dry Creek")
canalNodes <- rbind(canalNodes, c(39.369, -121.874))

# Oroville-Wyandotte ID
canalNames <- c(canalNames, "Oroville-Wyandotte ID")
canalNodes <- rbind(canalNodes, c(39.574, -121.474))

# Palermo Canal
canalNames <- c(canalNames, "Palermo Canal")
canalNodes <- rbind(canalNodes, c(39.499, -121.49))

# Miocine and Wilenor Canals
canalNames <- c(canalNames, "Miocine and Wilenor Canals")
canalNodes <- rbind(canalNodes, c(39.603, -121.481))

# Tarr Ditch
canalNames <- c(canalNames, "Tarr Ditch")
canalNodes <- rbind(canalNodes, c(39.128, -121.303))

# Little Chico Creek
canalNames <- c(canalNames, "Little Chico Creek")
canalNodes <- rbind(canalNodes, c(39.734, -121.794))

# Stony Creek
canalNames <- c(canalNames, "Stony Creek")
canalNodes <- rbind(canalNodes, c(39.807, -122.353))

# Clear Creek
canalNames <- c(canalNames, "Clear Creek")
canalNodes <- rbind(canalNodes, c(40.495, -122.431))

# Cottonwood Creek
canalNames <- c(canalNames, "Cottonwood Creek")
canalNodes <- rbind(canalNodes, c(39.538, -121.679))

# Whiskeytown and Shasta
canalNames <- c(canalNames, "Whiskeytown and Shasta")
canalNodes <- rbind(canalNodes, c(40.631, -122.477))


# Make a list of unique diversions based on name --------------------------
unique_names = c()
for (i in 1:dim(divSpec[[2]])[1]) {
  if (divSpec[[2]][i,2] == 0){
    pr <- TRUE
    for (j in 1:length(canalNames)) {
      nstr <- nchar(canalNames[j])
      if (substr(divSpec[[7]][i], 1, nstr) == canalNames[j]){
        unique_names <- c(unique_names, canalNames[j])
        pr <- FALSE
        break
      }
    }
    if (pr)
      print(i)
  }
  else{
    unique_names <- c(unique_names, divSpec[[7]][i])
  }
}

unique_names <- levels(factor(unique_names))


# Make a list of element for each diversion name --------------------------
uDivCoords <- matrix(data = NA, nrow = length(unique_names), ncol = 2)
uDivElem <- vector(mode = "list", length = length(unique_names))
type <- vector(mode = "numeric", length = length(unique_names))

for (i in 1:dim(divSpec[[2]])[1]) {
  if (divSpec[[2]][i, 2] != 0){
    id <- which(unique_names == divSpec[[7]][i])
    uDivCoords[id,] <- as.numeric(slot(rivernodes_4326, "coords")[divSpec[[2]][i,2], c(2,1)])
    uDivElem[[id]] <- divSpec[[3]][i][[1]][,1]
    type[id] <- 1
  }
  else{
    # find which of the unique name containts part of this
    for (j in 1:length(unique_names)) {
      nstr <- nchar(unique_names[j])
      if (substr(divSpec[[7]][i], 1, nstr) == unique_names[j]){
        id_cn <- which(canalNames == unique_names[j])
        uDivCoords[j, ] <- canalNodes[id_cn,]
        uDivElem[[j]] <- c(uDivElem[[j]], divSpec[[3]][i][[1]][,1])
        type[j] <- 2
        break
      }
    }
  }
}



# Write diversion nodes as js file ---------------------------------------------
filename = "../js_scripts/C2vsimDiversions.js"
cat("var C2VsimDivs = [\n", file = filename)
for (i in 1:length(unique_names)) {
  xy <- uDivCoords[i,]
  
  tmp <- paste0("\t{ point: [", xy[1], ", ", xy[2], "], type: ", type[i], ", name: '", gsub("'", "", unique_names[i]), "', elids: [")
  
  tmp <- paste0(tmp, toString (uDivElem[[i]]), "]},\n")
  cat(tmp , sep = "",  file = filename, append = TRUE)
}
cat("];", file = filename, append = TRUE)


# Write all mesh polygons -------------------------------------------------
c2vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")
c2vsim_mesh_4326 <-spTransform(c2vsim_mesh, CRS("+init=epsg:4326"))

filename = "../js_scripts/c2vsimMesh.js"
cat("var c2vsimMesh = [\n", file = filename)

for (i in 1:length(c2vsim_mesh_4326)) {
  coords <- shapefile.coords(c2vsim_mesh_4326, i)
  coords <- coords[-dim(coords)[1],]
  tmp <- paste0("\t{ id:", c2vsim_mesh_4326$IE[i], ", Polygon: [")
  for (j in 1:dim(coords)[1]) {
    tmp <- paste0(tmp, "[", coords[j,2], ", ", coords[j,1], "]")
    if (j < dim(coords)[1])
      tmp <- paste0(tmp, ",")
  }
  tmp <- paste(tmp, "]},\n")
  cat(tmp, file = filename, append = TRUE)
  
}
cat("];", file = filename, append = TRUE)


# Write only the diversion polygons ---------------------------------------

# make a unique list of polygons that receive diverted water
div_polys = c()
for (i in 1:length(uDivElem)) {
  div_polys <- c(div_polys, uDivElem[[i]])
}

div_polys <- unique(div_polys)

filename = "../js_scripts/polys_with_diversions.js"
cat("var div_polys = [\n", file = filename)

for (i in 1:length(div_polys)) {
  if (div_polys[i] == 0)
    next
  coords <- shapefile.coords(c2vsim_mesh_4326, div_polys[i])
  coords <- coords[-dim(coords)[1],]
  tmp <- paste0("\t{ id:", c2vsim_mesh_4326$IE[div_polys[i]], ", Polygon: [")
  for (j in 1:dim(coords)[1]) {
    tmp <- paste0(tmp, "[", coords[j,2], ", ", coords[j,1], "]")
    if (j < dim(coords)[1])
      tmp <- paste0(tmp, ",")
  }
  tmp <- paste(tmp, "]},\n")
  cat(tmp, file = filename, append = TRUE)
  
}
cat("];", file = filename, append = TRUE)
