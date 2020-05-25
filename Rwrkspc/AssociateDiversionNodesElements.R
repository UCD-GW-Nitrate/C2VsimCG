library(rgdal)
library(stringr)
library(gwtools)
# Read the diversion specification data
divSpec <- gwtools::c2vsim.readDivSpec(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdivspec.dat")

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

# Make a unique list of diversions that are accosiated with river nodes ---------
diversionRiverNodes <- unique(divSpec$RDV[,2])
diversionRiverNodes <- diversionRiverNodes[-which(diversionRiverNodes == 0)]
uDivCoords <- matrix(data = NA, nrow = length(diversionRiverNodes), ncol = 2)
#uDivNames <- vector(mode = "list", length = length(diversionRiverNodes))
uDivNames <- c()
uDivElem <- vector(mode = "list", length = length(diversionRiverNodes))

for (i in 1:length(diversionRiverNodes)) {
  # assign the coordinates of the river node
  ind <- which(rivernodes_4326$IRV == diversionRiverNodes[i])
  uDivCoords[i,] <- as.numeric(slot(rivernodes_4326, "coords")[ind, c(2,1)])
  
  # find how many diversion exists from this river node
  ii <- which(divSpec$RDV[,2] == diversionRiverNodes[i])
  # Find a common name for those diversions
  divname <- divSpec$RDVnames[ii[1]]
  elemIds <- divSpec$RDVELEM[ii[1]][[1]][,1]
  if (length(ii)>1){
    for (j in 2:length(ii)) {
      elemIds <- c(elemIds, divSpec$RDVELEM[ii[j]][[1]][,1])
      divname <- paste0(divname, "</br>", divSpec$RDVnames[ii[j]])
      #divname <- extractCommonRoot(divname, divSpec[[7]][ii[j]])
    }
  }
  #uDivNames[[i]] <- divname
  uDivNames <- c(uDivNames, divname)
  uDivElem[[i]] <- unique(elemIds)
}

# Append the diversions that have no river node associated --------------------------
# List find the elements that receive diversion from the 0 river nodes
canalElem <- vector(mode = "list", length = length(canalNames))
for (i in 1:dim(divSpec$RDV)[1]){
  if (divSpec$RDV[i,2] == 0){
    not_found <- T
    for (j in 1:length(canalNames)){
      nstr <- nchar(canalNames[j])
      if (substr(divSpec$RDVnames[i], 1, nstr) == canalNames[j]){
        canalElem[[j]] <- unique(c(canalElem[[j]], divSpec$RDVELEM[[i]][,1]))
        not_found <- F
        break
      }
    }
    if(not_found)
      print(i)
  }
}

AllNames <- c(uDivNames,canalNames)
AllCoords <- rbind(uDivCoords, canalNodes)
AllElem <- c(uDivElem, canalElem)
type <- c(rep(1,length(uDivNames)), rep(2,length(canalNames)))

save("AllNames", "AllCoords", "AllElem", "type",file = "AssociateDiversionNodes.RData")

# Write diversion nodes as js file ---------------------------------------------
filename = "../js_scripts/C2vsimDiversions.js"
cat("var C2VsimDivs = [\n", file = filename)
for (i in 1:length(AllNames)) {
  xy <- AllCoords[i,]
  
  tmp <- paste0("\t{ point: [", xy[1], ", ", xy[2], "], type: ", type[i], ", name: '", gsub("'", "", AllNames[i]), "', elids: [")
  
  tmp <- paste0(tmp, toString (AllElem[[i]]), "]},\n")
  cat(tmp , sep = "",  file = filename, append = TRUE)
}
cat("];", file = filename, append = TRUE)


# Write all mesh polygons -------------------------------------------------
c2vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")
c2vsim_mesh_4326 <-spTransform(c2vsim_mesh, CRS("+init=epsg:4326"))

filename = "../js_scripts/c2vsimMesh.js"
cat("var c2vsimMesh = [\n", file = filename)

for (i in 1:length(c2vsim_mesh_4326)) {
  clr_elem <- color_elem_price(elem_price[i]/1000000, 2.5, 1000)

  coords <- shapefile.coords(c2vsim_mesh_4326, i)
  coords <- coords[-dim(coords)[1],]

  tmp <- paste0("\t{ id:", c2vsim_mesh_4326$IE[i], ",color:'", clr_elem,  "', cost:", elem_price[i]/1000000,  ", Polygon: [")
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
  clr_elem <- color_elem_price(elem_price[div_polys[i]]/1000000, 2.5, 1000)
  
  coords <- shapefile.coords(c2vsim_mesh_4326, div_polys[i])
  coords <- coords[-dim(coords)[1],]
  tmp <- paste0("\t{ id:", c2vsim_mesh_4326$IE[div_polys[i]], ",color:'", clr_elem,  "', cost: ", elem_price[div_polys[i]]/1000000,  ", Polygon: [")
  for (j in 1:dim(coords)[1]) {
    tmp <- paste0(tmp, "[", coords[j,2], ", ", coords[j,1], "]")
    if (j < dim(coords)[1])
      tmp <- paste0(tmp, ",")
  }
  tmp <- paste(tmp, "]},\n")
  cat(tmp, file = filename, append = TRUE)
  
}
cat("];", file = filename, append = TRUE)



