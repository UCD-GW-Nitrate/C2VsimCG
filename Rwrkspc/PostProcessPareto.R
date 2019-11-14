## write Pareto solution as javascript variable
cat("var paretoPoints = [\n", file = "../js_scripts/ParetoSolution.js")
for (i in 1:dim(pS[[2]])[1]) {
  ln <- 
    #{ x: 1.00000, y: 5.07249, sy: new Date(1993,10,1), ey: new Date(1994,9,1) },
  cat("\t{ x: ", pS[[2]][i,2], ", y: ", -pS[[2]][i,1], ", id: ", i,  "},\n" , file = "../js_scripts/ParetoSolution.js", append = TRUE)
  
}

cat("];", file = "../js_scripts/ParetoSolution.js", append = TRUE)

## Pareto desicion variables
cat("var paretoVariables = [\n", file = "../js_scripts/ParetoVariables.js")
for (i in 1:dim(pS[[2]])[1]){
  ids <- which(pS[[1]][i,]==1)
  cat("\t[", ids[1], ", ", file = "../js_scripts/ParetoVariables.js", append = TRUE)
  cat(ids[-1],"],\n", sep = ",", file = "../js_scripts/ParetoVariables.js", append = TRUE)
}
cat("];", file = "../js_scripts/ParetoVariables.js", append = TRUE)

### write the polygon coordinates
## First we have to read few things
# The element ids that were used in the optimization
nDivPoints <- scan(file = "../RunC2Vsim/divElem.dat", skip = 0, n = 1)
divelem <- c()
for (i in 1:nDivPoints) {
  temp <- scan(file = "../RunC2Vsim/divElem.dat", skip = i, nlines = 1)
  divelem = c(divelem, temp[-1:-2])
}

# Read the element shapefile and convert it to WGS84 i.e. 4326
cvmesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")
cvmesh4326 <- spTransform(cvmesh, CRS("+init=epsg:4326"))
cvmesh_polygons <- slot(cvmesh4326,"polygons")
cvmesh_area <- slot(cvmesh,"polygons")


cat("var polygons = [\n", file = "../js_scripts/polygons.js")
cat("var polyarea = [\n", file = "../js_scripts/polyarea.js")
for (i in 1:length(divelem)){
  polycoords <- slot(slot(cvmesh_polygons[[divelem[i]]], "Polygons")[[1]],"coords")
  temp <- "["
  for (k in 1:dim(polycoords)[1]){
    temp <- paste0(temp,"[", polycoords[k,2],",", polycoords[k,1], "]")
    if (k != dim(polycoords)[1]){
      temp <- paste0(temp,",")
    }
  }
  temp <- paste0(temp,"],")
  cat("\t", temp, "\n", file = "../js_scripts/polygons.js", append = TRUE)
}
cat("];", file = "../js_scripts/polygons.js", append = TRUE)
cat("];", file = "../js_scripts/polyarea.js", append = TRUE)

# Read the shapefile with the diversion polygons
diversion_polygons <- readOGR(dsn = "../gis_data/diversionElements.shp")
diversion_polygons4326 <- spTransform(diversion_polygons, CRS("+init=epsg:4326"))

# write them as javascript variable
cat("var all_candidate_polys = [\n", file = "../js_scripts/all_candidate_polys.js")
for (i in 1:length(diversion_polygons$IE)) {
  if (is.na(diversion_polygons$DivNode[i])){
    next
  }
  poly <- slot(slot(slot(diversion_polygons4326, "polygons")[[i]], "Polygons")[[1]], "coords")
  cat("\t{ id: ", diversion_polygons$IE[i], ", DivND: ", diversion_polygons$DivNode[i], file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
  cat(", Polygon: [", file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
  for (j in 1:(dim(poly)[1]-1)) {
    cat("[", poly[j,2], ",", poly[j,1], "]", file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
    if (j != dim(poly)[1]-1)
      cat(",", file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
  }
  if (i != length(diversion_polygons$IE))
    cat("]", "},\n", file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
  else
    cat("]", "}\n", file = "../js_scripts/all_candidate_polys.js", sep = "", append = TRUE)
}
cat("];", file = "../js_scripts/all_candidate_polys.js", append = TRUE)





