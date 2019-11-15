library("lubridate")
library("pracma")
source('~/GitHub/nsgaii/Rnsgaii/nsgaii_io.R')

fl1 <- "d:\\giorgk\\Documents\\GitHub\\C2VsimCG\\OptimResults\\maxWTminArea\\paretoSolutions_44791.dat"
ps <- nsgaii.readParetoSolution(fl1)


# Read the element ids that were used in the optimization
nDivPoints <- scan(file = "../RunC2Vsim/divElem.dat", skip = 0, n = 1)
divelem <- c()
for (i in 1:nDivPoints) {
  temp <- scan(file = "../RunC2Vsim/divElem.dat", skip = i, nlines = 1)
  divelem = c(divelem, temp[-1:-2])
}

## write the Pareto Solutions in one variable with the following fields:
# id, x (objective 1), y (objective 2), idLand (an array with the polygon ids)
## write Pareto solution as javascript variable
js_file <- "../js_scripts/paretoSolutions_44791.js"
cat("var paretoPoints = [\n", file = js_file)
for (i in 1:dim(ps[[2]])[1]) {
  cat("\t{ id: ", i, ", x: ", ps[[2]][i,2], ", y: ", -ps[[2]][i,1], sep = "",  file = js_file, append = TRUE)
  cat(", idLand: [", sep = "", file = js_file, append = TRUE)
  bitid <- which(ps[[1]][i,]==1)
  if (length(which(ps[[1]][i,]==1)) != 0){
    for (j in 1:length(bitid)) {
      el_id <- divelem[bitid[j]]
      # We have to find in which row this element id is written in the all_candidate_polys variable
      # This is show in tempv calculated from a snippet bellow
      
      id <- which(tempv == el_id)
      cat(id, sep = "", file = js_file, append = TRUE)
      if (j !=length(bitid))
        cat(",", sep = "", file = js_file, append = TRUE)
    }
  }
  cat("]", sep = "", file = js_file, append = TRUE)
  cat("}", sep = "", file = js_file, append = TRUE)
  if (i!=dim(ps[[2]])[1])
    cat(",\n", sep = "", file = js_file, append = TRUE)
  else
    cat("\n", sep = "", file = js_file, append = TRUE)
}
cat("];", file = js_file, append = TRUE)

## Read the Diversion Time Series file
fileDTS <- "d:\\giorgk\\Documents\\GitHub\\C2VsimCG\\RunC2Vsim\\divTimeSeries.dat"
info <- scan(file = fileDTS, skip = 0, n = 2)
dts.data <- matrix(nrow = info[2], ncol = info[1])
dts.id <- vector(mode = "numeric",  length = info[1])

for (i in 1:info[1]) {
  temp <- scan(file = fileDTS, skip = i, n = info[2]+1)
  dts.id[i] <- temp[1]
  dts.data[,i] <- temp[-1]
}

## hand write the diversion coordinates. Normally we would get this information from the river file
dts.points <- matrix(nrow = 3, ncol = 2)
dts.points[1,] <- c(35.440852, -118.933162)
dts.points[2,] <- c(36.357291, -119.121089)
dts.points[3,] <- c(36.784392, -119.414387)

# Write diversion as javascript variable
js_DTS_file <- "../js_scripts/diversionTimeSeries.js"
cat("var divNodes = [\n", file = js_DTS_file)

sdateId = 528
nmonths = 528
for (i in 1:length(dts.id)) {
  cat("\t{ id: ", dts.id[i], ", point: [", dts.points[i,1], ",", dts.points[i,2], "],\n\t Values: [\n\t\t" , sep = "",  file = js_DTS_file, append = TRUE)
  sy <- 1965 # This is the year that the diversion started
  sm <- 10 # This is the month
  
  cumvalues <- cumsum(dts.data[sdateId:dim(dts.data)[1],i])
  for (j in 1:length(cumvalues)) {
    ndays <- as.integer(days_in_month(as.Date(paste0(sy,"-",sm,"-1"))))
    
    temp <- paste0("{x: new Date(", sy, ",", sm,",",ndays,")", ", y: ",  cumvalues[j],"}")
    cat(temp, sep = "",  file = js_DTS_file, append = TRUE)
    
    if (j !=length(cumvalues))
      cat(",\n\t\t", sep = "",  file = js_DTS_file, append = TRUE)
    
    sm = sm + 1
    if (sm >= 13){
      sm = 1
      sy = sy + 1
    }
  }
  
  
  cat("]}", file = js_DTS_file, append = TRUE)
  if (i != length(dts.id))
    cat(",\n", sep = "", file = js_DTS_file, append = TRUE)
  else
    cat("\n", sep = "", file = js_DTS_file, append = TRUE)
}
cat("];", file = js_DTS_file, append = TRUE)







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
tempv <- vector(mode = "numeric", length = 0)
cat("var all_candidate_polys = [\n", file = "../js_scripts/all_candidate_polys.js")
for (i in 1:length(diversion_polygons$IE)) {
  if (is.na(diversion_polygons$DivNode[i])){
    next
  }
  tempv <- c(tempv,diversion_polygons$IE[i])
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





