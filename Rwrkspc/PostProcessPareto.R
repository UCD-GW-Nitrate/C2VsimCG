library("lubridate")
library("pracma")
library("rgdal")
source('~/GitHub/nsgaii/Rnsgaii/nsgaii_io.R')
source("c2vsim_io.R")

fl1 <- "../OptimResults/maxWTminArea/paretoSolutions_58946.dat"
ps <- nsgaii.readParetoSolution(fl1)

# Read the element ids that were used in the optimization =================
# 
nDivPoints <- scan(file = "../RunC2Vsim/divElem.dat", skip = 0, n = 1)

divelem <- c()
divNodeId <- c()
divelem_perNode <- vector(mode = "list", length = 3)
for (i in 1:nDivPoints) {
  temp <- scan(file = "../RunC2Vsim/divElem.dat", skip = i, nlines = 1)
  divelem <-  c(divelem, temp[-1:-2])
  divNodeId <- c(divNodeId, vector(mode = "numeric", length = length(temp[-1:-2])) + temp[1])
  divelem_perNode[[i]] <- temp[-1:-2]
}


# Read the shapefile that contain the diversion polygons ============
# 
diversion_polygons <- readOGR(dsn = "../gis_data/diversionElements.shp")
diversion_polygons4326 <- spTransform(diversion_polygons, CRS("+init=epsg:4326"))


# Write them as javascript variable================
# 
# # This containts the element ids in the order they are printed in the following js file
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




# write the Pareto Solutions as one JS variable ===============
#  with the following fields:
# id, x (objective 1), y (objective 2), idLand (an array with the polygon ids)
## write Pareto solution as javascript variable
js_file <- "../js_scripts/paretoSolutions_58946.js"
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


# Read the Diversion Time Series file ===========================
# 
fileDTS <- "../RunC2Vsim/divTimeSeries.dat"
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



# Write diversion Timeseries as javascript variable ========================
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




# Execute the simulation for the pareto solution=================================
# to read  the groundwater storage change and the steam return flow First run the first
# code snippet to load the pareto solutions

# Run the base simulation so that we can compare the pareto solution
# Read the base simulation diversion data, print, run and read

# Read ...
divspec <- c2vsim.readDivSpec("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdivspec.dat")
divdata <- c2vsim.readDivData("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdiversions.dat")
# Write...
c2vsim.writeDivSpec("../RunC2Vsim/tempRspec.dat", DivSpec = divspec)
c2vsim.writeDivData("../RunC2Vsim/tempRdata.dat", data = divdata)

# run C2Vsim (after you make sure that the diversion files in the C2Vsim.in are the ones listed above)
proj_dir <- getwd()
setwd("../RunC2Vsim")
system("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/bin/Simulation3.02.exe CVsim.in")
system("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/bin/Budget3.02.exe CVBudget.in")

# Read the result 
GWBUDbase <- c2vsim.readGWBUD("Results/CVground.BUD", NtimeSteps = 528)
DivBUDbase <- c2vsim.readDiversionBUD(filename = "Results/CVdiverdtl.BUD", NtimeSteps = 528)
SWHYDbase <- c2vsim.readSWHYD(filename = "Results/CVSWhyd.out", Nskip = 5, maxChar = 7000, NtimeSteps = 528)
GWHYDbase <- c2vsim.readGWHYD(filename = "Results/CVGWhyd.out", NtimeSteps = 528)
setwd(proj_dir)


### Evaluate the Pareto Solutions
psGWBUD <- vector(mode = "list", length = dim(ps[[2]])[1])
psDVBUD <- vector(mode = "list", length = dim(ps[[2]])[1])
psSWHYD <- vector(mode = "list", length = dim(ps[[2]])[1])
psGWHYD <- vector(mode = "list", length = dim(ps[[2]])[1])
for (i in seq(1,dim(ps[[2]])[1]-1,1)) {
  # make a copy of divspec
  temp_divspec <- divspec
  temp_divdata <- divdata
  
  id_active <- which( ps[[1]][i,] == 1 )
  active_elements <- divelem[id_active]
  div_node_act_elem <- divNodeId[id_active]
  unique_div_nodes <- unique(div_node_act_elem)
  
  
  for(j in 1:length(unique_div_nodes)){
    # Find how many element recharged by this diversion node
    id_el <- which(div_node_act_elem == unique_div_nodes[j])
    if (is.na(id_el[1]))
      next
    
    temp_divspec[[1]][1] <- temp_divspec[[1]][1] + 1
    temp_divspec[[2]] <- rbind(temp_divspec[[2]], 
                               c(dim(temp_divspec[[2]])[1]+1, unique_div_nodes[j], 265,1, dim(divdata)[2]-1+j, 0.95,
                                 dim(divdata)[2]-1+j, 0.05,1,1,dim(divdata)[2]-1+j,0,1,1))
    
    a <- c()
    for (k in 1:length(id_el)) {
      ii <- which(diversion_polygons$IE == active_elements[id_el[k]])
      a <- c(a, slot(slot(diversion_polygons, "polygons")[[ii]],"area")/1000000)
    }
    m <- cbind(active_elements[id_el], a)
    temp_divspec[[3]][[length(divspec[[3]])+j]] <- m
    
    # Add the time series
    idiv <- which(dts.id == unique_div_nodes[j])
    temp_divdata <- cbind(temp_divdata, dts.data[,idiv])
  }
  
  # Write the files
  c2vsim.writeDivSpec("../RunC2Vsim/tempRspec.dat", DivSpec = temp_divspec)
  c2vsim.writeDivData("../RunC2Vsim/tempRdata.dat", data = temp_divdata, NCOLDV = dim(temp_divdata)[2]-1)
  
  setwd("../RunC2Vsim")
  system("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/bin/Simulation3.02.exe CVsim.in")
  system("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/bin/Budget3.02.exe CVBudget.in")
  
  # Read the result 
  GWBUD <- c2vsim.readGWBUD("Results/CVground.BUD", NtimeSteps = 528)
  DivBUD <- c2vsim.readDiversionBUD(filename = "Results/CVdiverdtl.BUD", NtimeSteps = 528, nExpectedDivs = 249)
  SWHYD <- c2vsim.readSWHYD(filename = "Results/CVSWhyd.out", Nskip = 5, maxChar = 7000, NtimeSteps = 528)
  GWHYD <- c2vsim.readGWHYD(filename = "Results/CVGWhyd.out", NtimeSteps = 528)
  psGWBUD[[i]] <- GWBUD
  psDVBUD[[i]] <- DivBUD
  psSWHYD[[i]] <- SWHYD
  psGWHYD[[i]] <- GWHYD
  
  setwd(proj_dir)
}

#=========SAVE ALL RUN RESULTS===========
save(GWBUDbase, DivBUDbase, SWHYDbase, GWHYDbase, psGWBUD,psDVBUD, psSWHYD, psGWHYD,
     file = "../OptimResults/maxWTminArea/ParetoSolutionsBUD_58946.RData")


# Load and Process Results ------------------------------------------------
load(file = "../OptimResults/maxWTminArea/ParetoSolutionsBUD_58946.RData")


# Compare pareto Groundwater storage and gain from stream against  --------

ES <- matrix(data = 0, nrow = length(GWBUDbase[[1]][[1]]), ncol = length(psGWBUD))
GFS <- matrix(data = 0, nrow = length(GWBUDbase[[1]][[1]]), ncol = length(psGWBUD))
baseBudAll <- c2vsim.cumBUD(GWBUDbase)
for (i in 1:length(psGWBUD)) {
  if (is.null( psGWBUD[[i]]))
    next
  
  scenarioAll <- c2vsim.cumBUD(psGWBUD[[i]])
  ES[,i] <- (scenarioAll$ES - baseBudAll$ES)/1000000
  GFS[,i] <- cumsum(scenarioAll$GFS - baseBudAll$GFS)/1000000
}


# Write Ending Storage  as javascript variables -------
js_ES_file <- "../js_scripts/ES_58946.js"
js_GFS_file <- "../js_scripts/GFS_58946.js"
myfnc.writeTSMatrix2JS(js_ES_file, data = ES, varName = "ES", sy = 1965, sm = 10)
myfnc.writeTSMatrix2JS(js_GFS_file, data = -GFS, varName = "GFS", sy = 1965, sm = 10)


# write river shapefile as geojson ----------------------------------------

river_shp <- readOGR(dsn = "../gis_data/C2Vsim_rivers.shp")
river_shp_4326 <- spTransform(river_shp, CRS("+init=epsg:4326"))
writeOGR(river_shp_4326, "../js_scripts/C2Vsim_rivers.geojson", layer = "C2Vsim_rivers", driver = "GeoJSON")


# Prepare data for PLOTS with PLOTLY --------------------------------------------------
library("plotly")

tm <- seq.Date(from = as.Date(paste0(1965,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")
div_df_plot <- data.frame(tm, apply(dts.data[529:1056,], 2, cumsum))
names(div_df_plot) <- c("Time", paste0("Dnd", as.character(dts.id)))

psDiv_df_plot <- vector(mode = "list", length = length(psDVBUD))
for (i in 1:length(psDVBUD)) {
  if (is.null(psDVBUD[[i]]))
    next
  
  psDiv_df_plot[[i]] <- data.frame("Time" = tm)
  psDiv_df_plot[[i]] <- cbind(psDiv_df_plot[[i]], apply(psDVBUD[[i]][[4]][,493:498]/1000, 2, cumsum))
  Reshape(rbind(paste0("DIV", as.character(dts.id)), paste0("DEF", as.character(dts.id))),1,6)
  
  names(psDiv_df_plot[[i]])[-1] <- Reshape(rbind(paste0("DIV", as.character(dts.id)), paste0("DEF", as.character(dts.id))),1,6)
  
}



# Make plots --------------------------------------------------------------
# Choose a pareto solution
ipar <- 15
tempdf <- merge(div_df_plot, psDiv_df_plot[[ipar]] )
plot_ly(tempdf, x = ~Time, y = ~Dnd1, type = 'scatter', mode = 'lines', name = "Node 1",
        line = list(color = '#1b9e77', width = 4, dash = 'dot')) %>%
  add_trace(y = ~Dnd421, type = 'scatter', mode = 'lines', name = "Node 421",
            line = list(color = '#d95f02', width = 4, dash = 'dot')) %>%
  add_trace(y = ~Dnd23, type = 'scatter', mode = 'lines', name = "Node 23",
            line = list(color = '#7570b3', width = 4, dash = 'dot')) %>% 
  add_trace(y = ~DIV1, type = 'scatter', mode = 'lines', name = "Node 1-act",
            line = list(color = '#1b9e77', width = 2, dash = 'solid')) %>% 
  add_trace(y = ~DIV421, type = 'scatter', mode = 'lines', name = "Node 421-act",
            line = list(color = '#d95f02', width = 2, dash = 'solid')) %>%
  add_trace(y = ~DIV23, type = 'scatter', mode = 'lines', name = "Node 23-act",
            line = list(color = '#7570b3', width = 2, dash = 'solid'))





tm <- seq.Date(from = as.Date(paste0(1965,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")
ipar <- 15
df <- data.frame("Time" = tm, "Base" = cumsum(CVBase$NSI), "Scen" = cumsum(c2vsim.cumBUD(psGWBUD[[ipar]])$NSI))

plot_ly(df, x = ~Time, y = ~Base, type = 'scatter', mode = 'lines', name = "Base",
        line = list(color = '#1b9e77', width = 2, dash = 'solid')) %>%
  add_trace(y = ~Scen, type = 'scatter', mode = 'lines', name = "Scenario",
            line = list(color = '#d95f02', width = 4, dash = 'solid'))





# Plot stream hydrographs for selected pareto solution --------------------
tm <- seq.Date(from = as.Date(paste0(1965,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")
ipar <- 15
iriv <- 1
iseq <- seq(1, 9,1)


df_plot <- data.frame("Time" = tm, "Base" = SWHYDbase[[2]][,iriv], "Scenario" = psSWHYD[[ipar]][[2]][,iriv])
plot_ly(df_plot, x = ~Time, y = ~Base, type = 'scatter', mode = 'lines', name = "Node 1",
        line = list(color = '#1b9e77', width = 4, dash = 'solid')) %>%
  add_trace(y = ~Scenario, type = 'scatter', mode = 'lines', name = "Scenario",
            line = list(color = '#d95f02', width = 2, dash = 'solid'))



df_plot <- data.frame("Time" = tm, "temp" = 
                        cumsum(psSWHYD[[ipar]][[2]][,iseq[1]] - SWHYDbase[[2]][,iseq[1]]))
names(df_plot)[2] = paste0("RIV", as.character(iseq[1]))

for (i in 2:length(iseq)) {
  df_plot <- cbind(df_plot, data.frame("temp" = cumsum(psSWHYD[[ipar]][[2]][, iseq[i]] - SWHYDbase[[2]][, iseq[i]])))
  names(df_plot)[i+1] = paste0("RIV", as.character(iseq[i]))
}



plot_ly(df_plot, x = ~Time, y = ~RIV1, type = 'scatter', mode = 'lines', name = "Node 1",
       line = list(color = '#1b9e77', width = 4, dash = 'solid')) %>%
  add_trace(y = ~RIV2, type = 'scatter', mode = 'lines', name = "Node 2",
            line = list(color = '#d95f02', width = 2, dash = 'solid'))


g <- ggplot(df_plot, aes_string(x="Time", y=paste0("RIV", as.character(iseq[1])))) + geom_line()
for (i in 2:length(iseq)){
  g <- g + geom_line(mapping = aes_string(x="Time", y=paste0("RIV", as.character(iseq[i]))))
}
g


##########=====#################================================
## OLD VERSION write Pareto solution as javascript variable -----------
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







