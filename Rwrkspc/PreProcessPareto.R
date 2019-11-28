library("plotly")
library("rgdal")

source("c2vsim_io.R")
SWHYD <- c2vsim.readSWHYD("../RunC2Vsim/Results/CVSWhydBase.out")

tm <- seq.Date(from = as.Date(paste0(1921,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")

riverNodeList <- c(1,23,421)
rm(DiversionTimeSeriesTAF)
DiversionTimeSeriesTAF <- data.frame("Time" = tm)
cumDTS <- data.frame("Time" = tm)
for (i in 1:length(riverNodeList)) {
  icol <- which(SWHYD[[1]] == riverNodeList[i])
  perc_95 <- as.numeric(quantile(SWHYD[[2]][,icol],0.95))
  dts <- SWHYD[[2]][,icol]
  id_above <- which(dts >= perc_95)
  id_below <- which(dts < perc_95)
  dts[id_above] <- dts[id_above] - perc_95
  dts[id_below] <-  0
  DiversionTimeSeriesTAF$TEMP <- dts/1000
  cumDTS$TEMP <- cumsum(dts/1000)
  names(DiversionTimeSeriesTAF)[names(DiversionTimeSeriesTAF) == "TEMP"] <- paste0("ND", as.character(riverNodeList[i]))
  names(cumDTS)[names(cumDTS) == "TEMP"] <- paste0("ND", as.character(riverNodeList[i]))
}

tempdf <- DiversionTimeSeriesTAF[528:1056,]
tempdf <- lapply(DiversionTimeSeriesTAF[528:1056,][-1],cumsum)
tempdf[-1] <- lapply(tempdf[-1], cumsum)

plot_ly(tempdf, x = ~Time, y = ~ND1, type = 'scatter', mode = 'lines', name = "KERN RIVER (1)",
        line = list(color = '#1b9e77', width = 4, dash = 'solid')) %>%
  add_trace(y = ~ND23, type = 'scatter', mode = 'lines', name = "KINGS RIVER (23)",
            line = list(color = '#d95f02', width = 4, dash = 'solid')) %>%
  add_trace(y = ~ND421, type = 'scatter', mode = 'lines', name = "KAWEAH RIVER (421)",
            line = list(color = '#7570b3', width = 4, dash = 'solid'))

# Print the diversion Time series into file -------------------------
dtsFile <- "../RunC2Vsim/divTimeSeries.dat"
cat(length(riverNodeList), " ", dim(DiversionTimeSeriesTAF)[1],"\n", sep = "",  file = dtsFile)
for (i in 1:length(riverNodeList)) {
  nm <- names((DiversionTimeSeriesTAF))[i+1]
  cat(substr(nm,3,nchar(nm))," ", file = dtsFile, append = TRUE)
  cat(DiversionTimeSeriesTAF[,i+1], sep = " ", file = dtsFile, append = TRUE)
  cat("\n", file = dtsFile, append = TRUE)
}

# Write the element info for each element. -------------
# At the moment each element is associated with an area which serves as objective function.
C2Vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")
a <- vector(mode = "numeric", length = length(C2Vsim_mesh))
for (j in 1:length(C2Vsim_mesh)) {
  a[j] <- slot(slot(C2Vsim_mesh, "polygons")[[j]], "area")/1000000
}



# Write element Info Distance ---------------------------------------------
C2Vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")
diversion_elements <- readOGR(dsn = "../gis_data/diversionElements.shp")
c2vsimRivers <- readOGR(dsn = "../gis_data/C2Vsim_riverNodes.shp")
a <- vector(mode = "numeric", length = length(C2Vsim_mesh))
for (i in 1:length(C2Vsim_mesh)) {
  id <- which(diversion_elements$IE == C2Vsim_mesh$IE[i])
  if (length(id) == 0){
    a[i] <- -9
    next
  }
  
  div_id <- diversion_elements$DivNode[id]
  if (is.na(div_id)){
    a[i] <- -9
    next
  }
  
  id_rivND <- which(c2vsimRivers$IRV == div_id)
  div_coord <- slot(c2vsimRivers, "coords")[id_rivND,]

  crd <- shapefile.coords(C2Vsim_mesh, i)
  cc <- calcCentroid(crd)
  a[i] <- sqrt(sum((cc - div_coord)^2))/1000
}

elemInfoFile <- "../RunC2Vsim/ElemInfoDist.dat"
cat(length(a), "\n", sep = "",  file = elemInfoFile)
write(x = t(cbind(C2Vsim_mesh$IE,a)), file = elemInfoFile, ncolumns = 2, sep = " ", append = TRUE)



# Write diversion node element correspondance -----------------------------
diversion_elements <- readOGR(dsn = "../gis_data/diversionElements.shp")
divElemFile <- "../RunC2Vsim/divElem.dat"
cat(length(riverNodeList), "\n", sep = "",  file = divElemFile)
for (i in 1:length(riverNodeList)){
  # find the elements that are recharged by the rivernode
  ielems <- which(diversion_elements$DivNode == riverNodeList[i]) 
  elem_ids <- diversion_elements$IE[ielems]
  cat(riverNodeList[i], " ", length(elem_ids), " ", file = divElemFile, append = TRUE)
  cat(elem_ids, sep = " ", file = divElemFile, append = TRUE)
  cat("\n", file = divElemFile, append = TRUE)
}


