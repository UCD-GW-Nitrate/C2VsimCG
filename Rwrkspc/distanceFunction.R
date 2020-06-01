library(rgdal)
library(gwtools)

rivers <- readOGR(dsn = "../gis_data", layer = "C2Vsim_rivers") 

{# Scenario 1
  scen1_poly <- readOGR(dsn = "../gis_data", layer = "MAY2020_scen1_elem")
  
  uniqueDivNodes <- unique(scen1_poly$divID)
  uniqueDivNodes <- uniqueDivNodes[-which(uniqueDivNodes == 0)]
  riv_per_div <- vector(mode = "list", length = length(uniqueDivNodes))
  riv_per_div[[1]] <- c(1,75)
  riv_per_div[[2]] <- c(2,3,4)
  riv_per_div[[3]] <- c(5)
  riv_per_div[[4]] <- c(6)
}

for (i in 1:length(scen1_poly)) {
  if (scen1_poly$divID[i] == 0){
    scen1_poly$divDist[i] = -1
    next
  }
  div_id <- scen1_poly$divID[i]
  n <- dim(scen1_poly@polygons[[i]]@Polygons[[1]]@coords)[1]-1
  xp <- mean(scen1_poly@polygons[[i]]@Polygons[[1]]@coords[1:n,1])
  yp <- mean(scen1_poly@polygons[[i]]@Polygons[[1]]@coords[1:n,2])
  
  min_dist <- 9999999999999
  
  for (j in 1:length(riv_per_div[[div_id]])) {
    riv_id <- riv_per_div[[div_id]][j]
    for(k in 1:length(rivers@lines[[riv_id]]@Lines)){
      crd <- rivers@lines[[riv_id]]@Lines[[k]]@coords
      nseg <- dim(crd)[1]-1
      for (ii in 1:nseg) {
        dst <- gwtools::distPointLineSeg(xp, yp, crd[ii,1], crd[ii,2], crd[ii+1,1], crd[ii+1,2])
        if (dst < min_dist){
          min_dist <-  dst
        }
      }
    }
  }
  scen1_poly$divDist[i] = min_dist
}

writeOGR(obj =scen1_poly, dsn = "../gis_data", layer = "MAY2020_scen1_elem_dst", driver = ogrDrivers()$name[17])
