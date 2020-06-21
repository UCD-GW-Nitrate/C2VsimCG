library(rgdal)
library(gwtools)

rivers <- readOGR(dsn = "../gis_data", layer = "C2Vsim_rivers") 
river_nodes <- readOGR(dsn = "../gis_data", layer = "C2Vsim_riverNodes")
river_nodes_coord <- matrix(data = NA, nrow = length(river_nodes), ncol = 2)
river_nodes_coord[,1] <- as.numeric(river_nodes@coords[,1]) 
river_nodes_coord[,2] <- as.numeric(river_nodes@coords[,2])
CVstrat <- gwtools::c2vsim.readStrat("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Preprocessor/CVstrat.dat")
CVmsh <- gwtools::c2vsim.readMesh("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Preprocessor/CVelement.dat")

{# Scenario 1
  scen1_poly <- readOGR(dsn = "../gis_data", layer = "MAY2020_scen1_elem")
  
  uniqueDivNodes <- unique(scen1_poly$divID)
  uniqueDivNodes <- uniqueDivNodes[-which(uniqueDivNodes == 0)]
  riv_per_div <- vector(mode = "list", length = length(uniqueDivNodes))
  riv_per_div[[1]] <- c(1)
  riv_per_div[[2]] <- c(2,3)
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
  mean_elevation <- mean(CVstrat$ELV[as.numeric(CVmsh[i, 2:(n+1)])])
  
  min_dist <- 9999999999999
  
  for (j in 1:length(riv_per_div[[div_id]])) {
    riv_id <- riv_per_div[[div_id]][j]
    for(k in 1:length(rivers@lines[[riv_id]]@Lines)){
      crd <- rivers@lines[[riv_id]]@Lines[[k]]@coords
      
      for (ii in 1:dim(crd)[1]) {
        tmp_dst <- sqrt((river_nodes_coord[,1] - crd[ii,1])^2 + (river_nodes_coord[,2] - crd[ii,2])^2)
        
        rnd_elev <- CVstrat$ELV[river_nodes$IGW[which(tmp_dst == min(tmp_dst))][1]]
        if (rnd_elev > mean_elevation){
          dst <- sqrt((crd[ii,1] - xp)^2 + (crd[ii,2] - yp)^2 )
          if (dst < min_dist){
            min_dist <-  dst
          }
        }
      }
      #nseg <- dim(crd)[1]-1
      #for (ii in 1:nseg) {
      #  dst <- gwtools::distPointLineSeg(xp, yp, crd[ii,1], crd[ii,2], crd[ii+1,1], crd[ii+1,2])
      #  if (dst < min_dist){
      #    min_dist <-  dst
      #  }
      #}
    }
  }
  scen1_poly$divDist[i] = min_dist
}

writeOGR(obj =scen1_poly, dsn = "../gis_data", layer = "MAY2020_scen1_elem_dst1", driver = ogrDrivers()$name[17])
