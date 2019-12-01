library("plotly")
library("rgdal")

source("c2vsim_io.R")
source("../../nsgaii/Rnsgaii/nsgaii_io.R")

c2vsim.path <- "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/"
gis.path <- "../gis_data/"


# Read the nodes
XY <- c2vsim.readNodes(paste0(c2vsim.path, "Preprocessor/CVnode.dat"))
# Read the mesh
MSH <- c2vsim.readMesh(paste0(c2vsim.path, "Preprocessor/CVelement.dat"))

### attempt to make plots
plot_ly( x = XY$X, y = XY$Y, 
         marker = list(size = 5,
                       color = 'rgba(255, 182, 193, .9)',
                       line = list(color = 'rgba(152, 0, 0, .8)',
                                   width = 0)))


## read outline of Central valley
cvoutline <- readOGR(dsn = paste0(gis.path, "C2Vsim_outline.shp"))
#cvoutline <- st_read(dsn = paste0(gis.path, "C2Vsim_outline.shp"))
cvoutline4326 <- spTransform(cvoutline, CRS("+init=epsg:4326"))

divpolys <- readOGR(dsn = paste0(gis.path, "C2Vsim_divPolys.shp"))
divpolys4326 <- spTransform(divpolys, CRS("+init=epsg:4326"))

plot_ly(cvoutline)

leaflet(cvoutline4326) %>% addTiles() %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE)) %>%
  addPolygons(data = divpolys4326, opacity = 0.5, fillOpacity = 0.5)

pS <- nsgaii.readParetoSolution("../OptimResults/maxWTminArea/paretoSolutions.dat")




plot_ly(x=pS[[2]][,2], y=-pS[[2]][,1], type='scatter',mode="markers",
marker = list(size = 10, symbol = "circle",
              color = 'rgba(255, 182, 193, .9)',
              line = list(color = 'rgba(152, 0, 0, .8)',
                          width = 2))) %>%
  layout(
    title = "Pareto Solutions",
    xaxis = list(title = "Recharge area"),
    yaxis = list(title = "Water level Rise")
  )


# How to create a shapefile from data
# https://gis.stackexchange.com/questions/214062/create-a-shapefile-from-dataframe-in-r-keeping-attribute-table


# Read hydraulic Conductivity ---------------------------------------------
CVparam <- c2vsim.readParam(paste0(c2vsim.path, "Simulation/CVparam.dat"))
CVstrat <- c2vsim.readStrat(paste0(c2vsim.path, "Preprocessor/CVstrat.dat"))
Zelev <- matrix(data = NA, nrow = dim(CVstrat)[1], ncol = 4)
Zelev[,1] <- CVstrat[,2]
for (i in 2:4) {
  Zelev[,i] <- Zelev[,1] - apply(CVstrat[,3:(i*2)], 1, sum) 
}

# calculate element barycenters
cc <- matrix(data = 0, nrow = dim(MSH)[1], ncol = 5)
for (i in 1:4) {
  id <- which(MSH[, i+1] != 0)
  cc[id,1] <- cc[id,1] + XY[MSH[id, i+1],2]
  cc[id,2] <- cc[id,2] + XY[MSH[id, i+1],3]
  
  for (j in 1:3) {
    cc[id,j+2] <- cc[id,j+2] + Zelev[MSH[id, i+1],j] + Zelev[MSH[id, i+1],j+1]
  }
}
cc[id,1:2] <- cc[id,1:2]/4
cc[id,3:5] <- cc[id,3:5]/8
id <- which(MSH[, 5] == 0)
cc[id,1:2]  <- cc[id,1:2]/3
cc[id,3:5] <- cc[id,3:5]/6

cbind(XY[,2:3], (Zelev[,1] + Zelev[,2])/2, CVparam[[1]]$PKH)

write.table(cbind(XY[,2:3], (Zelev[,1] + Zelev[,2])/2, CVparam[[1]]$PKH), file = "temp.dat", row.names = FALSE, col.names = FALSE ,append = FALSE)
write.table(cbind(XY[,2:3], (Zelev[,2] + Zelev[,3])/2, CVparam[[2]]$PKH), file = "temp.dat", row.names = FALSE, col.names = FALSE ,append = TRUE)
write.table(cbind(XY[,2:3], (Zelev[,3] + Zelev[,4])/2, CVparam[[3]]$PKH), file = "temp.dat", row.names = FALSE, col.names = FALSE ,append = TRUE)



