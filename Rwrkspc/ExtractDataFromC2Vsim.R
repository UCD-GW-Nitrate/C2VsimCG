source("c2vsim_io.R")
source("d:/giorgk/Documents/GitHub/nsgaii/Rnsgaii/nsgaii_io.R")

c2vsim.path <- "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/"
gis.path <- "../gis_data/"


# Read the nodes
XY <- c2vsim.Nodes(paste0(c2vsim.path, "Preprocessor/CVnode.dat"))
# Read the mesh
MSH <- c2vsim.Mesh(paste0(c2vsim.path, "Preprocessor/CVelement.dat"))

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
