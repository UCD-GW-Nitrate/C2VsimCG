library(rgdal)
library(plotly)
library(gwtools)
library(readxl)

# ------------------------------
# MAY20 input files
# In this run we will remove the delta contribution to run a local scenario
# See the PrepareInputFIles.R which is the parent script of this

# Read the element prices
IE_price <- read_xlsx('cv_land_price_022620.xlsx')
# Remove the missing prices
IE_price <- IE_price[-which(is.na(IE_price$pred_price_county_fe)),]

# ElemInfo File
# Runc first the code of the AssociateDiversionNodesElements.R script to create
# the AllNames and AllElem variables or load the saved data
load(file = "AssociateDiversionNodes.RData")
# @@@@@@@@@@@@@@@@@@@@@@@@@
c2vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_mesh.shp")

{# Without Friant Kern and Delta Scenario 1
  scen_name <- "scen1"
  index_AllNames <- c(76, 70, 71, 72, 73, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69)
  # Kern -> 1, Kings->2 Kaweah->3, Tule->4
  div_id <-         c( 3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  6)-2
  elem_mesh_db <- c2vsim_mesh@data
  elem_mesh_db$divID <- vector(mode = "numeric", length = dim(elem_mesh_db)[1])
  riverNodeList <- c(5, 29, 427, 19)
  riverNodeName <- c("Kern", "Kings", "Kaweah", "Tule")
  riverNodeColors <- c('#4daf4a','#984ea3','#ff7f00','#cccc29')
}

{# With Friant Kern but without Delta Scenario 2
  scen_name <- "scen2"
  index_AllNames <- c(77, 76, 70, 71, 72, 73, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69)
  # Friant->1, Kern -> 2, Kings->3 Kaweah->4, Tule->5
  div_id <-         c( 2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  6)-1
  elem_mesh_db <- c2vsim_mesh@data
  elem_mesh_db$divID <- vector(mode = "numeric", length = dim(elem_mesh_db)[1])
  riverNodeList <- c(55, 5, 29, 427, 19)
  riverNodeName <- c("Friant-Kern", "Kern", "Kings", "Kaweah", "Tule")
  riverNodeColors <- c('#377eb8','#4daf4a','#984ea3','#ff7f00','#cccc29')
}

{# With Delta and Friant Canal Scenario 3
  scen_name <- "scen3"
  index_AllNames <- c(79, 81, 82, 80, 47, 77, 76, 70, 71, 72, 73, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69)
  # Delta->1, Friant->2, Kern -> 3, Kings->4 Kaweah->5, Tule->6
  div_id <-         c( 1,  1,  1,  1,  1,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  6) 
  elem_mesh_db <- c2vsim_mesh@data
  elem_mesh_db$divID <- vector(mode = "numeric", length = dim(elem_mesh_db)[1])
  riverNodeList <- c(419, 55, 5, 29, 427, 19)
  riverNodeName <- c("Delta", "Friant-Kern", "Kern", "Kings", "Kaweah", "Tule")
  riverNodeColors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#cccc29')
  
}


for (i in seq(length(index_AllNames),1,-1)) {
  print(AllNames[index_AllNames[i]])
  ids <- AllElem[[index_AllNames[i]]]
  # remove the zero elements
  zero_el <- which(ids == 0)
  if (length(zero_el)){
    ids <- ids[-zero_el]
  }
  #### For Delta remove the elements that are located above Levis
  ## These are the elements with div_id[i] == 1 and ids <= 901
  if (scen_name == "scen3"){
    if (div_id[i] == 1){
      remove_ids <- which (ids <= 901)
      ids <- ids[-remove_ids]
    }
  }
  if (length(ids)!=0)
    elem_mesh_db$divID[ids] <- div_id[i]
}

elem_mesh_db$price <- 0
elem_mesh_db$price[IE_price$ie] <- IE_price$pred_price_county_fe


c2vsim_mesh@data <- elem_mesh_db
# write the shapefile
writeOGR(obj = c2vsim_mesh, dsn = "../gis_data", layer = paste0('MAY2020_', scen_name, '_elem'), driver = ogrDrivers()$name[17])

# ---- write the divElem file
divElemFile <- paste0("../OptimResults/inputFiles/divElem_MAY20_", scen_name,".dat")
cat(length(riverNodeList), "\n", sep = "",  file = divElemFile)
elem_unique <- c()
for (i in 1:length(riverNodeList)){
  ielems <- which(elem_mesh_db$divID == i)
  elem_ids <- c2vsim_mesh$IE[ielems]
  #has_price <- match(elem_ids,IE_price$ie)
  #elem_ids <- elem_ids[!is.na(has_price)]
  elem_unique <- c(elem_unique,elem_ids)
  cat(riverNodeList[i], " ", length(elem_ids), " ", file = divElemFile, append = TRUE)
  cat(elem_ids, sep = " ", file = divElemFile, append = TRUE)
  cat("\n", file = divElemFile, append = TRUE)
}

# Calculate the diversion time series
SWHYD <- gwtools::c2vsim.readSWHYD("../RunC2Vsim/CVSWhydBase.out")
tm <- seq.Date(from = as.Date(paste0(1921,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")
DiversionTimeSeriesTAF <- data.frame("Time" = tm)
cumDTS <- data.frame("Time" = tm)
delta_upp_lim <- 200000
for (i in 1:length(riverNodeList)){
  prcnts <- c(90, 95)
  # if (i == 1){prcnts <- c(98, 99)}
  
  icol <- which(SWHYD$NodeIds == riverNodeList[i])
  for (j in 1:length(prcnts)) {
    perc_tl <- as.numeric(quantile(SWHYD$SWHyd[,icol], prcnts[j]/100))
    dts <- SWHYD$SWHyd[,icol]
    id_above <- which(dts >= perc_tl)
    id_below <- which(dts < perc_tl)
    dts[id_above] <- dts[id_above] - perc_tl
    dts[id_below] <-  0
    
    # for the delta diversions we set an upper limit
    if (scen_name == "scen3"){
      if (i == 1){
        dts[id_above[which(dts[id_above] > delta_upp_lim)]] <- delta_upp_lim
      }
    }
    
    DiversionTimeSeriesTAF$TEMP <- dts/1000
    cumDTS$TEMP <- cumsum(dts/1000)
    names(DiversionTimeSeriesTAF)[names(DiversionTimeSeriesTAF) == "TEMP"] <- paste0("ND", as.character(riverNodeList[i]), "_", as.character(prcnts[j]))
    names(cumDTS)[names(cumDTS) == "TEMP"] <- paste0("ND", as.character(riverNodeList[i]), "_", as.character(prcnts[j]))
  }
}

prc_tl <- 90
if (prc_tl == 95) {iprc <- c(seq(2,17,2))} # 95
if (prc_tl == 90) {iprc <- c(seq(3,17,2))} # 90
tempdf <- DiversionTimeSeriesTAF[1:1056,-iprc]# 528:1056
cumtempdf <- apply(tempdf[,-1], 2, cumsum)

{# Write data for gnuplot
  write.table(data.frame(tempdf$Time[528:1056], apply(tempdf[528:1056,-1], 2, cumsum)/1000), 
              file = "../OptimResults/dts_local_90.data",append = F, row.names = F, col.names = F)
}

p <- plot_ly(tempdf)
for (i in 1:(dim(tempdf)[2]-1)) {
  p <- add_trace(p, x = ~Time, y = cumtempdf[,i]/1000,  type = 'scatter', mode = 'lines', name = riverNodeName[i],
                 line = list(color = riverNodeColors[i], width = 4, dash = 'solid'))
}
p %>%
  layout(title = paste0("Excess flow above ", prc_tl, "th percentile") , #and 1% from Delta"
         xaxis = list(title = "Time", showline = TRUE,linewidth = 1.0,zeroline = F,
                      tickfont = list(size = 14), titlefont = list(size = 16) ), 
         yaxis = list(title = "[MAF]", showline = TRUE,linewidth = 1.0,zeroline = F,
                      tickfont = list(size = 14), titlefont = list(size = 16) ),
         legend=list(x = 0.05, y = 0.95, font = list(size = 14))) #, type = "log"


# 850 x 530
# orca(p,"../OptimResults/DTS_local_95.png", width = 850*2, height = 530*2 )
  
{# plot stream flow hydrograph and the percentile thresholds
  id_riv <- 55
  p <- plot_ly()
  p <- add_trace(p, x = DiversionTimeSeriesTAF$Time, y=SWHYD$SWHyd[,id_riv]/1000000, type = 'scatter', mode = 'lines', name = 'Streamflow')
  p90 <- as.numeric(quantile(SWHYD$SWHyd[,id_riv],0.9))/1000000
  p95 <- as.numeric(quantile(SWHYD$SWHyd[,id_riv],0.95))/1000000
  p <- add_trace(p, x = c(DiversionTimeSeriesTAF$Time[1], DiversionTimeSeriesTAF$Time[1056]),
                 y = c(p90,p90), type = 'scatter', mode = 'lines', name = '90th percentile')
  p <- add_trace(p, x = c(DiversionTimeSeriesTAF$Time[1], DiversionTimeSeriesTAF$Time[1056]),
                 y = c(p95,p95), type = 'scatter', mode = 'lines', name = '95th percentile')
  p %>% layout(title = "Streamflow hydrograph for Delta" , 
            xaxis = list(title = "Time"), 
            yaxis = list(title = "[MAF]"),
            legend = list(x = 0.1,y = 0.95))
}

{# Write selected hydrographs as files
  write.table(data.frame(cumDTS$Time, SWHYD$SWHyd[,55]/1000000, cumDTS[,2:3]/1000), 
              file = "../OptimResults/FrianKernHyd.data",append = F, row.names = F, col.names = F)
}

# ---- write the diversion Time series file
dtsFile <- paste0('../OptimResults/inputFiles/DTS_MAY20_', scen_name, "_" , prc_tl ,'.dat')
if (scen_name == "scen3"){
  dtsFile <- paste0('../OptimResults/inputFiles/DTS_MAY20_', scen_name, "_", delta_upp_lim/1000, "K_" , prc_tl ,'.dat')
}


cat(length(riverNodeList), " ", dim(DiversionTimeSeriesTAF)[1],"\n", sep = "",  file = dtsFile)
for (i in 1:length(riverNodeList)) {
  cat(riverNodeList[i]," ", file = dtsFile, append = TRUE)
  cat(tempdf[,i+1], sep = " ", file = dtsFile, append = TRUE)
  cat("\n", file = dtsFile, append = TRUE)
}

# ---- write the Cost for each element
# Run the code in EconomicObjFnc.R to generate the elem_price
# The cost is in the IE_price data frame. However there are elements without cost.
# For those elements we will assign a very large value

{# Calculate element area
  area <- vector(mode = "numeric", length = length(c2vsim_mesh))
  for (i in 1:length(c2vsim_mesh)) {
    area[i] <- c2vsim_mesh@polygons[[i]]@area
  }
}

ie_no_price <- which(is.na(match(c2vsim_mesh$IE, IE_price$ie)))
tmp_price <- vector(mode = "logical", length = length(c2vsim_mesh$IE))
tmp_price[IE_price$ie] <- IE_price$pred_price_county_fe
tmp_price[ie_no_price] <- max(tmp_price)*3
elemInfoFile <- "../OptimResults/inputFiles/ElemCost_MAY20.dat"
cat(length(tmp_price), "\n", sep = "",  file = elemInfoFile)
write(x = t(cbind(c2vsim_mesh$IE, tmp_price/1000, area/(1000*1000))), file = elemInfoFile, ncolumns = 3, sep = " ", append = TRUE)
