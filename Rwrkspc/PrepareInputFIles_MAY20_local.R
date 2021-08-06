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
c2vsim_nodes <- readOGR(dsn = "../gis_data/C2Vsim_nodes.shp")
c2vsim_rivnodes <- readOGR(dsn = "../gis_data/C2Vsim_riverNodes.shp")
c2vsim_rivers <- readOGR(dsn = "../gis_data/C2Vsim_rivers.shp")
MSHfile <- c2vsim.readMesh(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Preprocessor/CVelement.dat")
{ # calculate element Barycenters
  bc_elem <- matrix(data = NA, nrow =  length(c2vsim_mesh), ncol = 2)
  for (i in 1:length(c2vsim_mesh)) {
    crds <- c2vsim_mesh@polygons[[i]]@Polygons[[1]]@coords
    n = dim(crds)[1]-1
    bc_elem[i,] <- apply(crds[1:n,],2,mean)
  }
}


{# Scen1 Local only
  # This is the list of river names where we divert water from
  riverNodeName <- c("Calaveras", "Stanislaus", "Tuolumne", 
                     "Merced", "Chowchilla", "Fresno",
                     "San Joaquin A", "Kings", "Kaweah",
                     "Tule", "Kern", "San Joaquin B")#), "Kern", "Kings", "Kaweah", "Tule")
  # This is the river ID where the water is extracted during the simulation
  riverNodeList <- c(166, 153, 141, 
                     120, 83, 72, 
                     61, 29, 427,
                     19, 5, 156)#5, 29, 427, 19)
  # This is the index of diversions in the AllNames list
  index_AllNames <- c(46, 86, 49, 51, 50, 85, 84, 
                      83, 53, 52, 54, 55, 
                      56, 60, 61, 62, 63, 64, 66, 67, 68,
                      69, 71, 72, 73, 70, 57, 58, 59) # 76, 70, 71, 72, 73, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69)
  # This is the index that associates the riverNodeList and index_AllNames.
  # It is possible that more than one index_AllNames may correspond to one riverNodeList
  div_id <- c(1, 2, 2, 3, 3, 3, 3,
             4, 4, 4, 5, 6, 
             7, 8, 8, 8, 8, 8, 9, 9, 9,
             10, 11, 11, 11, 11, 12, 12, 12)
  
  # This is the river id where the water is extracted in the original C2Vsim model
  div_riv_id_actual <- c(165, 146, 152, 142, 136, 138, 137,
                         116, 117, 119, 81, 70,
                         60, 24, 25, 28, 43, 37, 421, 422, 426,
                         18, 2, 3, 4, 7, 145, 144, 115)
  # this is the river id in c2vsim_rivers shapefiles
  riv_id <- c(25, 23, 21, 
              17, 11, 9, 
              7, 2, 5,
              6, 1, 24, 22, 20, 16)
}


{# make a unique association between diversion nodes and elements based on the nearest distance
  cnt <- 1
  div_df <- data.frame(matrix(data = NA, nrow = 1000, ncol = 10))
  # RiverName The name of the river where the water is diverted from
  # DivNdId River node ID (IRV) where the water is extracted from the river in C2Vsim simulation
  # DivNdElev The elevation of the above node
  # DivNdOpt The river node id where the water is extracted from the river during optimization.
  # ElemID: The element id that receives water from the diversion node
  # ElemElevAv : The average element elevation (avergae of the corner elements)
  # ElemElevMax : The elevation of the element corner with maximum elevation
  # DistAct: The actual distance between Elem barycenter and the DivNdId
  # DistNear: The distance between the ElemID barycenter and the closest stream node of the RiverName stream
  # NearNdElev: The elevation of the node that DistNear corresponds to
  div_df_names <- c("RiverName", "DivNdId", "DivNdElev", "DivNdOpt", "ElemID", 
                    "ElemElevAv", "ElemElevMax",  "DistAct", "DistNear", "NearNdElev") 
  names(div_df ) <- div_df_names
  for (i in 1:length(index_AllNames)) {
    # find the groundwater node id for the diversion node that the water is actually diverted from in the model
    igw <- c2vsim_rivnodes$IGW[which(c2vsim_rivnodes$IRV == div_riv_id_actual[i])]
    divNdActElev <- c2vsim_nodes$ELV[which(c2vsim_nodes$ID == igw)]
    elem_list <- AllElem[[index_AllNames[i]]]
    {# Corrections on element list
      if (i == 1){
        elem_list[1] <- 547
        elem_list[6] <- 585
      }
    }
    
    for (j in 1:length(elem_list)) {
      if (elem_list[j] == 0){
        next
      }
      if (c2vsim_mesh$ISGE[elem_list[j]] < 8){ # The diversions is northern of Stockton and not considered
        next
      }
      dist <- sqrt(sum((bc_elem[elem_list[j],] - as.numeric(c2vsim_nodes@coords[igw,]))^2))/1000
      
      el_nd_ind <- MSHfile[elem_list[j],2:5]
      el_nd_ind <- as.numeric(el_nd_ind[which(el_nd_ind != 0)])
      
      elem_elev_mean <- mean(c2vsim_nodes$ELV[el_nd_ind])
      elem_elev_max <- max(c2vsim_nodes$ELV[el_nd_ind])
      
      # Find the river node that is closest to the element
      if (div_id[i] == 12){
        rid = riv_id[12:15]
      }
      else{
        rid <- riv_id[div_id[i]]
      }
      dist_near <- 999999999
      near_elev <- 99999999
      for (ii in 1:length(rid)) {
        for (kk in 1:dim(c2vsim_rivers@lines[[rid[ii]]]@Lines[[1]]@coords)[1]) {
          temp_dst <- sqrt(sum((c2vsim_rivers@lines[[rid[ii]]]@Lines[[1]]@coords[kk,] - bc_elem[elem_list[j],])^2))/1000
          if (temp_dst < dist_near){
            dist_near <- temp_dst
            igw_near <- which(sqrt((c2vsim_nodes@coords[, 1] - c2vsim_rivers@lines[[rid[ii]]]@Lines[[1]]@coords[kk, 1])^2 + 
                        (c2vsim_nodes@coords[, 2] - c2vsim_rivers@lines[[rid[ii]]]@Lines[[1]]@coords[kk, 2])^2) < 0.1)
            near_elev <- c2vsim_nodes$ELV[igw_near]
          }
        }
      }
      
      # Check if the element has already a diversion node that can receive water from
      tmp_id <- which(div_df$ElemID == elem_list[j])
      if (length(tmp_id) == 0){# If it does not exist add it
        div_df$RiverName[cnt] <-  riverNodeName[div_id[i]] # River name
        div_df$DivNdId[cnt] <- div_riv_id_actual[i] #Diversion node id
        div_df$DivNdElev[cnt] <- divNdActElev # Diversion node elevation
        div_df$DivNdOpt[cnt] <- riverNodeList[div_id[i]] # Diversion node during optimization
        div_df$ElemID[cnt] <- elem_list[j] # Element ID
        div_df$ElemElevAv[cnt] <- elem_elev_mean # mean elevation
        div_df$ElemElevMax[cnt] <- elem_elev_max # max elevation
        div_df$DistAct[cnt] <- dist # actual distance
        div_df$DistNear[cnt] <- dist_near # distance to nearest node
        div_df$NearNdElev[cnt] <- near_elev # Elevation of the node
        cnt <- cnt + 1
      }
      else{
        print(c(i,j,tmp_id))
        if (div_df$DistAct[tmp_id] > dist){
          print(paste("Update node for element", div_df$Elem[tmp_id]))
          div_df$RiverName[tmp_id] <- riverNodeName[div_id[i]]
          div_df$DivNdId[tmp_id] <-  div_riv_id_actual[i]
          div_df$DivNdElev[tmp_id] <- divNdActElev
          div_df$DivNdOpt[tmp_id] <- riverNodeList[div_id[i]]
          div_df$DistAct[tmp_id] <- dist
          div_df$DistNear[tmp_id] <- dist_near
          div_df$NearNdElev[tmp_id] <- near_elev
        }
      }
    }
  }
  div_df <- div_df[-(cnt:1000),]
}

{# Write info to excel
  write.csv(div_df, file = "LocalDiversions_Jul20.csv",row.names = F, col.names = T)
}

{# Diversion elements from Friant Kern Canal and Delta
  # set ii 77 for Friant Kern Canal
  # set ii c(79, 47, 80) For Delta
  ii <- c(79, 47, 80)
  if (length(ii) == 1){
    riv_name <- "Friant-Kern Canal"
    divrividact <- 54
    fileName <- "FriantDiversions_Jul20.csv"
  }
  else{
    riv_name <- "Delta"
    divrividact <- 418
    fileName <- "DeltaDiversions_Jul20.csv"
  }
  igw <- c2vsim_rivnodes$IGW[which(c2vsim_rivnodes$IRV == divrividact)]
  divNdActElev <- c2vsim_nodes$ELV[which(c2vsim_nodes$ID == igw)]
  
  cnt <- 1
  CVwide_df <- data.frame(matrix(data = NA, nrow = 1000, ncol = 8))
  # RiverName The name of the river where the water is diverted from
  # DivNdId River node ID (IRV) where the water is extracted from the river in C2Vsim simulation
  # DivNdElev The elevation of the above node in ft
  # DivNdOpt The river node id where the water is extracted from the river during optimization.
  # ElemID: The element id that receives water from the diversion node
  # ElemElevAv : The average element elevation (average of the corner elements) in ft
  # ElemElevMax : The elevation of the element corner with maximum elevation in ft
  # DistAct: The actual distance between Elem barycenter and the DivNdId
  div_df_names <- c("RiverName", "DivNdId", "DivNdElev", "DivNdOpt", "ElemID", 
                    "ElemElevAv", "ElemElevMax",  "DistAct") 
  names(CVwide_df ) <- div_df_names
  
  for (i in ii) {
    print(i)
    elem_list <- AllElem[[i]]
    for (j in 1:length(elem_list)){
      if (elem_list[j] == 0){
        next
      }
      dist <- sqrt(sum((bc_elem[elem_list[j],] - as.numeric(c2vsim_nodes@coords[igw,]))^2))/1000
      
      el_nd_ind <- MSHfile[elem_list[j],2:5]
      el_nd_ind <- as.numeric(el_nd_ind[which(el_nd_ind != 0)])
      
      elem_elev_mean <- mean(c2vsim_nodes$ELV[el_nd_ind])
      elem_elev_max <- max(c2vsim_nodes$ELV[el_nd_ind])
      
      CVwide_df$RiverName[cnt] <-  riv_name # River name
      CVwide_df$DivNdId[cnt] <- divrividact #Diversion node id
      CVwide_df$DivNdElev[cnt] <- divNdActElev # Diversion node elevation
      CVwide_df$DivNdOpt[cnt] <- divrividact # Diversion node during optimization
      CVwide_df$ElemID[cnt] <- elem_list[j] # Element ID
      CVwide_df$ElemElevAv[cnt] <- elem_elev_mean # mean elevation
      CVwide_df$ElemElevMax[cnt] <- elem_elev_max # max elevation
      CVwide_df$DistAct[cnt] <- dist # actual distance
      cnt <- cnt + 1
    }
  }
  
  CVwide_df <- CVwide_df[-(cnt:1000),]
  
  {# Write to file
    write.csv(CVwide_df, file = fileName, row.names = F, col.names = T)
  }
  
}




{# Add Friant Kern Canal
  # Friant Kern Canal id in AllElem
  ii = 77
  Friant_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  names(Friant_df ) <- c("Elem", "NdAct", "Dist") 
  id_nd <- c2vsim_rivnodes$IGW[which(c2vsim_rivnodes$IRV == 54)]
  for (i in 1:length(AllElem[[ii]])) {
    dist <- sqrt(sum((bc_elem[AllElem[[ii]][i],]- as.numeric(c2vsim_nodes@coords[id_nd,]))^2))/1000
    Friant_df <- rbind(Friant_df, c(AllElem[[ii]][i],  54,  dist) )
  }
  names(Friant_df ) <- c("Elem", "NdAct", "Dist") 
}

{# Add Delta Canal
  # Delta canal and 
  Delta_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  names(Delta_df ) <- c("Elem", "NdAct", "Dist") 
  id_nd <- c2vsim_rivnodes$IGW[which(c2vsim_rivnodes$IRV == 417)]
  for (ii in c(79, 47, 80)) {# San Luis Canal, Delta to SR9, O'Neill Forebay
    for (i in 1:length(AllElem[[ii]])) {
      tmp_id <- which(Delta_df$Elem == AllElem[[ii]][i])
      if (length(tmp_id) == 0){
        dist <- sqrt(sum((bc_elem[AllElem[[ii]][i],]- as.numeric(c2vsim_nodes@coords[id_nd,]))^2))/1000
        Delta_df <- rbind(Delta_df, c(AllElem[[ii]][i],  417,  dist) )
      }
    }
  }
  names(Delta_df ) <- c("Elem", "NdAct", "Dist") 
}


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
  id_riv <- 419
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
