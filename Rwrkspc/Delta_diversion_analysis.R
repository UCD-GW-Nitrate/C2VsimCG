library(plotly)
library(gwtools)

# Read diversion data
divSpec <- gwtools::c2vsim.readDivSpec(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdivspec.dat")
divData <- gwtools::c2vsim.readDivData(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVdiversions.dat")

# first read the base hydrograph
swhydBase <- c2vsim.readSWHYD("../RunC2Vsim/CVSWhydBase.out")
# These are river the nodes ids associated with the node that the water is diverted for the CVP
deltaNodeIds <- c(204,417,418, 419)
# They also correspond to the column in the file above that all rivers have been printed
ii <- which(swhydBase$NodeIds == deltaNodeIds[3])

# set time sequence
c2vsimTime <- convertTime(swhydBase$Time)
nt <- length(c2vsimTime$y) 
tm <- seq.Date(from = as.Date(paste0(c2vsimTime$y[1],"/",c2vsimTime$m[1],"/1")),
               to = as.Date(paste0(c2vsimTime$y[nt],"/",c2vsimTime$m[nt],"/1")),
               by = "month")

# From these plots it seems that the flows from 217 + 417 ~= 418 
{
  p <- plot_ly()
  for (i in 1:length(deltaNodeIds)) {
    p <- add_trace(p, x = tm, y = swhydBase$SWHyd[, deltaNodeIds[i]],type = "scatter", mode = 'lines', name = deltaNodeIds[i])
  }
  p
}

{
  p <- plot_ly()
  #p <- add_trace(p, x = tm, y = swhydBase$SWHyd[, 204] + swhydBase$SWHyd[, 417] ,type = "scatter", mode = 'lines', name = "both")
  p <- add_trace(p, x = tm, y = swhydBase$SWHyd[, 418],type = "scatter", mode = 'lines', name = "418")
  p <- add_trace(p, x = tm, y = swhydBase$SWHyd[, 419],type = "scatter", mode = 'lines', name = "419")
  p
}

# Read the Delta Conditions
delta_cond <- readxl::read_xlsx("delta_conditions.xlsx")
delta_states <- unique(delta_cond$DeltaConditions)
delta_cond$DCnum <- 0
for (i in 1:length(delta_states)) {
  ii <- which(delta_cond$DeltaConditions == delta_states[i])
  delta_cond$DCnum[ii] <- i
}

{
  tm1 <- seq.Date(from = as.Date(paste0(2008,"/",4,"/1")),
                 to = as.Date(paste0(c2vsimTime$y[nt],"/",c2vsimTime$m[nt],"/1")),
                 by = "month")
  {
    p <- plot_ly()
    p <- add_trace(p, x = tm1, y = swhydBase$SWHyd[1039:nt, 418],type = "scatter", mode = 'lines',yaxis = "y2",
                   name = "418",line = list(width = 3.5))
    
    p <- add_trace(p, x = delta_cond$ST_Date, y = delta_cond$DCnum,type = "scatter", mode = 'lines',
                   fill = 'tozeroy', name="Delta conditions",fillcolor = 'rgba(255, 212, 96, 0.5)', line = list(width = 0.0))
    
    p %>% layout(yaxis2 = list(overlaying = "y",side = "right", title = "Delta streamflow"), 
                 yaxis = list(ticktext = delta_states, tickvals = c(1:5),  title = 'Excess conditions'),
                 legend = list(x = 0.5, y = 0.1), margin = list(l=50,r=50,b=50,t=50))
  }
}

{# plot existing diversions from the delta
  tm2 <- seq.Date(from = as.Date(paste0(1921,"/",10,"/1")),
                 to = as.Date(paste0(2009,"/",9,"/1")),
                 by = "month")
  ii <- which(divSpec$RDV[,2] == 418)
  {
    p <- plot_ly()
    for (i in 1:length(ii)) {
      p <- add_trace(p, x = tm2, y = cumsum(divData[,ii[i]+1])/1000, type = 'scatter', mode = "lines", name = divSpec$RDVnames[ii[i]] )
    }
    p %>% layout(legend = list(x = 0.1, y = 0.9), margin = list(l=50,r=30,b=30,t=30),
                 xaxis = list(title = "TIme"),
                 yaxis = list(title = "Cumulative Diversion [MAF]"),
                 title = "Existing Diversions from Delta")
  }
  
  
}
