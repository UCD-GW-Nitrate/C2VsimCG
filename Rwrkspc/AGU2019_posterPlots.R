library("gridExtra")
source("c2vsim_io.R")
source('~/GitHub/nsgaii/Rnsgaii/nsgaii_io.R')

# Plot Storage Change -----------------------------------------------------
MAF_KM3 <- 1.2334818375475
CVgwbud <- c2vsim.readGWBUD("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Results/CVground.BUD")
CVgwbud <- c2vsim.cumGWBUD(CVgwbud)
tm <- tm <- seq.Date(from = as.Date(paste0(1921,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")

ggplot(CVgwbud, aes(x = tm, y = ((ES - ES[1])/1000000)*MAF_KM3 ) ) +
  geom_line(size = 1.2, color = "red") + 
  xlim(as.Date(paste0(1921,"/",1,"/1")), as.Date(paste0(2009,"/",10,"/31"))) + 
  labs(title= "Groundwater storage with respect to 1921",
    x = "Time", 
    y = bquote('Storage [' ~km^3* ']'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
    size = 20, hjust = 0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12))


# Plot Urban Demand -------------------------------------------------------
LWBUD <- c2vsim.readLWBUD(filename = "../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Results/CVlandwater.BUD")
LWBUD <- c2vsim.cumLWBUD(LWBUD)
# 
TAF_KM3 <- 0.00123348
# urbanDem <- c2vsim.readDivData("../c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Simulation/CVurbandem.dat",skiplines = 94)
df <- data.frame("Time" = tm, "UrbanDemand" = LWBUD[[2]]$USR/1000*TAF_KM3, "AgDemand" =  LWBUD[[1]]$ASR/1000*TAF_KM3)

ggplot(df, aes(x = tm) ) +
  geom_line(aes(y = AgDemand/10, color = "Ag demand"), size = 1.2) + ##8DB359 
  geom_line(aes(y = UrbanDemand, color = "Urban demand"), size = 1.2) + #F29DAA
  scale_y_continuous(sec.axis = sec_axis(~.*10, name = bquote('Argicultural demand [' ~km^3* ']')))+
  xlim(as.Date(paste0(1921,"/",1,"/1")), as.Date(paste0(2009,"/",10,"/31"))) + 
  labs(title= "Central Valley water supply requirements",
       x = "Time", 
       y = bquote('Urban demand [' ~km^3* ']'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 20, hjust = 0.5) +
  scale_color_manual("", breaks = c("Ag demand", "Urban demand"),  values = c("#8DB359", "#F29DAA"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.15, 0.92),
        legend.text = element_text(size = 20),
        legend.key=element_blank(),
        legend.justification = "center")



# Streamflow hydrographs --------------------------------------------------
TAF_KM3 <- 0.00123348
tm <- tm <- seq.Date(from = as.Date(paste0(1921,"/",10,"/1")),to = as.Date(paste0(2009,"/",9,"/1")),by = "month")
swhyd <- c2vsim.readSWHYD("../RunC2Vsim/Results/CVSWhydBase.out")

df <- data.frame("Time" = tm, "Kern" =  (swhyd[[2]][,1]/1000)*TAF_KM3, 
                 "Kaweah" =  (swhyd[[2]][,421]/1000)*TAF_KM3, 
                 "Kings" =  (swhyd[[2]][,23]/1000)*TAF_KM3)

kern_95 <- as.numeric(quantile((swhyd[[2]][,1]/1000)*TAF_KM3,0.95))
kings_95 <- as.numeric(quantile((swhyd[[2]][,23]/1000)*TAF_KM3,0.95))
kaweah_95 <- as.numeric(quantile((swhyd[[2]][,421]/1000)*TAF_KM3,0.95))

df1 <- data.frame("kernX" = c(0, kern_95, kern_95), "kingsX" = c(0, kings_95, kings_95),
                  "kaweahX" = c(0, kaweah_95, kaweah_95), "Y"  =c(0.95, 0.95, 0))

ggplot(df, aes(x = tm) )+
  geom_line(aes(y = Kings, color = "Kings"), size = 1.1) +
  geom_line(aes(y = Kern, color = "Kern"), size = 1.1) +
  geom_line(aes(y = Kaweah, color = "Kaweah"), size = 1.1) +
  labs(title= "Stream hydrographs for 3 selected river nodes",
       x = "Time", 
       y = bquote('Stream flow [' ~km^3* '/month]'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 20, hjust = 0.5)+
  scale_color_manual("", breaks = c("Kern", "Kaweah", "Kings"),  values = c("#1b9e77", "#d95f02", "#7570b3"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.85, 0.9),
        legend.text = element_text(size = 20),
        legend.key=element_blank(),
        legend.justification = "center")


# Plot cumulative distributions for the streamflows -----------------------
ggplot(df) + 
  stat_ecdf(aes(Kern, color = "Kern"), geom = "step", size = 1.1) +
  stat_ecdf(aes(Kings, color = "Kings"), geom = "step", size = 1.1) +
  stat_ecdf(aes(Kaweah, color = "Kaweah"), geom = "step", size = 1.1) +
  geom_line(data = df1, aes(x = kingsX, y = Y), size = 1, color = "#7570b3", linetype = "dashed") +
  geom_line(data = df1, aes(x = kernX, y = Y), size = 1, color = "#d95f02", linetype = "dashed") +
  geom_line(data = df1, aes(x = kaweahX, y = Y), size = 1, color = "#1b9e77", linetype = "dashed") +
  labs(title= "Empirical Cumulative Distribution Function",
       x = bquote('Stream flow [' ~km^3* '/month]'),
       y = "Cumulative probability",  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 20, hjust = 0.5)+
  scale_color_manual("", breaks = c("Kern", "Kaweah", "Kings"),  values = c("#1b9e77", "#d95f02", "#7570b3"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.9, 0.83),
        legend.text = element_text(size = 20),
        legend.key=element_blank(),
        legend.justification = "center")


# Plot the stream hydrographs without and with diversion ------------------
load("../OptimResults/maxWTminArea/ParetoSolutionsBUD_58946.RData")
tm <- tm <- seq.Date(from = as.Date(paste0(1965, "/", 10, "/1")),to = as.Date(paste0(2009, "/", 9, "/1")),by = "month")
df <- data.frame("Time" = tm, "KernB" = SWHYDbase[[2]][,1]/1000*TAF_KM3, "KernA" = psSWHYD[[15]][[2]][,1]/1000*TAF_KM3,
                 "KingsB" = SWHYDbase[[2]][,23]/1000*TAF_KM3, "KingsA" = psSWHYD[[15]][[2]][,23]/1000*TAF_KM3,
                 "KaweahB" = SWHYDbase[[2]][,421]/1000*TAF_KM3, "KaweahA" = psSWHYD[[15]][[2]][,421]/1000*TAF_KM3)

g1 <- ggplot(df, aes(x = tm)) +
  geom_line(aes(y = KernB), color = 'red', size = 1.2)+
  geom_line(aes(y = KernA), size = 1.2, color = "#1b9e77") +
  labs(y = bquote('Kern [' ~km^3* '/month]'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 15, hjust = 0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 15),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 2.12 , color = "black"),
        axis.line.y = element_line(colour = "grey"))

g2 <- ggplot(df, aes(x = tm)) +
  geom_line(aes(y = KingsB), color = 'red', size = 1.2)+
  geom_line(aes(y = KingsA), size = 1.2, color = "#7570b3")+
  labs(y = bquote('Kings [' ~km^3* '/month]'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 15, hjust = 0.5) + 
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 15),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 2.12 , color = "black"),
        axis.line.y = element_line(colour = "grey"))

g3 <- ggplot(df, aes(x = tm)) +
  geom_line(aes(y = KaweahB), color = 'red', size = 1.2)+
  geom_line(aes(y = KaweahA), size = 1.2, color = "#d95f02")+
  labs(y = bquote('Kaweah [' ~km^3* '/month]'),
       x = "Time", #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 15, hjust = 0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 2.12 , color = "black"),
        axis.title.x = element_text(size = 15, hjust = 0.5, vjust = 1.12),
        axis.line = element_line(colour = "grey"))

grid.arrange(g1,g2,g3,nrow = 3)


# Plot cumulative diversion -----------------------------------------------

df <- data.frame("Time" = tm, "Kern" = cumsum(SWHYDbase[[2]][,1]/1000*TAF_KM3 - psSWHYD[[15]][[2]][,1]/1000*TAF_KM3), 
                 "Kings" = cumsum(SWHYDbase[[2]][,23]/1000*TAF_KM3 - psSWHYD[[15]][[2]][,23]/1000*TAF_KM3),
                 "Kaweah" = cumsum(SWHYDbase[[2]][,421]/1000*TAF_KM3 - psSWHYD[[15]][[2]][,421]/1000*TAF_KM3))

ggplot(df, aes(x = tm)) +
  geom_line(aes(y = Kings, color = "Kings"), size = 1.1) +
  geom_line(aes(y = Kern, color = "Kern"), size = 1.1) +
  geom_line(aes(y = Kaweah, color = "Kaweah"), size = 1.1) + 
  labs(title= "Cumulative diverted water",
       x = "Time", 
       y = bquote('Diversion amount [' ~km^3* ']'),  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
       size = 20, hjust = 0.5) + 
  scale_color_manual("", breaks = c("Kern", "Kaweah", "Kings"),  values = c("#1b9e77", "#d95f02", "#7570b3")) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.1, 0.9),
        legend.text = element_text(size = 20),
        legend.key=element_blank(),
        legend.justification = "center")


# Pareto hyper volume evolution -------------------------------------------
RandInitPop <- nsgaii.readParetoHistory(filename = "../OptimResults/maxWTminArea/paretoHistory_44715.dat", nGen = 400)
InformInitPop <- nsgaii.readParetoHistory(filename = "../OptimResults/maxWTminArea/paretoHistory_44791.dat", nGen = 400)
refPoint <- c(-200000/200, 8000/8)
randHV <- vector(mode = "numeric", length = length(RandInitPop))
infoHV <- vector(mode = "numeric", length = length(InformInitPop))

for (i in 1:length(RandInitPop)) {
  randHV[i] <- nsgaii.calculateHyperVolume(cbind(RandInitPop[[i]][,1]/200, RandInitPop[[i]][,2]/8), refPoint)
  m <- cbind(InformInitPop[[i]][,1]/200, InformInitPop[[i]][,2]/8)
  infoHV[i] <- nsgaii.calculateHyperVolume(m[-which(m[,1] == 0),], refPoint)
}  

df <- data.frame("Iteration" = seq(1,400,1), "RAND" = randHV/10000, "INFO" = infoHV/10000)

ggplot(df, aes(x = Iteration)) +
  geom_line(aes(y = RAND, color = "Random initial population"), size = 1.1) +
  geom_line(aes(y = INFO, color = "Informed random population"), size = 1.1) + 
  labs(title= "Hyper volume evolution",
       x = "Iteration", 
       y = "Hyper Volume",
       size = 20, hjust = 0.5) +
 scale_color_manual("", breaks = c("Random initial population", "Informed random population"),  values = c("#F8766D", "#00BFC4")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.75, 0.80),
        legend.text = element_text(size = 20),
        legend.key=element_blank(),
        legend.justification = "center")


# Plot final Pareto -------------------------------------------------------
ps <- nsgaii.readParetoSolution("../OptimResults/maxWTminArea/paretoSolutions_58946.dat")

df <- data.frame("Area" = ps[[2]][-dim(ps[[2]])[1],2], "WLR" = -ps[[2]][-dim(ps[[2]])[1],1]/1000)

ggplot(df) + 
  geom_point(aes(x = Area, y = WLR), size = 4, color = "#08519C")+
  labs(title= "Pareto Front",
       x = "Cost function (Area)", 
       y = "Environmental function",
       size = 20, hjust = 0.5) + 
  theme(panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2.12 , color = "black"),
        axis.line.y = element_line(colour = "grey"),
        axis.line.x = element_line(colour = "grey"))
  




