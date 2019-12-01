library("pracma")

#' c2vsim.Nodes Reads the Node coordinates of the C2Vsim
#'
#' @param filename is the name of the file. Usually is CVnode.dat
#' @param ND is the number of nodes. The default value corresponds to the coarse
#'   grig version
#' @param Nskip is the number of lins to skip before start reading the list of
#'   nodes, The default value corresponds to the coarse grig version
#'
#' @return a data frame with the following fields: ID, X and Y
#' @export
#'
#' @examples
#' XY <- c2vsim.Nodes("CVnode.dat")
c2vsim.readNodes <- function(filename, ND = 1393, Nskip = 80 ){
  XY <-   read.table(file = filename,
                     header = FALSE, sep = "", skip = Nskip, nrows = ND,
                     quote = "",fill = TRUE,
                     col.names = c("ID", "X", "Y"))
  return(XY)
}


#' c2vsim.Mesh Reads the Mesh element information
#'
#' @param filename is the name of the mesh element file. Usually is
#'   CVelement.dat
#' @param NE is the number of element. The default value corresponds to the
#'   coarse grid version
#' @param Nskip is the number of lins to skip before start reading the list of
#'   nodes, The default value corresponds to the coarse grig version
#'
#' @return a data frame with the following fields: ID, ND1...ND4. If the
#'   elements are triangles then the ND4 index is zero
#' @export
#'
#' @examples
#' MSH <- c2vsim.Mesh(CVelement.dat")
c2vsim.readMesh <- function(filename, NE = 1392, Nskip = 93){
  M <-   read.table(file = filename,
                     header = FALSE, sep = "", skip = Nskip, nrows = NE,
                     quote = "",fill = TRUE,
                    col.names = c("ID", "ND1", "ND2", "ND3", "ND4"))
  return(M)
}


#' c2vsim.readGWBUD Reads the groundwater budget file
#'
#' This file reads the groundwater budget file for the IWFM 302 version and
#' might not work for any modeli
#'
#' @param filename is the filename
#' @param Nsub the number of subregions printed in the file
#' @param Nskip is the number of line between the subregion tables. In fact this
#'   number corresponds to the number of lines from the first up to the second
#'   dashed line. Between the subregions there is an extra empty line that is
#'   added to the Nskip. Therefore before reading files make sure that this is
#'   tha case as this can change other versions
#' @param NtimeSteps is the number of time steps used in the simulation. By
#'   default is set to 1056
#'
#' @return Returns a list of size Nsub data frames with the budget time series.
#'   The columns in the data frame are the following: DP: Deep Percolation BS:
#'   Beginning Storage ES: Ending Storage NDP: Net Deep Percolation GFS: Gain
#'   From Stream R: Recharge GFL: Gain From Lake BI: Boundary Inflow S:
#'   Subsidence SI: Subsurface Irrigation TDO: Tile Drain Outflow P: Pumping
#'   NSI: Net Subsurface Inflow D: Discrepancy CS: Cumulative Subsidence
#'
#' @export
#'
#' @examples
#' We assume that the simulation in the file we are reading starts from Oct-1965.
#' Therefore the number ot time steps is 528. To read the file use
#' GWB <- c2vsim.readGWBUD(file, 21,8, 528)
#'
#' Then you can access the data as for example to plot the Ending storage of the third subregion:
#' plot(GWB[[3]]$ES)
c2vsim.readGWBUD <- function(filename, Nsub = 21, Nskip = 8, NtimeSteps = 1056){
  GWBList <- vector(mode = "list", length = Nsub)
  for (i in 1:Nsub) {
    GWB <- read.table(file =  filename, 
                      header = FALSE, sep = "", skip = Nskip + (i-1)*(NtimeSteps+Nskip+1) , nrows = NtimeSteps,
                      quote = "",fill = TRUE,
                      col.names = c("Time","DP", "BS", "ES", "NDP", "GFS", "R", "GFL", "BI", "S", "SI", "TDO", "P", "NSI", "D", "CS"))
    GWBList[[i]] <- GWB
  }
  return(GWBList)
}

c2vsim.readLWBUD <- function(filename, Nsub = 21, NtimeSteps = 1056, maxchar = 2000){
  out <- vector(mode = "list", length = Nsub)
  
  alllines <- readLines(filename)
  iline <- 1
  for (i in 1:Nsub) {
    ag_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 7))
    ur_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 6))
    ie_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 2))
    
    
    cnt <- 0
    while (TRUE){
      if (strcmp(substr(alllines[iline], 1, 2),"--")){
        cnt <-  cnt + 1
        if (cnt == 2)
          break
      }
      iline = iline + 1
    }
    iline = iline + 1
    for (j in 1:NtimeSteps) {
      t <- strsplit(substr(alllines[iline], 1, maxchar)[[1]], split = " ")[[1]]
      t <- as.numeric((t[which(t != "")])[-1])
      ag_df <- rbind(ag_df, t[1:7])
      ur_df <- rbind(ur_df, t[8:13])
      ie_df <- rbind(ie_df, t[14:15])
      iline <- iline + 1
    }
    
    x <- c("Area", "PCUAW", "ASR", "P", "D", "S", "RU")
    colnames(ag_df) <- x
    x <- c("Area", "USR", "P", "D", "S", "RU")
    colnames(ur_df) <- x
    x <- c("Import", "Export")
    colnames(ie_df) <- x
    
    temp <- vector(mode = "list", length = 3)
    temp[[1]] <- ag_df
    temp[[2]] <- ur_df
    temp[[3]] <- ie_df
    
    out[[i]] <- temp
  }
  return(out)
}


c2vsim.readGWHYD <- function(filename, Nskip = 4, NtimeSteps = 1056, maxChar = 40000){
  out <- vector(mode = "list", length = 2)
  alllines <- readLines(filename)
  
  t <- strsplit(substr(alllines[Nskip+1], 1, maxChar)[[1]], split = " ")[[1]]
  LayerIds <- as.numeric((t[which(t != "")])[-1:-2])
  
  t <- strsplit(substr(alllines[Nskip+2], 1, maxChar)[[1]], split = " ")[[1]]
  NodeIds <- as.numeric((t[which(t != "")])[-1:-2])
  
  M <- matrix(data = NA, nrow = NtimeSteps, ncol = length(NodeIds))
  for (i in 1:NtimeSteps) {
    t <- strsplit(substr(alllines[Nskip+3+i], 1, maxChar)[[1]], split = " ")[[1]]
    M[i,] <- as.numeric(t[which(t != "")][-1])
    
  }
  
  out[[1]] <- NodeIds
  out[[2]] <- M
  return(out)
}


c2vsim.readSWHYD <- function(filename, Nskip = 5, NtimeSteps = 1056, maxChar = 7000){
  out <- vector(mode = "list", length = 2)
  alllines <- readLines(filename)
  t <- strsplit(substr(alllines[Nskip+1], 1, maxChar)[[1]], split = " ")[[1]]
  NodeIds <- as.numeric((t[which(t != "")])[-1:-2])
  M <- matrix(data = NA, nrow = NtimeSteps, ncol = length(NodeIds))
  for (i in 1:NtimeSteps) {
    t <- strsplit(substr(alllines[Nskip+1+i], 1, maxChar)[[1]], split = " ")[[1]]
    M[i,] <- as.numeric(t[which(t != "")][-1])
    
  }
  
  out[[1]] <- NodeIds
  out[[2]] <- M
  return(out)
}

c2vsim.readDiversionBUD <- function(filename, Nsub = 21, NtimeSteps = 1056, maxChar = 2000, nExpectedDivs = 246){
  out <- vector(mode = "list", length = 4)
  surfDelMat <- matrix(data = NA, nrow = NtimeSteps, ncol = nExpectedDivs*2)
  divMat <- matrix(data = NA, nrow = NtimeSteps, ncol = nExpectedDivs*2)
  SWDrn <- matrix(data = NA, nrow = 1, ncol = nExpectedDivs)
  Divrn <- matrix(data = NA, nrow = 1, ncol = nExpectedDivs)
  
  alllines <- readLines(filename)
  iline <- 1
  while (TRUE){
    if (strcmp(substr(alllines[iline], 1, 2),"--"))
      break
    iline = iline + 1
  }
  iline = iline + 1
  for (i in 1:Nsub) {
    # Read the diversion node ids
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    divNode <- as.numeric((t[which(t != "")])[-1:-2])
    # Read the river nodes
    iline <- iline + 1
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    rivNode <- as.numeric(t[which(t != "")][-1:-3])
    # Read whether is SURFACE WATER DELIVERIES or DIVERSIONS
    iline <- iline + 1
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    t <- t[which(t != "")]
    idp <- which(t == "(+)")
    idm <- which(t == "(-)")
    
    idpmat <- Reshape(rbind(divNode[idp]*2-1,divNode[idp]*2), 1, length(idp)*2)
    idpt <- Reshape(rbind(idp*2-1,idp*2),1,length(idp)*2)
    
    idmmat <- Reshape(rbind(divNode[idm]*2-1,divNode[idm]*2), 1, length(idm)*2)
    idmt <- Reshape(rbind(idm*2-1,idm*2),1,length(idm)*2)
    SWDrn[divNode[idp]] <- rivNode[idp]
    Divrn[divNode[idm]] <- rivNode[idm]
    
    iline <- iline + 1
    for (j in 1:NtimeSteps) {
      iline <- iline + 1
      t <- substr(alllines[iline], 1, maxChar)
      t <- chartr(old = "(", new = " ", t)
      t <- chartr(old = ")", new = " ", t)
      t <- strsplit(t[[1]], split = " ")[[1]]
      t <- as.numeric(t[which(t != "")][-1])
      # put in the matrix the Surface water deliviries
      surfDelMat[j,idpmat] <- t[idpt]
      # and the diversion
      divMat[j,idmmat] <- t[idmt]
    }
    if (i != Nsub){
      while (TRUE){
        if (strcmp(substr(alllines[iline], 1, 2),"--"))
          break
        iline = iline + 1
      }
      iline = iline + 1
    }
  }
  out[[1]] <- SWDrn
  out[[2]] <- Divrn
  out[[3]] <- surfDelMat
  out[[4]] <- divMat
  return(out)
}


#' c2vsim.cumBud summarizes the groundwater budget for the entire domain or the
#' specified ids
#'
#' @param GWB This is a List of data frames with the budget terms for each
#'   subregion. You can obtain this as the output of the c2vsim.readGWBUD
#'   function
#' @param ids the ids of the subregion you want to summarize
#'
#' @return A data frame with the budget numbers corresponding to the entire area
#' @export
#'
#' @examples
#' GWBUDALL <- c2vsim.cumBUD(GWB, ids = c(5,10,21))
c2vsim.cumGWBUD <- function(GWB, ids = NA){
  if (is.na(ids)){
    ids <- 1:length(GWB)
  }
  GWBALL <- GWB[[1]]
  GWBALL[,-1] = 0
  for (i in ids) {
    GWBALL[,-1] = GWBALL[,-1] + GWB[[i]][,-1] 
  }
  return (GWBALL)
}

c2vsim.cumLWBUD <- function(LWB, ids = NA){
  if (is.na(ids)){
    ids <- 1:length(LWB)
  }
  AGALL <- LWB[[1]][[1]]
  AGALL[,] = 0
  URALL <- LWB[[1]][[2]]
  URALL[,] = 0
  IEALL <- LWB[[1]][[3]]
  IEALL[,] = 0
  for (i in ids) {
    AGALL = AGALL + LWB[[i]][[1]] 
    URALL = URALL + LWB[[i]][[2]]
    IEALL = IEALL + LWB[[i]][[3]]
  }
  out <- vector(mode = "list", length = 3)
  out[[1]] <- AGALL
  out[[2]] <- URALL
  out[[3]] <- IEALL
  
  return(out)
}


#' c2vsim.readDivSpec reads the Diversion specification file.
#'
#' The way this function work is to read the entire file into the memory and
#' then process it line by line. However this maybe problematic for large files
#' see
#' (https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line)
#'
#'
#'
#' @param filename is the name of the file to read
#'
#' @return a List with the following content: 1 -> a vector with the values
#'   NRDV, NDIVS, FACTX, FACTY 2 -> a matrix NRDV x 14 with the Surface Water
#'   Diversion Specifications 3 -> a List with length NRDV with matrices Nelem x
#'   2 that describe the Recharge zone for each diversion point 4 -> a matrix
#'   NDIVS x 6 with the Bypass Configuration Specifications 5 -> a list of NDIVS
#'   matrices with the rating tables for those bypasses that are defined. 6 -> a
#'   List with length NDIVS with matrices Nelem x 2 that describe the Seepage
#'   locations for bypass canals
#' @export
#'
#' @examples
c2vsim.readDivSpec <- function(filename){
  DIVSPEC = vector(mode = "list", length = 7)
  headers <- vector(mode = "numeric", length = 4)
  con <- file(filename, open = "r")
  fileLines <- readLines(con)
  close(con)
  section <- 1
  irdv <- 0
  i <- 1
  while (TRUE) {
    if (i == length(fileLines))
      break
    if (strcmp(substr(fileLines[i],1,1),"C")){
      i <- i + 1
      next
    }
    
    if (section == 1){
      NRDV <- scan(text = fileLines[i], n=1, quiet = TRUE)
      headers[1] <- NRDV
      section  <-  section + 1
      RDV <- matrix(data = NA,nrow = NRDV, ncol = 14)
      RDVELEM <- vector(mode = "list", length = NRDV)
      RDVnames <- c()
      i <- i + 1
    }
    else if (section == 2){
      div_name <- strsplit (substr(fileLines[i-2],4,1000), "\t")[[1]][1]
      RDVnames <- c(RDVnames, div_name)
      irdv <-  irdv + 1
      RDV[irdv,] <- scan(text = fileLines[i], n=14, quiet = TRUE)
      i <- i + 1
      if (irdv == NRDV){
        section <-  section + 1
        irdv <- 0
      }
    }
    else if (section == 3){
      #print(i)
      irdv <-  irdv + 1
      temp <- scan(text = fileLines[i], n = 4, quiet = TRUE)
      i <- i + 1
      elem <- matrix(data = NA, nrow = max(c(temp[2],1)), ncol = 2)
      elem[1,] <- temp[3:4] 
      if (temp[2] >= 2){
        for (j in 2:temp[2]) {
          elem[j,] <- scan(text = fileLines[i], n = 2, quiet = TRUE)
          i <- i + 1
        }
      }
      RDVELEM[[irdv]] <- elem
      if (irdv == NRDV){
        section <-  section + 1
        irdv <- 0
      }
    }
    else if (section == 4){
      NDIVS <- scan(text = fileLines[i], n = 1, quiet = TRUE)
      headers[2] <- NDIVS
      FACTX <- scan(text = fileLines[i+1], n = 1, quiet = TRUE)
      headers[3] <- FACTX
      #TUNITX <- scan(text = fileLines[i+2], n = 1, quiet = TRUE)
      FACTY <- scan(text = fileLines[i+3], n = 1, quiet = TRUE)
      headers[4] <- FACTY
      #TUNITY <- scan(text = fileLines[i+4], n = 1, quiet = TRUE)
      BYPS <- matrix(data = NA, nrow = NDIVS, ncol = 6)
      BYPSRT <- vector(mode = "list", length = NDIVS)
      BYPSELEM <- vector(mode = "list", length = NDIVS)
      i <- i + 5
      section <-  section + 1
    }
    else if (section == 5){
      irdv <-  irdv + 1
      BYPS[irdv,] <- scan(text = fileLines[i], n = 6, quiet = TRUE)
      i <- i + 1
      if (BYPS[irdv,4] < 0){
        RT <- matrix(data = NA, nrow = abs(BYPS[irdv,4]), ncol = 2)
        for (j in 1:abs(BYPS[irdv,4])) {
          RT[j,] <- scan(text = fileLines[i], n = 2, quiet = TRUE)
          i <- i + 1
        }
        BYPSRT[[irdv]] <- RT
      }
      if (irdv == NDIVS){
        section <-  section + 1
        irdv <- 0
      }
    }
    else if (section == 6){
      irdv <-  irdv + 1
      temp <- scan(text = fileLines[i], n = 4, quiet = TRUE)
      i <- i + 1
      elem <- matrix(data = NA, nrow = max(c(temp[2],1)), ncol = 2)
      elem[1,] <- temp[3:4] 
      if (temp[2] >= 2){
        for (j in 2:temp[2]) {
          elem[j,] <- scan(text = fileLines[i], n = 2, quiet = TRUE)
          i <- i + 1
        }
      }
      BYPSELEM[[irdv]] <- elem
      if (irdv == NDIVS){
        section <-  section + 1
        irdv <- 0
        break
      }
    }
  }
  DIVSPEC[[1]] <- headers
  DIVSPEC[[2]] <- RDV
  DIVSPEC[[3]] <- RDVELEM
  DIVSPEC[[4]] <- BYPS
  DIVSPEC[[5]] <- BYPSRT
  DIVSPEC[[6]] <- BYPSELEM
  DIVSPEC[[7]] <- RDVnames
  return(DIVSPEC)
}


#' c2vsim.readDivData Reads the Diversion Data file The default values
#' correspond to the Coars version of C2Vsim input file
#'
#' @param filename The name of the file
#' @param skiplines the number of files to skip up to the data
#' @param NtimeSteps the number of steps of the diversions timeseries
#'
#' @return a table with size NtimeSteps x the number of time series in file
#' @export
#'
#' @examples
c2vsim.readDivData <- function(filename, skiplines = 376, NtimeSteps = 1056){
  
  DivData <- read.table(file =  filename, 
                      header = FALSE, sep = "", skip = skiplines , nrows = NtimeSteps,
                      quote = "",fill = TRUE
                      )
  return(DivData)
}

c2vsim.writeDivSpec <- function(filename, DivSpec){
  con <- file(filename, open = "w")
  write(DivSpec[[1]][1], file = con)
  
  write.table(DivSpec[[2]], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  for (i in 1:length(DivSpec[[3]])) {
    write(c(DivSpec[[2]][i,1], dim(DivSpec[[3]][[i]])[1], DivSpec[[3]][[i]][1,]), file = con, sep = " ")
    if (dim(DivSpec[[3]][[i]])[1] > 2)
      write.table(DivSpec[[3]][[i]][-1,], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    else if (dim(DivSpec[[3]][[i]])[1] == 2)
      write(DivSpec[[3]][[i]][2,], file = con, sep = " ")
  }
  
  write(DivSpec[[1]][2], file = con)
  write(DivSpec[[1]][3], file = con)
  write("1min", file = con)
  write(DivSpec[[1]][4], file = con)
  write("1min", file = con)
  write("C", file = con)
  
  for (i in 1:dim(DivSpec[[4]])[1]){
    write(DivSpec[[4]][i,], file = con, sep = " ", ncolumns = 6, append = TRUE)
    if (DivSpec[[4]][i,4]<0){
      for (j in 1:dim(DivSpec[[5]][[i]])[1]) {
        temp <- sprintf("%.2f", DivSpec[[5]][[i]][j,])
        write(paste(temp[1], temp[2]), file = con)
      }
      #write.table(DivSpec[[5]][[i]], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  for (i in 1:length(DivSpec[[6]])) {
    write(c(DivSpec[[4]][i,1], dim(DivSpec[[6]][[i]])[1], DivSpec[[6]][[i]][1,]), file = con, sep = " ")
    if (dim(DivSpec[[6]][[i]])[1] > 2)
      write.table(DivSpec[[6]][[i]][-1,], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    else if (dim(DivSpec[[6]][[i]])[1] == 2)
      write(DivSpec[[6]][[i]][2,], file = con, sep = " ")
  }
  
  close(con)
}

c2vsim.writeDivData <- function(filename, data, NCOLDV = 265, FACTDV = 43560000.0, NSPDV = 1, NFQDV = 0, DSSFL= "", sep = " "){
  con <- file(filename, open = "w")
  write(NCOLDV,file = con)
  write(FACTDV,file = con)
  write(NSPDV,file = con)
  write(NFQDV,file = con)
  write(DSSFL,file = con)
  write.table(data, file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  #write(data, file = con, ncolumns = dim(data)[2])
  close(con)
}

# At the moment this reads only the OPTION 2 (for Aquifer Parameter Definition)
c2vsim.readParam <- function(filename, nSkip = 267, Nnodes = 1393, Nlay = 3){
  PRM <- vector(mode = "list", length = Nlay)
  df <- data.frame(matrix(nrow = Nnodes, ncol = 10))
  colnames(df) <- c("PKH", "PS", "PN", "PV", "PL", "SCE", "SCI",  "DC", "DCMIN", "HC")
  for (i in 1:Nlay) {
    PRM[[i]] <- df
  }
  
  allLines <- readLines(filename)
  index <- seq(from = 1, to = 1393*Nlay, by = Nlay) + nSkip
  ii <- 0
  for (i in index) {
    ii <- ii + 1
    PRM[[1]][ii,] <- scan(text = allLines[i], n = 11, quiet = TRUE)[-1]
    for (j in 2:Nlay) {
      PRM[[j]][ii,] <- scan(text = allLines[i+j-1], n = 10, quiet = TRUE)
    }
  }
  
  return(PRM)
}

c2vsim.readStrat <- function(filename, nSkip = 92, Nnodes = 1393, Nlay = 3){
  matlabels <- c("ID", "ELV")
  for (i in 1:Nlay) {
    matlabels <- c(matlabels, paste0("L", i, "1"), paste0("L", i, "2"))
  }
  matlabels <- c(matlabels, "dummy")
  ELEV <-   read.table(file = filename,
                    header = FALSE, sep = "", skip = nSkip, nrows = Nnodes,
                    quote = "",fill = TRUE,
                    col.names = matlabels)
  return(ELEV[,-dim(ELEV)[2]])
  
}
