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
c2vsim.readMesh <- function(filename, NE = 1392, Nskip = 93, Ncols = 5){
  if (Ncols == 5){
    M <-   read.table(file = filename,
                       header = FALSE, sep = "", skip = Nskip, nrows = NE,
                       quote = "",fill = TRUE,
                      col.names = c("ID", "ND1", "ND2", "ND3", "ND4"))
  }
  else if (Ncols == 6) {
    M <-   read.table(file = filename,
                      header = FALSE, sep = "", skip = Nskip, nrows = NE,
                      quote = "",fill = TRUE,
                      col.names = c("ID", "ND1", "ND2", "ND3", "ND4", "SUB"))
  }
  return(M)
}


#' c2vsim.readGWBUD Reads the groundwater budget file
#'
#' This file reads the groundwater budget file for the IWFM 302 version and
#' might not work for any modeli
#'
#' @param filename is the filename
#' @param Nsub the number of subregions printed in the file
#' @param Nskip is the number of lines between the subregion tables. In fact this
#'   number corresponds to the number of lines from the first up to the second
#'   dashed line. Between the subregions there is an extra empty line that is
#'   added to the Nskip. Therefore before reading files make sure that this is
#'   tha case as this can change other versions
#' @param NtimeSteps is the number of time steps used in the simulation. By
#'   default is set to 1056
#'   
#' @param CG set this to true when reading from coarse grid version. 
#'   set it False when reading from fine grid
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
c2vsim.readGWBUD <- function(filename, Nsub = 21, Nskip = 8, NtimeSteps = 1056, CG = TRUE){
  if (length(Nskip) == 1){
    Nskip <- rep(Nskip, Nsub)
  }
  if (length(Nskip) != Nsub){
    stop("length(Nskip) != Nsub or length(Nskip) != 1")
  }
  if (CG){
    fieldnames <- c("Time","DP", "BS", "ES", "NDP", "GFS", "R", "GFL", "BI", "S", "SI", "TDO", "P", "NSI", "D", "CS")
  }
  else{
    fieldnames <- c("Time","P", "BS", "ES", "DP", "GFS", "R", "GFL", "BI", "S", "SI", "TDO", "PMP", "ORZ", "NSI", "D", "CS")
  }
  GWBList <- vector(mode = "list", length = Nsub)
  for (i in 1:Nsub) {
    GWB <- read.table(file =  filename, 
                      header = FALSE, sep = "", skip = sum(Nskip[1:i]) + NtimeSteps*(i-1) , nrows = NtimeSteps,
                      quote = "",fill = TRUE,
                      col.names = fieldnames)
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
    tm <- c()
    
    
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
      tm <- c(tm,t[1])
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
    
    out[[i]] <- list(Time = tm, AG = ag_df, UR = ur_df, IE = ie_df)
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
  
  out <- list(LayerIds = LayerIds, NodeIds = NodeIds, GWHyd = M)
  return(out)
}


c2vsim.readSWHYD <- function(filename, Nskip = 5, NtimeSteps = 1056, maxChar = 7000){
  alllines <- readLines(filename)
  t <- strsplit(substr(alllines[Nskip+1], 1, maxChar)[[1]], split = " ")[[1]]
  NodeIds <- as.numeric((t[which(t != "")])[-1:-2])
  M <- matrix(data = NA, nrow = NtimeSteps, ncol = length(NodeIds))
  tm <- c()
  for (i in 1:NtimeSteps) {
    t <- strsplit(substr(alllines[Nskip+1+i], 1, maxChar)[[1]], split = " ")[[1]]
    tm <- c(tm,t[1])
    M[i,] <- as.numeric(t[which(t != "")][-1])
    
  }

  out <- list(Time = tm, NodeIds = NodeIds, SWHyd = M)
  return(out)
}

c2vsim.readDiversionBUD <- function(filename, Nsub = 21, NtimeSteps = 1056, maxChar = 2000, nExpectedDivs = 246){
  out <- vector(mode = "list", length = Nsub)
  alllines <- readLines(filename)
  iline <- 1
  for (i in 1:Nsub) {
    # loop until we find the first dashed line
    while (TRUE){
      if (strcmp(substr(alllines[iline], 1, 2),"--"))
        break
      iline = iline + 1
    }
    # Read the diversion ids
    iline = iline + 1
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    divNode <- as.numeric((t[which(t != "")])[-1:-2])
    # Read the river node ids
    iline <- iline + 1
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    rivNode <- as.numeric(t[which(t != "")][-1:-3])
    # Separate the SURFACE WATER DELIVERIES from the DIVERSIONS
    iline <- iline + 1
    t <- strsplit(substr(alllines[iline], 1, maxChar)[[1]], split = " ")[[1]]
    t <- t[which(t != "")]
    idp <- which(t == "(+)")
    idm <- which(t == "(-)")
    
    SWD <- matrix(data = NA, nrow = NtimeSteps, ncol = length(idp))
    DIV <- matrix(data = NA, nrow = NtimeSteps, ncol = length(idm))
    SWDshrtg <- matrix(data = NA, nrow = NtimeSteps, ncol = length(idp))
    DIVshrtg <- matrix(data = NA, nrow = NtimeSteps, ncol = length(idm))
    tm <- c()
    
    # read the timeseries
    iline <- iline + 1
    for (j in 1:NtimeSteps){
      iline <- iline + 1
      t <- substr(alllines[iline], 1, maxChar)
      t <- chartr(old = "(", new = " ", t)
      t <- chartr(old = ")", new = " ", t)
      t <- strsplit(t[[1]], split = " ")[[1]]
      tm <- c(tm, t[1])
      t <- as.numeric(t[which(t != "")][-1])
      t <- pracma::Reshape(t,2,length(t)/2)
      SWD[j,] <- t[1,idp]
      DIV[j,] <- t[1,idm]
      SWDshrtg[j,] <- t[2,idp]
      DIVshrtg[j,] <- t[2,idm]
    }
    out[[i]] <- list(Time = tm, 
                     SWD = list(DivId = divNode[idp], RivId = rivNode[idp], SWD = SWD, SHRT = SWDshrtg),
                     DIV = list(DivId = divNode[idm], RivId = rivNode[idm], SWD = DIV, SHRT = DIVshrtg))
  }
  return(out)
  
  
  
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
  ALL <- LWB[[1]]
  ALL$AG[] <- 0
  ALL$UR[] <- 0
  ALL$IE[] <- 0
  for (i in ids) {
    ALL$AG = ALL$AG + LWB[[i]]$AG
    ALL$UR = ALL$UR + LWB[[i]]$UR
    ALL$IE = ALL$IE + LWB[[i]]$IE
  }
  return(ALL)
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
  # DIVSPEC[[1]] <- headers
  # DIVSPEC[[2]] <- RDV
  # DIVSPEC[[3]] <- RDVELEM
  # DIVSPEC[[4]] <- BYPS
  # DIVSPEC[[5]] <- BYPSRT
  # DIVSPEC[[6]] <- BYPSELEM
  # DIVSPEC[[7]] <- RDVnames
  DIVSPEC <- list(RDVnames = RDVnames,
                  headers = headers,
                  RDV = RDV,
                  RDVELEM = RDVELEM,
                  BYPS = BYPS,
                  BYPSRT = BYPSRT,
                  BYPSELEM = BYPSELEM)
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
  write(DivSpec$headers[1], file = con)
  
  write.table(DivSpec$RDV, file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  for (i in 1:length(DivSpec$RDVELEM)) {
    write(c(DivSpec$RDV[i,1], dim(DivSpec$RDVELEM[[i]])[1], DivSpec$RDVELEM[[i]][1,]), file = con, sep = " ")
    if (dim(DivSpec$RDVELEM[[i]])[1] > 2)
      write.table(DivSpec$RDVELEM[[i]][-1,], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    else if (dim(DivSpec$RDVELEM[[i]])[1] == 2)
      write(DivSpec$RDVELEM[[i]][2,], file = con, sep = " ")
  }
  
  write(DivSpec$headers[2], file = con)
  write(DivSpe$headers[3], file = con)
  write("1min", file = con)
  write(DivSpec$headers[4], file = con)
  write("1min", file = con)
  write("C", file = con)
  
  for (i in 1:dim(DivSpec$BYPS)[1]){
    write(DivSpec$BYPS[i,], file = con, sep = " ", ncolumns = 6, append = TRUE)
    if (DivSpec$BYPS[i,4]<0){
      for (j in 1:dim(DivSpec$BYPSRT[[i]])[1]) {
        temp <- sprintf("%.2f", DivSpec$BYPSRT[[i]][j,])
        write(paste(temp[1], temp[2]), file = con)
      }
      #write.table(DivSpec[[5]][[i]], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  for (i in 1:length(DivSpec$BYPSELEM)) {
    write(c(DivSpec$BYPS[i,1], dim(DivSpec$BYPSELEM[[i]])[1], DivSpec$BYPSELEM[[i]][1,]), file = con, sep = " ")
    if (dim(DivSpec$BYPSELEM[[i]])[1] > 2)
      write.table(DivSpec$BYPSELEM[[i]][-1,], file = con, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    else if (dim(DivSpec$BYPSELEM[[i]])[1] == 2)
      write(DivSpec$BYPSELEM[[i]][2,], file = con, sep = " ")
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


#' c2vsim.readHeadALL
#' Reads the heads the groundwater head values for all layers all time steps
#' The default values correspond to the fine grid version.
#'
#' @param filename is the name of the file with the groundwater head values.
#' @param nNode is the total number of nodes in the simulation.
#' @param nLay is the Number of layers.
#' @param NtimeSteps is the number of monthly time steps.
#' @param nSkip is the number of lines to skip before the first time step values are printed
#'
#' @return a list of lists. The first item in the list is a vector with the node ids.
#' The second lement in the list is a data frame with the time stamps
#' The third length NtimeSteps where each element in the list is a matrix 
#' nNode x nLay.
#' @export
#'
#' @examples
c2vsim.readHeadALL <- function(filename, nNode = 30179, nLay = 4, NtimeSteps = 505, nSkip = 6, quiet = FALSE){
  # A list to hold the outputs
  out <- vector(mode = "list", length = 3)
  # A list to hold the head values for each time step
  Hall <- vector(mode = "list", length = NtimeSteps)
  # a Data frame to hold the dates
  timedf <- data.frame(matrix(data = NA, nrow = NtimeSteps, ncol = 3), row.names = NULL)
  colnames(timedf) <- c("Y", "M", "D")
  # Read all lines
  allLines <- readLines(filename)
  # allocate memory for head data for one time step
  H <- matrix(data = NA, nrow = nNode, ncol = nLay)
  # Specify the maximum number of characters to read from each line
  maxChar <- nNode*13 + 25
  
  # Read the node IDS
  tmp <- strsplit(substr(allLines[nSkip], 1, maxChar)[[1]], split = " ")[[1]]
  tmp <- tmp[which(tmp != "")]
  out[[1]] <- as.numeric(tmp[c(-1,-2)])
  
  iln <- nSkip+1
  itm <- 1
  while(T){
    tmp <- strsplit(substr(allLines[iln], 1, maxChar)[[1]], split = " ")[[1]]
    if (length(tmp) == 0){
      iln <- iln + 1
      next
    }
    tmp <- tmp[which(tmp != "")]
    if (!quiet){
      print(tmp[1])
    }
    timedf[itm,] <- as.numeric(strsplit( substr(tmp[1],1,10), "/")[[1]])[c(3,1,2)]
    H[,1] <- as.numeric(tmp[-1])
    for (j in 2:nLay) {
      iln <- iln + 1
      tmp <- strsplit(substr(allLines[iln], 1, maxChar)[[1]], split = " ")[[1]]
      tmp <- tmp[which(tmp != "")]
      H[,j] <- as.numeric(tmp)
    }
    Hall[[itm]] <- H
    H[,] = NA
    itm <- itm + 1  
    iln <- iln + 1
    if (itm > NtimeSteps || iln > length(allLines)){
      break
    }
  }
  out[[2]] <- timedf
  out[[3]] <- Hall
  return(out)
}

#' c2vsim.avHead 
#' Calculates Average head values for selected time period
#'
#' @param HeadList This is list of head. It is the third element of the output list from the 
#' c2vsim.readHeadALL function
#' @param ids the ids of the months to consider in the averaging calculations
#'
#' @return A matrix [nNodes x nLay] with the average Head values
#' @export
#'
#' @examples
c2vsim.avHead <- function(HeadList, ids){
  Hav <- HeadList[[1]]
  Hav[,] <- 0
  
  for (i in ids) {
    Hav <- Hav + HeadList[[i]]
  }
  return(Hav/length(ids))
}

c2vsim.calcVelocityField <- function(Head, XY, MSH, K, POR){
  
}

c2vsim.units.m2ft <- function(){
  return(3.280839895)
}

c2vsim.units.ft2m <- function(){
  return(1/c2vsim.units.m2ft())
}


#' c2vsim.read.LandUse
#' Reads the land use time series. For the IWFM version 3 it reads the CVland use. 
#' For the IWFM 15 it reads the files that have the _area suffix.
#' The default values correspond to the coarse grid case version 3
#'
#' @param filename The name of the file to read
#' @param NtimeSteps The number of time steps in the file
#' @param Nelem The number of elements
#' @param Ninfo The number of columns to read
#' @param Nskip The number of lines to skip before the timeseries
#' @param maxChar The number of characters to read for each line
#' 
#'
#' @return
#' @export
#'
#' @examples
c2vsim.read.LandUse <- function(filename, NtimeSteps = 88, Nelem = 1392, Ninfo = 5, colNames = NA, Nskip = 95, maxChar = 1000){
  if (is.na(colNames)){
    colNames = paste0("V", c(1:5))
  }
  # Prepare output data
  dataOut = vector(mode = "list", length = NtimeSteps)
  outALL = vector(mode = "list", length = 2)
  # Allocate for dates
  timedf <- data.frame(matrix(data = NA, nrow = NtimeSteps, ncol = 3), row.names = NULL)
  colnames(timedf) <- c("Y", "M", "D")
  
  # Read the entire file
  allLines <- readLines(filename)
  cnt_ln <- Nskip + 1
  
  for (i in 1:NtimeSteps) {
    
    tmp <- strsplit(substr(allLines[cnt_ln], 1, maxChar)[[1]], " ")[[1]]
    tmp <- tmp[which(tmp != "")]
    print(tmp)
    
    timedf[i,] <- as.numeric(strsplit( substr(tmp[1],1,10), "/")[[1]])[c(3,1,2)]
    df <- data.frame(matrix(data = NA, nrow = Nelem, ncol = Ninfo), row.names = NULL)
    colnames(df) <-colNames
    df[1,] <- as.numeric(tmp[-1])
    cnt_ln <- cnt_ln + 1
    for (j in 2:Nelem) {
      tmp <- strsplit(substr(allLines[cnt_ln], 1, maxChar)[[1]], " ")[[1]]
      tmp <- tmp[which(tmp != "")]
      df[j,] <- as.numeric(tmp)
      cnt_ln <- cnt_ln + 1
    }
    
    dataOut[[i]] <- df
  }
  outALL[[1]] <- timedf
  outALL[[2]] <- dataOut
  return(outALL)
}

