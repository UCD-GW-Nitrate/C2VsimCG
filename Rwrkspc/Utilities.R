# requires library("lubridate")
myfnc.writeTSMatrix2JS <- function(fileName, data, varName = "jsvar", sy = 2000, sm = 1){
  cat("var ", varName, " = [\n", file = fileName)
  isy <- sy
  ism <- sm
  
  for (i in 1:dim(data)[2]){
    sy <- isy
    sm <- ism
    cat("\t{ id: ", i, ",\n\t Values: [\n\t\t" , sep = "",  file = fileName, append = TRUE)
    
    for (j in 1:dim(data)[1]) {
      ndays <- as.integer(days_in_month(as.Date(paste0(sy,"-",sm,"-1"))))
      temp <- paste0("{x: new Date(", sy, ",", sm,",",ndays,")", ", y: ",  sprintf("%.2f", data[j,i]),"}")
      cat(temp, sep = "",  file = fileName, append = TRUE)
      
      if (j != dim(data)[1])
        cat(",\n\t\t", sep = "",  file = fileName, append = TRUE)
      
      sm = sm + 1
      if (sm >= 13){
        sm = 1
        sy = sy + 1
      }
    }
    cat("]}", file = fileName, append = TRUE)
    
    if (i != dim(data)[2])
      cat(",\n" , sep = "",  file = fileName, append = TRUE)
    else
      cat("\n" , sep = "",  file = fileName, append = TRUE)
    
  }
  cat("];", file = fileName, append = TRUE)
}

shapefile.coords <- function(sh, id){
  coords <- slot( slot( slot(sh, "polygons")[[id]], "Polygons" )[[1]] , "coords")
  return(coords)
}

calcCentroid <- function (coords){
  n <- dim(coords)[1]
  if (sum(coords[1,] - coords[n,]) == 0)
    coords <- coords[-n,]
  
  return(apply(coords,2,sum)/dim(coords)[1])
}


# library(stringr)
extractCommonRoot <- function(str1, str2){
  n <- min(c(str_length(str1), str_length(str2)))
  for (i in 1:n) {
    if (substr(str1,1,i) != substr(str2,1,i))
      break
  }
  return(substr(str1,1,i-1))
}

readHistogramLine <- function(s){
  s <- gsub("\"","", s)
  s <- substr(s,2,nchar(s)-2)
  s <- strsplit(s,",")
  out <- vector(mode = "list", length = length(s[[1]]))
  m <- matrix(data = NA, nrow = length(s[[1]]), ncol = 2)
  for (i in 1:length(s[[1]])) {
    t <- strsplit(s[[1]][i],":")
    m[i,] <- as.numeric(t[[1]])
  }
 
  return(m)
}

color_elem_price <- function(price, min_price, max_price){
  if(is.na(price)){
    return (paste0('#',paste(as.hexmode(c(0,0,0)),collapse = '')))
  }
  if (price > max_price){
    return (paste0('#',paste(as.hexmode(c(127,0,0)),collapse = '')))
  }
  
  if (price < min_price){
    return (paste0('#',paste(as.hexmode(c(255,247,236)),collapse = '')))
  }
  
  u = (price - min_price)/(max_price - min_price)
  r = round(255*(1-u) + 215*u)
  g = round(247*(1-u) + 48*u)
  b = round(236*(1-u) + 31*u)
  return (paste0('#',paste(as.hexmode(c(r,g,b)),collapse = '')))
}
