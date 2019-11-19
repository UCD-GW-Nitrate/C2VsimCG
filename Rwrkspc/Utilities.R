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