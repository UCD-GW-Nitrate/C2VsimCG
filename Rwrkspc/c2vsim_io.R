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
c2vsim.Nodes <- function(filename, ND = 1393, Nskip = 80 ){
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
c2vsim.Mesh <- function(filename, NE = 1392, Nskip = 93){
  M <-   read.table(file = filename,
                     header = FALSE, sep = "", skip = Nskip, nrows = NE,
                     quote = "",fill = TRUE,
                     col.names = c("ID", "ND1", "ND2", "ND3", "ND4"))
  return(M)
  
}