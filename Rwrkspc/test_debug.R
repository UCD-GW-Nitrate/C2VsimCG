elem_price <- vector(mode = "numeric", length = length(m))
for (iel in 1:length(m)) {
  print(paste("element:", iel))
  crop_el <- LU_elem_crop_percent[[iel]]
  ncrops <- dim(crop_el)[1]
  crop_el <- cbind(crop_el, vector(mode = "numeric", length = ncrops))
  for(icrop in 1:ncrops){
    davis_class <- ca_cdl_xwalk$davis_class[ which(ca_cdl_xwalk$class == crop_el[icrop])]
    if (isempty(davis_class)){
      print(paste("Not found", crop_el[icrop]))
      next
    }
    elem_with_class <- which(landPriceTable$davis_class_any == davis_class)
    if (isempty(elem_with_class)){
      print(paste(davis_class, "not found"))
      next
    }
    elem_with_class <- unique(landPriceTable$ie[elem_with_class])
    # Find the distance between the elements that contain price for that class with the current element
    dst_ml <- sqrt((elem_bc[iel,1] - elem_bc[elem_with_class,1])^2 + (elem_bc[iel,2] - elem_bc[elem_with_class,2])^2)*(1/1609.34)
    
    if (length(which(dst_ml<20)) > 1){
      # If there are more than one close elements then choose the one with the closet depth
      closest_elem <- sort(dst_ml,index.return=T)
      if (closest_elem$x[1] == 0){
        selected_elem <- elem_with_class[closest_elem$ix[1]]
      }
      else{
        close_picked_elem <- elem_with_class[closest_elem$ix[which(closest_elem$x<20)]]
        depth_elem_list <- Depth_per_elem[close_picked_elem]
        depth_current_el <- Depth_per_elem[iel]
        depth_sort <- sort(sqrt((depth_elem_list-depth_current_el)^2), index.return=T)
        selected_elem <- close_picked_elem[depth_sort$ix[1]]
      }
    }
    else{
      # if there is only one or none close elements choose the closest
      closest_elem <- sort(dst_ml,index.return=T)
      selected_elem <- elem_with_class[closest_elem$ix[1]]
      
    }
    iprice <- which(landPriceTable$davis_class_any == davis_class & landPriceTable$ie == selected_elem)
    crop_el[icrop,3] <- landPriceTable$mean_ie_lu_price[iprice]
  }
  
  elem_area <- slot(cv_mesh_poly[[iel]], "area")*(1/4046.86)
  crop_acres <- elem_area*crop_el[,2]/sum(crop_el[,2])
  elem_price[iel] <- sum(crop_acres*crop_el[,3])
}