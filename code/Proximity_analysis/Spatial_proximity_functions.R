filter <- dplyr::filter
select <- dplyr::select
pull <- dplyr::pull

## load functions 
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## calculate distance function
calc_d <- function(pt1, pt2) {
  d <- sqrt( (( pt2[2] - pt1[2] )^2 ) + (( pt2[1] - pt1[1] )^2 ) )
  return(d)
}

## calculate angle
getDir <- function(origin, pt2){
  pt2 <- pt2 - origin
  deg <- abs(atan(pt2[2]/pt2[1]) * (180/pi))
  if(pt2[1] > 0 & pt2[2] > 0){deg <- deg
  } else if(pt2[1] < 0 & pt2[2] > 0){ deg <- 180-deg
  } else if(pt2[1] < 0 & pt2[2] < 0){ deg <- 180 + deg
  } else if(pt2[1] > 0 & pt2[2] < 0){ deg <- 360 - deg}
  return(deg)
}

# assign angle
assign.angle <- function(degrees, range){
  require(data.table)
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  
  b=data.table(range=range)
  b[,merge:=range]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  return(Merge_a_b)
}

anchorPoint <- function(degrees, anchor){
  object <- (degrees + (360 - anchor)) 
  object <- ifelse(object > 360, object - 360, object)
  return(object)
}


get_neighbors <- function(cell_x_coords, cell_y_coords, r, cell_coords_df, cell_id){
  
  all_cell_coord <- data.frame(x = cell_coords_df$x_centroid, y = cell_coords_df$y_centroid)
  circ_coord <- circleFun(c(cell_x_coords, cell_y_coords), r, 1000)
  circ_coord <- circ_coord[!duplicated(circ_coord),]
  # get cells within radius
  tmp <- circ_coord[,c(1:2)]
  io <- splancs::inout(all_cell_coord, tmp)
  
  if(mean(io) == 0){message("no proximal cells found to ", cell_id); return(NULL)}
  extract_cells <- as.data.frame(cell_coords_df)[which(io == TRUE),] ## get cell within radius
  if(nrow(extract_cells) == 1){message("no proximal cells found to ", cell_id); return(NULL)}
  
  # calculate distance and degrees
  ## Calculate distance between cells
  tmp <- extract_cells[,c('x_centroid', 'y_centroid')]
  tmp <- rbind(tmp[which(rownames(tmp) %in% cell_id),], tmp)
  tmp <- tmp[!duplicated(tmp),]
  
  if(nrow(tmp) == 1){message("no proximal cells found to ", cell_id); return(NULL)}
  
  toCheck <- combn(rownames(tmp), 2, simplify = FALSE)
  names(toCheck) <-
    sapply(toCheck, paste, collapse = " - ")
  toCheck <- toCheck[grep(pattern = cell_id, toCheck)]
  
  if(length(toCheck) == 0){return(NULL)}
  
  ## calculate distance
  tmp.d <- sapply(toCheck, function(j){
    calc_d(tmp[j[1],c(1,2)], tmp[j[2],c(1,2)]) })
  names(tmp.d) <- gsub('\\.y_centroid', '', names(tmp.d))
  tmp.d <- as.data.frame(do.call(rbind, tmp.d))
  tmp.d$cell_id <- rownames(tmp.d)
  tmp.d <- tidyr::separate(tmp.d, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
  rownames(tmp.d) <- NULL
  
  ## calculate degrees
  tmp.degree <- sapply(toCheck, function(j){
    getDir(c(cell_x, cell_y), tmp[j[2],c(1,2)]) })
  names(tmp.degree) <- gsub('\\.y_centroid', '', names(tmp.degree))
  tmp.degree <- as.data.frame(do.call(rbind, tmp.degree))
  tmp.degree$cell_id <- rownames(tmp.degree)
  tmp.degree <- tidyr::separate(tmp.degree, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
  rownames(tmp.degree) <- NULL
  
  ## add degrees to distance table
  tmp.d$degree <- tmp.degree$V1[match(tmp.d$cell.b, tmp.degree$cell.b)]
  tmp.d$celltypeA <- extract_cells$final_CT[match(tmp.d$cell.a, rownames(extract_cells))]
  tmp.d$celltypeB <- extract_cells$final_CT[match(tmp.d$cell.b, rownames(extract_cells))]
  
  range <- seq(60, 360, 60) -30
  Merge_a_b <- assign.angle(degrees = tmp.d$degree, range = range)
  tmp.d$angle <- Merge_a_b$range[match(tmp.d$degree, Merge_a_b$degrees)]
  tmp.d <- tmp.d[!duplicated(tmp.d),]
  return(tmp.d)
}
