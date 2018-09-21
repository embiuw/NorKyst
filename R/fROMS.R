#' \code{readROMS} Read Norfjords data
#'
#' @param roms.dir Repository of Norkyst data (in NetCDF format)
#' @param date Text string matching date info in file name
#' @param var Name of variable to read (currently limited to one variable only)
#' @param xlim Specify longitude extent of subregion of grid. Note! This is only possible
#' if data are converted to xyz format
#' @param ylim Specify latitude extent of subregion of grid. Note! This is only possible
#' if data are converted to xyz format
#' @param to.xyz Logical, should data be converted to xyz points (TRUE)
#' or be kept in original array form (FALSE)
#'
#' @return Returns either 1) list containing latitude (matrix), longitude (matrix) and variable (matrix or array)
#' or 2) data frame of lat/lon/variable data
#' @details Reads data corresponding to one variable from Norkyst model output covering the Troms coast.
#' Expects input in NetCDF4 format, and uses the \code{\link{ncdf4}} package.
#' @family Norkyst model output
#' @seealso \code{\link{cropROMS}} to crop data in matrix & array form,
#'   \code{\link{imageROMS}} to plot data
#' @author Martin Biuw
#' @importFrom ncdf4 nc_open ncvar_get nc_close
#' @export


readROMS <- function(roms.dir='R:/Jofrid/Kaldfjord_160m_2014', date='20141101', var='u',
                     xlim=c(18.3, 18.8), ylim=c(69.6, 69.9), to.xyz=F, layer=35) {

  varname <- var

  filelist <- dir(roms.dir)
  filedates <- unlist(lapply(filelist, function(x) {
    substr(tail(unlist(strsplit(x, '_')), 1), 1, 8)
  }))

  which.file <- match(date, filedates)
  nc <- nc_open(paste(roms.dir, filelist[which.file], sep='/'))

  dim.nam <- paste(c('lon', 'lat'), var, sep='_')

  if(var=='temp' | var=='salt' | var=='h') {
    dim.nam <- paste(c('lon', 'lat'), 'rho', sep='_')
  }

  eval(parse(text=paste('lon <- ncvar_get(nc, nc$var$', dim.nam[1], ')', sep='')))
  eval(parse(text=paste('lat <- ncvar_get(nc, nc$var$', dim.nam[2], ')', sep='')))
  eval(parse(text=paste('var <- ncvar_get(nc, nc$var$', var, ')', sep='')))

  nc_close(nc)

  if(any(!is.na(xlim))) {
    lon.clip <- which(lon>=xlim[1] & lon<=xlim[2])
    lat.clip <- which(lat>=ylim[1] & lat<=ylim[2])
  }


  if(to.xyz) {
    if(length(dim(var))==3) {
      xyz <- data.frame(lon=lon[intersect(lon.clip, lat.clip)],
                      lat=lat[intersect(lon.clip, lat.clip)],
                      var=var[,,layer][intersect(lon.clip, lat.clip)])
    } else {
      xyz <- data.frame(lon=lon[intersect(lon.clip, lat.clip)],
                        lat=lat[intersect(lon.clip, lat.clip)],
                        var=var[,][intersect(lon.clip, lat.clip)])
    }
    names(xyz)[3] <- varname
    xyz
  } else {
    out.l <- list(lon=lon, lat=lat, var=var)
    names(out.l)[3] <- varname
    out.l
  }
}

#' \code{cropROMS} Spatial crop of Norfjords model output
#'
#' @param uc u current vector data (list with components lon, lat, and u)
#' @param vc v current vector data (list with components lon, lat, and v)
#' @param xdim Specify subregion (row indexes) of matrix in first dimension.
#' @param ydim Specify subregion (column indexes) of matrix in second dimension.
#'
#' @return Returns a cropped version of the dataset, as a list with sublists u and v, whcih are
#' themselves lists containing objects lat, lon and either u or v
#' (i.e. similar to unconverted output from readROMS)
#' @details Crops matrices and arrays of u and v variables, given specified row and column matrices.
#' Note! Function currently only functional for current vector (u & v) data combined.
#' @family Norkyst model output
#' @seealso \code{\link{readROMS}} to read Norfjords data,
#'   \code{\link{imageROMS}} to plot data
#' @author Martin Biuw
#' @export

cropROMS <- function(uc=u, vc=v, xdim=c(100:200), ydim=c(1:150)) {
  vdim=list(xdim=c(xdim, max(xdim)+1), ydim=ydim[-length(ydim)])

  uc$lon <- uc$lon[xdim, ydim]
  uc$lat <- uc$lat[xdim, ydim]
  uc$u <- uc$u[xdim, ydim,]

  vc$lon <- vc$lon[vdim$xdim, vdim$ydim]
  vc$lat <- vc$lat[vdim$xdim, vdim$ydim]
  vc$v <- vc$v[vdim$xdim, vdim$ydim,]
  list(u=uc, v=vc)
}


#' \code{imageROMS} Plot current vector data from Norfjords model
#'
#' @param u current vector data (list with components lon, lat, and u)
#' @param v current vector data (list with components lon, lat, and v)
#' @param arr.dens Plot arrows for every \code{arr.dens} grid cell
#' @param scale Scale the length of arrows to provide visually pleasing results
#' @param arr.col Colour of arrows
#'
#' @return Silent. Plots the desired current field
#' @details Creates an image map, where the colour intensity represents the current speed
#' and arrows represent direction. This version works only for data in original matrix/array format,
#' and currently only plots surface currents
#' @family Norkyst model output
#' @seealso \code{\link{readROMS}} to read Norfjords data,
#'   \code{\link{cropROMS}} to crop Norkyst data
#' @author Martin Biuw
#' @importFrom fields image.plot tim.colors
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export

imageROMS <- function(u=u, v=v, arr.dens=2, scale=100, arr.col='slategrey') {
  options(warn=-1)
  del.dims <- c(max(c(dim(u$u)[1], dim(v$v)[1])),
               max(c(dim(u$u)[2], dim(v$v)[2])))
  im.dims <- list(c(1:(del.dims[1]-1)),
                  c(1:(del.dims[2]-1)))

  arr.mids <- as.matrix(expand.grid(seq(1, del.dims[1]-1, by=arr.dens),
                   seq(1, del.dims[2]-1, by=arr.dens), KEEP.OUT.ATTRS = F))
  arr.starts <- t(apply(arr.mids, 1, function(x) {
    c(x[1]-(0.5*((scale*u$u[x[1],x[2],dim(u$u)[3]]))), x[2]-(0.5*(scale*v$v[x[1], x[2], dim(v$v)[3]])))
  }))
  arr.stops <- t(apply(arr.mids, 1, function(x) {
    c(x[1]+(0.5*((scale*u$u[x[1],x[2],dim(u$u)[3]]))), x[2]+(0.5*(scale*v$v[x[1], x[2], dim(v$v)[3]])))
  }))

  image.plot(im.dims[[1]],
             im.dims[[2]],
             sqrt(u$u[,-del.dims[2],1]^2+v$v[-del.dims[1],,1]^2),
             col=colorRampPalette(brewer.pal(3, 'Reds'))(256),
             axes=F, xlab='', ylab='')
              box()
  arrows(arr.starts[,1], arr.starts[,2], arr.stops[,1], arr.stops[,2], length=0.05, col=arr.col)
  options(warn=0)
}
