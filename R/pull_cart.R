#' Pull cartographic information
#'
#' @param iso3c Country iso3c code
#' @param year Year of prevalence outputs
#' @param vector_rasters List of vector species rasters.
#'   See \code{\link[malariaAtlas]{listRaster}} for available rasters. Specified
#'    rasters must have an extent >= to the extent of the country specified.
#'   Please note the different definitions of vector species metrics -
#'   especially relative abundance vs probability of occurrence.
#' @param prevalence_rasters List of prevalence rasters. See
#'   \code{\link[malariaAtlas]{listRaster}} for available rasters. Default is
#'   for both falciparum and vivax.
#' @param spatial_limits_rasters List of spatial distribution limits rasters.
#'   Default is for both falciparum and vivax. See
#'   \code{\link[malariaAtlas]{listRaster}} for available rasters.
#'
#' @return A raster stack
#' @export
pull_cart <- function(iso3c, year,
                      vector_rasters = list(
                        funestus = "Anopheles funestus",
                        arabiensis = "Anopheles arabiensis Patton, 1905",
                        gambiae = "Anopheles gambiae Giles, 1902"),
                      prevalence_rasters = list(
                        pfpr = "Plasmodium falciparum PR2 - 10 version 2020",
                        pvpr = "Plasmodium vivax PR1-99 version 2020"),
                      spatial_limits_rasters = list(
                        pf_limits = "Plasmodium falciparum Spatial Limits",
                        pv_limits = "Plasmodium vivax Spatial Limits")){

  pop <- get_pop(iso3c = iso3c, year = year)
  prev <- get_prev(prevalence_rasters = prevalence_rasters, year = year, pop = pop)
  vectors <- get_vectors(vector_rasters = vector_rasters, year = year, pop = pop)
  spatial_limits <- get_spatial_limits(spatial_limits_rasters = spatial_limits_rasters, pop = pop)

  # Create stack
  rasters <- terra::rast(c(
    list(pop),
    prev,
    vectors,
    spatial_limits
  ))

  return(rasters)
}

#' Get population raster.
#'
#' Downloads the unconstrained individual countries 2000-2020 UN adjusted
#' (1km resolution) from the \href{https://www.worldpop.org/}{WorldPop server}.
#'
#' @inheritParams pull_cart
#'
#' @return Population raster
#' @export
get_pop <- function(iso3c, year) {
  raster_meta <- jsonlite::fromJSON(paste0("https://www.worldpop.org/rest/data/pop/wpicuadj1km?iso3=", iso3c))

  raster_files <- raster_meta$data %>%
    filter(popyear == year) %>%
    dplyr::select(files) %>%
    unlist()

  raster_file <- raster_files[grepl(".tif", raster_files)]

  td <- tempdir()
  raster_address <- paste0(td, "/", iso3c, ".tif")
  df <- utils::download.file(url = raster_file, destfile = raster_address, mode = "wb")
  pop <- terra::rast(raster_address)
  names(pop) <- "pop"

  return(pop)
}

#' Get prevalence rasters.
#'
#' Downloads the specified \href{https://malariaatlas.org/}{malaria atlas project} prevalence rasters.
#'
#' @inheritParams pull_cart
#'
#' @return A List of prevalence rasters.
#' @export
get_prev <- function(prevalence_rasters, year, pop) {
  pop_extent <- matrix(terra::ext(pop), ncol = 2, byrow = TRUE)

  prev <- list()
  for(i in seq_along(prevalence_rasters)){
    prev[[i]] <- malariaAtlas::getRaster(surface = prevalence_rasters[[i]], year = year, extent = pop_extent)
    names(prev[[i]]) <- names(prevalence_rasters)[i]
  }

  prev <- lapply(prev, terra::rast)

  prev <- lapply(prev, function(x, y){
    terra::resample(x, y)
  }, y = pop)

  return(prev)
}

#' Get vector species.
#'
#' Downloads the specified \href{https://malariaatlas.org/}{malaria atlas project}
#'   vector species rasters.
#'
#' @inheritParams pull_cart
#'
#' @return A List of vector species rasters.
#' @export
get_vectors <- function(vector_rasters, year, pop) {
  pop_extent <- matrix(terra::ext(pop), ncol = 2, byrow = TRUE)

  vectors <- list()
  for(i in seq_along(vector_rasters)){
    vectors[[i]] <- malariaAtlas::getRaster(surface = vector_rasters[[i]], extent = pop_extent)
    names(vectors[[i]]) <- names(vector_rasters)[i]
  }

  vectors <- lapply(vectors, terra::rast)

  vectors <- lapply(vectors, function(x, y){
    if(!identical(res(x), res(y))){
      x <- terra::resample(x, y)
      values(x)[values(x) < 0] <- 0
    }
    return(x)
  }, y = pop)

  return(vectors)
}

#' Get spatial limits rasters.
#'
#' Downloads the specified \href{https://malariaatlas.org/}{malaria atlas project}
#'   spatial limits of transmission rasters.
#'
#' @inheritParams pull_cart
#'
#' @return A List of spatial limits rasters.
#' @export
get_spatial_limits <- function(spatial_limits_rasters, pop) {
  pop_extent <- matrix(terra::ext(pop), ncol = 2, byrow = TRUE)

  spatial_limits <- list()
  for(i in seq_along(spatial_limits_rasters)){
    spatial_limits[[i]] <- malariaAtlas::getRaster(surface = spatial_limits_rasters[[i]], extent = pop_extent)
    names(spatial_limits[[i]]) <- names(spatial_limits_rasters)[i]
  }

  spatial_limits <- lapply(spatial_limits, terra::rast)

  spatial_limits <- lapply(spatial_limits, function(x, y){
    if(all.equal(res(x), res(y))){
      # Use method near for classes
      x <- terra::resample(x, y, method = "near")
    }
    return(x)
  }, y = pop)

  return(spatial_limits)
}
