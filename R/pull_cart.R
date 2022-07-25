#' Pull cartographic information
#'
#' Creates a terra raster stack for specified country and datasets. Only years
#' for which there are population data will be returned. Note, this
#' function is fragile and may fail if specified rasters are not available for
#' a given country.
#'
#' @param iso3c Country iso3c code
#' @param years Years
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
#' @param intervention_rasters List of intervention coverage rasters.
#'   Default is for both ITN, IRS and treatment. See
#'   \code{\link[malariaAtlas]{listRaster}} for available rasters.
#' @param urban_density_threshold Threshold for urban classification (people per square km).
#' Default is 1500 people per square km. Based on urban classification by Eurostat,
#' detailed in \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8065116/}{Khavari et al}.
#' @param quiet If TRUE, suppress status messages (if any), and the progress bar.
#'
#'
#' @return A raster stack
#' @export
pull_cart <- function(iso3c,
                      years,
                      vector_rasters = list(
                        funestus = "Anopheles funestus",
                        arabiensis = "Anopheles arabiensis Patton, 1905",
                        gambiae = "Anopheles gambiae Giles, 1902"
                      ),
                      prevalence_rasters = list(
                        pfpr = "Plasmodium falciparum PR2 - 10 version 2020",
                        pvpr = "Plasmodium vivax PR1-99 version 2020"
                      ),
                      spatial_limits_rasters = list(
                        pf_limits = "Plasmodium falciparum Spatial Limits",
                        pv_limits = "Plasmodium vivax Spatial Limits"
                      ),
                      intervention_rasters = list(
                        itn_use = "Insecticide treated bednet (ITN) use version 2020",
                        irs_cov = "Indoor Residual Spraying (IRS) coverage version 2020",
                        tx = "Effective treatment with an Antimalarial drug version 2020"
                      ),
                      urban_density_threshold = 1500,
                      quiet = FALSE){

  # Temporal
  years <- years[pop_available(iso3c = iso3c, years = years)]

  raster_list <- malariaAtlas::listRaster(printed = FALSE) |>
    dplyr::select(.data$title, .data$min_raster_year, .data$max_raster_year)

  prevalence_rasters_available <- lapply(prevalence_rasters, raster_available,
                                         years = years, raster_list = raster_list)
  intervention_rasters_available <- lapply(intervention_rasters, raster_available,
                                           years = years, raster_list = raster_list)

  temporal <- list()
  for(i in seq_along(years)){
    # Population
    pop <- get_pop(iso3c = iso3c, year = years[i], quiet = quiet)
    ur <- make_urban_rural(pop = pop, urban_density_threshold = urban_density_threshold)
    # Spatial limits
    if(i == 1){
      spatial_limits <- list()
      for(j in seq_along(spatial_limits_rasters)){
        spatial_limit <- get_map(raster = spatial_limits_rasters[[j]], pop = pop, method = "near", quiet = quiet)
        names(spatial_limit) <- names(spatial_limits_rasters)[j]
        spatial_limits[[j]] <- spatial_limit
      }
      spatial_limits <- terra::rast(spatial_limits)
    }
    temporal[[i]] <- c(pop, ur, spatial_limits)
    # Prevalence
    for(j in seq_along(prevalence_rasters)){
      if(prevalence_rasters_available[[j]][i]){
        prev <- get_map(raster = prevalence_rasters[[j]], pop = pop, year = years[i], quiet = quiet)
        names(prev) <- names(prevalence_rasters)[j]
        temporal[[i]] <- c(temporal[[i]], prev)
      }
    }
    # Interventions
    for(j in seq_along(intervention_rasters)){
      if(intervention_rasters_available[[j]][i]){
        int <- get_map(raster = intervention_rasters[[j]], pop = pop, year = years[i], quiet = quiet)
        names(int) <- names(intervention_rasters)[j]
        temporal[[i]] <- c(temporal[[i]], int)
      }
    }
  }
  names(temporal) <- years

  # Static
  vectors <- list()
  if(!is.null(vector_rasters)){
    for(i in seq_along(vector_rasters)){
      vector <- get_map(raster = vector_rasters[[i]], pop = pop, quiet = quiet)
      if(!is.null(vector)){
        names(vector) <- names(vector_rasters)[i]
        vectors <- c(vectors, vector)
      }
    }
    vectors <- terra::rast(vectors)
  }

  # Combined
  rasters <- list(temporal = temporal, vectors = vectors)

  return(rasters)
}

#' Get population raster.
#'
#' Downloads the unconstrained individual countries 2000-2020 UN adjusted
#' (1km resolution) from the \href{https://www.worldpop.org/}{WorldPop server}.
#'
#' @inheritParams pull_cart
#' @param year Year
#' @param quiet If TRUE, suppress status messages (if any), and the progress bar.
#'
#' @return Population raster
#' @export
get_pop <- function(iso3c, year, quiet = FALSE) {
  raster_meta <- jsonlite::fromJSON(paste0("https://www.worldpop.org/rest/data/pop/wpicuadj1km?iso3=", iso3c))

  raster_files <- raster_meta$data %>%
    dplyr::filter(.data$popyear == year) %>%
    dplyr::select(.data$files) %>%
    unlist()

  raster_file <- raster_files[grepl(".tif", raster_files)]

  td <- tempdir()
  raster_address <- paste0(td, "/", iso3c, "_", year, ".tif")
  df <- utils::download.file(url = raster_file, destfile = raster_address, mode = "wb", quiet = quiet)
  pop <- terra::rast(raster_address)
  names(pop) <- "pop"

  return(pop)
}

#' Get urban rural raster
#'
#' Recodes the population raster into urban and rural cells, based on a threshold
#' population density.
#'
#' @param pop Population raster (at 1km squared)
#' @inheritParams pull_cart
#'
#' @return Categorical raster
make_urban_rural <- function(pop, urban_density_threshold){
  ur <- pop
  terra::values(ur)[terra::values(pop) < urban_density_threshold] <- 0
  terra::values(ur)[terra::values(pop) >= urban_density_threshold] <- 1
  levels(ur) <- c("rural", "urban")
  names(ur) <- "urban_rural"
  return(ur)
}

#' Get Malaria Atlas Project (MAP) rasters.
#'
#' Downloads the specified \href{https://malariaatlas.org/}{malaria atlas project} rasters.
#'
#' @param raster Raster
#' @param pop Population raster
#' @param year Year
#' @param method See \code{\link[terra]{rast}} for more information
#' @param quiet If TRUE, suppress status messages (if any), and the progress bar
#' @param resample resample to match resolution and extent of pop
#'
#' @return A raster.
#' @export
get_map <- function(raster, pop, year = NA, method = "bilinear", quiet = FALSE, resample = TRUE) {
  pop_extent <- get_extent(pop)

  getRaster_safe <- purrr::possibly(malariaAtlas::getRaster, otherwise = NULL)
  if(quiet){
    out <- suppressMessages(getRaster_safe(surface = raster, year = year, extent = pop_extent))
  } else {
    out <- getRaster_safe(surface = raster, year = year, extent = pop_extent)
  }

  if(!is.null(out)){
    out <- terra::rast(out)

    if(resample){
      if(!identical(terra::res(out), terra::res(pop))){
        out <- terra::resample(out, pop, method = method)
      }
    }
  }

  return(out)
}

#' Extract and format raster extent
#'
#' @param x raster
#'
#' @return Extent as matrix
get_extent <- function(x){
  matrix(terra::ext(x), ncol = 2, byrow = TRUE)
}

#' Check if population data are available for a given country and given years
#'
#' @inheritParams pull_cart
#'
#' @return Boolean vector
pop_available <- function(iso3c, years){
  raster_meta <- jsonlite::fromJSON(paste0("https://www.worldpop.org/rest/data/pop/wpicuadj1km?iso3=", iso3c))
  years %in% raster_meta$data$popyear
}

#' Check if MAP raster data are available for a given country and given years
#'
#' @param raster_name Name of raster
#' @param raster_list Raster metadata from Malaria Atlas Project
#' @inheritParams pull_cart
#'
#' @return Boolean vector
raster_available <- function(raster_name, years, raster_list){
  raster_list[raster_list$title == raster_name, "min_raster_year"] <= years &
    raster_list[raster_list$title == raster_name, "max_raster_year"] >= years
}



