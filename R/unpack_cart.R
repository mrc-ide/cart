#' Extract information from rasters
#'
#' @param iso3c_sf A simple feature shape file to extract for
#' @param stack The raster stack from \code{\link{pull_cart}}
#'
#' @return Tibble with list columns of raw extracted values
#' @export
unpack_cart <- function(iso3c_sf, stack){
  sitesv <- methods::as(iso3c_sf, "SpatVector")
  raw_values <- terra::extract(x = stack, y = sitesv) %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), list)) %>%
    dplyr::select(-(.data$ID)) %>%
    dplyr::ungroup()
  sf_tibble <- tibble::as_tibble(sf::st_drop_geometry(iso3c_sf))
  out <- dplyr::bind_cols(sf_tibble, raw_values)
  return(out)
}
