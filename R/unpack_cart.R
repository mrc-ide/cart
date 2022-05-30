unpack_cart <- function(iso3c_sf, stack){
  sitesv <- as(iso3c_sf, "SpatVector")
  raw_values <- terra::extract(x = stack, y = sitesv) %>%
    group_by(ID) %>%
    summarise(across(everything(), list)) %>%
    select(-ID) %>%
    ungroup()
  sf_tibble <- tibble::as_tibble(sf::st_drop_geometry(iso3c_sf))
  out <- dplyr::bind_cols(sf_tibble, raw_values)
  return(out)
}
