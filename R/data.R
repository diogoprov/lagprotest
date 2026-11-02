#' Seasonal hepatic and climatic series of *Dendropsophus minutus*
#'
#' Monthly means of four hepatic modules and of climate for the treefrog
#' *Dendropsophus minutus* at a site in southeastern Brazil (20.8992° S,
#' 49.3502° W), sampled independently each month from April 2016 to March
#' 2017. Each module is a set of variables measured on the same individuals,
#' averaged by month, so every element is a 12-row series covering one
#' annual cycle.
#'
#' @format A named list of five numeric matrices, each with 12 rows (months,
#'   named `"2016-04"` to `"2017-03"`):
#'   \describe{
#'     \item{`Somatic`}{Scaled mass index (`mean_smi`), hepatosomatic index
#'       (`mean_ihs`), body mass (`mean_weight`) and snout–vent length
#'       (`mean_size`).}
#'     \item{`Histochemistry`}{Mean area of hemosiderin (`hemosiderin`),
#'       glycogen (`mean_gli`), melanin (`mean_mel`) and lipofuscin
#'       (`mean_lipo`) per image.}
#'     \item{`Morphometry`}{Hepatocyte area and volume (`mean_cellA`,
#'       `mean_cellV`), nuclear volume and area (`mean_nuV`, `mean_nuA`) and
#'       nucleus–cytoplasm ratio (`mean_ncr`).}
#'     \item{`Density`}{Volumetric density of hepatocytes, sinusoids,
#'       melanomacrophage centers (`cmm`) and portal tissue, as centered
#'       log-ratios.}
#'     \item{`Climate`}{Mean daily minimum and maximum temperature (`Tmin`,
#'       `Tmax`, °C), Palmer drought severity index (`Drought`), mean daily
#'       rainfall (`Rainfall`, mm), relative humidity (`RH`, %), UV-B
#'       irradiance (`UVB`, J m\eqn{^{-2}} per day) and photoperiod with
#'       civil twilight (`photoperiod`, h).}
#'   }
#'
#' @source Hepatic data: Franco-Belussi, L., Moraes, L., de Oliveira, C.,
#'   Fernandes, C. E. and Provete, D. B. (2026). Data and code for: Seasonal
#'   hepatic plasticity follows a temporal response hierarchy in a
#'   Neotropical frog. Zenodo. \doi{10.5281/zenodo.21613832}. Licensed under
#'   CC BY 4.0.
#'
#'   Climate: temperature, rainfall and humidity from BR-DWGD (Xavier et al.
#'   2022, \doi{10.1002/joc.7731}); drought from TerraClimate (Abatzoglou et
#'   al. 2018, \doi{10.1038/sdata.2017.191}); UV-B from glUV (Beckmann et al.
#'   2014, \doi{10.1111/2041-210X.12168}); photoperiod computed for the site
#'   coordinates. Cite these sources when using the corresponding variables.
#'
#' @examples
#' str(dminutus_liver)
#' fit <- lagged_protest(dminutus_liver$Climate,
#'                       dminutus_liver$Histochemistry, lags = 0:3,
#'                       nperm = 199, seed = 1)
#' fit$optimal_lag
"dminutus_liver"
