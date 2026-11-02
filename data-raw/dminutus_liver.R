# Builds data/dminutus_liver.rda: monthly means of four hepatic modules and
# of climate for Dendropsophus minutus, April 2016 to March 2017.
#
# Source: monthly_fixed.rds, produced by rerunning the data chunks of
# Data_analysis_workflow_v13.qmd from Franco-Belussi et al. (2026,
# doi:10.5281/zenodo.21613832; CC BY 4.0), with photoperiod recomputed for
# the correct coordinates (the workflow passed longitude and latitude in
# the wrong order to LightLogR::photoperiod()). Climate: BR-DWGD (Tmin,
# Tmax, Rainfall, RH; Xavier et al. 2022), TerraClimate (Drought, PDSI;
# Abatzoglou et al. 2018), glUV (UVB; Beckmann et al. 2014), photoperiod
# with civil twilight.
months <- format(seq(as.Date("2016-04-01"), by = "month", length.out = 12),
                 "%Y-%m")
m <- readRDS(file.path("~/Library/CloudStorage/Dropbox/Anlises_MS D.minutus",
                       "Current Zoology/reanalysis/monthly_fixed.rds"))
dminutus_liver <- lapply(m, function(M) {
  M <- as.matrix(as.data.frame(M))
  storage.mode(M) <- "double"
  rownames(M) <- months
  M
})
stopifnot(identical(names(dminutus_liver),
                    c("Somatic", "Histochemistry", "Morphometry", "Density",
                      "Climate")),
          all(vapply(dminutus_liver, nrow, integer(1)) == 12))
usethis::use_data(dminutus_liver, overwrite = TRUE)
