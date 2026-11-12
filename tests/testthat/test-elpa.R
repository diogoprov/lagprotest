# Small synthetic archive with the structure of the PELD-ELPA DwC-A: two
# habitats, two sites each, three months, three species; month zero-padded
# in the event and measurement tables but not in the occurrence table.
write_fake_elpa <- function(dir, months = 1:3) {
  base <- "BR:PELD-ELPA:FURG:patos:ichthyofauna:2000"
  sites <- data.frame(site = c("A", "B", "C", "D"),
                      habitat = c("estuary", "estuary", "sandy beach",
                                  "sandy beach"))
  ev <- merge(data.frame(month = months), sites)
  ev$eventID <- sprintf("%s:%02d:%s", base, ev$month, ev$site)
  ev$okey <- sprintf("%s:%d:%s", base, ev$month, ev$site)
  write_tsv <- function(x, f) {
    utils::write.table(x, file.path(dir, f), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "")
  }
  write_tsv(data.frame(id = ev$eventID, eventID = ev$eventID,
                       habitat = ev$habitat), "event.txt")
  spp <- c("Sp one", "Sp two", "Sp three")
  occ <- merge(ev, data.frame(j = 1:3, sp = spp))
  occ$occurrenceID <- sprintf("%s:%05d", occ$okey, occ$j)
  # Counts: 10 * month for sp one, 5 for sp two, 1 for sp three at site A
  # only; nothing of sp three elsewhere.
  occ$n <- ifelse(occ$j == 1, 10 * occ$month, ifelse(occ$j == 2, 5, 1))
  occ <- occ[!(occ$j == 3 & occ$site != "A"), ]
  write_tsv(data.frame(id = occ$okey, occurrenceID = occ$occurrenceID,
                       scientificName = occ$sp), "occurrence.txt")
  cnt <- data.frame(id = occ$eventID,
                    occurrenceID = sub(":([0-9]):", ":0\\1:",
                                       occ$occurrenceID),
                    measurementType = "individualCount",
                    measurementValue = occ$n)
  env <- function(type, val) {
    data.frame(id = ev$eventID, occurrenceID = NA, measurementType = type,
               measurementValue = val)
  }
  temp <- 20 + ev$month + (ev$site %in% c("B", "D"))
  sal <- as.character(ev$month + (ev$site %in% c("B", "D")))
  sal[ev$site == "B" & ev$month == 2] <- "ND"
  secchi <- ifelse(ev$habitat == "estuary" & ev$month == 2, NA, 0.5)
  write_tsv(rbind(cnt, env("Water temperature", temp),
                  env("Salinity", sal), env("Secchi disk", secchi)),
            "extendedmeasurementorfact.txt")
  invisible(dir)
}

fake_dir <- function(...) {
  dir <- tempfile("elpa")
  dir.create(dir)
  write_fake_elpa(dir, ...)
}

test_that("elpa_blocks() builds monthly blocks from the archive", {
  dir <- fake_dir()
  el <- elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 3),
                    species = "all", hellinger = FALSE, path = dir)
  expect_named(el, c("Environment", "Fish"))
  expect_identical(rownames(el$Fish), c("2000-01", "2000-02", "2000-03"))
  expect_identical(colnames(el$Environment), c("temp", "sal", "secchi"))
  # Catch per haul (count / 5), averaged over sites A and B.
  expect_equal(unname(el$Fish[, "Sp one"]), 2 * (1:3))
  expect_equal(unname(el$Fish[, "Sp three"]), rep(0.1, 3))
  expect_equal(unname(el$Environment[, "temp"]), 20.5 + 1:3)
  # "ND" salinity is missing, not zero: month 2 uses site A only.
  expect_equal(unname(el$Environment[, "sal"]), c(1.5, 2, 3.5))
  # Secchi missing at every estuary site in month 2: interpolated.
  expect_equal(unname(el$Environment[, "secchi"]), rep(0.5, 3))
})

test_that("elpa_blocks() selects species and applies Hellinger", {
  dir <- fake_dir()
  el <- elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 3),
                    path = dir)
  expect_identical(colnames(el$Fish), c("Sp one", "Sp two"))
  raw <- elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 3),
                     species = "all", hellinger = FALSE, path = dir)
  expect_equal(el$Fish, sqrt(raw$Fish / rowSums(raw$Fish))[, 1:2])
  beach <- elpa_blocks("sandy beach", from = c(2000, 1), to = c(2000, 3),
                       species = "all", path = dir)
  expect_false("Sp three" %in% colnames(beach$Fish))
  expect_error(elpa_blocks("sandy beach", from = c(2000, 1),
                           to = c(2000, 3), species = "Sp three",
                           path = dir), "not found")
})

test_that("elpa_blocks() reads a zipped archive", {
  skip_if(Sys.which("zip") == "")
  dir <- fake_dir()
  zipfile <- tempfile(fileext = ".zip")
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  utils::zip(zipfile, list.files(), flags = "-q")
  setwd(old)
  expect_equal(elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 3),
                           path = zipfile),
               elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 3),
                           path = dir))
})

test_that("elpa_blocks() validates the period and the path", {
  dir <- fake_dir(months = c(1, 3))
  expect_error(elpa_blocks(from = c(2000, 1), to = c(2000, 3), path = dir),
               "No sampling in 1 month")
  expect_error(elpa_blocks(from = c(2000, 3), to = c(2000, 1), path = dir),
               "later than")
  expect_error(elpa_blocks(from = c(2000, 13), path = dir), "c\\(year")
  expect_error(elpa_blocks(path = file.path(dir, "nope")), "does not exist")
})

test_that("elpa_download() reuses a cached archive", {
  dir <- tempfile("cache")
  dir.create(dir)
  dest <- file.path(dir, "dwca-composition_abundance_fish_elpa-v1.16.zip")
  file.create(dest)
  expect_identical(elpa_download(dir), dest)
  expect_match(elpa_url(), "^https://ipt\\.sibbr\\.gov\\.br/")
})

test_that("elpa_download() fetches the archive when it is not cached", {
  dir <- tempfile("cache")
  dir.create(dir)
  local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      expect_identical(url, elpa_url())
      file.create(destfile)
      0L
    },
    .package = "utils"
  )
  expect_message(dest <- elpa_download(dir), "CC BY-NC")
  expect_true(file.exists(dest))
})

test_that("elpa_blocks() handles a one-month period", {
  dir <- fake_dir()
  el <- elpa_blocks("estuary", from = c(2000, 1), to = c(2000, 1),
                    species = "all", hellinger = FALSE, path = dir)
  expect_equal(dim(el$Fish), c(1L, 3L))
  expect_equal(dim(el$Environment), c(1L, 3L))
  expect_identical(rownames(el$Fish), "2000-01")
  # Secchi was not recorded in the estuary in month 2: nothing to fill from.
  expect_error(elpa_blocks("estuary", from = c(2000, 2), to = c(2000, 2),
                           path = dir), "Too few observed values")
})

test_that("elpa_blocks() downloads the archive when no path is given", {
  dir <- fake_dir()
  local_mocked_bindings(elpa_download = function(...) dir)
  expect_equal(elpa_blocks(from = c(2000, 1), to = c(2000, 3)),
               elpa_blocks(from = c(2000, 1), to = c(2000, 3), path = dir))
})
