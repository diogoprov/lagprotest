#' Patos Lagoon juvenile fish and environment blocks (PELD-ELPA)
#'
#' Builds two monthly blocks, environment and fish, for one habitat of the
#' Patos Lagoon estuary, southern Brazil, from the PELD-ELPA beach-seine
#' survey (Darwin Core Archive v1.16). The data are licensed CC BY-NC 4.0
#' and are not shipped with \pkg{lagprotest}: pass a local copy of the
#' archive in `path`, or leave `path = NULL` to download it from the SiBBr
#' IPT into the session's temporary directory.
#'
#' Each sampling event is one site in one month, with five seine hauls.
#' Counts are summed by species within an event and divided by five (catch
#' per haul), and both catch and environment are averaged over the sites of
#' the chosen habitat in each month. Months with no record of an
#' environmental variable are filled by linear interpolation; a month with
#' no fish sampling at all is an error, since the blocks must be regular
#' series. No transformation, detrending or deseasonalization is applied to
#' the environment; see the example for the steps used in the lagprotest
#' preprint, which also added freshwater discharge from the Brazilian
#' National Water Agency (not included here).
#'
#' @param habitat `"estuary"` or `"sandy beach"`.
#' @param from,to Integer vectors `c(year, month)` delimiting the period.
#'   The default covers 1996-08 to 2016-12, the longest stretch without
#'   unsampled months.
#' @param species Species (scientific names) to keep in the fish block.
#'   `NULL` (default) keeps the most abundant species that together make up
#'   90% of the total catch over all sampling events in the chosen habitat
#'   and period; `"all"` keeps
#'   every species.
#' @param hellinger Logical. If `TRUE` (default), the fish block holds
#'   Hellinger-transformed catches, computed on all species before
#'   `species` are selected; if `FALSE`, mean catch per haul.
#' @param path Path to the archive (`.zip`) or to a directory with its
#'   extracted files. `NULL` downloads version 1.16 from the SiBBr IPT
#'   (about 3 MB), once per session.
#'
#' @return A named list of two numeric matrices with one row per month (row
#'   names `"YYYY-MM"`):
#'   \describe{
#'     \item{`Environment`}{Water temperature (`temp`, °C), salinity (`sal`)
#'       and Secchi depth (`secchi`).}
#'     \item{`Fish`}{One column per species in `species`.}
#'   }
#'
#' @source Vieira, J. P., Garcia, A. M. and Lemos, V. M. PELD-ELPA:
#'   Species composition and abundance patterns of fish assemblages at
#'   shallow waters of Patos Lagoon estuary. Version 1.16. Sistema de
#'   Informação sobre a Biodiversidade Brasileira (SiBBr), via GBIF,
#'   \doi{10.15468/kci8zb}. Licensed CC BY-NC 4.0: cite the dataset and do
#'   not use it commercially. The program is described in Lemos et al.
#'   (2022), Earth System Science Data 14: 1015–1041,
#'   \doi{10.5194/essd-14-1015-2022}.
#'
#' @examplesIf interactive()
#' el <- elpa_blocks("estuary")
#' detrend <- function(m) apply(m, 2, function(x) residuals(lm(x ~ seq_along(x))))
#' X <- detrend(log1p(el$Environment))
#' Y <- detrend(el$Fish)
#' fit <- lagged_protest(X, Y, lags = 0:6, period = 12, nperm = 999,
#'                       seed = 1)
#' fit
#' @export
elpa_blocks <- function(habitat = c("estuary", "sandy beach"),
                        from = c(1996, 8), to = c(2016, 12),
                        species = NULL, hellinger = TRUE, path = NULL) {
  habitat <- match.arg(habitat)
  check_ym <- function(x) {
    is.numeric(x) && length(x) == 2 && !anyNA(x) && x[2] %in% 1:12
  }
  if (!check_ym(from) || !check_ym(to)) {
    stop("`from` and `to` must be c(year, month).", call. = FALSE)
  }
  ym <- function(y, m) y * 12 + m - 1
  if (ym(from[1], from[2]) > ym(to[1], to[2])) {
    stop("`from` must not be later than `to`.", call. = FALSE)
  }
  stopifnot(is.null(species) || is.character(species),
            is.logical(hellinger), length(hellinger) == 1, !is.na(hellinger))

  if (is.null(path)) path <- elpa_download()
  read <- elpa_reader(path)
  ev <- read("event.txt")
  mf <- read("extendedmeasurementorfact.txt")
  oc <- read("occurrence.txt")

  # Event identifiers pad the month with a zero in some tables only.
  norm <- function(x) {
    sub("(ichthyofauna:[0-9]{4}):0?([0-9]{1,2})(:|$)", "\\1:\\2\\3", x)
  }
  key <- norm(ev$eventID)
  part <- function(i) vapply(strsplit(key, ":", fixed = TRUE),
                             function(z) z[i], character(1))
  ev <- data.frame(key = key, t = ym(as.integer(part(6)),
                                     as.integer(part(7))),
                   habitat = ev$habitat)
  t0 <- ym(from[1], from[2])
  t1 <- ym(to[1], to[2])
  ev <- ev[ev$habitat == habitat & ev$t >= t0 & ev$t <= t1, ]
  months <- t0:t1
  missing <- setdiff(months, ev$t)
  if (length(missing)) {
    stop("No sampling in ", length(missing), " month(s) of the period (",
         paste(ym_label(missing[seq_len(min(3, length(missing)))]),
               collapse = ", "), if (length(missing) > 3) ", ...",
         "); choose `from` and `to` within a continuous stretch.",
         call. = FALSE)
  }

  # Environment: one value per event and variable, averaged by month.
  env_vars <- c(temp = "Water temperature", sal = "Salinity",
                secchi = "Secchi disk")
  ef <- mf[mf$measurementType %in% env_vars, ]
  ef_key <- norm(ef$id)
  env <- vapply(env_vars, function(v) {
    s <- ef$measurementType == v
    val <- suppressWarnings(as.numeric(ef$measurementValue[s]))
    by_event <- val[match(ev$key, ef_key[s])]
    m <- tapply(by_event, factor(ev$t, levels = months), mean, na.rm = TRUE)
    m[is.nan(m)] <- NA
    m <- as.numeric(m)
    if (anyNA(m)) fill_linear(m) else m
  }, numeric(length(months)))
  if (!is.matrix(env)) env <- matrix(env, nrow = 1)
  colnames(env) <- names(env_vars)

  # Fish: catch per haul by event and species, averaged by month.
  cf <- mf[mf$measurementType == "individualCount", ]
  okey <- norm(cf$occurrenceID)
  sp <- oc$scientificName[match(okey, norm(oc$occurrenceID))]
  ekey <- sub(":[0-9]+$", "", okey)
  keep <- ekey %in% ev$key & !is.na(sp)
  n <- suppressWarnings(as.numeric(cf$measurementValue[keep]))
  catch <- tapply(n, list(factor(ekey[keep], levels = ev$key),
                          factor(sp[keep])), sum)
  catch[is.na(catch)] <- 0
  catch <- catch / 5
  fish <- apply(catch, 2, function(x) {
    tapply(x, factor(ev$t, levels = months), mean)
  })
  if (!is.matrix(fish)) fish <- matrix(fish, nrow = 1,
                                       dimnames = list(NULL, colnames(catch)))

  total <- sort(colSums(catch), decreasing = TRUE)
  if (is.null(species)) {
    species <- names(total)[seq_len(which(cumsum(total) / sum(total) >=
                                            0.9)[1])]
  } else if (identical(species, "all")) {
    species <- names(total)
  } else {
    absent <- setdiff(species, colnames(fish))
    if (length(absent)) {
      stop("Species not found in the selected habitat and period: ",
           paste(absent, collapse = ", "), call. = FALSE)
    }
  }
  if (hellinger) {
    rs <- rowSums(fish)
    fish <- sqrt(fish / ifelse(rs > 0, rs, 1))
  }
  fish <- fish[, species, drop = FALSE]

  rn <- ym_label(months)
  dimnames(env) <- list(rn, colnames(env))
  dimnames(fish) <- list(rn, species)
  list(Environment = env, Fish = fish)
}

ym_label <- function(t) sprintf("%d-%02d", t %/% 12, t %% 12 + 1)

elpa_url <- function() {
  "https://ipt.sibbr.gov.br/peld/archive.do?r=composition_abundance_fish_elpa&v=1.16"
}

elpa_download <- function(dir = tempdir()) {
  dest <- file.path(dir, "dwca-composition_abundance_fish_elpa-v1.16.zip")
  if (!file.exists(dest)) {
    message("Downloading the PELD-ELPA archive (CC BY-NC 4.0) from the ",
            "SiBBr IPT.")
    old <- options(timeout = max(300, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    utils::download.file(elpa_url(), dest, mode = "wb", quiet = TRUE)
  }
  dest
}

elpa_reader <- function(path) {
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  is_dir <- dir.exists(path)
  function(f) {
    con <- if (is_dir) file.path(path, f) else unz(path, f)
    utils::read.delim(con, quote = "", colClasses = "character",
                      na.strings = "", comment.char = "",
                      encoding = "UTF-8")
  }
}
