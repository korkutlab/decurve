#' Set a distinct number for each experiment
#'
#' Each experiment has a cell line names followed by blank rows. Use this information to set a distinct number for each experiment.
#'
#' @param tb A tibble of dose effect data for many experiments with a cell_line column. In this column, each experiment has a cell line name in the first row followed by blank rows.
#' @return \item{A tibble}{It has a distinct number for each expt.}
#' @examples
#' library(decurve)
#' tb <- add_expt_num(tb)
#' @export
add_expt_num <- function(tb) {
  tb_all <- dplyr::tibble()
  expt <- 0
  for (i in 1:nrow(tb)) {
    if (!is.na(tb[i, 1])) {
      expt <- expt + 1
    }
    tb_all <- tb_all %>%
      dplyr::bind_rows(dplyr::mutate(tb[i, ], expt = expt))
  }

  tb_all <- tb_all %>%
    select(expt, everything())

  return(tb_all)
}


#' Tidy dose effect data of a single experiment
#'
#' tidy data to have attributes: cell_line, drug1, drug2, dose1, dose2, effect
#'
#' @param tb A tibble of dose effect data. Each element is of string type.
#'
#' Columns: cell_line, drug1, drug2, effect1_..., effect12_..., effect2_....
#'
#' Rows:
#'
#' first row: cell line name, drug1 name, drug2 name, backgrounds.
#'
#' other rows: doses and effects.
#'
#' @return \item{A tidied tibble}{Columns: cell_line, drug1, drug2, dose1, dose2, effect}
#' @examples
#' library(decurve)
#' tb <- tidy_dose_effect(tb)
#' @export
tidy_dose_effect <- function(tb) {
  # get cell line
  cell_line <- as.character(tb$cell_line[1])
  tb <- tb %>%
    dplyr::select(-cell_line)

  # get drug names and background
  drug1 <- as.character(tb$drug1[1])
  drug2 <- as.character(tb$drug2[1])

  background <- tb[1, grepl("effect", colnames(tb))] %>%
    as.double() %>%
    mean(na.rm = TRUE)
  tb <- tb[-1, ]

  # get doses and effects
  tb <- tb %>%
    dplyr::rename(dose1 = drug1,
                  dose2 = drug2) %>%
    dplyr::select(dose1, dose2, everything())

  # step1: effect0 (no drug); remove background; remove NA
  effect0 <- tb[1, -c(1:2)] %>%
    as.double()
  effect0 <- effect0 - background
  effect0 <- effect0[!is.na(effect0)]

  # normalize effect0
  norm_const <- mean(effect0)
  effect0 <- effect0/norm_const
  tb <- tb[-1, ]

  # step2: conver doses and effects to double; remove background; normalize effects
  mt <- as.matrix(tb)
  mt <- apply(mt, 2, as.numeric)
  mt[, -c(1:2)] <- mt[, -c(1:2)] - background
  mt[, -c(1:2)] <- mt[, -c(1:2)]/norm_const

  # wide to long (separate conditions and replicates to different rows)
  tb <- as_tibble(mt)
  tb <- tb %>%
    tidyr::gather("type", "effect", -dose1, -dose2) %>%
    dplyr::filter(!is.na(effect))

  # correct doses when single agents are used
  tb$dose2[grepl("effect1_", tb$type)] <- 0
  tb$dose1[grepl("effect2_", tb$type)] <- 0
  tb <- tb %>%
    dplyr::select(-type)

  # add effect0 (no drug)
  tb <- dplyr::tibble(dose1 = 0,
                      dose2 = 0,
                      effect = effect0) %>%
    dplyr::bind_rows(tb)

  # prepare output
  tb <- tb %>%
    dplyr::mutate(cell_line = cell_line,
                  drug1 = drug1,
                  drug2 = drug2) %>%
    dplyr::select(cell_line, drug1, drug2, everything())

  return(tb)
}


# GR
# formula 3 in the paper:
# GRcalculator: an online tool for calculating and mining dose–response data
# Clark, et al. BMC Cancer volume 17, Article number: 698 (2017)
calc_GR <- function(effect, T_assay, T_double) {
  if (effect < 0) {
    GR <- -1
  } else {
    GR <- 2^( (log2(effect) + T_assay/T_double) / (T_assay/T_double) ) - 1
  }

  return(GR)
}

#' Calculate growth rates for a single experiment
#'
#' Given assay duration and doubling time, convert effect to GR for each dose.
#'
#' GR(c) = 2^( (log2( x(c)/x(0) ) + T/Td) / (T/Td) ) - 1
#'       = 2^( (log2( effect(c) ) + T/Td) / (T/Td) ) - 1
#'
#' x's are the number of cells at the end of the assay
#' in an untreated (or vehicle-treated) control well (x(0)) and
#' in a drug-treated well (x(c)).
#'
#' T is assay duration.
#'
#' Td is doubling time of untreated cells.
#'
#' @param tb A tibble of dose effect data.
#'
#' Columns: cell_line (string), drug1 (string), drug2 (string), dose1 (double), dose2 (double), effect (double)
#'
#' @param assay_duration Assay duration.
#' @param doubling_time Doubling time of untreated cells.
#' @return \item{A tibble}{Columns: cell_line, drug1, drug2, dose1, dose2, GR}
#' @examples
#' library(decurve)
#' tb <- effect2gr(tb, doubling_time = 44, assay_duration = 120)
#' @export
effect2GR <- function(tb, assay_duration, doubling_time) {
  tb <- tb %>%
    dplyr::mutate(GR = purrr::map_dbl(effect,
                                      calc_GR,
                                      T_assay = assay_duration,
                                      T_double = doubling_time)) %>%
    dplyr::select(-effect)

  return(tb)
}
