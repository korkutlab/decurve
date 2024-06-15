#' Calculate Bliss scores of combinations for a single experiment
#'
#' Bliss = mean effect of comb drugs  / (mean effect of drug1 * mean effect of drug2)
#'
#' @param tb A tibble of dose effect data.
#'
#' Columns: cell_line (string), drug1 (string), drug2 (string), dose1 (double), dose2 (double), effect (double)
#'
#' @return \item{A tibble}{Columns: cell_line, drug1, drug2, dose1, dose2, effect1 (mean), effect2 (mean), effect12 (mean), bliss}
#' @examples
#' library(decurve)
#' tb <- calc_bliss(tb)
#' @export
calc_bliss <- function(tb) {
  tb_comb_dose <- tb %>%
    dplyr::filter(dose1 != 0 & dose2 != 0) %>%
    dplyr::distinct(dose1, dose2, .keep_all = TRUE) %>%
    dplyr::select(-effect)

  tb_all <- tibble()
  for (i in 1:nrow(tb_comb_dose)) {
    d1 <- tb_comb_dose$dose1[i]
    d2 <- tb_comb_dose$dose2[i]
    effect1 <- tb %>%
      dplyr::filter(dose1 == d1 & dose2 == 0) %>%
      .$effect %>%
      mean()
    effect2 <- tb %>%
      dplyr::filter(dose1 == 0 & dose2 == d2) %>%
      .$effect %>%
      mean()
    effect12 <- tb %>%
      dplyr::filter(dose1 == d1 & dose2 == d2) %>%
      .$effect %>%
      mean()
    bliss <- effect12 / (effect1 * effect2)
    tb_all <- tb_all %>%
      dplyr::bind_rows(dplyr::tibble(effect1 = effect1,
                                     effect2 = effect2,
                                     effect12 = effect12,
                                     bliss = bliss))
  }

  tb_bliss <- tb_comb_dose %>%
    dplyr::bind_cols(tb_all)

  return(tb_bliss)
}

