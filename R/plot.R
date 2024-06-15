se <- function(x) sqrt(var(x)/length(x))

#' Dose effect curve
#'
#' Plot dose effect curve of a single experiment
#'
#' @import ggpubr
#'
#' @param tb A tibble of dose effect data.
#'
#' Columns: cell_line (string), drug1 (string), drug2 (string), dose1 (double), dose2 (double), effect (double)
#'
#' @return \item{A ggplot object}{returned ggplot object by ggline}
#' @examples
#' library(decurve)
#' p <- plot_dose_effect(tb)
#' print(p)
#' @export
plot_dose_effect <- function(tb, xlog10 = FALSE) {
  # cell line
  cell_line <- tb$cell_line[1]

  # drugs
  drug1 <- tb$drug1[1]
  drug2 <- tb$drug2[1]

  # sort doses and convert to string
  tb <- tb %>%
    dplyr::arrange(dose1, dose2)

  # drug12
  tb12 <- tb %>%
    dplyr::filter((dose1 == 0 & dose2 == 0) |
                    (dose1 > 0 & dose2 > 0) ) %>%
    dplyr::mutate(drugs = paste0(drug1, "+", drug2)) %>%
    dplyr::mutate(doses = paste0(dose1, "\n", dose2))

  # dose map to label (doses)
  tb_dose_map <- tb12 %>%
    dplyr::select(dose1, dose2, doses) %>%
    dplyr::distinct()

  tb12 <- tb12 %>%
    dplyr::select(drugs, doses, effect)

  # drug1
  tb1 <- tb %>%
    dplyr::filter(dose2 == 0) %>%
    dplyr::mutate(drugs = drug1) %>%
    dplyr::left_join(tb_dose_map, by = "dose1") %>%
    dplyr::select(drugs, doses, effect)

  # drug2
  tb2 <- tb %>%
    dplyr::filter(dose1 == 0) %>%
    dplyr::mutate(drugs = drug2) %>%
    dplyr::left_join(tb_dose_map, by = "dose2") %>%
    dplyr::select(drugs, doses, effect)

  # add space before single drugs so that they are ordered before the comb
  tb1$drugs <- paste0(" ", tb1$drugs)
  tb2$drugs <- paste0(" ", tb2$drugs)

  # merge singles and comb
  tb <- tb1 %>%
    dplyr::bind_rows(tb2) %>%
    dplyr::bind_rows(tb12)

  # simplify if two doses are same
  if (sum(!sapply(str_split(tb$doses, "\n"), function(x) x[1] == x[2])) == 0) {
    tb$doses <- sapply(str_split(tb$doses, "\n"), function(x) x[1])
  }

  # plot
  if (xlog10) {
    if (sum(grepl("\n", tb$doses)) > 0) {
      stop("To plot x-axis in log scale, dose1 and dose2 need to be equal.")
    }

    # calculate mean and se
    tb_mean <- tb %>%
      mutate(id = paste0(drugs, "|", doses)) %>%
      select(id, effect) %>%
      group_by(id) %>%
      summarise(effect_mean = mean(effect))

    tb_se <- tb %>%
      mutate(id = paste0(drugs, "|", doses)) %>%
      select(id, effect) %>%
      group_by(id) %>%
      summarise(effect_se = se(effect))

    tb <- tb %>%
      distinct(drugs, doses) %>%
      mutate(id = paste0(drugs, "|", doses)) %>%
      left_join(tb_mean, by = "id") %>%
      left_join(tb_se, by = "id") %>%
      select(-id)

    # recast doses to numeric
    tb <- tb %>%
      mutate(doses = as.double(doses))

    # change doses 0 to a small number (dummy)
    d <- sort(unique(tb$doses))
    d <- d[d > 0]
    dummy_dose <- max(d)/10
    while (dummy_dose >= min(d)) {
      dummy_dose <- dummy_dose/10
    }
    tb <- tb %>%
      mutate(doses = ifelse(doses == 0, dummy_dose, doses))

    p <- ggplot(tb, aes(x=doses, y=effect_mean,
                        ymin=effect_mean-effect_se,
                        ymax=effect_mean+effect_se,
                        group=drugs, color=drugs)) +
      geom_line() +
      geom_point()+
      geom_errorbar(width = 0.05)

    p <- p +
      labs(title=cell_line, x="doses", y = "effect")+
      theme_classic() +
      scale_color_manual(values=c("red", "blue", "green")) +
      #ylim(min(min(tb$effect_mean), 0), max(max(tb$effect_mean), 1)) +
      scale_y_continuous(limits = c(min(min(tb$effect_mean - tb$effect_se), 0),
                                    max(max(tb$effect_mean + tb$effect_se), 1)),
                         breaks = seq(0, 1, 0.25),
                         labels = seq(0, 1, 0.25) %>%
                           as.character() %>%
                           stringr::str_replace("\\.0+$", "") %>%
                           stringr::str_replace("(\\.[1-9]+)0+$", "\\1")
                         ) +
      scale_x_continuous(trans = "log10", # scales::pseudo_log_trans(sigma = 0.01, base = 10), # not working well
                         breaks = c(min(tb$doses), max(tb$doses)/10, max(tb$doses)),
                         labels = c(0, max(tb$doses)/10, max(tb$doses)) %>%
                           as.character() %>%
                           stringr::str_replace("\\.0+$", "") %>%
                           stringr::str_replace("(\\.[1-9]+)0+$", "\\1")
                         )

  } else {
    p <- ggpubr::ggline(tb, "doses", "effect",
                        color = "drugs",
                        add = "mean_se", # import ggpubr in this function, otherwise "mean_se" is not known by the system
                        palette = c("red", "blue", "green"),
                        title = cell_line) +
      # ylim(min(min(tb$effect), 0), max(max(tb$effect), 1)) +
      scale_y_continuous(limits = c(min(min(tb$effect), 0), max(max(tb$effect), 1)),
                         breaks = seq(0, 1, 0.25),
                         labels = seq(0, 1, 0.25) %>%
                           as.character() %>%
                           stringr::str_replace("\\.0+$", "") %>%
                           stringr::str_replace("(\\.[1-9]+)0+$", "\\1")
                         )
  }

  return(p)
}


#' Dose GR curve
#'
#' Plot dose GR curve of a single experiment
#'
#' @import ggpubr
#'
#' @param tb A tibble of dose GR data.
#'
#' Columns: cell_line (string), drug1 (string), drug2 (string), dose1 (double), dose2 (double), GR (double)
#'
#' @return \item{A ggplot object}{returned ggplot object by ggline}
#' @examples
#' library(decurve)
#' p <- plot_dose_GR(tb)
#' print(p)
#' @export
plot_dose_GR <- function(tb) {
  # cell line
  cell_line <- tb$cell_line[1]

  # drugs
  drug1 <- tb$drug1[1]
  drug2 <- tb$drug2[1]

  # sort doses and convert to string
  tb <- tb %>%
    dplyr::arrange(dose1, dose2)

  # drug12
  tb12 <- tb %>%
    dplyr::filter((dose1 == 0 & dose2 == 0) |
                    (dose1 > 0 & dose2 > 0) ) %>%
    dplyr::mutate(drugs = paste0(drug1, "+", drug2)) %>%
    dplyr::mutate(doses = paste0(dose1, "\n", dose2))

  # dose map to label (doses)
  tb_dose_map <- tb12 %>%
    dplyr::select(dose1, dose2, doses) %>%
    dplyr::distinct()

  tb12 <- tb12 %>%
    dplyr::select(drugs, doses, GR)

  # drug1
  tb1 <- tb %>%
    dplyr::filter(dose2 == 0) %>%
    dplyr::mutate(drugs = drug1) %>%
    dplyr::left_join(tb_dose_map, by = "dose1") %>%
    dplyr::select(drugs, doses, GR)

  # drug2
  tb2 <- tb %>%
    dplyr::filter(dose1 == 0) %>%
    dplyr::mutate(drugs = drug2) %>%
    dplyr::left_join(tb_dose_map, by = "dose2") %>%
    dplyr::select(drugs, doses, GR)

  # add space before single drugs so that they are ordered before the comb
  tb1$drugs <- paste0(" ", tb1$drugs)
  tb2$drugs <- paste0(" ", tb2$drugs)

  # merge singles and comb
  tb <- tb1 %>%
    dplyr::bind_rows(tb2) %>%
    dplyr::bind_rows(tb12)

  # plot
  p <- ggpubr::ggline(tb, "doses", "GR",
                      color = "drugs",
                      add = "mean_se", # import ggpubr in this function, otherwise "mean_se" is not known by the system
                      palette = c("red", "blue", "green"),
                      title = cell_line) +
    ylim(min(min(tb$GR), -1), max(max(tb$GR), 1))

  return(p)
}
