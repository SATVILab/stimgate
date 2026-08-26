plot_cyt_comp_points <- function(data_tidy_cyt_comp_to_manual) {
  plot_tbl <- data_tidy_cyt_comp_to_manual |>
    dplyr::mutate(
      cyt = factor(
        cyt,
        levels = c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6")
      ),
      pop = factor(
        pop,
        levels = c(
          "CD4 T cells", "CD8 T cells",
          "TCRgd T cells", "B cells", "NK cells"
        )
      )
    )

  plot_tbl <- plot_tbl |>
    dplyr::mutate(diff = freq_bs_auto - freq_bs_man) |>
    dplyr::arrange(pop, cyt, desc(abs(diff)))

  pop_vec <- unique(as.character(plot_tbl$pop))[c(2, 3, 1, 4, 5)]
  cyt_vec <- unique(as.character(plot_tbl$cyt))[c(1, 2, 6, 4, 5, 3)]
  p_points_list <- purrr::map(pop_vec, function(pop) {
    plot_tbl_pop <- plot_tbl |>
      dplyr::filter(pop == .env$pop)
    purrr::map(cyt_vec, function(cyt) {
      plot_tbl_pop_cyt <- plot_tbl_pop |>
        dplyr::filter(cyt == .env$cyt)
      p <- ggplot(
        plot_tbl_pop_cyt |>
          dplyr::mutate(
            freq_bs_man = pmax(freq_bs_man, 0.001),
            freq_bs_auto = pmax(freq_bs_auto, 0.001)
          ),
        aes(x = freq_bs_man, y = freq_bs_auto, col = stim)
      ) +
        cowplot::theme_cowplot(font_size = 9, font_family = "Helvetica") +
        geom_vline(xintercept = 0.01, linetype = "dotted", size = 0.4) +
        geom_hline(yintercept = 0.01, linetype = "dotted", size = 0.4) +
        cowplot::background_grid(major = "xy") +
        geom_abline(intercept = 0, slope = 1) +
        geom_point(size = 0.5) +
        labs(x = "Manual", y = "Auto") +
        scale_colour_manual(
          values = col_from_stim
        ) +
        theme(strip.background = element_blank()) +
        # labs(title = paste0(pop, "\n", cyt)) +
        theme(legend.position = "none") +
        theme(axis.title = element_blank()) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.8)) +
        scale_y_continuous(
          trans = "log10",
          labels = function(x) {
            scales::label_number()(x) |>
              stringr::str_remove("0$")
          }
        ) +
        scale_x_continuous(
          trans = "log10",
          labels = function(x) {
            scales::label_number()(x) |>
              stringr::str_remove("0$")
          }
        )
      UtilsGGSV::axis_limits(
        p = p,
        limits_expand = list(0.001),
        limits_equal = TRUE
      ) +
        coord_equal() +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    })
  }) |>
    purrr::flatten()

  p_points_grid <- cowplot::plot_grid(
    plotlist = p_points_list,
    ncol = length(unique(plot_tbl$cyt)),
    align = "hv",
    axis = "tblr",
  label_fontfamily = "Helvetica"
  )

  p_points_cyt <- ggplot(
    tibble::tibble(
      x = seq(0.1, 0.9, length.out = 6),
      y = 0.5,
      label = cyt_vec
    ),
    aes(x, y, label = label)
  ) +
    theme_void() +
    geom_text(hjust = 0.5, vjust = 0.5)

  p_points_pop <- ggplot(
    tibble::tibble(
      x = 0.5,
      y = seq(0.85, 0.25, length.out = 5),
      label = gsub(" cells", "", pop_vec)
    ),
    aes(x, y, label = label)
  ) +
    theme_void() +
    geom_text(angle = -90, vjust = 0.5, hjust = 0.5)

  p_points_stim <- ggplot(
    tibble::tibble(
      x = c(0, 0.15, 0.45, 0.825),
      y = 0.5,
      stim = names(col_from_stim)[-5],
      label = long_from_short_stim[names(col_from_stim)[-5]]
    )
  ) +
    theme_void() +
    geom_text(
      aes(x = x + 0.01, y = y, label = label),
      hjust = "left", size = 3
    ) +
    geom_point(aes(x = x, y = y, col = stim)) +
    scale_colour_manual(
      values = col_from_stim
    ) +
    theme(legend.position = "none") +
    lims(x = c(0, 1), y = c(0, 1))

  cowplot::ggdraw() +
    cowplot::draw_plot(
      plot = p_points_grid,
      x = 0.03,
      width = 0.93,
      y = 0.06,
      height = 0.905
    ) +
    cowplot::draw_plot(
      plot = p_points_cyt,
      x = 0.09,
      width = 0.855,
      y = 0.97,
      height = 0.03
    ) +
    cowplot::draw_plot(
      plot = p_points_pop,
      x = 0.98,
      width = 0.05,
      y = 0.14,
      height = 0.82,
      hjust = 0.5
    ) +
    cowplot::draw_text(
      text = c(
        "Automated gating frequency (%)",
        "Manual gating frequency (%)"
      ),
      x = c(0.015, 0.53),
      y = c(0.5375, 0.0585),
      size = 12,
      angle = c(90, 0)
    ) +
    cowplot::draw_plot(
      plot = p_points_stim,
      x = 0.0625,
      width = 0.925,
      y = 0,
      height = 0.035
    )
}

.get_corr_tbl_cyt <- function(data_tidy_cyt_comp_to_manual) {
  corr_tbl_cyt <- data_tidy_cyt_comp_to_manual |>
    dplyr::group_by(pop, cyt, stim) |>
    dplyr::filter(
      quantile(freq_stim_man, 0.75) > 3 * max(0.01, median(freq_uns_man)) &&
        quantile(freq_bs_man, 0.75) > 0.02
    ) |>
    dplyr::summarise(
      pcc = cor(freq_bs_auto, freq_bs_man),
      ccc = .calc_ccc(freq_bs_auto, freq_bs_man),
      .groups = "drop"
    )

  corr_tbl_cyt <- corr_tbl_cyt |>
    dplyr::mutate(
      pop = gsub(" T cells", "", pop),
      pop = gsub(" cells", "", pop),
    )

  corr_tbl_cyt <- corr_tbl_cyt |>
    dplyr::mutate(
      cyt = factor(cyt, levels = c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6")),
      pop = factor(pop, levels = c("CD4", "CD8", "TCRgd", "B", "NK"))
    )
}

.label_cytokine <- function(cyt) {
  purrr::map(
    cyt,
    function(cyt_ind) {
      switch(cyt_ind,
        "IFNg" = bquote(paste(plain(paste("IFN")), gamma)),
        cyt_ind
      )
    }
  )
}

.label_pop <- function(pop) {
  purrr::map(
    pop,
    function(pop_ind) {
      switch(pop_ind,
        "TCRgd" = bquote(paste(plain(paste("TCR")), gamma, delta)),
        pop_ind
      )
    }
  )
}

plot_p_cyt_corr_tbl_cyt <- function(corr_tbl_cyt) {
  col_vec <- RColorBrewer::brewer.pal(11, name = "RdBu") |> rev()
  p_pcc <- ggplot(
    corr_tbl_cyt |>
      dplyr::filter(stim != "p4") |>
      dplyr::mutate(stim = factor(
        .data$stim,
        levels = c("p1", "mtb", "ebv")
      )),
    aes(x = cyt, y = pop)
  ) +
    cowplot::theme_cowplot(line_size = 1, font_family = "Helvetica") +
    geom_raster(aes(fill = pcc)) +
    geom_text(aes(label = round(pcc, 2)), size = 2.25) +
    facet_wrap(
      ~stim,
      ncol = 3, scales = "fixed",
      labeller = labeller(
        stim = long_from_short_stim
      )
    ) +
    scale_fill_gradientn(
      colours = col_vec,
      values = seq(0, 1, length.out = length(col_vec)),
      na.value = "gray75",
      name = "Pearson\ncorrelation"
    ) +
    scale_y_discrete(breaks = levels(corr_tbl_cyt$pop), labels = .label_pop) +
    scale_x_discrete(breaks = levels(corr_tbl_cyt$cyt), labels = .label_cytokine) +
    labs(x = "Cytokine", y = "Population") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size = 9.5),
      legend.title = element_text(size = 10)
    )  +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    )

  p_ccc <- ggplot(
    corr_tbl_cyt |>
      dplyr::filter(stim != "p4") |>
      dplyr::mutate(stim = factor(
        .data$stim,
        levels = c("p1", "mtb", "ebv")
      )),
    aes(x = cyt, y = pop)
  ) +
    cowplot::theme_cowplot(line_size = 1, font_family = "Helvetica") +
    geom_raster(aes(fill = ccc)) +
    geom_text(aes(label = round(ccc, 2)), size = 2.25) +
    facet_wrap(
      ~stim,
      ncol = 3, scales = "fixed",
      labeller = labeller(
        stim = long_from_short_stim
      )
    ) +
    scale_fill_gradientn(
      colours = col_vec,
      values = seq(0, 1, length.out = length(col_vec)),
      na.value = "gray75",
      name = "Concordance\ncorrelation"
    ) +
    scale_y_discrete(breaks = levels(corr_tbl_cyt$pop), labels = .label_pop) +
    scale_x_discrete(breaks = levels(corr_tbl_cyt$cyt), labels = .label_cytokine) +
    labs(x = "Cytokine", y = "Population") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size = 9.5),
      legend.title = element_text(size = 10)
    ) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    )

  cowplot::plot_grid(
    p_pcc,
    p_ccc,
    align = "hv",
    axis = "ltbr",
    ncol = 1,
    labels = c("b", "c"),
  label_fontfamily = "Helvetica"
  )  +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
}

plot_fig_supp_cyt_auto_vs_manual <- function(p_cyt_corr,
                                             p_cyt_point) {
  # #
  cowplot::plot_grid(
    p_cyt_point,
    ggplot() +
      theme_void(),
    p_cyt_corr,
    labels = c("a", ""),
    ncol = 1,
    rel_heights = c(1, 0.03, 1),
    rel_widths = c(1, 1, 0.98),
    vjust = 1.5,
    label_fontfamily = "Helvetica"
  ) +
    theme(
      plot.background = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white")
    )
  
}

save_fig_supp_cyt_auto_vs_manual <- function(p,
                                             path_dir_save) {
  UtilsGGSV::ggsave2(
    file.path(path_dir_save, "p_supp_cyt"),
    plot = p,
    width = 18,
    height = 25,
    units = "cm"
  )
  c(
    "png" = paste0(file.path(path_dir_save, "p_supp_cyt"), ".png"),
    "pdf" = paste0(file.path(path_dir_save, "p_supp_cyt"), ".pdf")
  )
}