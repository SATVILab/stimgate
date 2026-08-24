.comp_against_manual_cyt_format_manual <- function(fn, pop, cyt) {
  # load
  manual_tbl <- suppressMessages(suppressWarnings(readr::read_csv(
    projr::projr_path_get(
      "raw-data-small",
      "comparison_data",
      "acscytof",
      fn
    )
  )))

  # remove empty columns
  manual_tbl <- manual_tbl |>
    dplyr::select(-(...58:...62))

  # get info
  manual_tbl <- manual_tbl |>
    dplyr::bind_cols(
      stringr::str_split(manual_tbl$`...1`, "_") |>
        purrr::map_df(function(x) {
          tibble::tibble(
            SampleID = paste0(x[1], "_", x[2]),
            stim = x[5]
          )
        })
    ) |>
    dplyr::select(SampleID, stim, everything()) |>
    dplyr::select(-`...1`)

  # remove "_FreqofParent" from column name
  colnames(manual_tbl) <- stringr::str_remove(
    colnames(manual_tbl),
    " FreqofParent"
  )

  # pivot longer
  manual_tbl <- manual_tbl |>
    tidyr::pivot_longer(
      -c(SampleID, stim),
      names_to = "pop",
      values_to = "freq"
    )

  # split pop into pop and cyt
  manual_tbl <- manual_tbl |>
    tidyr::separate(
      col = pop,
      sep = "/",
      into = c("pop", "cyt")
    )

  # remove Perf
  manual_tbl <- manual_tbl |>
    dplyr::filter(!cyt == "Perf")

  # rename mtbaux as mtb
  manual_tbl <- manual_tbl |>
    dplyr::mutate(stim = ifelse(stim == "mtbaux", "mtb", stim))

  # rename freq as freq_stim
  manual_tbl <- manual_tbl |>
    dplyr::rename(freq_stim = freq)

  # subtract unstim
  manual_tbl <- manual_tbl |>
    dplyr::group_by(SampleID, pop, cyt) |>
    dplyr::mutate(freq_uns = freq_stim[stim == "uns"]) |>
    dplyr::ungroup() |>
    dplyr::mutate(freq_bs = freq_stim - freq_uns) |>
    dplyr::select(
      SampleID,
      stim,
      pop,
      cyt,
      freq_stim,
      freq_uns,
      freq_bs
    ) |>
    dplyr::rename(
      freq_bs_man = freq_bs,
      freq_stim_man = freq_stim,
      freq_uns_man = freq_uns
    )

  # set freq_bs <0 to 0
  manual_tbl <- manual_tbl |>
    dplyr::mutate(freq_bs_man = pmax(freq_bs_man, 0))

  # rename TCRgd+ to TCRgd T cells
  manual_tbl <- manual_tbl |>
    dplyr::mutate(
      pop = ifelse(
        pop == "TCRgd+",
        "TCRgd T cells",
        pop
      )
    )

  # filter pop, if not NULL
  if (!is.null(pop)) {
    manual_tbl <- manual_tbl |> dplyr::filter(pop %in% .env$pop)
  }

  # filter cyt, if not NULL
  if (!is.null(cyt)) {
    manual_tbl <- manual_tbl |> dplyr::filter(cyt %in% .env$cyt)
  }

  manual_tbl
}
