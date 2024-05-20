# old determination of the # of clusters

dens_mat <- dens_tbl |>
  dplyr::select(-c(batch_sh, stim, ind)) |>
  as.matrix()
dist_mat <- dist(
  dens_mat,
  method = "manhattan"
)
hclust_obj <- hclust(dist_mat, method = "ward.D2")
height_vec <- hclust_obj$height
height_vec_diff <- diff(height_vec)
max_height_ind <- max(which(height_vec_diff < 0.5 * max(height_vec_diff)))
height <- height_vec[max_height_ind + 1]
cut_vec <- cutree(hclust_obj, h = height)
max_cluster <- min(6, length(gs) / 3) |>
  floor() |>
  max(1)



n_clus <- x |>
  stats::setNames(1:max_cluster) |>
  purrr::map_dbl("gap")

if (length(unique(cut_vec)) > max_cluster) {
  cut_vec <- cutree(hclust_obj, k = max_cluster)
}

n_grp <- length(unique(cut_vec))


# old calculating cp_join

cp_tbl <- prop_bs_by_cp_tbl |>
  tidyr::pivot_longer(
    names_to = "n_grp", values_to = "grp_level",
    -c(ind:grp)
  ) |>
  dplyr::mutate(n_grp_grp_level = paste0(
    n_grp, "_", grp_level
  )) |>
  dplyr::mutate(cp_grp_level = cp_grp_lab_vec[
    n_grp_grp_level
  ]) |>
  # dplyr::select(n_grp_grp_level, cp_grp_level) |>
  # dplyr::filter(!(is.na(cp_grp_level) | cp_grp_level == Inf)) |>
  dplyr::group_by(ind) |>
  dplyr::mutate(cp_grp_level = min(cp_grp_level, na.rm = TRUE)) |>
  # dplyr::slice(1) |>
  dplyr::rename(cp_join = cp_grp_level) |>
  # dplyr::mutate(
  # cp_join_der = cp_grp_level[cp_grp_level ==
  # min(cp_grp_level, na.rm = TRUE)][1]
  # ) |>
  dplyr::ungroup()


# old getting cp_tg

cp_tbl |>
  dplyr::group_by(ind) |>
  dplyr::mutate(
    cp_join_tg = min(cp[
      cp >= cp_join & cp >= cp_tg_ctrl &
        prop_bs_cp_diff_sd <= 2
    ]),
    cp_join_tg = ifelse(
      is.na(cp_join_tg), cp_join_lse, cp_join_tg
    ),
    cp_join_tg_orig = pmin(cp_join_tg, cp_orig_quant_min)
  ) |>
  # dplyr::select(ind, cp, prop_bs_cp_diff_sd,
  #               cp_orig_proper_min ,cp_join,
  #               cp_tg_ctrl:cp_join_lse_orig,
  #               cp_join_tg, cp_join_tg_orig) |>
  # dplyr::filter(cp >= cp_tg_ctrl) |>
  # dplyr::select(grp, ind, cp_join, cp_join_ind, cp_orig_quant_min) |>
  dplyr::ungroup()
