# Audit remaining .getCpTg() call sites
# Internal inventory of all remaining `.getCpTg()` call sites, their activation
# switches, downstream consumers, and whether they can affect current package
# output.
.get_cp_tg_call_audit <- function() {
  tibble::tribble(
    ~callSite,
    ~activation,
    ~downstreamConsumer,
    ~observableOutput,
    ~classification,
    ~intermediateArtifacts,
    ".gateBatchAll(): tgType = 'tolCtrl'",
    "chnlSettings$tolCtrl is non-NULL",
    "gateUse = 'ctrl' rows",
    "No effect on gateTblInit.rds or gateTbl.rds in current default flow",
    "backwards-compatible but optional",
    paste0(
      "intermediateData/<stage>/<chnl>/ind/*/tolCtrl/",
      "{cpTgIndInit,cpTgList,cpTgPrejoin*}.rds"
    ),
    ".gateBatchAll(): tgType = 'tolClust'",
    "chnlSettings$tolClust is non-NULL",
    "gateUse = 'tgClust' rows passed as gateTblCtrl into .getCpCluster()",
    paste(
      "Current .getCpCluster() Q15/Q60/Q85 local-FDR",
      "path does not read gateTblCtrl"
    ),
    "dead for current output",
    paste0(
      "intermediateData/<stage>/<chnl>/ind/*/tolClust/",
      "{cpTgIndInit,cpTgList,cpTgPrejoin*}.rds"
    )
  )
}

.get_cp_tg_migration_note_157 <- function() {
  paste(
    "Issue #157 migration note:",
    "the obsolete default init tgClust tailgate path has been removed.",
    "Current local-FDR cluster quantile gating does not consume gateTblCtrl,",
    "so the legacy tgClust plumbing is dead for current outputs.",
    "Single-positive branches have been removed per issue #196."
  )
}
