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
    paste(
      "gateUse = 'ctrl' rows used only by legacy",
      "single-positive Adj control path"
    ),
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
    ),
    ".gateBatchSingle(): tgType = 'tg'",
    paste(
      "Legacy branch entered only when chnlSettings$gateTbl is non-NULL",
      "and gateName contains 'tg'"
    ),
    "Produces gateUse = 'gate' single-positive thresholds",
    "Can affect output only for legacy/manual single-positive workflows",
    "backwards-compatible but optional",
    paste0(
      "intermediateData/<stage>/<chnl>/ind/*/tg/",
      "{cpTgIndInit,cpTgList,cpTgPrejoin*}.rds"
    ),
    ".gateBatchSingle(): tgType = 'Adj'",
    paste(
      "Legacy branch entered only when chnlSettings$gateTbl is non-NULL",
      "and gateName contains 'Adj'"
    ),
    "Produces gateUse = 'ctrl' single-positive control thresholds",
    "Can affect output only for legacy/manual single-positive workflows",
    "backwards-compatible but optional",
    paste0(
      "intermediateData/<stage>/<chnl>/ind/*/Adj/",
      "{cpTgIndInit,cpTgList,cpTgPrejoin*}.rds"
    ),
    ".gateBatchSingle(): tgType = 'Clust'",
    paste(
      "Legacy branch entered only when chnlSettings$gateTbl is non-NULL",
      "and gateName contains 'Clust'"
    ),
    "Produces gateUse = 'tgClust' single-positive control thresholds",
    "Can affect output only for legacy/manual single-positive workflows",
    "backwards-compatible but optional",
    paste0(
      "intermediateData/<stage>/<chnl>/ind/*/Clust/",
      "{cpTgIndInit,cpTgList,cpTgPrejoin*}.rds"
    )
  )
}

.get_cp_tg_migration_note_157 <- function() {
  paste(
    "Issue #157 migration note:",
    paste(
      "the only .getCpTg() path still on by default is",
      "the init tgClust plumbing,"
    ),
    "but current cluster quantile gating does not consume gateTblCtrl.",
    "Remaining migration work is therefore the optional legacy single-positive",
    paste(
      "tg/Adj/Clust branches plus cleanup of",
      "tgClust dead plumbing and associated"
    ),
    "intermediate artifact expectations."
  )
}
