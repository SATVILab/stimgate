.setDev <- function() {
  Sys.setenv("PROJR_PROFILE" = "dev")
}
.unsetDev <- function() {
  profileVec <- projr::projr_profile_get()
  if ("dev" %in% profileVec) {
    profileVec <- profileVec[profileVec != "dev"]
    profileVec <- profileVec[profileVec != "default"]
    profileVec <- paste0(profileVec, collapse = ",")
    Sys.setenv("PROJR_PROFILE" = profileVec)
  }
}
