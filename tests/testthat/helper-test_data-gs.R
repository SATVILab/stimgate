# Convert a flowSet to a GatingSet and save it to disk.
# This function creates a GatingSet from a flowSet by extracting individual frames,
# then saves it to the specified cache directory for use by StimGate.
#
# Arguments:
#   fs - A flowSet object containing flow cytometry data
#   dir_cache - Directory path where the GatingSet should be saved
#
# Returns:
#   Character string with the path to the saved GatingSet directory
get_gatingset <- function(fs,
                          dir_cache) {
  # Extract individual flowFrames from the flowSet
  frames_list <- lapply(seq_along(fs), function(i) fs[[i]])
  # Create new flowSet from frames
  fs2 <- flowCore::flowSet(frames = frames_list)
  # Convert flowSet to GatingSet
  gs <- flowWorkspace::GatingSet(fs2)
  # Define save path
  path_save <- file.path(dir_cache, "gs")
  # Remove existing GatingSet if present
  if (dir.exists(path_save)) {
    unlink(path_save, recursive = TRUE)
  }
  # Save GatingSet to disk
  flowWorkspace::save_gs(
    gs,
    path = path_save
  )
  path_save
}

# Create a complete GatingSet from scratch for a specific simulation scenario.
# This is a convenience wrapper that loads the base flowSet, generates simulated
# channel data for the chosen scenario, and creates a GatingSet ready for analysis.
#
# Arguments:
#   dir_cache - Directory path where the GatingSet should be cached
#   scenario - Either "easy" or "poor_separation", determining simulation parameters
#
# Returns:
#   A list containing:
#     - chnl_list: the simulated channel list
#     - batch_list: batch groupings from the first channel
#     - marker: vector of marker/channel names
#     - path_gs: path to the saved GatingSet
get_gatingset_from_scratch <- function (dir_cache, scenario) {
  # Load base flowSet
  fs <- get_fs()
  # Generate simulated data based on scenario
  chnl_list <- switch(scenario,
    "easy" = get_chnl_list_easy(fs = fs),
    "poor_separation" = get_chnl_list_poor_separation(fs = fs),
    stop("Unknown scenario: ", scenario)
  ) 
  # Extract the final flowSet from the last channel in the list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  # Clean up cache directory if it exists
  if (dir.exists(dir_cache)) {
    unlink(dir_cache, recursive = TRUE)
  }
  dir.create(dir_cache, recursive = TRUE)
  # Create and save GatingSet
  path_gs <- get_gatingset(
    fs = fs_gate,
    dir_cache = dir_cache
  )
  # Return all components needed for downstream analysis
  list(
    chnl_list = chnl_list,
    batch_list = chnl_list[[1]]$batch_list,
    marker = names(chnl_list),
    path_gs = path_gs
  )
}
