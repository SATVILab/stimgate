# Summary of Changes: Remove .debug Parameters and Add stage Parameter

## Overview

This PR removes all `.debug` parameters from function definitions and calls throughout the stimgate R package, replacing them with environment variable-based debugging. It also adds a `stage` parameter system to track algorithm stages and enable intermediate data saving.

## Key Changes

### 1. Removed `.debug` Parameters

**Affected Files:**
- `R/gate.R` - Main entry point `stimgate_gate()` function
- `R/stats.R` - `.get_stats()` internal function
- `R/cp_cluster.R` - `.get_cp_cluster()` and related functions
- `R/cp_cluster-helper.R` - Multiple helper functions (22 function signatures updated)
- `R/stats-helper-overall.R` - Statistics helper functions
- `R/gate_marker-helper.R` - Gate marker helper functions
- `R/gate_batch.R` - Batch gating functions
- `R/gate_batch-helper.R` - Batch helper functions

**Changes Made:**
- Removed `.debug = FALSE` from all function signatures
- Removed all `.debug = .debug` parameter passing
- **Preserved all `.debug("message")` function calls** - these are the actual debug logging functions

**New Debug Mechanism:**
- Debugging is now controlled entirely by the `STIMGATE_DEBUG` environment variable
- Set `STIMGATE_DEBUG="TRUE"` (or "yes", "y", "1") to enable debug output
- The `stimgate_gate()` function automatically initializes the environment variable if not set
- All debug logging via `.debug()` function checks the environment variable automatically

### 2. Added `stage` Parameter System

**Valid stage values:** `"init"`, `"cyt_pos"`, `"single"`

**Affected Files:**
- `R/gate_batch.R` - Added `stage` parameter to `.gate_batch()`
- `R/gate_batch-helper.R` - Added `stage` to `.gate_batch_all()` and `.gate_batch_single()`
- `R/gate_marker-helper.R` - Added `stage` to `.gate_marker_pre_adj_gates_gate()` and fixed typo (`stagae` → `stage`)
- `R/gate_marker.R` - Passes `stage` to helper functions
- `R/cyt_pos_gates.R` - Added `stage` to all functions in the cyt_pos gate chain
- `R/gate.R` - Passes `stage = "cyt_pos"` to `.gate_cyt_pos()`

**Stage Flow:**
1. **"init" stage**: `.gate_init()` → `.gate_marker()` → `.gate_marker_pre_adj_gates_gate()` → `.gate_batch()` → `.get_cp_uns_loc()`
2. **"single" stage**: `.gate_single()` → `.gate_marker()` → `.gate_marker_pre_adj_gates_gate()` → `.gate_batch()` → `.get_cp_uns_loc()`
3. **"cyt_pos" stage**: `.gate_cyt_pos()` → `.get_cyt_pos_gates_gate_name()` → `.get_cyt_pos_gates_ind()` → `.get_cp_pos_gates_chnl()` → `.int_save_nm()`

**Purpose:**
- Enables tracking of which algorithm stage is executing
- Allows intermediate data saving organized by stage
- Used with `.int_save()` and `.int_save_nm()` functions
- Intermediate saving controlled by `STIMGATE_INTERMEDIATE` environment variable

### 3. Updated Documentation

**Documentation Changes:**
- Removed `.debug` parameter from `stimgate_gate()` function documentation
- Created new documentation files:
  - `man/dot-debug.Rd` - Documents the `.debug()` logging function
  - `man/dot-debug_file_create.Rd` - Documents debug file creation
  - `man/dot-debug_file_get_path.Rd` - Documents debug file path retrieval
  - `man/stimgate_debug_copy.Rd` - Documents copying debug file to working directory
  - `man/stimgate_debug_print.Rd` - Documents printing debug output
- Deleted `man/dot-debug_msg.Rd` (obsolete function)

### 4. Updated Copilot Instructions

**File:** `.github/copilot-instructions.md`

**Changes:**
- **Guideline #2**: Updated to describe using `.debug()` function with `STIMGATE_DEBUG` environment variable (instead of `.debug` parameter)
- **Guideline #3**: Added guidance on `stage` parameter usage and intermediate data saving
- **Required before commit**: Reordered to run `devtools::document()` first, then `styler::style_pkg()`, then tests

### 5. Updated Tests

**File:** `tests/testthat/test-run.R`

**Changes:**
- Updated test that previously used `debug = TRUE` parameter
- Now sets `STIMGATE_DEBUG="TRUE"` environment variable instead
- Uses `on.exit(Sys.unsetenv("STIMGATE_DEBUG"))` for cleanup

### 6. Code Formatting

- Ran `styler::style_pkg()` to ensure consistent code formatting
- 14 files reformatted with minor whitespace changes
- All trailing whitespace removed

## Validation

✅ **All R source files parse successfully** (29 files checked)
✅ **All test files parse successfully** (10 files checked)  
✅ **Documentation updated successfully** via `roxygen2::roxygenise()`
✅ **Code formatted successfully** via `styler::style_pkg()`

## Migration Guide

### For Users

**Before:**
```r
stimgate_gate(
  .data = gs,
  path_project = path_project,
  batch_list = batch_list,
  marker = marker,
  .debug = TRUE
)
```

**After:**
```r
# Set environment variable before calling
Sys.setenv("STIMGATE_DEBUG" = "TRUE")

stimgate_gate(
  .data = gs,
  path_project = path_project,
  batch_list = batch_list,
  marker = marker
)
```

### For Developers

**Before:**
```r
.my_function <- function(param1, param2, .debug = FALSE) {
  .debug("Processing data") # This won't work without proper setup
  # ... function body ...
  result <- .helper_function(param1, .debug = .debug)
  result
}
```

**After:**
```r
.my_function <- function(param1, param2) {
  .debug("Processing data") # Automatically checks STIMGATE_DEBUG
  # ... function body ...
  result <- .helper_function(param1)
  result
}
```

**For stage parameter:**
```r
.my_function <- function(param1, param2, stage, path_project) {
  # Save intermediate data
  .int_save(ind, stage, path_project, 
    important_result, 
    intermediate_data
  )
  # ... rest of function ...
}
```

## Benefits

1. **Cleaner API**: No need to pass `.debug` through multiple function levels
2. **Better Control**: Debug mode can be enabled globally via environment variable
3. **Easier Testing**: Tests can enable/disable debug mode without function signature changes
4. **Stage Tracking**: Algorithm stages are now explicit and traceable
5. **Intermediate Data**: Organized saving of intermediate results by stage and sample index
6. **Maintainability**: Reduced parameter clutter in internal functions

## Files Changed

**Total: 33 files modified**

- 22 R source files
- 1 copilot instructions file
- 8 documentation files (man/)
- 2 test files
- 2 old reference scripts (data-raw/)

**Lines Changed:**
- 198 insertions
- 159 deletions
- Net: +39 lines (mostly new documentation)
