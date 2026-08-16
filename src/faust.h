/*
 * faust.h — thin redirect to the stimgate package-local FAUST header.
 *
 * The vendored FAUST source files include "faust.h" by name.  This file
 * provides that name within the stimgate src/ directory and delegates to
 * stimgate_faust.h, which contains only the declarations actually required
 * by the vendored files.
 */
#include "stimgate_faust.h"
