# Upload bodenmiller_bcr_xl_fs.rds to GitHub release
# Run this script once manually to upload the test data

# Check if gh package is available
if (!requireNamespace("gh", quietly = TRUE)) {
  stop("The 'gh' package is required. Install with: install.packages('gh')")
}

if (!requireNamespace("httr", quietly = TRUE)) {
  stop("The 'httr' package is required. Install with: install.packages('httr')")
}

# GitHub repository
repo <- "SATVILab/stimgate"
tag <- "test_data"
release_name <- "Test Data"

fs <- .get_fs()
rds_path <- file.path(tempdir(), "bodenmiller_bcr_xl_fs.rds")
saveRDS(fs, rds_path)

# Check if release exists, create if not
release_info <- tryCatch(
  {
    gh::gh("GET /repos/{owner}/{repo}/releases/tags/{tag}",
      owner = strsplit(repo, "/")[[1]][1],
      repo = strsplit(repo, "/")[[1]][2],
      tag = tag
    )
  },
  error = function(e) {
    message("Creating release 'test_data'")
    gh::gh("POST /repos/{owner}/{repo}/releases",
      owner = strsplit(repo, "/")[[1]][1],
      repo = strsplit(repo, "/")[[1]][2],
      tag_name = tag,
      name = release_name,
      body = "Test data for stimgate package examples and tests",
      draft = FALSE,
      prerelease = FALSE
    )
  }
)

if (!is.null(release_info)) {
  message("Release 'test_data' exists")
}

# Get upload URL
upload_url <- gsub("\\{\\?name,label\\}", "", release_info$upload_url)

# Upload the asset using httr
message("Uploading bodenmiller_bcr_xl_fs.rds...")
response <- httr::POST(
  url = paste0(upload_url, "?name=", basename(rds_path)),
  httr::add_headers(
    "Authorization" = paste("token", gh::gh_token()),
    "Content-Type" = "application/octet-stream"
  ),
  body = httr::upload_file(rds_path)
)

if (httr::status_code(response) %in% c(200, 201)) {
  message("Upload complete!")
} else {
  stop("Upload failed with status: ", httr::status_code(response))
}
