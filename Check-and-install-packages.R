# Check and install the R packages required by scripts 01-11.
#
# Standard use from the project root:
#   source("Check-and-install-packages.R")
#
# Check only (do not install):
#   "D:/Program Files/R/R-4.4.1/bin/Rscript.exe" \
#     Check-and-install-packages.R --check-only

minimum_r_version <- "4.3.1"

cran_packages <- c(
  # Data input and manipulation
  "readxl", "dplyr", "stringr", "reshape2",
  # Normative modeling
  "gamlss", "gamlss.add", "gamlss.dist",
  # Parallel processing
  "foreach", "doParallel",
  # Classification and prediction
  "caret", "glmnet", "nnet", "pROC", "kernlab",
  # Prognostic analysis and harmonization
  "survival", "mgcv",
  # Visualization
  "ggplot2"
)

# ComBatFamily is distributed from its official GitHub repository rather than
# CRAN. It is used only by scripts/10_ComBatGAM_sensitivity.R.
github_packages <- c(
  ComBatFamily = "andy1764/ComBatFamily"
)

package_is_installed <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

package_status <- function() {
  packages <- c(cran_packages, names(github_packages))
  sources <- c(
    rep("CRAN", length(cran_packages)),
    rep("GitHub", length(github_packages))
  )
  installed <- vapply(packages, package_is_installed, logical(1))
  versions <- vapply(
    packages,
    function(package) {
      if (package_is_installed(package)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    },
    character(1)
  )

  data.frame(
    package = packages,
    source = sources,
    installed = installed,
    version = versions,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

install_required_packages <- function(install_missing = TRUE) {
  if (getRversion() < minimum_r_version) {
    warning(
      "R ", minimum_r_version, " or later is recommended; current version is ",
      as.character(getRversion()), ".",
      call. = FALSE
    )
  }

  repositories <- getOption("repos")
  if (is.null(repositories) ||
      length(repositories) == 0L ||
      !("CRAN" %in% names(repositories)) ||
      is.na(repositories[["CRAN"]]) ||
      identical(unname(repositories[["CRAN"]]), "@CRAN@")) {
    repositories["CRAN"] <- "https://cloud.r-project.org"
    options(repos = repositories)
  }

  initial_status <- package_status()
  cat("\nRequired-package status:\n")
  print(initial_status, row.names = FALSE)

  missing_cran <- cran_packages[!vapply(
    cran_packages, package_is_installed, logical(1)
  )]
  missing_github <- names(github_packages)[!vapply(
    names(github_packages), package_is_installed, logical(1)
  )]

  if (!install_missing) {
    all_available <- length(c(missing_cran, missing_github)) == 0L
    if (all_available) {
      message("\nAll required R packages are installed.")
    } else {
      message(
        "\nCheck-only mode: missing packages: ",
        paste(c(missing_cran, missing_github), collapse = ", ")
      )
    }
    return(invisible(initial_status))
  }

  if (length(missing_cran) > 0L) {
    message("\nInstalling missing CRAN packages: ", paste(missing_cran, collapse = ", "))
    utils::install.packages(missing_cran, dependencies = TRUE)
  }

  if (length(missing_github) > 0L) {
    if (!package_is_installed("remotes")) {
      message("\nInstalling the CRAN package 'remotes' for GitHub installation.")
      utils::install.packages("remotes", dependencies = TRUE)
    }
    if (!package_is_installed("remotes")) {
      stop(
        "The 'remotes' package could not be installed; ComBatFamily cannot be installed.",
        call. = FALSE
      )
    }

    for (package in missing_github) {
      repository <- unname(github_packages[[package]])
      message("\nInstalling ", package, " from GitHub: ", repository)
      remotes::install_github(
        repository,
        dependencies = TRUE,
        upgrade = "never",
        build_vignettes = FALSE
      )
    }
  }

  final_status <- package_status()
  still_missing <- final_status$package[!final_status$installed]

  cat("\nFinal required-package status:\n")
  print(final_status, row.names = FALSE)

  if (length(still_missing) > 0L) {
    stop(
      "Installation was incomplete. Missing packages: ",
      paste(still_missing, collapse = ", "),
      ". Review the installation messages above; on Windows, source packages may require Rtools.",
      call. = FALSE
    )
  }

  message("\nAll required R packages are installed.")
  invisible(final_status)
}

check_only <- "--check-only" %in% commandArgs(trailingOnly = TRUE) ||
  identical(Sys.getenv("NCOMMS_CHECK_ONLY"), "1")

install_required_packages(install_missing = !check_only)
