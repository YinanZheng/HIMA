## Welcome message and citations when package is loaded

.onAttach <- function(libname, pkgname) {
  pkg_version <- packageVersion("HIMA")

  packageStartupMessage(
    "Welcome to HIMA (version ", pkg_version, ")!\n",
    "Thank you for using the HIMA package for high-dimensional mediation analysis.\n",
    "To access full functionality, please ensure you are using the latest version.\n\n",
    "Getting Started:\n",
    "- Explore the tutorial at: https://YinanZheng.github.io/HIMA/articles/hima-vignette.html\n",
    "- For additional documentation, use ?hima.\n\n",
    "If you use HIMA in your work, please run citation(\"HIMA\") to see how to cite the package."
  )
}
