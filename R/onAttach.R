## Welcome message when package is loaded

.onAttach <- function(libname, pkgname) {
  pkg_version <- packageVersion("HIMA")

  packageStartupMessage(rep("*", 85))
  
  packageStartupMessage("HIMA version ", pkg_version, "\n",
                        "To access full functionality of HIMA, please make sure this version is current.")
  
  citation <- paste0("\nCitation:\n",
                     "  1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J,\n",
                     "     Just A, Colicino E, Vokonas P, Zhao L, Lv J, Baccarelli A, Hou L, Liu L.\n",
                     "     Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.\n",
                     "     Bioinformatics. 2016;32(20):3150-3154.\n",
                     "     PMID: 27357171; PMCID: PMC5048064.\n",
                     "\n",
                     "  2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L.\n",
                     "     Mediation Analysis for Survival Data with High-Dimensional Mediators.\n",
                     "     Bioinformatics. 2021;37(21):3815-3821.\n",
                     "     PMID: 34343267. PMCID: PMC8570823.")
  
  packageStartupMessage(citation)
  
  packageStartupMessage(rep("*", 85))
}
