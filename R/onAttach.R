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
                     "     Bioinformatics. 2016;32(20):3150-4.\n",
                     "     PMID: 27357171; PMCID: PMC5048064.\n",
                     "\n",
                     "  2. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L.\n",
                     "     HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data.\n",
                     "     BMC Bioinformatics. 2022;23(1):296.\n",
                     "     PMID: 35879655 PMCID: PMC9310002.\n",
                     "\n",
                     "  3. Zhang H, Zheng Y, Hou L, Zheng C, Liu L.\n",
                     "     Mediation Analysis for Survival Data with High-Dimensional Mediators.\n",
                     "     Bioinformatics. 2021;37(21):3815-21.\n",
                     "     PMID: 34343267. PMCID: PMC8570823.\n",
                     "\n",
                     "  4. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L.\n", 
                     "     Mediation effect selection in high-dimensional and compositional microbiome data.\n",
                     "     Stat Med. 2021;40(4):885-96.\n",
                     "     PMID: 33205470. PMCID: PMC7855955.\n")
  
  packageStartupMessage(citation)
  
  packageStartupMessage(rep("*", 85))
}
