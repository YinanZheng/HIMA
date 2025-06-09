test_that("hima_survival runs when d is 1", {
  data(SurvivalData)
  pheno <- SurvivalData$PhenoData
  mediator <- SurvivalData$Mediator
  expect_silent(
    hima_survival(
      X = pheno$Treatment,
      OT = pheno$Time,
      status = pheno$Status,
      M = mediator,
      COV = pheno[, c("Sex", "Age")],
      topN = 1,
      scale = FALSE,
      FDRcut = 1,
      verbose = FALSE,
      parallel = FALSE
    )
  )
})
