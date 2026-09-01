test_that("hima_dblasso works on continuous outcome data", {
  data(ContinuousOutcome)
  pheno_data <- ContinuousOutcome$PhenoData
  mediator_data <- ContinuousOutcome$Mediator
  
  fit <- hima_dblasso(
    X = pheno_data$Treatment,
    Y = pheno_data$Outcome,
    M = mediator_data,
    COV = pheno_data[, c("Sex", "Age")],
    scale = FALSE,
    FDRcut = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(fit) || is.null(fit))
  if (is.data.frame(fit)) {
    expect_true(all(c("Index", "alpha_hat", "beta_hat", "IDE", "pmax") %in% colnames(fit)))
    expect_true(nrow(fit) > 0)
  }
})

test_that("hima wrapper works with DBlasso penalty", {
  data(ContinuousOutcome)
  pheno_data <- ContinuousOutcome$PhenoData
  mediator_data <- ContinuousOutcome$Mediator
  
  fit <- hima(
    formula = Outcome ~ Treatment + Sex + Age,
    data.pheno = pheno_data,
    data.M = mediator_data,
    penalty = "DBlasso",
    scale = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "hima")
})

test_that("hima works with MCP penalty", {
  data(ContinuousOutcome)
  pheno_data <- ContinuousOutcome$PhenoData
  mediator_data <- ContinuousOutcome$Mediator
  
  fit <- hima(
    formula = Outcome ~ Treatment + Sex + Age,
    data.pheno = pheno_data,
    data.M = mediator_data,
    penalty = "MCP",
    scale = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "hima")
})

test_that("hima_survival works on SurvivalData", {
  data(SurvivalData)
  pheno_data <- SurvivalData$PhenoData
  mediator_data <- SurvivalData$Mediator
  
  fit <- hima_survival(
    X = pheno_data$Treatment,
    OT = pheno_data$Time,
    status = pheno_data$Status,
    M = mediator_data,
    COV = pheno_data[, c("Sex", "Age")],
    scale = FALSE,
    FDRcut = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(fit) || is.null(fit))
})

test_that("hima_microbiome works on MicrobiomeData", {
  data(MicrobiomeData)
  pheno_data <- MicrobiomeData$PhenoData
  mediator_data <- MicrobiomeData$Mediator
  
  fit <- hima_microbiome(
    X = pheno_data$Treatment,
    Y = pheno_data$Outcome,
    OTU = mediator_data,
    COV = pheno_data[, c("Sex", "Age")],
    FDRcut = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(fit) || is.null(fit))
})

test_that("hima_efficient works on ContinuousOutcome", {
  data(ContinuousOutcome)
  pheno_data <- ContinuousOutcome$PhenoData
  mediator_data <- ContinuousOutcome$Mediator
  
  fit <- hima_efficient(
    X = pheno_data$Treatment,
    Y = pheno_data$Outcome,
    M = mediator_data,
    COV = pheno_data[, c("Sex", "Age")],
    scale = FALSE,
    FDRcut = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(fit) || is.null(fit))
})


