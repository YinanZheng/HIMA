## Re-submission
This is a re-submission of HIMA following its archival due to the archival of its upstream dependency hdi (which in turn depended on scalreg).

In version 2.3.4:
* Completely removed the dependency on hdi (and transitive dependency on scalreg).
* Implemented an internal, native de-biased LASSO inference routine (dlasso_proj) using glmnet.
* Added comprehensive unit test coverage with 	estthat.

## Test environments
* local Windows 10 x64, R 4.4.2

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.
