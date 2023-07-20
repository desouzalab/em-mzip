source("../dzip.R")
source("../llmgzip.R")
source("../rmgzip.R")
source("../EMmgzip.R")
source("../NRmgzip.R")

source("mgzip_estimation.R")

mgzip_estimation("Scenario 21")
mgzip_estimation("Scenario 22")
mgzip_estimation("Scenario 23")

mgzip_estimation("Scenario 25", init = c(TRUE, "k-means"))
mgzip_estimation("Scenario 26")

mgzip_estimation("Scenario 28")
