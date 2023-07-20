source("mzip_summary.R")
source("../../helper-functions/v_measure.R")

scen11 <- mzip_summary("Scenario 11")
vmat11 <- mzip_summary_Z("Scenario 11")

scen12 <- mzip_summary("Scenario 12")
vmat12 <- mzip_summary_Z("Scenario 12")

scen13 <- mzip_summary("Scenario 13")
vmat13 <- mzip_summary_Z("Scenario 13")

scen14 <- mzip_summary("Scenario 14")
vmat14 <- mzip_summary_Z("Scenario 14")

scen15 <- mzip_summary("Scenario 15")
vmat15 <- mzip_summary_Z("Scenario 15")

scen16 <- mzip_summary("Scenario 16")
vmat16 <- mzip_summary_Z("Scenario 16")

save(scen11, scen12, scen13, scen14, scen15, scen16,
	 vmat11, vmat12, vmat13, vmat14, vmat15, vmat16,
	 file = "summary_scen_11-16.RData")
cat("Summaries saved\n")
