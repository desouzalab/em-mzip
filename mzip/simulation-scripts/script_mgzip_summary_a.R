source("mgzip_summary.R")
source("../../helper-functions/v_measure.R")

main_dir <- "~/scratch/Data/"

scen21 <- mgzip_summary(paste0(main_dir, "Scenario 21"))
vmat21 <- mgzip_summary_Z(paste0(main_dir, "Scenario 21"))

scen22 <- mgzip_summary(paste0(main_dir, "Scenario 22"))
vmat22 <- mgzip_summary_Z(paste0(main_dir, "Scenario 22"))

scen23 <- mgzip_summary(paste0(main_dir, "Scenario 23"))
vmat23 <- mgzip_summary_Z(paste0(main_dir, "Scenario 23"))

scen25 <- mgzip_summary(paste0(main_dir, "Scenario 25"))
vmat25 <- mgzip_summary_Z(paste0(main_dir, "Scenario 25"))

scen26 <- mgzip_summary(paste0(main_dir, "Scenario 26"))
vmat26 <- mgzip_summary_Z(paste0(main_dir, "Scenario 26"))

scen28 <- mgzip_summary(paste0(main_dir, "Scenario 28"))
vmat28 <- mgzip_summary_Z(paste0(main_dir, "Scenario 28"))


# Estimates for beta0 in Scenario 22, Case 1 (G = 12)

case_dir <- paste0(main_dir, "Scenario 22/Case 1")
source(paste0(case_dir, "/_Parameters.txt"))

n <- list.files(case_dir)
n <- regmatches(n, regexec("Data ([0-9]+) Y.dat", n))
n <- as.numeric(sapply(n, `[`, 2))
n <- sort(n[!is.na(n)])

b0_22_1 <- t(sapply(n, function(x) {
	scan(paste0(case_dir, "/EM ", x, " beta0.dat"), quiet = TRUE)
}))

save(scen21, scen22, scen23, scen25, scen26, scen28,
     vmat21, vmat22, vmat23, vmat25, vmat26, vmat28,
     b0_22_1,
     file = "summary_scen_21-28.RData")
cat("Summaries saved\n")
