#------------------------------------------------------------------------------#
# Main Module: runs a pre-specified simulation from start to end               #
# based on auxiliary modules and saves results on the local machine            #
#------------------------------------------------------------------------------#


# IMPORTING PACKAGES -----------------------------------------------------------

if (!require("batchtools")) {install.packages("batchtools"); library(batchtools)}
if (!require("psych")) {install.packages("psych"); library(psych)}
if (!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}
if (!require("ineq")) {install.packages("ineq"); library(ineq)}
if (!require("purrr")) {install.packages("purrr"); library(purrr)}
if (!require("mlr")) {install.packages("mlr"); library(mlr)}
if (!require("data.table")) {install.packages("data.table"); library(data.table)}
if (!require("BBmisc")) {install.packages("BBmisc"); library(BBmisc)}
if (!require("semTools")) {install.packages("semTools"); library(semTools)}
if (!require("matrixcalc")) {install.packages("matrixcalc"); library(matrixcalc)}
if (!require("xgboost")) {install.packages("xgboost"); library(xgboost)}


source("xgb-functions.R") 
source("CD_approach.R")
source("experiment_functions.R")
modxgb <- readRDS(file = "tunedxgb.rds") 


# SETTING CONDITIONS -----------------------------------------------------------

REPLS = 500 # 1 # 250
N     = c(200, 1000)
K     = c(1, 5)
VPF   = c(4, 10)
RHO   = c(0, 0.5)
PF    = c("low","high")
SF    = c("none", "medium")
normality    = c("normal", "non_normal")


cond <- CJ(N = N, vpf = VPF, k = K,
           rho = RHO, pf = PF,
           sf = SF,  sk = SK,
           kr = KR) 


cond <- cond[ !((cond$k == 1 & cond$rho == 0.5) | (cond$k == 1 & cond$sf == "medium")) ]


n_cond <- nrow(cond)
n_jobs <- nrow(cond) * REPLS
cond.id <- rep(1:nrow(cond), each = REPLS) 
simulation_design <- list(data_simulation = cond)


# SETTING UP THE EXPERIMENT  ------------------------------------------------------

reg <- makeExperimentRegistry("experiment", seed = 95558)
reg$packages <- c("psych","mvtnorm","ineq","purrr","mlr",
                  "data.table", "BBmisc", "semTools", "matrixcalc") 

reg$cluster.functions =  makeClusterFunctionsSocket() # Parallelization 



batchExport(export = list(load_mat = load_mat, EFA.Comp.Data = EFA.Comp.Data,
                          GenData = GenData, Factor.Analysis = Factor.Analysis,
                          modxgb = modxgb,
                          makeRLearner.classif.xgboost.earlystop = makeRLearner.classif.xgboost.earlystop,
                          trainLearner.classif.xgboost.earlystop = trainLearner.classif.xgboost.earlystop,
                          predictLearner.classif.xgboost.earlystop = predictLearner.classif.xgboost.earlystop,
                          createDMatrixFromTask = createDMatrixFromTask, fac_definition = fac_definition, fac_estimation_v3 = fac_estimation_v3, save_smt = save_smt, check_fit = check_fit))


addProblem("data_simulation", fun = simulation) 
addAlgorithm("analyze", fun = analyze)
addExperiments(prob.designs = simulation_design, repls = REPLS)


# SUBMITING THE EXPERIMENT -----------------------------------------------------

submitJobs(reg = reg, resources = list(walltime = 60))
waitForJobs()


# GETTING THE RESULTS ----------------------------------------------------------

capture.output(getStatus(reg = reg), file = "Status.txt")
capture.output(getErrorMessages(missing.as.error = FALSE, reg = reg), file = "Error_log.txt")


job.ids = findDone(reg = reg)
job.pars = unwrap(getJobPars(job.ids, reg = reg)) # Parameters per REP
job.pars = job.pars[, -c("problem", "algorithm")] 

raw_res = unwrap(reduceResultsDataTable(job.ids, fun = NULL, reg = reg)) 
write.table(raw_res, file = "raw_res.csv", sep = ";", dec = ",")


if(length(findDone(reg = reg)$job.id) != length(findJobs(reg = reg)$job.id)) {
    cond.id <- cond.id[-findJobs(reg = reg)$job.id[-findDone(reg = reg)$job.id]]
}
    
complete_res =  cbind(cond.id, rjoin(job.pars, raw_res))
write.table(complete_res, file = "complete_res.csv", sep = ";", dec = ",")



