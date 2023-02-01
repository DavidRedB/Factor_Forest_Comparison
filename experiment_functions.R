
################################################################################
# This script is and adaptation from the code findable on:
# Goretzko y Bühner, 2020 ->  https://osf.io/mvrau/
# Goretzko y Bühner, 2022 ->  https://psyarxiv.com/u8eca/
# Auerswald y Moshagen, 2019 -> https://osf.io/gqma2/?view_only=d03efba1fd0f4c849a87db82e6705668
################################################################################


if (!require("lavaan")) {install.packages("lavaan"); library(lavaan)}
if (!require("data.table")) {install.packages("data.table"); library(data.table)}




# SMT APPROACH -----------------------------------------------------------------
fac_definition <- function(k = 8, dat) {
    # Input maximum number of factors, data
    # Output: definition of viable efa models as strings in a list for lavaan
    
    p <- ncol(dat)  # number of variables
    N <- nrow(dat)  # number of observations
    
    # Model definition 
    cx <- ' =~ '
    for(i in 1:p){
        cx <- paste(cx, sprintf('X%1.0d',i),sep='')
        if(i < p){ cx <- paste(cx,' + ',sep='') }
    }
    
    
    model_names <- list()
    for(i in 1:k) { 
        
        # Degrees of freedom check
        mod_pars <- p*i + (i*(i + 1) / 2) + p - i^2 # Loading matrix + Factor Var-Cov + Error - Lavaan Computational constraint (Ver Lorenzo-Seva et al, 2011)
        df <- (p*(p + 1) / 2) - mod_pars
        
        if(df <= 0) {break}
        
        
        # Model definition
        model_name <- ''
        if(i > 1) {
            for(j in 2:i) {
                model_name <- paste0(model_name, '+', 'efa("efa")*f', j)
            }
        }
        
        model_name <- paste0('efa("efa")*f1', model_name, cx)
        model_names[i] <- model_name
        
    }
    
    return(model_names)
    
}    


fac_estimation_v3 <- function(dat, models, use.estimator = "WLSMV", use.max_iters = 1000, use.max_time = 12) {
    # Input models' definition to fit in a list
    # Output: 
    #    comp_table: table with relevant info about the procedure
    #    pred: best model based on Satorra-Bentler or NA if none of the models were adequate
    #    pred_rmsea: best model based on RMSEA or NA if none of the models were adequate
    
    m <- length(models)
    comp_table <- data.frame(SBchisq = double(m), chisq.pvalue = double(m), RMSEA = double(m), RMSEA.pvalue = double(m), convergence = logical(m), TERMINATED = double(m)) 
    pred.chisq <- NA
    pred.rmsea <- NA
    
    
    
    if(m == 0) { # Error when no models are to be estimated
        comp_table[1, ] <- NA
        comp_table['TERMINATED'] <- -1
        return(list(comp_table = comp_table[1, ], pred.chisq = pred.chisq, pred.rmsea = pred.rmsea)) 
    }
    
    for(i in 1:m) { 
        
        
        fit <- check_fit(dat, models[[i]], use.estimator, use.max_iters, use.max_time)
        if(typeof(fit) == "S4") { # Checks if it is a lavaan model
            
            convergence <- inspect(fit, "converged")
            if(convergence == FALSE) {
                comp_table[i, 'convergence'] <- convergence
                next
            }
            
            
            measures <- fitmeasures(fit)
            comp_table[i, 'SBchisq'] <- measures["chisq.scaled"]
            comp_table[i, 'chisq.pvalue'] <- measures["pvalue.scaled"]
            comp_table[i, 'RMSEA'] <- measures["rmsea"]
            comp_table[i, 'RMSEA.pvalue'] <- measures["rmsea.pvalue"]
            comp_table[i, 'convergence'] <- convergence
            
            
            if(is.na(comp_table[i, 'chisq.pvalue']) == FALSE & comp_table[i, 'chisq.pvalue'] > 0.05){
                if(is.na(pred.chisq)){
                    pred.chisq <- i
                }
            }
            
            if(is.na(comp_table[i, 'RMSEA.pvalue']) == FALSE & comp_table[i, 'RMSEA.pvalue'] > 0.05){
                if(is.na(pred.rmsea)){
                    pred.rmsea <- i
                }
            }
            
            if(!is.na(pred.chisq) & !is.na(pred.rmsea)){
                comp_table[['TERMINATED']] <- 0
                return(list(comp_table = comp_table[1:max(pred.chisq, pred.rmsea), ], 
                            pred.chisq = pred.chisq, pred.rmsea = pred.rmsea ))
            }
            
        } else { 
            if(fit == "TIMEOUT") { # Could not estimate the model due to timeout
                comp_table[i, ] <- NA
                comp_table[['TERMINATED']] <- -50 
                return(list(comp_table = comp_table[1:i, ], 
                            pred.chisq = pred.chisq, pred.rmsea = pred.rmsea ))
            } else { # Unknown error
                comp_table[i, ] <- NA
                comp_table[['TERMINATED']] <- -100
                return(list(comp_table = comp_table[1:i, ], 
                            pred.chisq = pred.chisq, pred.rmsea = pred.rmsea ))
                
            }
        }
    } 
    comp_table['TERMINATED'] <- -25
    return(list(comp_table = comp_table, 
                pred.chisq = pred.chisq, pred.rmsea = pred.rmsea)) 
}


save_smt <- function(comp_table, job.id = 1) {
    # Saves smt data to an external file
    comp_table <-cbind(job_id = rep(job.id, dim(comp_table)[1]), comp_table)
    
    if(!file.exists('SMT_table.csv')) {
        write.table(comp_table, 'SMT_table.csv', sep = ";", dec = ',',
                    append = F, col.names = T, row.names = F)
    } else{
        write.table(comp_table, 'SMT_table.csv', sep = ";", dec = ",",
                    append = TRUE, col.names = F, row.names = F)
    }
}


check_fit <- function(dat, model, use.estimator, use.max_iters, use.max_time){
    
    setTimeLimit(elapsed = use.max_time)
    fit <- tryCatch(cfa(model=model,
                        data=dat,
                        estimator = use.estimator, rotation = "none", 
                        control=list(iter.max= use.max_iters)),
                    
                    error = function(e) {
                        
                        if(e[[1]] == "reached elapsed time limit") {
                            return("TIMEOUT")
                        } else {
                            return("ERROR")
                        }
                    })
    
    setTimeLimit(elapsed = Inf)
    return(fit)
}



# EXPERIMENT FUNCTIONS ---------------------------------------------------------


load_mat <- function(k, p, vpf, pf, sf){
    
    switch(pf,
           "low" = {
               pf_low = 0.35 
               pf_upper = 0.5
           },
           "medium" = {
               pf_low = 0.5
               pf_upper = 0.65
           },
           "high" = {
               pf_low = 0.65
               pf_upper = 0.8
           })
    switch(sf,
           "none" = {
               sf_low = 0 
               sf_upper = 0
           },
           "low" = {
               sf_low = 0
               sf_upper = 0.1
           },
           "medium" = {
               sf_low = 0.1
               sf_upper = 0.2
           })
    
    x <- runif(p, pf_low, pf_upper) 
    y <- runif(p*(k-1) , sf_low, sf_upper)
    
    i <- 1:(p)
    
    j <- rep(1:k, each=vpf) 
    
    L <- matrix(NA, p, k) 
    L[cbind(i, j)] <- x 
    L[is.na(L)] <- y
    L
}


simulation <- function(data, job,  # mandatory arguments # 
                           N,vpf, k, rho, pf, sf,  
                          normality # argument for the continuous distribution
){
    
    if(normality == "normal") {sk<-0; kr <- 0} else {sk <- 3; kr <- 8}
    p = vpf*k 
    L <- load_mat(k, p, vpf, pf, sf)
    fcor <- matrix(rho, k,k) + diag(1-rho,k,k)
    Sigma <- L%*%fcor%*%t(L) + diag(diag(diag(p)-L%*%fcor%*%t(L)))
    
    # simulate data with specific distribution depending on simulation condition
    dat <- tryCatch(data.frame(semTools::mvrnonnorm(n = N,
                                                    mu       = rep(0,p),
                                                    Sigma    = Sigma,
                                                    skewness = rep(sk,p),
                                                    kurtosis = rep(kr,p))),
                    error = function(e) {return(-500)})  # In case sigma is not  positive definite. 
    
    return(list(dat = dat, L = L, p = p, k = k, N = N, sk = sk, kr = kr))
}


analyze <- function(data, job,
                        instance, save = TRUE, ...){
    
    res <- data.frame(pa_solution = 0, guttman = 0,
                      xgb = 0, ekc = 0, cd = 0, SB = 0, RMSEA = 0) # DF to store results
    
    # Checking if simulate gave an error
    if(is.null(dim(instance$dat))){
        res[1, ] <- rep(-500, dim(res)[2])
        return(res)
    }
    
    # Starting the procedure
    dat_cor <- cor(instance$dat, use = "pairwise.complete.obs")
    eigval <- eigen(dat_cor)$values
    vareig <- cumsum(eigval)/instance$p
    
    # EKC - empirical kaiser criterion -----------------------------------------
    
    lref <- rep(0,instance$p)
    for (i in 1:instance$p) {
        lref[i] <- max(((1 + sqrt(instance$p/instance$N))^2) * (instance$p-sum(lref))/(instance$p-i+1),1)
    }
    ekc <- which(eigval<=lref)[1]-1
    
    # calculate eigenvalue based features
    
    eiggreater1 <- sum(eigval > 1)  
    releig1 <- eigval[1]/instance$p
    releig2 <- sum(eigval[1:2])/instance$p
    releig3 <- sum(eigval[1:3])/instance$p
    eiggreater07 <- sum(eigval > 0.7)
    sdeigval <- sd(eigval)
    var50 <- min(which(vareig > 0.50))
    var75 <- min(which(vareig > 0.75))
    
    # calculate matrix norm features
    
    onenorm <- norm(dat_cor,"O")
    frobnorm <- norm(dat_cor,"F")
    maxnorm <- norm(dat_cor-diag(instance$p),"M")
    avgcor <- sum(abs(dat_cor-diag(instance$p)))/(instance$p*(instance$p-1))
    specnorm <- sqrt(eigen(t(dat_cor)%*%dat_cor)$values[1])
    
    smlcor <- sum(dat_cor <= 0.1)
    avgcom <- mean(smc(dat_cor))
    det <- det(dat_cor)
    
    # "inequality" features 
    
    KMO <- KMO(dat_cor)$MSA
    Gini <- ineq(lower.tri(dat_cor), type = "Gini")
    Kolm <- ineq(lower.tri(dat_cor), type = "Kolm")
    
    # PA - parallel analysis ---------------------------------------------------
    
    pa <- fa.parallel(instance$dat,fa="fa",plot = FALSE, cor = "cor")
    pa_solution <- pa$nfact
    fa_eigval <- pa$fa.values
    
    # adding -1000 to shorter eigenvalue vectors
    
    eigval[(length(eigval)+1):80] <- -1000
    fa_eigval[(length(fa_eigval)+1):80] <- -1000
    names(eigval) <- paste("eigval", 1:80, sep = "")
    names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
    
    
    # CD - comparison data approach (with maximum number of factors = 8) -------
    
    cd <- EFA.Comp.Data(Data = instance$dat, F.Max = 8)
    
    
    
    # XGB ----------------------------------------------------------------------
    
    # Input for trained model
    
    features <- cbind(data.frame(k = instance$k,N=instance$N,p=instance$p,eiggreater1,
                                 releig1,releig2,releig3,eiggreater07,sdeigval,var50,
                                 var75,onenorm,frobnorm,maxnorm, avgcor,specnorm, smlcor,
                                 avgcom,det, KMO, Gini, Kolm, pa_solution, ekc, cd),
                      t(eigval), t(fa_eigval))
    
    print(features)
    # predict number of factors based on observed features
    
    pred_xgb <- predict(modxgb, newdata=features)
    
    
    # SMT approach ------------------------------------------------------------
    
    models <- fac_definition(dat = instance$dat)
    smt_res <- fac_estimation_v3(dat = instance$dat, models = models)
    if(save == TRUE) {save_smt(smt_res$comp_table, job.id = job$id)} # Saving smt_results
    
    SB <- smt_res$pred.chisq
    #RMSEA <- smt_res$pred.rmsea
    RMSEA <- NA
    # store results of parallel analysis, Kaiser-Guttman rule, factor forest,
    # empirical Kaiser criterion and comparison data approach)
    
    res[1, ] <- c(pa_solution, eiggreater1, as.numeric(pred_xgb$data$response),ekc, cd, SB, RMSEA)
    
    return(res)
}





