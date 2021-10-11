# For Section 5.1 Evaluation and Comparison with Average Treatment Effect
# and Section 5.2 Evaluation of Multiple Constraints 

#Load R Packages

library(foreach)
library(doParallel)
library(policytree)
library(randomForestSRC)
#install.packages('devtools')
#devtools::install_github("grf-labs/policytree", subdir = "r-package/policytree") #install.packages('policytree')

#for multiple covariates
#simulation
#pre-setting
repnum = 200

info = matrix(0, nrow=repnum, ncol=42)
# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl, cores=cores)

info <- foreach(ith=1:repnum, .combine='rbind', .packages=c("policytree",'randomForestSRC')) %dopar%
{ 
    set.seed(2357)
    seeds <- ceiling(runif(repnum, 10000, 1e+09))
    set.seed(seeds[ith])
    cat("Reps", ith, "\n") 
     

    n.list = c(200, 500, 1000) # sample sizes 
    K = length(n.list)
    size_heu = size_hat = effect_hat = cond_flag_hat = size_mc = effect_mc = pcd_mc = pos_mc = cond_flag_mc = matrix(0, nrow=1, ncol=K)

    # Set Scenario
    
    scen_id = 3
    del_id = 1

    p = 10 # dimension of X 
    pi = 0.5 #propensity score 
    
    for(k in 1:length(n.list)){

        n = n.list[k]

        x = matrix(runif(n*p, -2, 2), nrow=n, ncol=p)
        a = rbinom(n, 1, pi)
        
        if(scen_id == 1){
            y = x[, 1] + 2 * x[, 2] + a * (x[, 1]) + rnorm(n, 0, 1)
            if(del_id == 1){
                delta = 0.7
                cut = -0.6
            }
            if(del_id == 2){ 
                delta = 1
                cut = 0
            }    
            if(del_id == 3){ 
                delta = 1.3
                cut = 0.6  
            }   
        }
        if(scen_id == 2){
            y = x[, 1] + 2 * x[, 2] + a * (x[, 1] * x[, 2]) + rnorm(n, 0, 1)
            if(del_id == 1){
                delta = 0.7
                cut = -0.42
            }
            if(del_id == 2){ 
                delta = 1
                cut = 0
            }    
            if(del_id == 3){ 
                delta = 1.3
                cut = 0.28 
            }             
        }        
        if(scen_id == 3){
            y = x[, 1] + 2 * x[, 2] + a * (x[, 1] + x[, 2]) + rnorm(n, 0, 1)
            if(del_id == 1){
                delta = 0.7
                cut = -1.18
            }
            if(del_id == 2){ 
                delta = 1
                cut = -0.55
            }    
            if(del_id == 3){ 
                delta = 1.3
                cut = -0.05  
            }           
        }                 
         
        df = cbind(x, a, y) # dataset

        ###random forest fit the ATE
        Y1 = y[a==1]
        X1 = x[a==1,]
        data1 = data.frame(cbind(Y1,X1))
        fit1 = rfsrc(Y1~.,data=data1)
        Y0 = y[a==0]
        X0 = x[a==0,]
        data0 = data.frame(cbind(Y0,X0))
        fit0 = rfsrc(Y0~.,data=data0)

        mu11 = predict(fit1)
        mu10 = predict(fit1,data0)
        mu01 = predict(fit0,data1)
        mu00 = predict(fit0)
        ## $predicted.oob: out of bag prediction
        tau1 = mu11$predicted.oob - mu01$predicted
        tau0 = mu10$predicted - mu00$predicted.oob

        tau = numeric(length(y))
        tau[a==1] = tau1
        tau[a==0] = tau0 
        
        #calculate the rank and reward
        r_sorted_list = sort(tau - delta, decreasing=TRUE)
        r_rank_list = rank(tau - delta) 
        cumu_r_sorted_list = cumsum(r_sorted_list) / seq(1, n, 1)
        
        size_heu[k] = sum(cumu_r_sorted_list>0) / n
        
        Gamma_mat = matrix(0, nrow=n, ncol=2)
        
        # outside subgroup
        Gamma_mat[,1] = 0 # outside subgroup
        
        # in subgroup
        # Use sign of cumulative r as Reward 1
        Gamma_mat[,2] = 2 * (cumu_r_sorted_list[n - r_rank_list + 1] >= 0) - 1  
        
        # Use value of cumulative r as Reward 2
        #Gamma_mat[,2] = cumu_r_sorted_list[n - r_rank_list + 1]
        
        # Use penalty as Reward 3
        # lambda = 0.1 penalty size
        #Gamma_mat[,2] = Gamma_mat[,2] + lambda * (tau < 0) * tau
        
        #find the ODR for subgroup
        tree <- policy_tree(x, Gamma_mat, depth = 2) 
        
        size_hat[k] = sum((predict(tree, x) == 2)) / n
        effect_hat[k] = sum((predict(tree, x) == 2) * y * ((a == 1) / pi - (a == 0) / (1 - pi))) / sum((predict(tree, x) == 2))
        cond_flag_hat[k] = effect_hat[k] >= delta
        
        
        NN = 50000

        X = matrix(runif(NN*p, -2, 2), nrow=NN, ncol=p) 
        
        est_rule = (predict(tree, X) == 2)
        size_mc[k] = sum(est_rule) / NN
        
        if(scen_id == 1){
            effect_mc[k] = sum(est_rule *(X[, 1])) / sum(est_rule)
            pcd_mc[k] = mean(est_rule == (X[, 1] >= cut)) 
            pos_mc[k] = sum((est_rule == 1) * (X[, 1] >= 0)) / sum(est_rule == 1)
        }
        if(scen_id == 2){
            effect_mc[k] = sum(est_rule *(X[, 1] * X[, 2])) / sum(est_rule)
            pcd_mc[k] = mean(est_rule == (X[, 1] * X[, 2] >= cut))      
            pos_mc[k] = sum((est_rule == 1) * (X[, 1] * X[, 2] >= 0)) / sum(est_rule == 1)
        }        
        if(scen_id == 3){
            effect_mc[k] = sum(est_rule *(X[, 1] + X[, 2])) / sum(est_rule)
            pcd_mc[k] = mean(est_rule == (X[, 1] + X[, 2] >= cut))          
            pos_mc[k] = sum((est_rule == 1) * (X[, 1] + X[, 2] >= 0)) / sum(est_rule == 1)
        }         
         
        cond_flag_mc[k] = effect_mc[k] >= delta 
         
        
        cat("n=",n,"size_pct=",size_mc[k],"effect=",effect_mc[k],'satisfies?',pcd_mc[k],"\n")         
      
        
      
    }
    cat("Repe",ith,"\n") 
    c(size_mc,effect_mc,pcd_mc,pos_mc,cond_flag_mc,size_hat,effect_hat,cond_flag_hat,size_heu) 
    
}
stopImplicitCluster()
stopCluster(cl)


save(info,file="Res.RData")
