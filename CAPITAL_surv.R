# For Section 5.3 Evaluation of Survival Data

#Load R Packages

library(foreach)
library(doParallel)
library(policytree)
library(survival)
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

info <- foreach(ith=1:repnum, .combine='rbind', .packages=c("policytree",'randomForestSRC','survival')) %dopar%
{ 
    set.seed(2357)
    seeds <- ceiling(runif(repnum, 10000, 1e+09))
    set.seed(seeds[ith])
    cat("Reps", ith, "\n") 
     

    n.list = c(500, 1000) # sample sizes 
    K = length(n.list)
    size_heu = size_hat = size_mc = effect_mc = pcd_mc = cond_flag_mc = matrix(0, nrow=1, ncol=K)
    
    # Set Scenario
    
    noise_id = 3
    cencor_id = 2
        
    p = 10 # dimension of X 
    pi = 0.5 #propensity score 
    
    for(k in 1:length(n.list)){

        n = n.list[k]

        x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)
        a = rbinom(n, 1, pi)
    
        if(noise_id == 1){
            y = 0.1 * x[, 1] + 0.2 * x[, 2] + a * (x[, 1]) + rnorm(n, 0, 1)
            if(cencor_id == 1){
                L = 12
                delta = 1.068  
            }
            if(cencor_id == 2){ 
                delta = 0.864
                L = 7  
            }            
        }
        if(noise_id == 2){
            y = 0.1 * x[, 1] + 0.2 * x[, 2] + a * (x[, 1]) + rlogis(n)
            if(cencor_id == 1){
                L = 21
                delta = 1.337 
            }
            if(cencor_id == 2){ 
                delta = 0.874
                L = 10  
            }            
        }        
        if(noise_id == 3){
            y = 0.1 * x[, 1] + 0.2 * x[, 2] + a * (x[, 1]) + log(-log(1-runif(n)))
            if(cencor_id == 1){
                L = 8
                delta = 0.730
            }
            if(cencor_id == 2){ 
                delta = 0.544
                L = 4
            }            
        }         
        
        t = exp(y)
        c = runif(n, 0, L)
        time = apply(cbind(t,c), 1, min)

        status = 1*(t <= c)
        
        df = cbind(x, a, time, status) # dataset

        ###random forest fit the Mean Survival
        time1 = time[a==1]
        status1 = status[a==1]
        X1 = x[a==1,]
        data1 = data.frame(cbind(time1,status1,X1))
        fit1 = rfsrc(Surv(time1,status1)~., data1)
        time0 = time[a==0]
        status0 = status[a==0]
        X0 = x[a==0,]
        data0 = data.frame(cbind(time0,status0,X0))
        fit0 = rfsrc(Surv(time0,status0)~., data0)
        #plot.survival(fit0)
        
        
        mu11 = predict(fit1)
        mu10 = predict(fit1,data0)
        mu01 = predict(fit0,data1)
        mu00 = predict(fit0)

        time_intl11 = (c(mu11$time.interest[2:(length(mu11$time.interest))],L) - mu11$time.interest)
        time_intl10 = (c(mu10$time.interest[2:(length(mu10$time.interest))],L) - mu10$time.interest) 
        time_intl01 = (c(mu01$time.interest[2:(length(mu01$time.interest))],L) - mu01$time.interest) 
        time_intl00 = (c(mu00$time.interest[2:(length(mu00$time.interest))],L) - mu00$time.interest) 


        aoc11 = function(xxx){ 
            sum(xxx * time_intl11)
        }

        aoc10 = function(xxx){ 
            sum(xxx * time_intl10)
        }

        aoc01 = function(xxx){ 
            sum(xxx * time_intl01)
        }

        aoc00 =function(xxx){ 
            sum(xxx * time_intl00)
        }
 
 
        ## $survival.oob: out of bag prediction
        tau1 = apply(mu11$survival.oob, 1, aoc11) - apply(mu01$survival, 1, aoc01)
            #mu11$predicted.oob - mu01$predicted
        tau0 = apply(mu10$survival, 1, aoc10) - apply(mu00$survival.oob, 1, aoc00)
            #mu10$predicted - mu00$predicted.oob

        tau = numeric(length(y))
        tau[a==1] = tau1
        tau[a==0] = tau0 
        
        #calculate the rank and reward
        r_sorted_list = sort(tau - delta, decreasing=TRUE)
        r_rank_list = rank(tau - delta) 
        cumu_r_sorted_list = cumsum(r_sorted_list) / seq(1, n, 1)
        
        size_heu[k] = sum(cumu_r_sorted_list>0) / n
        
        abs_cumu_r_sorted_list = abs(cumu_r_sorted_list)
        
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
        
        NN = 10000  
        NNN = 1000  

        X = matrix(runif(NN*p, -1, 1), nrow=NN, ncol=p) 

        Y1 = exp(0.1 * X[, 1] + 0.2 * X[, 2] + (X[, 1]))

        Y0 = exp(0.1 * X[, 1] + 0.2 * X[, 2])

        CX_gen_E1 = function (x){
            epsilon = exp(rnorm(NNN, 0, 1))  # L = 8
            mean(apply(cbind(x * epsilon, L),1,min))
        }
        CX_gen_E2 = function (x){
            epsilon = exp(rlogis(NNN))
            mean(apply(cbind(x * epsilon, L),1,min))
        }
        CX_gen_E3 = function (x){
            epsilon = exp(log(-log(1-runif(NNN))))  # L = 8
            mean(apply(cbind(x * epsilon, L),1,min))
        } 
        if(noise_id == 1){
            CX_surv = apply(as.matrix(Y1), 1, CX_gen_E1) - apply(as.matrix(Y0), 1, CX_gen_E1)
        }
        if(noise_id == 2){
            CX_surv = apply(as.matrix(Y1), 1, CX_gen_E2) - apply(as.matrix(Y0), 1, CX_gen_E2)
        }
        if(noise_id == 3){
            CX_surv = apply(as.matrix(Y1), 1, CX_gen_E3) - apply(as.matrix(Y0), 1, CX_gen_E3)
        }

        est_rule = (predict(tree, X) == 2)
        size_mc[k] = sum(est_rule) / NN
        
        effect_mc[k] = sum(est_rule * CX_surv) / sum(est_rule) 

        cond_flag_mc[k] = effect_mc[k] >= delta
        
        pcd_mc[k] = mean(est_rule == (CX_surv>=0))
         
        
        cat("n=",n,"size_pct=",size_mc[k],"effect=",effect_mc[k],'satisfies?',pcd_mc[k],"\n")         
      
    }
    cat("Repe",ith,"\n") 
    c(size_mc,effect_mc,pcd_mc,cond_flag_mc,size_hat,size_heu) 
    
}
stopImplicitCluster()
stopCluster(cl) 

save(info,file="Res.RData")
