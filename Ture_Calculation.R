# Calculate the true for Section 5

# For Section 5.1 Evaluation and Comparison with Average Treatment Effect
# and Section 5.2 Evaluation of Multiple Constraints 

# Scenario 1

p = 10 # dimension of X

n = 500000 

x = matrix(runif(n*p, -2, 2), nrow=n, ncol=p)
# Kernel Density Plot
CX =  x[, 1] 
CX_den <- density(CX) # returns the density data 

plot(CX_den, main='') 

sum((CX > 0) * CX) / sum(CX > 0) # \delta = 1, expected selected sample ratio : 50% S1
mean(CX > 0)  

# Scenario 2 
# Kernel Density Plot
CX = (x[, 1] * x[, 2])
CX_den <- density(CX) # returns the density data

plot(CX_den, main='') 

sum((CX > 0) * CX) / sum(CX > 0) # \delta = 1.0 , expected selected sample ratio : 50% S2 
mean(CX > 0)  


# Scenario 3  
# Kernel Density Plot
CX = (x[, 1] + x[, 2])
CX_den <- density(CX) # returns the density data 

plot(CX_den, main='') 

sum((CX > -0.55) * CX) / sum(CX > -0.55) # \delta = 1.0, expected selected sample ratio : 63% S2 
mean(CX > -0.55)  
 

# For Section 5.3 Evaluation of Survival Data

# Scenario 4


set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 7

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
    epsilon = exp(rnorm(NN, 0, 1))   # L = 12
#     epsilon = exp(rlogis(NN))        # L = 21
#     epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv1 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)


set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 12

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
    epsilon = exp(rnorm(NN, 0, 1))   # L = 12
#     epsilon = exp(rlogis(NN))        # L = 21
#     epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv2 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)
         


set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 10

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
#     epsilon = exp(rnorm(NN, 0, 1))   # L = 12
    epsilon = exp(rlogis(NN))        # L = 21
#     epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv3 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)
         

set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 21

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
#     epsilon = exp(rnorm(NN, 0, 1))   # L = 12
    epsilon = exp(rlogis(NN))        # L = 21
#     epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv4 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)
         

set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 4

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
#     epsilon = exp(rnorm(NN, 0, 1))   # L = 12
#     epsilon = exp(rlogis(NN))        # L = 21
    epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv5 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)
         

set.seed(2333)
# calculate the true:
p = 10 # dimension of X

n = 10000  
NN = 1000  

x = matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

L = 8

Y1 = exp(0.1 * x[, 1] + 0.2 * x[, 2] + (x[, 1]))

Y0 = exp(0.1 * x[, 1] + 0.2 * x[, 2])

CX_gen = function (x){
    
#     epsilon = exp(rnorm(NN, 0, 1))   # L = 12
#     epsilon = exp(rlogis(NN))        # L = 21
    epsilon = exp(log(-log(1-runif(NN))))  # L = 8
    mean(apply(cbind(x * epsilon, L),1,min))
    
}
# runif(n, 0, 7) 25% # runif(n, 0, 12) #15% #normal
# runif(n, 0, 10) 25% # runif(n, 0, 21) #15% #logsitic
# runif(n, 0, 4) 25% # runif(n, 0, 8) #15% #extreme
 
CX_surv6 = apply(as.matrix(Y1), 1, CX_gen) - apply(as.matrix(Y0), 1, CX_gen)
         


CX_den <- density(CX_surv1) # returns the density data
CX_den2 <- density(CX_surv2) # returns the density data
CX_den3 <- density(CX_surv3) # returns the density data
CX_den4 <- density(CX_surv4) # returns the density data
CX_den5 <- density(CX_surv5) # returns the density data
CX_den6 <- density(CX_surv6) # returns the density data
pdf("S2_den.pdf") 

plot(CX_den, main='', xlab='Restricted Mean Survival Time', col = 'blue', lty=1, xlim=c(-2,3), ylim = c(0,1.1))
lines(CX_den2, col = 'blue', lty=2)
lines(CX_den3, col = 'green', lty=1)
lines(CX_den4, col = 'green', lty=2)
lines(CX_den5, col = 'red', lty=1)
lines(CX_den6, col = 'red', lty=2)

legend("topright", legend=c('Normal with Censoring 25%','Normal with Censoring 15%', 
                            'Log with Censoring 25%','Log with Censoring 15%',
                           'Extreme with Censoring 25%','Extreme with Censoring 15%'), 
       lty = c(1,2,1,2,1,2), col = c("blue","blue",'green','green', "red","red"))

dev.off()