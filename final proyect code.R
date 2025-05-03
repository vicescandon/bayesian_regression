#### FINAL PROYECT
rmvnorm<-function(n,mu,Sigma) 
{ # samples from the multivariate normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(t(E%*%chol(Sigma))+c(mu))
}

library(readr)
setwd("C:/Users/victo/OneDrive/Documents/UNIFI/Classes/Bayesian Data Analysis/final proyect/DATA BASE")
base_final <- read_csv("base_final.csv")

#### 10%
library(tidyverse)
library(dplyr)
library(MASS)
library(mvtnorm)
library(patchwork)

data<-base_final

#### DATA CLEANING
#eliminate obs from 2019-2021 and NA in variables necessary
data2 <- data %>%
  filter(!is.na(treat_prop) &
           !is.na(nteacher) &
           !is.na(abandono) &
           !is.na(nteacher2) &
           !is.na(teachratio) &
           !is.na(prdoclice) &
           !is.na(csize) &
           marginal == "Muy Alto" &
           year <= 2018)


#abandono en el tiempo

# Convert the binary variable to a factor for better labeling
data2$treatment <- factor(data2$treat_prop, levels = c(0, 1), labels = c("Control", "Treated"), ordered = TRUE)

# Create the plot
ggplot(data2, aes(x = year, y = abandono, color = treatment)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_x_continuous(breaks = unique(data$year)) +
  labs(title = "Average Outcome Over Years by Treatment Group",
       x = "Year",
       y = "Average Outcome",
       color = "Group") +
  theme_minimal()

#another useful plot
boxplot(abandono ~ treat_prop,data=data2, xlab = "Treatment", ylab = "Outcome", outline= FALSE)

#reg abandono en funcion de tratamiento (sin controles)
reg0<-lm(data=data2,formula = abandono~treat_prop)
summary(reg0)

#reg with controls
form<-data2$abandono~data2$treat_prop+
  data2$nteacher+
  data2$nteacher2+
  data2$teachratio+
  data2$prdoclice+
  data2$csize


#reg abandono en funcion de tratamiento (con controles)
reg1<-lm(data=data2,formula = form)
summary(reg1)

#abandono frequency
hist(x=data2$abandono, xlab = "Dropout", main = "Outcome distribution")
hist(x=data2$treat_prop)

#priors and models

#PRIOR CHECKING :::::::::::::::::::::::::::::::::
#data to use
datab <-data2 %>% 
  dplyr::select(abandono, treat_prop, nteacher, nteacher2, teachratio, prdoclice, csize)

X <- as.matrix(datab %>% dplyr::select(-abandono))
X<-as.data.frame(X)

#intercept
X<- X %>%
  add_column(intercept=rep(1,nrow(X)),.before = "treat_prop")
X<-as.matrix(X)

XtX<-t(X)%*%X



#### prior checking ####
y <- datab$abandono

n_c<-length(y)

s_c<-sd(y)

m_c<-mean(y)

N<-200

#prior sets
beta <-reg1$coefficients
sigma0<-(summary(reg1)$sigma)^2
Sigma0<- solve(t(X)%*%X)*sigma0*n_c
nu<-1

####nreplicate data from posterior####

#prepare vectors
s2_rep<-rep(0,N)
beta_rep<-matrix(0,nrow=N, ncol=length(beta))
yrep <- vector("list", N)

#loop for generation of r.v.
for (j in 1:N) {
  s2_rep[j]=1/rgamma(1,nu/2,(nu*sigma0)/2) #(200x1)
  beta_rep[j,]=mvrnorm(n = 1, mu = beta, Sigma= Sigma0)#(200x1)
  
  yrep[[j]]=rnorm(n=n_c,mean = X%*%beta_rep[j,],sd=sqrt(s2_rep[j]))  #(6x1) NOTE IT IS JUST A VECTOR
  #this is ONE yrep vector for each S draws of theta
}

#plot densities
# Convert yrep list into a data frame with a replication identifier
yrep_df <- bind_rows(lapply(seq_along(yrep), function(i) {
  data.frame(Predicted = yrep[[i]], Replicate = as.factor(i))
}))

# Plotting each replication's density
prior_densities<-ggplot(yrep_df, aes(x = Predicted, group = Replicate, color = Replicate)) +
  geom_density(size = 0.3, show.legend = FALSE) +
  labs(title = "Density of Each Replication of Prior Predictive Simulations",
       x = "Predicted Values", 
       y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ # Improve x-axis label readability
  xlim(min=-150, max=150) # Limit x-axis to 80% of max to compact the view
  

#means
mean_rep <- sapply(yrep, mean, na.rm = TRUE)
mean_rep_list <- lapply(yrep, mean, na.rm = TRUE)
mean_df <- data.frame(Replication = 1:length(mean_rep), Mean = mean_rep)
prior_means<-ggplot(mean_df, aes(x = Replication, y = Mean)) +
  geom_line() +  # Line connecting the means
  geom_point() +  # Points for each mean
  labs(title = "Means of Each Replication",
       x = "Replication Number",
       y = "Mean") +
  theme_minimal()

#sd
sd_rep_list <- lapply(yrep, sd, na.rm = TRUE)
sd_rep <- unlist(sd_rep_list)
sd_df <- data.frame(Replication = 1:length(sd_rep), StdDev = sd_rep)
prior_sd<-ggplot(sd_df, aes(x = Replication, y = StdDev)) +
  geom_line() +
  geom_point() +
  labs(title = "Standard Deviations of Each Replication",
       x = "Replication Number",
       y = "Standard Deviation") +
  theme_minimal()

#plots for prior predictive check
prior_densities + prior_means+ prior_sd +
  plot_layout(ncol = 1)



#POSTERIOR ESTIMATION USING UNIT INFORMATION PRIOR:::::::::::::::::::
n<-length(y)
S<-5000

# unit information prior setting
beta.0 <-reg1$coefficients #SPECOFY PARAMETERS, THIS ONE AS TE FREQ ESTIMATE
p<-length(beta.0)

nu.0<-1  ; sigma2.0<-sum(reg1$res^2)/(n-p) #FREQ EST AGAIN
Sigma.0<- solve(t(X)%*%X)*sigma2.0*n

### some convenient quantities
iSigma.0<-solve(Sigma.0) #INVERSE OF SIGMA ZERO
XtX<-t(X)%*%X #X TRANSPOSE X
# Sigma.0 is the covariance matrix of the unit information prior

### store mcmc samples in these objects
beta.post<-matrix(nrow=S,ncol=p) #S ROW ARE THE S SIMULATIONS
sigma2.post<-rep(NA,S)

### starting value
sigma2<- var(residuals(reg1))


### MCMC algorithm
for( scan in 1:S) {
  
  #update beta
  V.beta<- solve(  iSigma.0 + XtX/sigma2 )
  E.beta<- V.beta%*%( iSigma.0%*%beta.0 + t(X)%*%y/sigma2 )
  beta<-t(rmvnorm(1, E.beta,V.beta) )
  
  #update sigma2
  nu.n<- nu.0+n
  ss.n<-nu.0*sigma2.0 + sum(  (y-X%*%beta)^2 )
  sigma2<-1/rgamma(1,nu.n/2, ss.n/2)
  
  #save results of this scan
  beta.post[scan,]<-beta
  sigma2.post[scan]<-sigma2
}

#correlation in mcmc
par(mfrow=c(1,1))
plot(beta.post)
#HIST POSTERIORS
par(mfrow=c(3,3))
hist(beta.post[,1],main = "Intercept", xlab = "Posterior support")
abline(v=mean(beta.post[,1]),col = "blue")
abline(v=reg1$coefficients[1],col="red")
hist(beta.post[,2], main = "Treatment", xlab = "Posterior support")
abline(v=mean(beta.post[,2]),col = "blue")
abline(v=reg1$coefficients[2],col="red")
hist(beta.post[,3], main = "nteacher", xlab = "Posterior support")
abline(v=mean(beta.post[,3]),col = "blue")
abline(v=reg1$coefficients[3],col="red")
hist(beta.post[,4], main = "nteacher square", xlab = "Posterior support")
abline(v=mean(beta.post[,4]),col = "blue")
abline(v=reg1$coefficients[4],col="red")
hist(beta.post[,5], main = "teacher ratio", xlab = "Posterior support")
abline(v=mean(beta.post[,5]),col = "blue")
abline(v=reg1$coefficients[5],col="red")
hist(beta.post[,6], main = "teacherbcs", xlab = "Posterior support")
abline(v=mean(beta.post[,6]),col = "blue")
abline(v=reg1$coefficients[6],col="red")
hist(beta.post[,7], main = "csize", xlab = "Posterior support" )
abline(v=mean(beta.post[,7]),col = "blue")
abline(v=reg1$coefficients[7],col="red")
hist(sigma2.post, main = "sigmasq", xlab = "Posterior support")
abline(v=mean(sigma2.post),col="blue")
abline(v=sigma2.0,col="red")

# credible intervals
# Number of betas (assuming 7 in your case)
num_betas <- ncol(beta.post)

# Initialize matrices to store credible intervals
beta_ci <- matrix(NA, nrow = num_betas, ncol = 2)
sigma2_ci <- numeric(2)

# Calculate 95% credible intervals for each beta
for (i in 1:num_betas) {
  beta_ci[i, ] <- quantile(beta.post[, i], probs = c(0.025, 0.975))
}

# Calculate 95% credible interval for sigma^2
sigma2_ci <- quantile(sigma2.post, probs = c(0.025, 0.975))

# Print results
colnames(beta_ci) <- c("2.5%", "97.5%")
rownames(beta_ci) <- paste0("beta_", 1:num_betas)
print("95% Credible Intervals for Betas:")
print(beta_ci)

print("95% Credible Interval for Sigma^2:")
print(sigma2_ci)


###### POSTERIOR PREDICTIVE CHECK
#store in 
N<-1500
yrep <- vector("list", N)
betas<-matrix(NA,nrow = N, ncol = 7)
s2<-rep(NA, N)

Ty_mean<-mean(datab$abandono)

# generate yrep
set.seed(1)
for (j in 1:N) {
  betas[j,]= beta.post[sample(nrow(beta.post), size =  1), ] 
  s2[j]=sample(sigma2.post, size = 1)    
  
  yrep[[j]]=rnorm(n=length(y),mean = X%*%betas[j,],sd=sqrt(s2[j]))  
  #this is ONE yrep vector for each S draws of theta
}


# T STAT
#means
mean_yrep <- lapply(yrep, mean, na.rm = TRUE)
mean_yrep<-as.numeric(unlist(mean_yrep))
#tstat
sum(mean_yrep>mean(datab$abandono))/N

#max's
max_yrep<-lapply(yrep,max)
max_yrep<-abs(as.numeric(unlist(max_yrep)))
sum(abs(max_yrep)>abs(max(datab$abandono)))


#plot means of generated y distributions
#means
par(mfrow=c(1,3))
hist(x=data2$abandono, xlab = "Dropout", main = "Outcome distribution")
hist(mean_yrep,main="yrep mean distribution")
abline(v=mean(mean_yrep),col="blue")
abline(v=mean(datab$abandono),col="red")
#max value
hist(max_yrep,xlim=c(0,110), main = "yrep max distribution")
abline(v=abs(max(datab$abandono)),col="red")

sum(datab$abandono>70)

