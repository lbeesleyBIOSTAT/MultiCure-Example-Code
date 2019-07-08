

####################
### Load in Code ###
####################

library(devtools)
install_github("lbeesleyBIOSTAT/MultiCure")
library(MultiCure)

#Set working directory here

######################
### View Vignettes ###
######################

browseVignettes('MultiCure')

#####################
### Simulate Data ###
#####################

### The following code simulates datasets under a multistate cure model with 1) no covariate missingness or unequal censoring, 2) covariate missingness, and 3) unequal censoring
NONE = SimulateMultiCure(type = 'NoMissingness') 
COV = SimulateMultiCure(type = 'CovariateMissingness')
CENS = SimulateMultiCure(type = 'UnequalCensoring')

######################
### Visualize Data ###
######################

### The following is a RShiny function for visualizing the observed information for recurrence and death. We note that there is not unequal censoring in this simulated dataset.Note: We have a very small probability of experiencing a recurrence after about time 50 (although they may rarely happen), so it seems generally reasonable to assume subjects still at risk for recurrence and death after time 50 are cured

datWIDE = data.frame( Y_R = NONE$Y_R, Y_D = NONE$Y_D, delta_R = NONE$delta_R, delta_D = NONE$delta_D)
VISUALIZEDATA(datWIDE)


### Data with unequal censoring
datWIDE = data.frame( Y_R = CENS$Y_R, Y_D = CENS$Y_D, delta_R = CENS$delta_R, delta_D = CENS$delta_D)
VISUALIZEDATA(datWIDE)





####################################################
####################################################
### Fitting Multistate Cure Model to "Complete" Data 
### with Weibull Baseline Hazards
####################################################
####################################################

### Setting up 
Cov = data.frame(X1 = NONE$X1,X2 = NONE$X2)
VARS = names(Cov)
TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
datWIDE = data.frame( Y_R = NONE$Y_R, Y_D = NONE$Y_D, delta_R = NONE$delta_R , delta_D = NONE$delta_D, G = NONE$G)

### Parameter Estimation
fit = MultiCure(iternum = 50, datWIDE  = datWIDE, Cov = Cov, ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib')

beta = fit[[1]]
alpha = fit[[2]]
scale = fit[[3]]
shape = fit[[4]]
beta_save = fit[[5]]
alpha_save = fit[[6]]
scale_save = fit[[7]]
shape_save = fit[[8]]

### Variance Estimation
OUT = VarianceEM(fit,iternum=20, bootnum=50, datWIDE, Cov, ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib')







####################################################################
####################################################################
### Fitting Multistate Cure Model to Data with Covariate Missingness 
### with Weibull Baseline Hazards
####################################################################
####################################################################

Cov = data.frame(X1 = COV$X1, X2 = COV$X2)
VARS = names(Cov)
TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
datWIDE = data.frame( Y_R = COV$Y_R, Y_D = COV$Y_D, delta_R = COV$delta_R , delta_D = COV$delta_D, G = COV$G)

### Parameter Estimation (This may take a moment)
ITERNUM = 200
fit = MultiCure(iternum = ITERNUM, datWIDE  = datWIDE, Cov = Cov, COVIMPUTEFUNCTION  = COVIMPUTEFUNCTION_Example, COVIMPUTEINITIALIZE = COVIMPUTEINITIALIZE_Example, IMPNUM = 10,
	 ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib')
beta_save = fit[[5]]
alpha_save = fit[[6]]
scale_save = fit[[7]]
shape_save = fit[[8]]

### Theta-Hat can be estimated as the mean of the last few iterations
beta = apply(beta_save[,(ITERNUM-10): ITERNUM], 1, mean)
alpha =  apply(alpha_save[,(ITERNUM-10): ITERNUM], 1, mean)
scale =  apply(scale_save[,(ITERNUM-10): ITERNUM], 1, mean)
shape =  apply(shape_save[,(ITERNUM-10): ITERNUM], 1, mean)

### Variance Estimation
OUT = VarianceMCEM(fit,var_method = 'default', datWIDE = datWIDE,  ASSUME = 'SameHazard', 
          TransCov = TransCov, BASELINE = 'weib', COVIMPUTEFUNCTION = COVIMPUTEFUNCTION_Example, 
          COVIMPUTEINITIALIZE = COVIMPUTEINITIALIZE_Example, POSTITER = 5)	



################################################################
################################################################
### Fitting Multistate Cure Model to Data with Unequal Censoring 
### with Weibull Baseline Hazards
################################################################
################################################################

Cov = data.frame(X1 = CENS$X1, X2 = CENS$X2)
VARS = names(Cov)
TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
datWIDE = data.frame( Y_R = CENS$Y_R, Y_D = CENS$Y_D, delta_R = CENS$delta_R, delta_D = CENS$delta_D, G = CENS$G)

### Parameter Estimation (This may take a moment)
ITERNUM = 200
fit = MultiCure(iternum = ITERNUM, datWIDE  = datWIDE, Cov = Cov, IMPNUM = 10,
	 ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib', 
	 UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION)

beta_save = fit[[5]]
alpha_save = fit[[6]]
scale_save = fit[[7]]
shape_save = fit[[8]]
### Theta-Hat can be estimating as the mean of the last few iterations
beta = apply(beta_save[,(ITERNUM-10): ITERNUM], 1, mean)
alpha =  apply(alpha_save[,(ITERNUM-10): ITERNUM], 1, mean)
scale =  apply(scale_save[,(ITERNUM-10): ITERNUM], 1, mean)
shape =  apply(shape_save[,(ITERNUM-10): ITERNUM], 1, mean)

### Variance Estimation
OUT = VarianceMCEM(fit,var_method = 'default', datWIDE = datWIDE,  ASSUME = 'SameHazard', 
          TransCov = TransCov, BASELINE = 'weib', UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION, 
          POSTITER = 5)


####################################################
####################################################
### Fitting Multistate Cure Model to "Complete" Data 
### with Cox Baseline Hazards
####################################################
####################################################

Cov = data.frame(X1 = NONE$X1,X2 = NONE$X2)
VARS = names(Cov)
TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
datWIDE = data.frame( Y_R = NONE$Y_R, Y_D = NONE$Y_D, delta_R = NONE$delta_R , delta_D = NONE$delta_D, G = NONE$G)

### Parameter Estimation
fit = MultiCure(iternum = 50, datWIDE  = datWIDE, Cov = Cov, ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'cox')
beta = fit[[1]]
alpha = fit[[2]]
beta_save = fit[[3]]
alpha_save = fit[[4]]
p_save = fit[[5]]

### Variance Estimation
OUT = VarianceEM(fit,iternum=20, bootnum=50, datWIDE, Cov, ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'cox')




###############################################
###############################################
### Examine Convergence of Parameter Values ###
###############################################
###############################################

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

par(mfrow=c(1,3))

### Beta for X1
matplot(t(beta_save[c(1,3,5,7),]), type = 'l',lty = 1, lwd = 2, col =  cbPalette[c(6,7,4,8)], main = 'Beta Parameters for X1', xlab = 'Iterations', ylab = 'Estimated Value')
legend(x='topright', fill = cbPalette[c(6,7,4,8)], legend = c( 'Transition 1->3', 'Transition 2->4', 'Transition 1->4', 'Transition 3->4'), cex = 0.8)
#Note: Some parameters may be equal by assumption

### Beta for X2
matplot(t(beta_save[c(2,4,6,8),]), type = 'l',lty = 1, lwd = 2, col =  cbPalette[c(6,7,4,8)], main = 'Beta Parameters for X2', xlab = 'Iterations', ylab = 'Estimated Value')
legend(x='bottomright', fill = cbPalette[c(6,7,4,8)], legend = c( 'Transition 1->3', 'Transition 2->4', 'Transition 1->4', 'Transition 3->4'), cex = 0.8)

### Alpha
matplot(t(alpha_save), type = 'l',lty = 1, lwd = 2, col =  cbPalette[c(6,7,4)], main = 'Alpha Parameters', xlab = 'Iterations', ylab = 'Estimated Value')
legend(x='topright', fill = cbPalette[c(6,7,4,8)], legend = c( 'Intercept', 'X1', 'X2'), cex = 0.8)




###############################################
###############################################
### Predicted State Occupancy Probabilities ### 
###############################################
###############################################


### These are RShiny functions that will plot the predicted state occupancy probabilities based on the multistate cure model fits. The probabilities can be estimated for a single set of covariates or for several sets of covariates simultaneously. The plots can be downloaded to the working directory


### Weibull Baseline Hazards
STATEOCCUPANCYWEIB(times = seq(0,max(datWIDE$Y_D),1), TransCov, newCov = data.frame(X1 = c(0,0.5), X2 = c(0,0.5)), 
	beta, alpha, scale, shape) 


### Cox Baseline Hazards, EM Algorithm
Haz = BaselineHazard_NOIMP(datWIDE, Cov, beta, alpha, TransCov, ASSUME = 'SameHazard', p=p_save[,length(p_save[1,])])	
STATEOCCUPANCYCOX_NOIMP(times = seq(0,max(datWIDE$Y_D),1), TransCov, newCov = data.frame(X1 = c(0,0.5), X2 = c(0,0.5)), 
	beta, alpha, Haz_13 = Haz[[1]], Haz_24 = Haz[[2]], Haz_14 = Haz[[3]], Haz_34 = Haz[[4]]) 





