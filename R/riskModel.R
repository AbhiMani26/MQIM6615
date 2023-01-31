library(dplyr)
library(fpp2)
library(tidyquant)
library(tidyverse)
library(quantmod)
library(dplyr)
library(readxl)
library(rugarch)
library(fGarch)
library(BurStFin)
library(matrixcalc)

calculateRiskModel <-function(){
  #loads data
  factor_return <- read.table("factors_2022.csv", header=TRUE,sep=",",as.is=TRUE)
  stock_returns <- read.table("universe_2022.csv", header=TRUE,sep=",",as.is=TRUE,)
  factor_return <- head(factor_return,250)
  stock_returns <- head(stock_returns,250)


  #removes date from the data of stock returns and factor returns
  factor_return <- factor_return[-1]
  stock_returns <- stock_returns[-1]
  stock_matrix <- do.call("cbind",stock_returns)
  factor_matrix <- do.call("cbind",factor_return)

  #One stock has missing data, fixes the issue
  x = which(colSums(is.na(stock_returns))>0)
  my_list <- c()
  k=0
  num_rows <- nrow(stock_returns)
  for(i in 1:num_rows){
    if(is.na(stock_returns[i,x])){
      k=i
      break
    }
    my_list <- append(my_list,stock_returns[i,x])
  }

  mean <- mean(my_list)
  sd <- sd(my_list)

  for(j in k+1:num_rows){
    stock_returns[j,x] <- rnorm(1,mean,sd)
    if(j == num_rows){
      break
    }
  }

  # Intializes matrix to hold stock betas and stock specific risk variance
  num_factors <- ncol(factor_return)
  num_stocks <- ncol(stock_returns)
  beta_matrix  <- matrix(0,num_stocks, num_factors)
  stockrisk_matrix <- matrix(0, num_stocks, num_stocks)

  #defines a function to calculate stock specific risk from the residuals of regression
  idiosyncratic_risk <- function(residuals){
    garchmodel1 <- garchFit( ~ garch(1,1), data=coredata(residuals), trace=FALSE)
    return(volatility(garchmodel1))
  }

  #Runs multivariate regression between stock returns and various factors
  residual_vector <- c()
  for(j in 1:ncol(stock_matrix)){
    result <- lm(stock_matrix[,j]~factor_matrix[,1] + factor_matrix[,2] + factor_matrix[,3]+ factor_matrix[,4] + factor_matrix[,5]+ factor_matrix[,6] + factor_matrix[,7]+ factor_matrix[,9] + factor_matrix[,10]+ factor_matrix[,11] + factor_matrix[,12]+ factor_matrix[,13] + factor_matrix[,14] + factor_matrix[,15])
    result_coff <- summary(result)
    for(i in 1:15){
      beta_matrix[j,i] = result_coff$coefficients[i+1]
    }
    residual_vector <- c(residual_vector,idiosyncratic_risk(result_coff$residuals))
  }

  #creates stock specific risk matrix
  for(i in 1:num_stocks){
    stockrisk_matrix[i,i] = residual_vector[i]
  }

  #creates the final risk model
  factor_covariance_matrix <- var.shrink.eqcor(factor_return)
  risk_model <- (beta_matrix %*% factor_covariance_matrix %*% t(beta_matrix))
  risk_model <- risk_model + stockrisk_matrix
  risk_model <- round(risk_model,8)
  if(is.positive.semi.definite(risk_model, tol=1e-8)){
    print("Risk Model is Positive Semi definite")
  }
  #since all eigen values are positive, hence we can say that matrix is semi positive definite
  return(risk_model)
}
