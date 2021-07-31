#' Simulates the scenarios of Section 3
#'
#'
#' @param scenario an integer that represents the scenario we explore in the paper, eg. scenario=1, simulates the scenario of Section 3.4.1, scenario=2, simulates the scenario of Section 3.4.2 etc
#' @param model an integer that represents the model from which we simulate the data, eg model=1, model=2 etc 
#' @param predictive_amplitude a number that captures the strength of the predictive signal, e.g. the parameter theta_pred of the paper
#' @param sample_size the number of exmamples (e.g. patients)
#' @param num_features the number of variables (features, covarietes)
#' 
#' @return the simulated data 
generate_scenarios_predictive = function(scenario, model, sample_size, num_features, predictive_amplitude) {
  
  
  if(scenario == 'S1')
  {
    # Sample T:
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y
    if(model=='M1')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(11,20,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(1,10,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M3')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,10,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M4')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(1,5,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M5')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = 0
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  
  
  if(scenario == 'S2')
  {
    # Sample t
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y 
    if(model=='M1')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = (apply(exp(matrix(runif(sample_size*length(nonzero_prog),0.50,1.50), sample_size, length(nonzero_prog))*X[,nonzero_prog]),1,sum))
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude !=0){
        pred_part = predictive_amplitude*t*(apply(exp(kronecker(t(runif(length(nonzero_pred),0.50,1.50)), matrix(1, sample_size, 1))*X[,nonzero_pred]),1,sum))
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] * X[,nonzero_prog[2]] + X[,nonzero_prog[3]]*X[,nonzero_prog[4]] + X[,nonzero_prog[5]]*X[,nonzero_prog[6]]+ X[,nonzero_prog[7]]*X[,nonzero_prog[8]] + X[,nonzero_prog[9]]*X[,nonzero_prog[10]] 
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0 ){
        pred_part = X[,nonzero_pred[1]] * X[,nonzero_pred[2]] + X[,nonzero_pred[3]]*X[,nonzero_pred[4]] + X[,nonzero_pred[5]]*X[,nonzero_pred[6]]+ X[,nonzero_pred[7]]*X[,nonzero_pred[8]] + X[,nonzero_pred[9]]*X[,nonzero_pred[10]]
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M3')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] + X[,nonzero_prog[2]]^2 + X[,nonzero_prog[3]]^3 + X[,nonzero_prog[4]]*X[,nonzero_prog[5]] + X[,nonzero_prog[6]]*X[,nonzero_prog[7]] + X[,nonzero_prog[8]]*X[,nonzero_prog[9]]*X[,nonzero_prog[10]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0){
        pred_part = X[,nonzero_pred[1]] + X[,nonzero_pred[2]]^2 + X[,nonzero_pred[3]]^3 + X[,nonzero_pred[4]]*X[,nonzero_pred[5]] + X[,nonzero_pred[6]]*X[,nonzero_pred[7]] + X[,nonzero_pred[8]]*X[,nonzero_pred[9]]*X[,nonzero_pred[10]]
        pred_part = predictive_amplitude*t*pred_part
        mu =  0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  
  
  if(scenario == 'S3')
  {
    # Sample t
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y
    if(model=='M1')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,5,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0){
        pred_part = (X[,nonzero_pred[1]]>0) + (X[,nonzero_pred[2]]<0) + (X[,nonzero_pred[3]]<0.675) + (X[,nonzero_pred[4]]>0.675) + (X[,nonzero_pred[5]]<0.675)
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,5,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      if(predictive_amplitude != 0){
        pred_part = abs(X[,nonzero_pred[1]])* (X[,nonzero_pred[2]]>0.675) + abs(X[,nonzero_pred[3]])*(X[,nonzero_pred[4]]<0)*(X[,nonzero_pred[5]]>0)
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  
  if(scenario == 'S4')
  {
    
    if(model=='M1')
    { 

      # Sample t
      prob_t = 0.50 
      t = rbinom(sample_size,1,prob_t)
      # Sample Xs
      rho = 0.50
      sigma_z =  toeplitz(rho^(0:(num_features - 1)))
      X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
      # Generate y
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(11,20,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M2')
    {     

      # Sample t
      prob_t = 0.50 
      t = rbinom(sample_size,1,prob_t)
      # Sample Xs
      rho = 0.50
      sigma_z =  toeplitz(rho^(0:(num_features - 1)))
      X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
      # Generate y
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] * X[,nonzero_prog[2]] + X[,nonzero_prog[3]]*X[,nonzero_prog[4]] + X[,nonzero_prog[5]]*X[,nonzero_prog[6]]+ X[,nonzero_prog[7]]*X[,nonzero_prog[8]] + X[,nonzero_prog[9]]*X[,nonzero_prog[10]] 
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      if(predictive_amplitude != 0 ){
        pred_part = X[,nonzero_pred[1]] * X[,nonzero_pred[2]] + X[,nonzero_pred[3]]*X[,nonzero_pred[4]] + X[,nonzero_pred[5]]*X[,nonzero_pred[6]]+ X[,nonzero_pred[7]]*X[,nonzero_pred[8]] + X[,nonzero_pred[9]]*X[,nonzero_pred[10]] 
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  output = c()
  output$t = t
  output$X = X
  output$y = y[,1]
  
  output$predictive = nonzero_pred

  
  return(output)
}
