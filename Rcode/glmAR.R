# GLM Anchor regression hands on for 1 dimensional case (p=1)



#x = x values (i.e. call it x in data matrix when using "optim")



# create GLM objective function for likelihood optimization
GLM.objective <- function(b){
  
  # eta computation
  eta <- b*x
  
  # construction case specific objective
  if (family == "Gaussian") {
    return(sum((Y-eta)^2)) # obvious
  } 
  if (family == "Gaussian loglink") {
    return(sum((Y-exp(eta))^2)) # obvious
  } 
  if (family == "Poisson") {
    return(-sum(Y*eta-exp(eta))) # over using log-likelihood and avoid term without eta
  } 
  
} 

# create AR objective function for likelihood optimization
AR.objective <- function(b) {
  
  # eta computation
  eta <- b*x
  
  # construction case specific objective
  if (family == "Gaussian") {
    return(sum((Y-eta)^2)) # obvious
  } 
  if (family == "Gaussian loglink") {
    return(sum((Y-exp(eta))^2)) # obvious
  } 
  if (family == "Poisson") {
    return(-sum(Y*eta-exp(eta))) # over using log-likelihood and avoid term without eta
  } 
  
} 

# testdata:
last_14 = data.frame(rbind(
  c(3460,  14,    0),
  c(3558,  17,    1),
  c(3802,  21,    2),
  c(3988,  22,    3),
  c(4262,  28,    4),
  c(4615,  36,    5),
  c(4720,  40,    6),
  c(5404,  47,    7),
  c(5819,  54,    8),
  c(6440,  63,    9),
  c(7126,  85,   10),
  c(7905, 108,   11),
  c(8733, 118,   12),
  c(9867, 200,   13)))
names(last_14) = c('World', 'US', 'days')

x <- last_14$days
Y <- last_14$US

# Gaussian
glm(Y~x-1, family=gaussian())
family = "Gaussian"
optimize(GLM.objective, interval = c(-20,20))

# Poisson
glm(Y~x-1, family=poisson())
family = "Poisson"
optimize(GLM.objective, interval = c(-20,20))


