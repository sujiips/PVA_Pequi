install.packages("hierfstat")  
require(hierfstat)             
 
Dir <- "/Users/sujii/Model_pequi_2"  
setwd(Dir)    

source('functions_pequi.R')
parms = list(
  Start <- 1
  , rep  <- 50
  , nloci <- 10
  , time <- 500
  , xDim <- 500
  , yDim <- 400
  , dPollen <- 1500    
  , dSeed <- 5     
  , avSeed <- c(10)
  , maxAge <- c(200)
  , adultAge <- 8
  , maxFathers <- 10
  , germ <- c(0.2)     
  , selection <- c(0.68)  
  , output_year <- c(1, 50, 75, 100, 200,250, 300, 350, 400,450, 500)
)

begin(parms)



