#power simulation based on logistic regression
logisticpowersimulation <- function(pd, pa, b, h, samp.size, analysis = 'unknown', inflation = 
                                      c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- creatediploidpopulation(pd, pa, b, h, analysis, pop.size)
  
  #inflation vector contains the different inflation factors that will be used
  data <- numeric(length(inflation))
  upperncp <- numeric(length(inflation))
  lowerncp <- numeric(length(inflation))
  
  #subsets the population into cases and controls based on phenotype
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  #loops through the different inflation factors to run chi-squared test
  for(t in inflation){
    logistic.ps <- numeric(n.sim)
    case.inflation <-t
    
    #repeats for the desired number of simulations
    for(i in 1:n.sim){
      num.casesample <- case.inflation*pd*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      sample <- rbind(casesample, controlsample)
      logistic <- glm(sample$pop.phenotypes ~ sample$pop.genotypes,family = 'binomial'( link = 'logit'))
      logistic.ps[i] <- summary(logistic)$coefficients[2,4]
      #stores chi-squared statistic and p-value
      
    }
    data[match(t,inflation)] <- mean(logistic.ps <= alpha)
    print(mean(logistic.ps))
  }
  #plots the data points for ncp
  plot(inflation*pd, data, 
       xlab = "Case Fraction", ylab = "power", las = 1, bty = "n", pch =19)
}