#power simulation based on Cochran-Armitage Trend Test
install.packages("DescTools")
armitagepowersimulation <- function(pd, pa, b, h, samp.size, analysis = 'unknown', 
                                    inflation = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, 
                                    pop.size = 1000000, alpha = 0.00000005, df = 2){
  pop <- creatediploidpopulation(pd, pa, b, h, analysis, pop.size)
  
  powerdata <- numeric(length(inflation))
  pvalues <- numeric(length(inflation))
  
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  for(t in inflation){
    armitage.stats <- numeric(n.sim)
    armitage.ps <- numeric(n.sim)
    case.inflation <-t
    
    for(i in 1:n.sim){
      num.casesample <- case.inflation*pd*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 0,]),
                    nrow(casesample[casesample$pop.genotypes == 1,]), 
                    nrow(casesample[casesample$pop.genotypes == 2,]) 
      )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 0,]),
                       nrow(controlsample[controlsample$pop.genotypes == 1,]),
                       nrow(controlsample[controlsample$pop.genotypes == 2,]) 
      )
      
      contingencytable <- matrix(c(caselist, controllist), ncol = 3, byrow = TRUE)
      #print(contingencytable)
      armitage.fit <- CochranArmitageTest(contingencytable, alternative = "one.sided")
      
      armitage.stats[i] <- armitage.fit$statistic
      armitage.ps[i] <- armitage.fit$p.value
    }
    #data[match(t,inflation)] <- mean(armitage.stats)
    powerdata[match(t,inflation)] <- mean(armitage.ps <= alpha)
    pvalues[match(t,inflation)] <- mean(armitage.ps)
  }
  #plots the data points for ncp
  plot(inflation*pd, powerdata, 
       xlab = "Case Fraction", ylab = "Power of Armitage trend test", las = 1, bty = "n", pch =19)
}