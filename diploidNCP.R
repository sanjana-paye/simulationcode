creatediploidpopulation <- function(pd, pa, b, h, analysis, pop.size = 1000000){
  #creates vector of genotypes
  pop.genotypes <- rbinom(pop.size, 2, pa)
  #creates empty vector for phenotypes
  pop.phenotypes <- numeric(length(pop.genotypes))
  
  if (analysis == 'dominant'){
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd))/((-1+pa)*(-1-pa+2*h*pa))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd)/((-1+pa)*(-1-pa+2*h*pa))))  
    pop.genotypes[pop.genotypes == 2] <- 1
    
  }else if (analysis == 'recessive'){
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd))/((-1+pa)*(-1-pa+2*h*pa))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd)/((-1+pa)*(-1-pa+2*h*pa))))  
    pop.genotypes[pop.genotypes == 1] <- 0
    pop.genotypes[pop.genotypes== 2] <-1
    
  }else{
    #assigns phenotypes based on the three genotypes using probabilities from chi-squared tables
    pop.phenotypes[pop.genotypes == 2] <- rbinom(sum(pop.genotypes == 2), 1, b)
    pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, 
                                                 (b*h+((1-h)*(-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd))/((-1+pa)*(-1-pa+2*h*pa))))
    pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, 
                                                 ((-2*b*h*pa-b*pa^2+2*b*h*pa^2+pd)/((-1+pa)*(-1-pa+2*h*pa))))  
  }
  
  #creates and returns a data frame with the population genotypes and phenotypes
  pop <- data.frame(pop.genotypes, pop.phenotypes)
  return(pop)
}

diploidcalculation <- function(pd, pa, b, h, samp.size, analysis, x, df = 2){
  inflation <- seq(0.001, x, length.out = 1000)
  
  data <- vector()
  
  for (t in inflation){
    
    #calculate non-centrality parameter
    if (analysis == 'dominant'){
      ncp <- (t*(pa-1)*(pa*(-2*h*(pa-1)+pa)^2)*((b-pd)^2)*(t*pd-1))/(pd*(pd-1+(2*h*(pa-1)-pa)*pa*(b-1-b*t+t*pd))*(-2-b*(t-1)*(2*h*(pa-1)-pa)*(pa-1)-pa+2*pd+pa*(pa+t*pd-t*pa*pd)+2*h*((t-1)*pd+(pa-2)*pa*(-1+t*pd))))
      
    }else if (analysis == 'recessive'){
      ncp <- (t*(pa^2)*((b-pd)^2)*(t*pd-1))/(pd*(-1+b-b*t+t*pd)*(1-pd+(pa^2)*(-1+b-b*t+t*pd)))
      
    }else{
      ncp <- -(t*pa*((b-pd)^2)*(-1+t*pd)*(2*(h^2)*(-1+pa)*(1+2*(-1+pa)*pa)*(1+b*(-1+t)-t*pd)+pa*(-1+pd+(pa^2)*(1+b*(-1+t)-t*pd))+h*pa*(-b*(-1+t)*((1-2*pa)^2)+(-1+t)*pd+4*(-1+pa)*pa*(-1+t*pd))))/((pd*(-1+b-b*t+t*pd)*(-1-2*h*(-1+pa)*pa+(pa^2)-b*(-1+t)*(h+2*h*(-1+pa)*pa-(pa^2))+pd-t*(pa^2)*pd+h*(-1+t+2*t*(-1+pa)*pa)*pd)*(-1+pd+(2*h*(-1+pa)-pa)*pa*(-1+b-b*t+t*pd))))
      
    }
    
    
    data[match(t, inflation)] <- ncp + (df/samp.size)
  }
  plot(inflation*pd, data,  xlab = "Case Fraction", ylab = expression(lambda), type = 'l', las = 1, bty = "n", pch =19, xlim = c(0, 1), ylim = c(0, max(data)*1.1))
  
} 

#functions which calculate the ncp at different inflation factors using noncentrality parameter
#the function plots line on graph and prints case fraction at which power is greatest
ncpcalculation <- function(pd, pa, samp.size, df, b, h, analysis, x){
  inflation <- seq(0, x, length.out = 1000)
  
  data <- vector()
  
  for (t in inflation){
    
    #calculate non-centrality parameter
    if (analysis == 'dominant'){
      ncp <- (t*(pa-1)*(pa*(-2*h*(pa-1)+pa)^2)*((b-pd)^2)*(t*pd-1))/(pd*(pd-1+(2*h*(pa-1)-pa)*pa*(b-1-b*t+t*pd))*(-2-b*(t-1)*(2*h*(pa-1)-pa)*(pa-1)-pa+2*pd+pa*(pa+t*pd-t*pa*pd)+2*h*((t-1)*pd+(pa-2)*pa*(-1+t*pd))))
      
    }else if (analysis == 'recessive'){
      ncp <- (t*(pa^2)*((b-pd)^2)*(t*pd-1))/(pd*(-1+b-b*t+t*pd)*(1-pd+(pa^2)*(-1+b-b*t+t*pd)))
      
    }else{
      ncp <- -(t*pa*((b-pd)^2)*(-1+t*pd)*(2*(h^2)*(-1+pa)*(1+2*(-1+pa)*pa)*(1+b*(-1+t)-t*pd)+pa*(-1+pd+(pa^2)*(1+b*(-1+t)-t*pd))+h*pa*(-b*(-1+t)*((1-2*pa)^2)+(-1+t)*pd+4*(-1+pa)*pa*(-1+t*pd))))/((pd*(-1+b-b*t+t*pd)*(-1-2*h*(-1+pa)*pa+(pa^2)-b*(-1+t)*(h+2*h*(-1+pa)*pa-(pa^2))+pd-t*(pa^2)*pd+h*(-1+t+2*t*(-1+pa)*pa)*pd)*(-1+pd+(2*h*(-1+pa)-pa)*pa*(-1+b-b*t+t*pd))))
      
    }
    
    
    data[match(t, inflation)] <- ncp + (df/samp.size)
  }
  lines(inflation*pd, data, type = "l") 
}

diploidncpsimulation <- function(pd, pa, b, h, samp.size, analysis, inflation = 
                                   c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, 
                                 pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- creatediploidpopulation(pd, pa, b, h, analysis, pop.size)
  
  #inflation vector contains the different inflation factors that will be used
  ncp <- numeric(length(inflation))
  upperncp <- numeric(length(inflation))
  lowerncp <- numeric(length(inflation))
  
  #subsets the population into cases and controls based on phenotype
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  #loops through the different inflation factors to run chi-squared test
  for(t in inflation){
    chisq.stats <- numeric(n.sim)
    chisq.ps <- numeric(n.sim)
    case.inflation <-t
    
    #repeats for the desired number of simulations
    for(i in 1:n.sim){
      num.casesample <- case.inflation*pd*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 2,]), 
                    nrow(casesample[casesample$pop.genotypes == 1,]),
                    nrow(casesample[casesample$pop.genotypes == 0,]) )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 2,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 1,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 0,]) ) 
      
      sample <- rbind(casesample, controlsample)
      chisq.fit <- chisq.test(sample$pop.genotypes, sample$pop.phenotypes)
      
      #stores chi-squared statistic and p-value
      chisq.stats[i] <- chisq.fit$statistic
      chisq.ps[i] <- chisq.fit$p.value
    }
    
    #average value of non-centrality parameter
    ncp[match(t,inflation)] <- mean(chisq.stats)/samp.size
    upperncp[match(t,inflation)] <- mean(chisq.stats)/samp.size + 
      (2*sd(chisq.stats/samp.size)/sqrt(length(chisq.stats)))
    lowerncp[match(t,inflation)] <- mean(chisq.stats)/samp.size - 
      (2*sd(chisq.stats/samp.size)/sqrt(length(chisq.stats)))
    
  }
  #plots the data points for ncp
  plot(inflation*pd, ncp, 
       xlab = "Case Fraction", ylab = "Noncentrality Parameter", las = 1, bty = "n", pch =19, xlim = c(0, 1), ylim = c(0,max(ncp)*1.1))
  for (t in inflation){
    x2 <- c(t*pd, t*pd)
    y2 <- c((upperncp[match(t,inflation)]), (lowerncp[match(t,inflation)]))
    print(x2)
    print(y2)
    lines(x2, y2)
  }
  
  #overlays ncp values from function
  ncpcalculation(pd, pa, samp.size, df, b, h, analysis, 1/pd)
  
}

ncppowersimulation <- function(pd, pa, b, h, samp.size, analysis, 
                               inflation = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3),
                               n.sim = 1000, pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- creatediploidpopulation(pd, pa, b, h, analysis, pop.size)
  
  #inflation vector contains the different inflation factors that will be used
  power <- numeric(length(inflation))
  upperpower <- numeric(length(inflation))
  lowerpower <- numeric(length(inflation))
  
  
  #subsets the population into cases and controls based on phenotype
  cases <- pop[pop$pop.phenotypes == 1,]
  controls <- pop[pop$pop.phenotypes == 0,]
  
  #loops through the different inflation factors to get the chi-squared value
  for(t in inflation){
    chisq.stats <- numeric(n.sim)
    chisq.ps <- numeric(n.sim)
    case.inflation <-t
    
    #repeats for the desired number of simulations
    for(i in 1:n.sim){
      num.casesample <- case.inflation*pd*samp.size
      num.controlsample <- samp.size - num.casesample
      casesample <- (cases[sample(nrow(cases), num.casesample),])
      controlsample <- (controls[sample(nrow(controls),num.controlsample),])
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 2,]), 
                    nrow(casesample[casesample$pop.genotypes == 1,]), 
                    nrow(casesample[casesample$pop.genotypes == 0,]) )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 2,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 1,]), 
                       nrow(controlsample[controlsample$pop.genotypes == 0,]) ) 
      
      sample <- rbind(casesample, controlsample)
      chisq.fit <- chisq.test(sample$pop.genotypes, sample$pop.phenotypes)
      
      #stores chi-squared statistic and p-value
      chisq.stats[i] <- chisq.fit$statistic
      chisq.ps[i] <- chisq.fit$p.value
    }
    
    #approximation for power
    power[match(t, inflation)] <- mean(chisq.ps <= alpha)
    
    #confidence intervals
    upperpower[match(t,inflation)] <- mean(chisq.ps <= alpha) + 
      (2*sqrt((power[match(t, inflation)])*(1-(power[match(t, inflation)])))/sqrt(n.sim))
    lowerpower[match(t,inflation)] <- mean(chisq.ps <= alpha) - 
      (2*sqrt((power[match(t, inflation)])*(1-(power[match(t, inflation)])))/sqrt(n.sim))
    
  }
  
  #plots data points for power
  plot(inflation*pd, power, xlab = "Case Fraction", ylab = "Power", las = 1, bty = "n", pch =19, xlim = c(0, 1))
  for (t in inflation){
    x2 <- c(t*pd, t*pd)
    y2 <- c((upperpower[match(t,inflation)]), (lowerpower[match(t,inflation)]))
    lines(x2, y2)
  }
  
  #overlays values of power calculated using calculated ncp and prints case frac
  #tion at which power is greatest
  powercalculation(pd, pa, samp.size, analysis, alpha, df, b, h, max(inflation))
}

#function which calculates power at different inflation factors using noncentrality parameter
#function plots line on graph and prints case fraction at which power is greatest
powercalculation <- function(pd, pa, samp.size, analysis, alpha, df, b, h, x){
  
  inflation <- seq(0.0001, x, length.out = 1000)
  
  data <- vector()
  
  for (t in inflation){
    
    #calculate non-centrality parameter
    if (analysis == 'dominant'){
      ncp.base <- (t*(pa-1)*(pa*(-2*h*(pa-1)+pa)^2)*((b-pd)^2)*(t*pd-1))/
        (pd*(pd-1+(2*h*(pa-1)-pa)*pa*(b-1-b*t+t*pd))*(-2-b*(t-1)*(2*h*(pa-1)-pa)*
                                                        (pa-1)-pa+2*pd+pa*(pa+t*pd-t*pa*pd)+2*h*((t-1)*pd+(pa-2)*pa*(-1+t*pd))))
      
    }else if (analysis == 'recessive'){
      ncp.base <- (t*(pa^2)*((b-pd)^2)*(t*pd-1))/(pd*(-1+b-b*t+t*pd)*(1-pd+(pa^2)
                                                                      *(-1+b-b*t+t*pd)))
      
    }else{
      ncp.base <- -(t*pa*((b-pd)^2)*(-1+t*pd)*(2*(h^2)*(-1+pa)*(1+2*(-1+pa)*pa)*
                                                 (1+b*(-1+t)-t*pd)+pa*(-1+pd+(pa^2)*(1+b*(-1+t)-t*pd))+h*pa*(-b*(-1+t)*((1-2*pa)^2)+
                                                                                                               (-1+t)*pd+4*(-1+pa)*pa*(-1+t*pd))))/((pd*(-1+b-b*t+t*pd)*(-1-2*h*(-1+pa)*pa+(pa^2)-b*(-1+t)*
                                                                                                                                                                           (h+2*h*(-1+pa)*pa-(pa^2))+pd-t*(pa^2)*pd+h*(-1+t+2*t*(-1+pa)*pa)*pd)*(-1+pd+(2*h*(-1+pa)-pa)*pa*(-1+b-b*t+t*pd))))
      
    }
    
    ncp <- ncp.base * samp.size
    
    #calculate critical value of chi-squared
    criticalvalue <- qchisq(1-alpha, df)
    
    #calculate power
    power <- pchisq(criticalvalue, df, ncp, lower.tail=FALSE, log.p = FALSE)
    data[match(t, inflation)] <- power
  }
  lines(inflation*pd, data, xlab = "Case Fraction", ylab = "Power", type = 'l', las = 1, bty = "n", pch =19)
  
}

