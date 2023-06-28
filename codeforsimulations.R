#functions creates a population with genotypes and phenotypes based on the frequencies 
#in the contingency tables in either the diploid or haploid case
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

createhaploidpopulation <- function(pd, pa, b, pop.size = 1000000){
  #creates vector of genotypes
  pop.genotypes <- rbinom(pop.size, 1, pa)
  #creates empty vector for phenotypes
  pop.phenotypes <- numeric(length(pop.genotypes))
  
  
  #assigns phenotypes based on the three genotypes using probabilities from chi-squared tables
  pop.phenotypes[pop.genotypes == 1] <- rbinom(sum(pop.genotypes == 1), 1, b)
  pop.phenotypes[pop.genotypes == 0] <- rbinom(sum(pop.genotypes == 0), 1, (pd-pa*b)/(1-pa))
  
  
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


haploidncpsimulation<- function(pd, pa, b, samp.size, analysis, inflation = 
                                  c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), n.sim = 200, 
                                pop.size = 1000000, alpha = 0.00000005, df = 1, allele = "risk"){
  #creates population
  pop <- createhaploidpopulation(pd, pa, b, pop.size)
  
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
      
      caselist <- c(nrow(casesample[casesample$pop.genotypes == 1,]),
                    nrow(casesample[casesample$pop.genotypes == 0,]) )
      controllist <- c(nrow(controlsample[controlsample$pop.genotypes == 1,]), 
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
  plot(inflation*pd, ncp, main = "NCP as a function of the case fraction", 
       xlab = "Case Fraction", ylab = "Noncentrality Parameter")
  for (t in inflation){
    x2 <- c(t*pd, t*pd)
    y2 <- c((upperncp[match(t,inflation)]), (lowerncp[match(t,inflation)]))
    print(x2)
    print(y2)
    lines(x2, y2)
  }
  
  #overlays ncp values from function
  #haploidncpcalculation(pd, pa, samp.size, df, b, h, allele, max(inflation))
  
}

ncppowersimulation <- function(pd, pa, b, h, samp.size, analysis, 
                               inflation = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3),
                               n.sim = 1000, pop.size = 1000000, alpha = 0.00000005, df = 2){
  #creates population
  pop <- createpopulation(pd, pa, b, h, analysis, pop.size)
  
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

haploidncpcalculation <- function(d, a, samp.size, df, b, h, allele, x){
  inflation <- seq(0, x, length.out = 1000)
  
  data <- vector()
  
  for (t in inflation){
    
    #calculate non-centrality parameter
    if (allele == 'risk'){
      ncp <- ((b-d)^2*a*t*(d*t-1))/(d*(1+b*(t-1)-d*t)*(-1+d+a+a*b*(t-1)-d*a*t))
      
    }else if (allele == 'protective'){
      ncp <- (a*t*(-1+b+d)^2*(-1+t*d))/((b*(-1+t)+t*(-1+d))*(1+a*b*(-1+t)+a*t*(-1+d)-d)*d)
      
    }
    
    data[match(t, inflation)] <- ncp 
  }
  plot(inflation*pd, data, xlab = "Case Fraction", ylab = expression(lambda), type = "l", las = 1, bty = "n", pch =19, xlim = c(0, 1))
}

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

