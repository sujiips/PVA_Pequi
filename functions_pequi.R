begin  <-  function( parms ){  
  for(g in 1:length(germ)){
    for(s in 1:length(selection)){ 
      for(m in 1:length(maxAge)){
        for(av in 1:length(avSeed)){  
  fsRate <- function(age){
    ifelse(age==1, 0.2, 
           ifelse(age < 3, 0.52,
                  ifelse(age < 7, 0.8,
                     ifelse(age >= maxAge[m], 0, 0.925))))
  }
  
  f_reproduction <- function( ind1 = NULL, ind2 = NULL, seedInd){
    theta <- runif(seedInd*fathers, 0, 2*pi)
    disp <- 1:dSeed
    probDisp <- (0.6942213/81.4806772)*((disp/81.4806772)^(0.6942213-1))*exp(-(disp/81.4806772)^0.6942213)
    x <- as.numeric(dat[ind1,1]) + round(cos(theta)*sample(disp, seedInd*fathers, prob=probDisp, replace = T))
    y <- as.numeric(dat[ind1,2]) + round(sin(theta)*sample(disp, seedInd*fathers, prob=probDisp, replace = T))
    genotypes <- matrix(unlist(c(
      sample(dat[ind1,3:4], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,3:4], seedInd, replace=T))),
      sample(dat[ind1,5:6], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,5:6], seedInd, replace=T))),
      sample(dat[ind1,7:8], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,7:8], seedInd, replace=T))),
      sample(dat[ind1,9:10], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,9:10], seedInd, replace=T))),
      sample(dat[ind1,11:12], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,11:12], seedInd, replace=T))),
      sample(dat[ind1,13:14], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,13:14], seedInd, replace=T))),
      sample(dat[ind1,15:16], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,15:16], seedInd, replace=T))),
      sample(dat[ind1,17:18], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,17:18], seedInd, replace=T))),
      sample(dat[ind1,19:20], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,19:20], seedInd, replace=T))),
      sample(dat[ind1,21:22], seedInd*fathers, replace=T), as.vector(t(sample(dat[ind2,21:22], seedInd, replace=T)))
    )), nrow=seedInd*fathers, byrow=F)
    genotypes <- data.frame(genotypes)
    colnames(genotypes) <- c("l1_1","l1_2","l2_1","l2_2","l3_1","l3_2","l4_1","l4_2","l5_1","l5_2",  
                             "l6_1","l6_2","l7_1","l7_2","l8_1","l8_2","l9_1","l9_2","l10_1","l10_2")
    age <- rep(0, seedInd*fathers)
    a1 <- genotypes[,seq( 1,(nloci*2),2 )]
    a2 <- genotypes[,seq(2,(nloci*2),2 )]
    homo <- a1==a2
    ho <- apply(homo, 1, sum)
    Ho <- 1-(ho/nloci)
    Fis <- ifelse(He == 0, 1,  1-(Ho/He))
    sRate <- ifelse(Fis <= 0, germ[g], germ[g]*(1-selection[s]*Fis))   ## Probability of germination - function of inbreeding    
    parent1 <- ind1
    parent2 <- rep(ind2, each=seedInd)
    old <- c(ind1, ind2)
    x1 <- adult$x
    y1 <- adult$y
    f_sameSpot <- function(o){
      which( x1[o] == x & y1[o] == y )}
    sameSpot <- unique(unlist(lapply(old, f_sameSpot)))
    age[sameSpot] <- NA
    new_ind <<- cbind(x, y, genotypes, age, parent1, parent2, sRate)
  }
  
  basic.stats1 <- function (data, diploid = TRUE, digits = 4) {
    loc.names <- names(data)[-1]
    if (length(table(data[, 1])) < 2) 
      data[dim(data)[1] + 1, 1] <- 2
    p <- pop.freq(data, diploid)
    n <- t(ind.count(data))
    if (diploid) {
      dum <- getal.b(data[, -1])
      Ho <- dum[, , 1] == dum[, , 2]
      sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x, 
                                                           data[, 1], mean, na.rm = TRUE))))
      mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    }
    else {
      sHo <- NA
      mHo <- NA
    }
    sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
    sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
    if (diploid) {
      Hs <- (1 - sp2 - sHo/2/n)
      Hs <- n/(n - 1) * Hs
      Fis  <- ifelse(Hs==0, Fis <- 1, Fis <- 1 - sHo/Hs)
    }
    else {
      Hs <- n/(n - 1) * (1 - sp2)
      Fis <- NA
    }
    np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
    mn <- apply(n, 1, fun <- function(x) {
      np <- sum(!is.na(x))
      np/sum(1/x[!is.na(x)])
    })
    msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
    mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
    mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
    if (diploid) {
      mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
      Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
      mFis <- ifelse(Ht==0, mFis <- 1, mFis  <- 1 - mHo/mHs)
    }
    else {
      mHs <- mn/(mn - 1) * (1 - msp2)
      Ht <- 1 - mp2 + mHs/mn/np
      mFis <- NA
    }
    Dst <- Ht - mHs
    Dstp <- np/(np - 1) * Dst
    Htp = mHs + Dstp
    Fst = Dst/Ht
    Fstp = Dstp/Htp
    Dest <- Dstp/(1 - mHs)
    res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst, 
                            Fstp, mFis, Dest))
    names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp", "Fst", 
                    "Fstp", "Fis", "Dest")
    if (diploid) {
      rownames(sHo) <- loc.names
      rownames(Fis) <- loc.names
    }
    overall <- apply(res, 2, mean, na.rm = TRUE)
    overall[7] <- overall[4]/overall[3]
    overall[8] <- overall[6]/overall[5]
    overall[9] <- 1 - overall[1]/overall[2]
    overall[10] <- overall[6]/(1 - overall[2])
    names(overall) <- names(res)
    all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
                                                      digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
                    Fis = round(Fis, digits), perloc = round(res, digits), 
                    overall = round(overall, digits))
    class(all.res) <- "bas.stats"
    all.res
  }
  
  ## Start simulations
  for (r in Start:rep){
    dat <- read.table( 'datain.txt', header = TRUE, sep = '\t')
    
    sRate <- fsRate(dat$age)  
    dat <- cbind(dat, sRate)
    
    dir.create("/Users/sujii/Model_pequi_2/out1")  ##escreva o diretorio aqui, mas mantenha o /out1
    setwd("/Users/sujii/Model_pequi_2/out1")   ##escreva o diretorio aqui, mas mantenha o /out1
    
    for(t in 1:time){
      print(t)
      ## Define new sRate 
      dat$sRate <- fsRate(dat$age)  
      
      ### Reproduction
      ## Expected heterozygosity - HWE
      l1 <- c(dat[,3], dat[,4]) 
      l2 <- c(dat[,5], dat[,6])
      l3 <- c(dat[,7], dat[,8])
      l4 <- c(dat[,9], dat[,10])
      l5 <- c(dat[,11], dat[,12])
      l6 <- c(dat[,13], dat[,14])
      l7 <- c(dat[,15], dat[,16])
      l8 <- c(dat[,17], dat[,18])
      l9 <- c(dat[,19], dat[,20])
      l10 <- c(dat[,21], dat[,22])
      
      He1 <- 1-sum((table(l1)/length(l1))^2)
      He2 <- 1-sum((table(l2)/length(l2))^2)
      He3 <- 1-sum((table(l3)/length(l3))^2)
      He4 <- 1-sum((table(l4)/length(l4))^2)
      He5 <- 1-sum((table(l5)/length(l5))^2)
      He6 <- 1-sum((table(l6)/length(l6))^2)
      He7 <- 1-sum((table(l7)/length(l7))^2)
      He8 <- 1-sum((table(l8)/length(l8))^2)
      He9 <- 1-sum((table(l9)/length(l9))^2)
      He10 <- 1-sum((table(l10)/length(l10))^2)
      
      He <- mean(He1, He2, He3, He4, He5, He6, He7, He8, He9, He10)
      
      
      ## Find reproductive pairs 
      adult <- subset(dat, dat$age >= adultAge)
      
      if(length(adult[,1]) > 0){
        for(parent in 1:length(adult[,1])){
          pairs <- expand.grid(parent, 1:length(adult[,1]))
          neigh <- which(apply(pairs, 1, function(p){
            (adult$x[p[1]] - adult$x[p[2]])^2 + (adult$y[p[1]]-adult$y[p[2]])^2 <= dPollen^2}))
          
          if(length(neigh) < maxFathers){
            seedInd <- round(avSeed[av]/length(neigh))
            sexPair <- neigh
            fathers <- length(sexPair)
          } else {
            seedInd <- round(avSeed[av]/maxFathers)
            sexPair <- sample(neigh, maxFathers, replace=T)
            fathers <- maxFathers
          }
          f_reproduction(parent, sexPair, seedInd)
          dat <- rbind(dat, new_ind) 
        }
      }
      
      ### Checking if the individuals will survive to the next time
      rand <- runif(length(dat[,1]), 0, 1)
      sRate_death <- which( rand > dat$sRate )
      dat$age[sRate_death] <- NA
      
      out <- which(dat$x < 1 | dat$x > xDim | dat$y < 1 | dat$y > yDim)
      dat$age[out] <- NA
      
      dat <- dat[!is.na(dat$age),]
      
      ## Individuals that survived get 1 year older
      dat$age <- dat$age + 1
      
      
      ### Save output
      if(t %in% output_year){
        mypath <- file.path(getwd(), paste("dataout_", t, ".txt", sep=""))
        write.table( dat, file=mypath, row.names= FALSE, sep = '\t' )
      }
    }
    
    ## Estimate genetic parameters
    setwd("/Users/sujii/Model_pequi_2/out1")  ##escreva o diretorio aqui, mas mantenha o /out1
    directory <- "/Users/sujii/Model_pequi_2/out1"   ##escreva o diretorio aqui, mas mantenha o /out1
    
    filename1 <- list.files(directory, full.names=F)
    if(length(filename1)==0){break}
    else if(length(filename1)!=0){
      file.n <- (length(filename1))
      filename <- filename1[1:file.n]
      n1 <- gsub("dataout_", "", filename)
      n <- as.integer(gsub(".txt", "", n1)) 
      npop <- (1:length(filename))
      
      dir.create("/Users/sujii/Model_pequi_2/out1/fstat")   ##escreva o diretorio aqui, mas mantenha o /out1/fstat
      
      
      for(f in 1:length(filename)){
        setwd("/Users/sujii/Model_pequi_2/out1")    ##escreva o diretorio aqui, mas mantenha o /out1
        DataTrees<-read.table(filename[f], header=TRUE, sep="\t");
        DataTrees <- subset(DataTrees, age > 1)
        DataTrees<-DataTrees[-c(1:2,23:26)]
        
        #Create a column with the required format (combine columns by locus) as defined in the adegenet documentation 
        for (i in seq(1,19,2)){
          DataTrees[,i]<-as.character(paste(DataTrees[,i],DataTrees[,i+1],sep= ""));
        }
        DataTrees2<-DataTrees[,-seq(2,20,2)];   #Remove the now useless second column of each locus
        colnames(DataTrees2)<- gsub("_1","",colnames(DataTrees2));  #Remove the 1 and 2 designations as they are combined
        pop <- rep(n[f], length(DataTrees2[,1]))  # Add pop column
        Data_fstat <- cbind(pop, DataTrees2)
        
        setwd("/Users/sujii/Model_pequi_2/out1/fstat")   ##escreva o diretorio aqui, mas mantenha o /out1/fstat
        write.table(Data_fstat, file = filename[f], quote=F, row.names=F)
      }
      
      ## Make dataset
      all <- function(id) {
        files <- list.files("/Users/sujii/Model_pequi_2/out1/fstat", full.names=F)  ##escreva o diretorio aqui, mas mantenha o /out1/fstat
        all_data <-data.frame()
        for (i in 1:length(files)){
          all_data <- rbind(all_data, read.table(files[i], head=T))
        }
        alldata <- all_data[order(all_data[,1]),]
      }
      
      data <- all(npop)
      
      pop <- length(npop)
      pop.names <- unique(data$pop)
      
      alleles <- allele.count(data, diploid=T)
      
      # Allelic richness per locus
      allelic.r <- allelic.richness(data, min.n=NULL, diploid=T)
      mean.rich <- round(apply(allelic.r$Ar, 2, mean), digits=1)
      
      # Mean number of alleles per locus
      a <- nb.alleles(data)
      mean.a <- round(apply(a, 2, mean), digits=1)
      
      basic <- basic.stats1(data, diploid=T, digits=3)
      
      ## number of individuals per locus
      n.ind <- basic$n.ind.samp
      mean.ind <- round(apply(n.ind, 2, mean), digits=1) 
      
      # Fis per pop
      fis <- data.frame(basic$Fis)
      mean.fis <- round(apply(fis, 2, mean), digits=3)
      
      ## Ho per pop
      Ho <- basic$Ho
      mean.Ho <- round(apply(Ho, 2, mean), digits=3)
      sd.Ho <- round(apply(Ho, 2, sd), digits=3)
      
      ## Hs per pop
      Hs <- basic$Hs
      mean.Hs <- round(apply(Hs, 2, mean), digits=3)
      sd.Hs <- round(apply(Ho, 2, sd), digits=3)
      
      
      ## Summary table
      basic.table <- cbind(pop.names, mean.ind, mean.a, mean.rich, mean.Ho, sd.Ho, mean.Hs, sd.Hs, mean.fis)
      colnames(basic.table) <- c("pop", "ind", "A", "R", "Ho", "Ho(SD)", "Hs", "Hs(SD)", "fis")  
      write.table(basic.table, file="/Users/sujii/Model_pequi_2/out1/fstat/basic.txt", quote=F, row.names=F, sep="\t")   ##escreva o diretorio aqui, mas mantenha o /out1/...
      
    }
    setwd("/Users/sujii/Model_pequi_2")  ##escreva o diretorio aqui
    N <- paste(g,s,m,av, sep="")
    name <- paste("s", N, r, sep="_")
    file.rename("out1", name) 
  }
}
}
}
}
}