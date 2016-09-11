# This file contain all of the functions needed for the general 
# analysis of all Six data sets.

## The Following function returns a summing table of the data for a given cariable
## ( i.e. number of samples, Mean, SD,...)
aggregate.data <- function(variable.name,d1){
  variable.mean <- aggregate(d1[,paste(variable.name)],by=list(d1$lab,d1$strain),mean,na.rm=T)
  colnames(variable.mean) <- c("lab","strain","mean")
  variable.sd <- aggregate(d1[,paste(variable.name)],by=list(d1$lab,d1$strain),sd,na.rm=T)
  colnames(variable.sd) <- c("lab","strain","sd")
  variable.n <- aggregate(d1[,paste(variable.name)],by=list(d1$lab,d1$strain),function(x){return(sum(!is.na(x)))})
  colnames(variable.n) <- c("lab","strain","n")
  
  
  variable <- merge(variable.mean,variable.sd,by=(intersect(names(variable.mean[,c("lab"   ,"strain")]), names(variable.sd[,c("lab","strain")]))))
  variable <- merge(variable,variable.n,by=(intersect(names(variable[,c("lab","strain")]),   names(variable.n[,c("lab","strain")]))))
  
  return(variable)
}


#### ANOVA table - (fixed & mixed) + Interaction ####
bulid.anova.table <- function(d,y){
  fixed.anova <- aov(y ~ strain*lab,data=d)
  tbl <- as.data.frame(summary(fixed.anova)[[1]])
  fit.lme <<- lme(y ~ strain, random=~1|lab/strain,data=d,method="REML",na.action="na.omit")
  tbl[1,6] <- anova(fit.lme)["strain","F-value"]
  tbl[1,7] <- anova(fit.lme)["strain","p-value"]
  colnames(tbl) <- c("Df","Sum Sq","Mean Sq","F.FLM","p.FLM","F.RLM","p.RLM")
  
  s2.int <- as.numeric(VarCorr(fit.lme)[4,1])
  return(list(tbl,s2.int))
}

### pairwise multi-lab comparisons ####
## Generate a table for multi-Lab t-tests ( Both fixed effect model & random lab model  based t-tests)
pairwise.multi.comparisons.function <- function(d,y,variable) {
  S <- length(unique(d$strain))
  L <- length(unique(d$lab))
  n <- mean(variable$n)
  
  tbl <- data.frame(matrix(0,ncol=11,nrow=dim(combn(levels(d$strain),2))[2]))
  colnames(tbl) <- c("strain1","strain2","diff","Std.Error.FLM","df.FLM","t.FLM","p.FLM","Std.Error.RLM","df.RLM","t.RLM","p.RLM")   
  
  #fit.lme <- lme(y ~ strain, random=~1|lab/strain,data=d,method="REML",na.action='na.omit')
  anova.fit.lme <- anova(fit.lme)
  
  
  pairwise <- summary(glht(fit.lme, linfct=mcp(strain="Tukey"),na.action='na.omit'))
  temp <- cbind(pairwise$test$coefficients,pairwise$test$sigma,pairwise$test$tstat)
  strains <- strsplit(rownames(temp), " - ")
  for (i in 1:length(strains)){
    tbl[i,"strain1"] <- strains[[i]][1]
    tbl[i,"strain2"] <- strains[[i]][2]
  }
  tbl[,c("diff","Std.Error.RLM","t.RLM")] <- temp
  tbl[,"df.RLM"] <- rep(anova.fit.lme["strain","denDF"],dim(tbl)[1])
  tbl[,"p.RLM"] <- pt(abs(tbl[,"t.RLM"]),df=tbl[,"df.RLM"],lower.tail=F)*2
  
  
  anova.tbl <-as.data.frame(summary(aov(y ~ strain*lab,data=d))[[1]])
  tbl[,"df.FLM"] <- anova.tbl["Residuals","Df"]
  
  tem=aggregate(variable$n, by = list(variable$strain), sum)
  colnames(tem)=c('strain','1/n')
  tem[,'1/n']=1/tem[,'1/n']
  tbl=cbind(tbl,rep(0,dim(tbl)[1]))
  colnames(tbl)[12]='1/n'
  for (ss in 1:(dim(tem)[1])){
    tbl[tbl$strain1==(tem$strain[ss]),'1/n']=tbl[tbl$strain1==(tem$strain[ss]),'1/n']+
      rep(tem[ss,'1/n'],sum(tbl$strain1==(tem$strain[ss])))
    tbl[tbl$strain2==(tem$strain[ss]),'1/n']=tbl[tbl$strain2==(tem$strain[ss]),'1/n']+
      rep(tem[ss,'1/n'],sum(tbl$strain2==tem$strain[ss]))
  }
  
  tbl[,"Std.Error.FLM"] <- sqrt(anova.tbl["Residuals","Mean Sq"]*tbl[,"1/n"])
  
  tbl[,"t.FLM"] <- abs(tbl[,"diff"])/(tbl[,"Std.Error.FLM"])
  tbl[,"p.FLM"] <- pt(abs(tbl[,"t.FLM"]),df=tbl[,"df.FLM"],lower.tail=F)*2 
  
  return(tbl)
}
#### In lab Analysis ####
# Computing degrees of freedom according to Satterwait 
compute.df.sater <- function(s2int,df.int,s2pooled,n1,n2){
  inv.v <- (((1/n1+1/n2)*s2pooled)^2*1/(n1+n2-2)+(2*s2int)^2*1/df.int)/
    (((1/n1+1/n2)*s2pooled+2*s2int)^2)
  df <- 1/inv.v
  return(df)
}


# Performing a pairwise comparison in a lab (Given the selected pair of genotypes and a lab ) - Immitating the real practical situation of a single researcer int thier single lab:
# the regular t-test Vs. GXL adjusted t-test
pair.comparison.function <- function(d,y,lab,strain1,strain2,anova.str){
  anova.tbl=anova.str[[1]]
  s2int=anova.str[[2]]
  S <- length(unique(d$strain))
  L <- length(unique(d$lab))
  
  g1 <- y[(d$lab==lab)&(d$strain==strain1)]
  n1 <- sum(!is.na(g1))
  g2 <- y[(d$lab==lab)&(d$strain==strain2)]
  n2 <- sum(!is.na(g2))
  
  t.tbl <- t.test(g1,g2,var.equal=T,alternative = "two.sided")
  diff <- unname(t.tbl$estimate["mean of x"] - t.tbl$estimate["mean of y"])
  s2pooled <- ((n1-1)*sd(g1,na.rm=T)^2+(n2-1)*sd(g2,na.rm=T)^2)/(n1+n2-2)
  t.FLM <- ifelse(s2pooled==0,0,abs(diff)/sqrt(s2pooled*(1/n1+1/n2)))
  df.FLM <- unname(t.tbl$parameter)
  p.FLM <- pt(t.FLM,df=df.FLM,lower.tail=F)*2
  
  df.int <- anova.tbl[3,"Df"]
  df.RLM <- compute.df.sater(s2int,df.int,s2pooled,n1,n2)
  t.RLM <- abs(diff)/sqrt((s2pooled*(1/n1+1/n2)+2*s2int))
  p.RLM <- pt(t.RLM,df=df.RLM,lower.tail=F)*2
  
  return(list(diff=diff,s2pooled=s2pooled,df.FLM=df.FLM,t.FLM=t.FLM,p.FLM=p.FLM,
              df.RLM=df.RLM,t.RLM=t.RLM,p.RLM=p.RLM))
}


## This function takes as input the results of the Fixed model (via anova.tbl) 
## and the Mixed Model (via pairwise) for measure y to create a table containing all pairwise comparisons in each lab separately.
## Using the "pair.comparison.function" function,it returns the regular t-test and GXL adjusted test for each pair in a lab. 

within.lab.function <- function (d,y,anova.str){
  anova.tbl=anova.str[[1]]
  s2.int=anova.str[[2]]
  lab <- levels(d$lab)
  strain <- levels(d$strain)
  strains <- t(combn(strain,2))
  colnames(strains) <- c("strain1","strain2")
  strain2 <- strains[,1] # To be consistent with the other tables
  strain1 <- strains[,2]
  tbl <- data.frame(matrix(0,ncol=11,nrow=dim(combn(levels(d$strain),2))[2]*length(lab)))
  colnames(tbl) <- c("lab","strain1","strain2","diff","s2pooled","df.FLM",
                     "t.FLM","p.FLM","df.RLM","t.RLM","p.RLM")
  
  tbl[,"lab"] <-  rep(lab,each=length(strains[,"strain1"]))
  
  tbl[,"strain1"] <- rep(strain1,times=length(lab))
  tbl[,"strain2"] <- rep(strain2,times=length(lab))
  for (i in 1:dim(tbl)[1]){
    pair <- pair.comparison.function(d,y,tbl[i,"lab"],tbl[i,"strain1"],tbl[i,"strain2"],anova.str)
    tbl[i,"diff"] <- pair$diff
    tbl[i,"s2pooled"] <- pair$s2pooled
    tbl[i,"df.FLM"] <- pair$df.FLM
    tbl[i,"t.FLM"] <- pair$t.FLM 
    tbl[i,"p.FLM"] <- pair$p.FLM 
    tbl[i,"df.RLM"] <- pair$df.RLM
    tbl[i,"t.RLM"] <- pair$t.RLM
    tbl[i,"p.RLM"] <- pair$p.RLM
    
  }
  
  return(tbl) 
  
}