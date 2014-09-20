library(ecodist)
library(car)
library(AICcmodavg)
library(MuMIn)
library(gplots)
library(Hmisc)

setwd("H:/PhD/Landscape_Genetics/Data")
MP_individuals <- read.table("MP_Genotyped_Ringed_Individuals_5.csv", sep=",", header=TRUE)
MP_Habitat <- read.table("MP_Populations_Revised_Habitat_5.csv", sep=",", header=TRUE)
MP<-merge(x=MP_individuals, y=MP_Habitat, by="Pop")

MP.3.bc <- as.matrix(distance(MP_Habitat$Clim_3, "bray-curtis"))
MP.12.bc <- as.matrix(distance(MP_Habitat$Clim_12, "bray-curtis"))
MP_Path <- as.matrix(read.table("MP_Cost_Path_Matrix_Clean.csv", sep=",", header=FALSE))
MP_Path2 <- as.matrix(read.table("MP_Coast_Watershed_Path_Matrix_10km.csv", sep=",", header=FALSE))
MP_Path3 <- as.matrix(read.table("MP_Distance_Matrix.csv", sep=",", header=FALSE))
MP_Path4 <- as.matrix(read.table("MP_Cost_Path_Matrix_New.csv", sep=",", header=FALSE))
MP_Sample <- as.matrix(read.table("MP_Sample_Size_Matrix.csv", sep=",", header=FALSE))
MP_Fst <- as.matrix(read.table("MP_Fst_Matrix.csv", sep=",", header=FALSE))
MP_Fst[MP_Fst <= 0] <- 0

MakeNullModel<- function(varin, varout)
{
  varin[col(varin) == row(varin) | lower.tri(varin)] <- NA
  varout <- c(varin, recursive=TRUE)
  varout<- varout[!is.na(varout)]
}

MP.Fst.Null<- MakeNullModel(MP_Fst, MP.fst.Null)
MP.3.bc.list.Null<- MakeNullModel(MP.3.bc, MP.3.bc.List.Null)
MP.12.bc.list.Null<- MakeNullModel(MP.12.bc, MP.12.bc.List.Null)
MP.Path.list.Null<- MakeNullModel(MP_Path, MP.Path.list.Null)
MP.Path2.list.Null<- MakeNullModel(MP_Path2, MP.Path2.list.Null)
MP.Path3.list.Null<- MakeNullModel(MP_Path3, MP.Path3.list.Null)
MP.Path4.list.Null<- MakeNullModel(MP_Path4, MP.Path4.list.Null)
MP.Sample.list.Null<- MakeNullModel(MP_Sample, MP.Sample.list.Null)
MP.Full <- cbind(MP.Fst.Null,MP.Path.list.Null,MP.12.bc.list.Null,MP.3.bc.list.Null)
MP.Full <- as.data.frame(MP.Full)
attach(MP.Full)
summary(MP.Full)

MP.Model <- list()
MP.Model[[1]]<- lm(MP.Fst.Null~ 1)
MP.Model[[2]]<- lm(MP.Fst.Null~ MP.Path.list.Null)
MP.Model[[3]]<- lm(MP.Fst.Null~ MP.Path.list.Null+MP.12.bc.list.Null)
MP.Model[[4]]<- lm(MP.Fst.Null~ MP.Path.list.Null+MP.3.bc.list.Null)
MP.Model[[5]]<- lm(MP.Fst.Null~ MP.12.bc.list.Null)
MP.Model[[6]]<- lm(MP.Fst.Null~ MP.3.bc.list.Null)
MP.Model[[7]]<- lm(MP.Fst.Null~ MP.Path.list.Null*MP.12.bc.list.Null)
MP.Model[[8]]<- lm(MP.Fst.Null~ MP.Path.list.Null*MP.3.bc.list.Null)
#MP.Model[[7]]<- lm(MP.Fst.Null~ MP.Path2.list.Null)
#MP.Model[[8]]<- lm(MP.Fst.Null~ MP.Path2.list.Null+MP.12.bc.list.Null)
#MP.Model[[9]]<- lm(MP.Fst.Null~ MP.Path2.list.Null+MP.3.bc.list.Null)
#MP.Model[[10]]<- lm(MP.Fst.Null~ MP.Path3.list.Null)
#MP.Model[[11]]<- lm(MP.Fst.Null~ MP.Path3.list.Null+MP.12.bc.list.Null)
#MP.Model[[12]]<- lm(MP.Fst.Null~ MP.Path3.list.Null+MP.3.bc.list.Null)
MP.Modnames <- paste("model", 1:length(MP.Model), sep = " ")
MP.model_table<-aictab(cand.set = MP.Model, modnames = MP.Modnames, sort = TRUE)
MP.model_table

min(MP.Path.list.Null)
max(MP.Path.list.Null)
min(MP.Path2.list.Null)/1000
max(MP.Path2.list.Null)/1000
max(MP.Path3.list.Null)/1000

# Normality of Residuals
# qq plot for studentized resid
qqPlot(MP.Model[[2]], main="QQ Plot")
# distribution of studentized residuals
library(MASS)
sresid <- studres(MP.Model[[2]]) 
hist(sresid, freq=FALSE, 
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)

plot(MP.Model[[2]])
plot(MP_Path3,MP_Path4)
summary(MP.Model[[2]])
summary(MP.Model[[12]])
summary(MP.Model[[4]])


par(mar=c(5,5,2,2))  
plot(MP.Path.list.Null,MP.Fst.Null,xlim=c(0,80000000000),xlab="Cost-weighted dispersal distance between populations", ylab= expression(Genetic~Dissimilarity~(italic(F[ST]))),las=1,xaxt ="n",cex.lab=1.3, pch=20)
#axis(1, at=c(0,8000000000,16000000000,24000000000,32000000000,40000000000,48000000000,56000000000,64000000000,72000000000,80000000000),labels=c("","","","","","","","","","",""),tick = TRUE)
#axis(1, at=c(443758848,76667502592), labels=c("Low Cost", "High Cost"),tick = TRUE)
axis(1, at=c(6000000000,75000000000), labels=c("Low Cost", "High Cost"),tick = TRUE)
newx<-seq(min(MP.Path.list.Null),max(MP.Path.list.Null),1000000)
prd<-predict(MP.Model[[2]],newdata=data.frame(MP.Path.list.Null=newx),interval = c("confidence"),level = 0.95,type="response")
lines(newx,prd[,1],col="black",lty=1)
lines(newx,prd[,2],col="black",lty=2)
lines(newx,prd[,3],col="black",lty=2)
text(8000000000,0.3,expression(italic(r) == 0.287),cex=1)
text(9000000000,0.27,expression(italic(P) == 0.0108),cex=1)
#axis(1, at=MP.Path.list.Null,labels=round(MP.Path.list.Null,digits=2),col.axis="blue", las=2, cex.axis=0.7, tck=-.01)
#minor.tick(nx=3,tick.ratio=1)
#text(0.4,2.2,"Y = 0.7293(X) + 0.7617",cex=1.3)

min(MP.Path.list.Null)
dev.new()
library(gplots)
library(Hmisc)

vif(MP.Model[[4]]) # variance inflation factors 
sqrt(vif(MP.Model[[4]])) > 2 # problem?

library(gvlma)
gvmodel <- gvlma(MP.Model[[2]]) 
summary(gvmodel)

Avg.MP<-model.avg(MP.Model, subset = delta < 2)
summary(Avg.MP)

### get AIC for null model, here we should also define the link function. 
AICnull<- AIC(glm(MP.Fst.Null~ MP.Path.list.Null+MP.3.bc.list.Null))
NullModel<- glm(MP.Fst.Null~ MP.Path.list.Null+MP.3.bc.list.Null)
TValuesNullModel<- data.frame(summary(NullModel)$coefficient[,3])

scrambling<- function(variable, varout)
{
  mtemp <- variable
  nameVector<- names(mtemp)
  NewOrder  <- sample(c(1:ncol(mtemp)), ncol(mtemp) ,replace=FALSE)        ### these are the lines that do the scrambling
  NewMatrix <- mtemp[NewOrder,NewOrder]                                    ### these are the lines that do the scrambling
  names(NewMatrix)<- nameVector[NewOrder]                                  ### place the correct column names
  NewMatrix[col(NewMatrix) == row(NewMatrix) | upper.tri(NewMatrix)] <- NA ### set the diagonal of the matrix and the duplicates to NA(upper triangle)
  ListM <- c(NewMatrix, recursive=TRUE)                                    ### convert the matrix to list
  varout<- ListM[!is.na(ListM)]								 ### remove the NA from the list, these are always at the same location
}

threshold<- 10000  	### initial number of runs, I have put this at 10.000
counter <- array(NA,threshold)
counter2 <- array(NA,threshold)	### to count
for(i in 1:threshold){
  # select the variables from the following list by using ## in front of the ones not to be used
  MP_Path_temp<- scrambling(MP_Path,MP_Path_temp)
  #MP.3.bc.temp<- scrambling(MP.3.bc,MP.3.bc.temp)
  #tempind3<- scrambling(IndependentTemp3,tempind3)  ### scramble the variable
  #tempind4<- scrambling(IndependentTemp4,tempind4)
  
  # get the AIC for the new model:
  ScramModel<- glm(MP.Fst.Null~ MP_Path_temp+MP.3.bc.list.Null)
  AICtemp<- AIC(ScramModel) 
  TValuesScramModel<- data.frame(summary(ScramModel)$coefficient[,3])
  counter[i]<- (AICtemp<AICnull)  ### add to the counter vector whether the new AIC is lower than the old
  counter2[i]<- (abs(TValuesScramModel[3,])>=abs(TValuesNullModel[5,]))  # 2 for tempind1, 3 for tempind2, 4 for tempind3, 5 for tempind4	
  par(mar=c(0,0,1,0))
  plot(0,0, xlim=c(-1,-2), ylim=c(-1,-2), xaxt="n", bty="n", yaxt="n", main=i, xlab="", ylab="") # built in counter on screen, to see how far it is.
}
pvalue_indep<- sum(counter, na.rm=TRUE)/threshold    # line for p-value based on AIC
pvalue_indep2<- sum(counter2, na.rm=TRUE)/threshold  # line for p-value based on t-value
pvalue_indep                                         # p-value AIC
pvalue_indep2                                        # p=value based on t-values
