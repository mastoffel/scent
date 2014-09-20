# This is the code submitted with Kurvers et al. Contrasting context-dependence of familiarity and kinship in animal social networks
# The code is written by Rudy M. Jonker. Contact: mrjonker@gmail.com
# in the code below the code is shown for the female foraging association analysis and the mate choice analysis.


### empty memory
rm(list=ls())
### define folder where files are for females 
setwd("D:/Dropbox/Geese/Testfolder/")

# load files
	DependentTemp<-    read.table("Females_Ass_Strength_M.CSV", sep=",", header=TRUE)
	IndependentTemp1<- read.table("Females-groupmatrix.CSV", sep=",", header=TRUE)
	IndependentTemp2<- read.table("Females_Relatedness_M.CSV", sep=",", header=TRUE)
	IndependentTemp3<- read.table("Females_Distance_Dominance_M.CSV", sep=",", header=TRUE)
	IndependentTemp4<- read.table("Females_Distance_Boldness_M.CSV", sep=",", header=TRUE)

### define folder where files are for males 
setwd("D:/Dropbox/Geese/2.Males/")
	DependentTemp<-    read.table("Males_Ass_Strength_M.CSV", sep=",", header=TRUE)
	IndependentTemp1<- read.table("Males-groupmatrix.CSV", sep=",", header=TRUE)
	IndependentTemp2<- read.table("Males_Relatedness_M.CSV", sep=",", header=TRUE)
	IndependentTemp3<- read.table("Males_Distance_Dominance_M.CSV", sep=",", header=TRUE)
	IndependentTemp4<- read.table("Males_Distance_Boldness_M.CSV", sep=",", header=TRUE)


#### Make Full Null Model
# first function to make column from matrix and remove duplicates
MakeNullModel<- function(varin, varout)
	{
	varin[col(varin) == row(varin) | upper.tri(varin)] <- NA
	varout <- c(varin, recursive=TRUE)
	varout<- varout[!is.na(varout)]
	}

## The Null model will have to be updated when terms are dropped.


# run function for all variables for null model
DepListNull<- MakeNullModel(DependentTemp, DepListNull)
Indep1ListNull<- MakeNullModel(IndependentTemp1, Indep1ListNull)
Indep2ListNull<- MakeNullModel(IndependentTemp2, Indep2ListNull)
Indep3ListNull<- MakeNullModel(IndependentTemp3, Indep3ListNull)
Indep4ListNull<- MakeNullModel(IndependentTemp4, Indep4ListNull)

##################################################################
##################################################################
# Null Model first.
##################################################################
##################################################################

### get AIC for null model, here we should also define the link function. 
AICnull<- AIC(glm(DepListNull~ Indep1ListNull+Indep2ListNull+Indep3ListNull+Indep4ListNull))
NullModel<- glm(DepListNull~ Indep1ListNull+Indep2ListNull+Indep3ListNull+Indep4ListNull)
TValuesNullModel<- data.frame(summary(NullModel)$coefficient[,3])




#Plotting the histogram of the residuals of the Null model
hist(resid(glm(DepListNull~ Indep1ListNull+Indep2ListNull+ Indep3ListNull+Indep4ListNull)), breaks=20)

#############################
### now do the scrambling ###
#############################

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

############# Actual running part per variable ###############
#########################
## PROCEDURE TO FOLLOW ##
#########################
# each variable needs a separate run. In each run, you first have to build the tempind file, which is a scrambled version of the matrix.
# This is done by commenting the correct lines below.
# Important is to also include the tempind object in the list of fixed effects in the model below.
# Same procedure to follow for the pairs.




# set initial number of runs
threshold<- 10000		### initial number of runs, I have put this at 10.000
counter <- array(NA,threshold)
counter2 <- array(NA,threshold)	### to count
for(i in 1:threshold){
# select the variables from the following list by using ## in front of the ones not to be used
#tempind1<- scrambling(IndependentTemp1,tempind1)
#tempind2<- scrambling(IndependentTemp2,tempind2)
tempind3<- scrambling(IndependentTemp3,tempind3)  ### scramble the variable
#tempind4<- scrambling(IndependentTemp4,tempind4)

# get the AIC for the new model:
ScramModel<- glm(DepListNull~ tempind3+ Indep2ListNull +Indep1ListNull +Indep4ListNull)
AICtemp<- AIC(ScramModel) 
TValuesScramModel<- data.frame(summary(ScramModel)$coefficient[,3])
counter[i]<- (AICtemp<AICnull)  ### add to the counter vector whether the new AIC is lower than the old
counter2[i]<- (abs(TValuesScramModel[2,])>=abs(TValuesNullModel[5,]))  # 2 for tempind1, 3 for tempind2, 4 for tempind3, 5 for tempind4	
par(mar=c(0,0,1,0))
plot(0,0, xlim=c(-1,-2), ylim=c(-1,-2), xaxt="n", bty="n", yaxt="n", main=i, xlab="", ylab="") # built in counter on screen, to see how far it is.
}
pvalue_indep<- sum(counter, na.rm=TRUE)/threshold    # line for p-value based on AIC
pvalue_indep2<- sum(counter2, na.rm=TRUE)/threshold  # line for p-value based on t-value
pvalue_indep                                         # p-value AIC
pvalue_indep2                                        # p=value based on t-values


################
### PAIRS
################
setwd("D:/Dropbox/Geese/Testfolder/pairs/")
DependentTemp<- read.table("Breeding_Pairs_M.CSV", sep=",", header=TRUE)
IndependentTemp1<- read.table("Breeding_group_M.CSV", sep=",", header=TRUE)
IndependentTemp2<- read.table("Breeding_Relatedness_M.CSV", sep=",", header=TRUE)
IndependentTemp3<- read.table("Breeding_Boldness_Distance_M.CSV", sep=",", header=TRUE)

#### Make Null Model
# first function to make column from matrix and remove duplicate. Because the males and females had to be scrambled separately, to prevent unrealistic same-sex pairs to occur, the columns that contain males and females are defined first.

MakeNullModel<- function(varin, varout)
{
	males<- c(1,3,9,10,11,12,13,17,19,21,22,27,28,29,31,34,35,37,38,39,41,42,44)
	females<- c(2,4,5,6,7,8,14,15,16,18,20,23,24,25,26,30,32,33,36,40,43)
	nameVector<- names(varin)
	MalesFemalesApart<- c(males,females)
	MalesFemalesApartMatrix<- (varin[MalesFemalesApart,MalesFemalesApart])

	MalesFemalesApartMatrix[col(MalesFemalesApartMatrix) == row(MalesFemalesApartMatrix) | upper.tri(MalesFemalesApartMatrix)] <- NA
	varout <- c(MalesFemalesApartMatrix, recursive=TRUE)
	varout<- varout[!is.na(varout)]
}
# run function for all variables for null model
DepListNull<- MakeNullModel(DependentTemp, DepListNull)
Indep1ListNull<- MakeNullModel(IndependentTemp1, Indep1ListNull)
Indep2ListNull<- MakeNullModel(IndependentTemp2, Indep2ListNull)
Indep3ListNull<- MakeNullModel(IndependentTemp3, Indep3ListNull)

### get AIC for null model, here we should also define the link function because the response variable is a binomial variable. 
AICNull<- AIC(glm(DepListNull~ Indep1ListNull+Indep2ListNull+Indep3ListNull 
	, 	family=binomial(link="logit")
))
NullModel<- glm(DepListNull~ Indep1ListNull+Indep2ListNull+Indep3ListNull, family=binomial(link="logit"))
TValuesNullModel<- data.frame(summary(NullModel)$coefficient[,3])



scramblingPairs<- function(variable, varout)
{
	males<- c(1,3,9,10,11,12,13,17,19,21,22,27,28,29,31,34,35,37,38,39,41,42,44)
	females<- c(2,4,5,6,7,8,14,15,16,18,20,23,24,25,26,30,32,33,36,40,43)
	mtemp <- variable
	nameVector<- names(mtemp)
	MalesFemalesApart<- c(males,females)
	MalesFemalesApartMatrix<- (mtemp[MalesFemalesApart,MalesFemalesApart])

	NewOrderMales<- sample(c(1:23), 23, replace=FALSE)
	NewOrderFemales<- sample(c(24:44), 21, replace=FALSE)

	NewMatrix <- MalesFemalesApartMatrix[c(NewOrderMales,NewOrderFemales),c(NewOrderMales,NewOrderFemales)]                                    ### these are the lines that do the scrambling
	#names(NewMatrix)<- nameVector[NewOrder]                                  ### place the correct column names
	NewMatrix[col(NewMatrix) == row(NewMatrix) | upper.tri(NewMatrix)] <- NA ### set the diagonal of the matrix and the duplicates to NA(upper triangle)
	ListM <- c(NewMatrix, recursive=TRUE)                                    ### convert the matrix to list
	varout<- ListM[!is.na(ListM)]								 ### remove the NA from the list, these are always at the same location
}


# set initial number of runs
threshold<- 1000				### initial number of runs
counter <- array(NA,threshold)	### to count
counter2 <- array(NA,threshold)
for(i in 1:threshold)
	{
	#tempind1<- scramblingPairs(IndependentTemp1,tempind1)  ### scramble the variable
	#tempind2<- scramblingPairs(IndependentTemp2,tempind2)  ### scramble the variable
	tempind3<- scramblingPairs(IndependentTemp3,tempind3)  ### scramble the variable
ScramModel<- glm(DepListNull~ tempind3+ Indep1ListNull +Indep2ListNull , family=binomial(link="logit"))
AICtemp<- AIC(ScramModel) 
TValuesScramModel<- data.frame(summary(ScramModel)$coefficient[,3])



counter[i]<- (AICtemp<AICNull) 	### add to the counter vector whether the new AIC is lower than the old
counter2[i]<- (abs(TValuesScramModel[2,])>=abs(TValuesNullModel[4,])) # 2 for tempind1, 3 for tempind2, 4 for tempind3 
par(mar=c(0,0,1,0))
plot(0,0, xlim=c(-1,-2), ylim=c(-1,-2), xaxt="n", bty="n", yaxt="n", main=i, xlab="", ylab="") # built in counter on screen, to see how far it is.
}
pvalue_indep<- sum(counter, na.rm=TRUE)/threshold          # line for p-value based on AIC
pvalue_indep                                               # line for p-value based on t-value
pvalue_indep2<- sum(counter2, na.rm=TRUE)/threshold        # p-value AIC
pvalue_indep2                                              # p-value based on t-values

