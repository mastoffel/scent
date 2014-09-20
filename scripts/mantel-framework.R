library(ecodist)
library(car)
library(AICcmodavg)
library(MuMIn)
library(gplots)
library(Hmisc)
library(sna)
source("ScentResults.R")

res <- ScentResults("mums","fa",1)

# subset results list output
relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]
factors <- res[[8]]

## all relevant similarity matrices
relatedness <- as.matrix(relatedness)
F1 <- as.matrix(relate.list[[1]])
F2 <- as.matrix(relate.list[[2]])
F3 <- as.matrix(relate.list[[3]])
F4 <- as.matrix(relate.list[[4]])

## make matrices symmetric
for (c in 1:(ncol(F1))) {
        for (r in 1:(nrow(F1))) {
                if (!is.na(F1[c,r])) {
                relatedness[r,c] <- relatedness[c,r]
                F1[r,c] <- F1[c,r]
                F2[r,c] <- F2[c,r]
                F3[r,c] <- F3[c,r]
                F4[r,c] <- F4[c,r]
                }
        }
}

## fill diagonals
diag(relatedness) <- 0
diag(F1) <- 0
diag(F2) <- 0
diag(F3) <- 0
diag(F4) <- 0

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
DepListNull<- MakeNullModel(relatedness, DepListNull)
Indep1ListNull<- MakeNullModel(F1, Indep1ListNull)
Indep2ListNull<- MakeNullModel(F2, Indep2ListNull)
Indep3ListNull<- MakeNullModel(F3, Indep3ListNull)
Indep4ListNull<- MakeNullModel(F4, Indep4ListNull)

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
        tempind1<- scrambling(F1,tempind1)
        #tempind2<- scrambling(F2,tempind2)
        #tempind3<- scrambling(F3,tempind3)  ### scramble the variable
        #tempind4<- scrambling(F4,tempind4)
        
        # get the AIC for the new model:
        ScramModel<- glm(DepListNull~ tempind1 + Indep2ListNull + Indep3ListNull  + Indep4ListNull)
        AICtemp<- AIC(ScramModel) 
        TValuesScramModel<- data.frame(summary(ScramModel)$coefficient[,3])
        counter[i]<- (AICtemp<AICnull)  ### add to the counter vector whether the new AIC is lower than the old
        counter2[i]<- (abs(TValuesScramModel[2,])>=abs(TValuesNullModel[2,]))  # 2 for tempind1, 3 for tempind2, 4 for tempind3, 5 for tempind4	
        par(mar=c(0,0,1,0))
        plot(0,0, xlim=c(-1,-2), ylim=c(-1,-2), xaxt="n", bty="n", yaxt="n", main=i, xlab="", ylab="") # built in counter on screen, to see how far it is.
}
pvalue_indep<- sum(counter, na.rm=TRUE)/threshold    # line for p-value based on AIC
pvalue_indep2<- sum(counter2, na.rm=TRUE)/threshold  # line for p-value based on t-value
pvalue_indep                                         # p-value AIC
pvalue_indep2                                        # p=value based on t-values

## 10000, F1 scrambled
#> pvalue_indep                                         # p-value AIC
#[1] 0.0069
#> pvalue_indep2                                        # p=value based on t-values
#[1] 0.0069

## get estimates and SE

model <- glm()