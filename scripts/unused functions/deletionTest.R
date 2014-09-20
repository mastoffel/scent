ParsModel <- function(df, test) 
## GLM deletion testing function to get the most parsimonous model, and output results
## data frame consists of dependent variable in the first column and independent variables in the other columns
## test is "Chisq" or "F"

stopdeleting = FALSE
dfBase <- df
while (stopdeleting == FALSE) {
# number of variables
numVars <- ncol(df)

# model
model1 <- glm(dfBase[, 1] ~., data = dfBase[, 2:numVars])
summary(model1)
model1ANOVA <- anova(model1, test = test)

# Remove independent variable with highest p value
del <- which.max(model1ANOVA[,5])
dfNew <- dfBase[, -del]

# get new model
model2 <- glm(df[, 1] ~., data = dfNew[, 2:(numVars-1)])
summary(model2)

# ANOVA between models
modelCompare <- anova(model1, model2, test = test)

# check stop criterion
        if (modelCompare[,5] <= 0.05) {
                stopdeleting == T
                bestmodel <- dfBase
                
        } else {
                dfBase <- dfNew
        }


}
