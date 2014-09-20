ParsModel <- function(df, test = "Chisq", family = gaussian) {
## GLM deletion testing function to get the most parsimonous model, and output results
## data frame format: dependent variable in the first column and independent variables in the other columns
## test is "Chisq" or "F"

stopdeleting <- FALSE
dfBase <- df
justOne <- FALSE
if (test == "Chisq") {indat <- 5}
if (test == "F") {indat <- 6}

## check case of just one independent variable and return glm if true
if (ncol(dfBase)==2) {
        model <- glm(dfBase[, 1] ~., data = subset(dfBase, select = 2), family = family)
        print(summary(model))
        return(model)
}
        
while (stopdeleting == FALSE) {
                
                # number of variables
                numVars <- ncol(dfBase)
                
                # model
                
                model1 <- glm(dfBase[, 1] ~., data = subset(dfBase, select = 2:numVars), family = family)
                summary(model1)
                model1ANOVA <- anova(model1, test = test)
                
                # Remove independent variable with highest p value
                del <- which.max(model1ANOVA[, indat])
                dfNew <- dfBase[, -del]
                
                # get new model
                model2 <- glm(dfBase[, 1] ~., data = subset(dfNew, select = 2:(numVars-1)), family = family)
                #summary(model2)
                
                # ANOVA between models
                modelCompare <- anova(model1, model2, test = test)
                
                # check stop criterion
                if (modelCompare[2, indat] <= 0.05) {
                        stopdeleting <- TRUE
                        bestmodel <- dfBase
                        ## check if just model with one variable left        
                } else if ((numVars-1)==2){
                        stopdeleting <- TRUE
                        bestmodel <- dfNew
                        justOne <- TRUE
                } else {
                        dfBase <- dfNew
                }
                
        }
        ## if model contains just one or no significant variable
        if (justOne == TRUE) {
                print(summary(model2))
                return(model2)
        } else {
                print(summary(model1))
                return(model1)
        }
}
