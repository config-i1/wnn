#' Weighted Nearest Neighbours
#'
#' Function prepares the data and produces weights, fitted and residuals for the
#' combined values based on Weighted Nearest Neighbours
#'
#' This function measures distance from each observation to the provided data,
#' calculates Gower's coefficient based on the provided explanatory variables
#' and then returns a weighted average based on it for each observation. This
#' is a non-parametric method for calculating averages based on a set of
#' explanatory variables
#'
#' @author Carlos Eduardo Rodriguez Calderon
#' @author Ivan Svetunkov
#' @keywords mdoels nonparametric  
#'
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.
#' @param data a data frame or a matrix, containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action	a function which indicates what should happen when the
#' data contain NAs. The default is set by the na.action setting of
#' \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#' factory-fresh default is \link[stats]{na.omit}. Another possible value
#' is NULL, no action. Value \link[stats]{na.exclude} can be useful.
#'
#' @return Function returns \code{model} - the final model of the class
#' "wnn", which contains:
#' \itemize{
#' \item fitted - fitted values,
#' \item residuals - residuals of the model,
#' \item call - how the model was called,
#' \item data - data used for the model construction,
#' \item weights - the matrix of weights for each observation. The columns
#' contain weights between the current obsrevation and all the others in sample.
#' }
#'
#' @seealso \code{\link[greybox]{alm}, \link[greybox]{lmCombine}}
#'
#' @examples
#'
#' ### An example with mtcars data and factors
#' mtcars2 <- within(mtcars, {
#'    vs <- factor(vs, labels = c("V", "S"))
#'    am <- factor(am, labels = c("automatic", "manual"))
#'    cyl  <- factor(cyl)
#'    gear <- factor(gear)
#'    carb <- factor(carb)
#' })
#' # The standard model with Log Normal distribution
#' ourModel <- wnn(mpg~., mtcars2[1:30,])
#'
#' # Produce predictions with the one sided interval (upper bound)
#' predict(ourModel, mtcars2[-c(1:30),])
#' plot(predict(ourModel, mtcars2[-c(1:30),]))
#'
#' @importFrom stats formula model.frame nobs residuals
#' @export wnn
wnn <- function(formula, data, subset, na.action){

    # Record the call of the function
    cl <- match.call();

    #### Form the necessary matrices ####
    # Call similar to lm in order to form appropriate data.frame
    mf <- match.call(expand.dots = FALSE);
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L);
    mf <- mf[c(1L, m)];
    mf$drop.unused.levels <- TRUE;
    mf[[1L]] <- quote(stats::model.frame);


    # If data is provided explicitly, check it
    if(exists("data",inherits=FALSE,mode="numeric") || exists("data",inherits=FALSE,mode="list")){
        if(!is.data.frame(data)){
            data <- as.data.frame(data);
        }
        else{
            dataOrders <- unlist(lapply(data,is.ordered));
            # If there is an ordered factor, remove the bloody ordering!
            if(any(dataOrders)){
                data[dataOrders] <- lapply(data[dataOrders],function(x) factor(x, levels=levels(x), ordered=FALSE));
            }
        }
        mf$data <- data;

        # If there are NaN values, remove the respective observations
        if(any(sapply(mf$data,is.nan))){
            warning("There are NaN values in the data. This might cause problems. Removing these observations.", call.=FALSE);
            NonNaNValues <- !apply(sapply(mf$data,is.nan),1,any);
            # If subset was not provided, change it
            if(is.null(mf$subset)){
                mf$subset <- NonNaNValues
            }
            else{
                mf$subset <- NonNaNValues & mf$subset;
            }
            dataContainsNaNs <- TRUE;
        }
        else{
            dataContainsNaNs <- FALSE;
        }
    }
    else{
        dataContainsNaNs <- FALSE;
    }

    # Prepare the data
    dataWork <- eval(mf, parent.frame());

    #### Produce fitted values and residuals ####
    # Record the dimensions of the new data and the original one
    nRows <- nrow(dataWork);
    nCols <- ncol(dataWork)-1;

    ## Prepare stuff for calculating weights
    # Define, where the categorical variables are used
    factors <- sapply(dataWork[,-1,drop=FALSE], is.factor);
    # Setup ranges of numerical variables
    ranges <- vector("numeric",nCols);
    rangesMatrix <- matrix(0,2,nCols);
    rangesMatrix[,!factors] <- apply(dataWork[,-1,drop=FALSE][,!factors,drop=FALSE], 2, range);
    ranges[] <- rangesMatrix[2,] - rangesMatrix[1,];

    ## Prepare the matrix of distances:
    # nrows - the nrows of the new data
    # ncols - the nrows of the original data
    distanceMatrix <- matrix(0,nRows,nRows);
    for(i in 1:nRows){
        for(j in 1:nRows){
            for(k in 1:nCols){
                distanceMatrix[i,j] <- distanceMatrix[i,j] + gower(dataWork[,-1,drop=FALSE][i,k],
                                                                   dataWork[,-1,drop=FALSE][j,k], ranges[k]);
            }
            # Produce weights based on the distances
            distanceMatrix[i,j] <- 1-distanceMatrix[i,j] / nCols;
        }
    }
    # Calculate weights based on the distance matrix
    weights <- distanceMatrix / apply(distanceMatrix,1,sum);
    # Fitted values
    yFitted <- weights %*% dataWork[,1];
    # Residuals
    errors <- dataWork[,1] - yFitted;

    return(structure(list(data=dataWork,call=cl, weights=weights,
                          fitted=yFitted, residuals=errors),class=c("wnn")));
}

#' @importFrom greybox actuals
#' @export
actuals.wnn <- function(object, ...){
    return(object$data[,1]);
}

#' @export
print.wnn <- function(x, ...){
    cat("Call:\n");
    print(x$call);
}

#' @importFrom stats nobs
#' @export
nobs.wnn <- function(object, ...){
    return(nrow(object$data));
}

#' @importFrom stats predict
#' @export
predict.wnn <- function(object, newdata=NULL, interval=c("none", "confidence", "prediction"),
                        level=0.95, side=c("both","upper","lower"), ...){
    # interval <- match.arg(interval);
    # side <- match.arg(side);
    #
    # if(side=="upper"){
    #     levelLow <- rep(0,length(level));
    #     levelUp <- level;
    # }
    # else if(side=="lower"){
    #     levelLow <- 1-level;
    #     levelUp <- rep(1,length(level));
    # }
    # else{
    #     levelLow <- (1 - level) / 2;
    #     levelUp <- (1 + level) / 2;
    # }

    # If the new data is not provided, do fitted values procedure
    if(is.null(newdata)){
        newdataExpanded <- object$data;
        newdataProvided <- FALSE;
    }
    else{
        newdataProvided <- TRUE;
        # If this is not a dataframe, make one
        if(!is.data.frame(newdata)){
            if(is.vector(newdata)){
                newdataNames <- names(newdata);
                newdata <- matrix(newdata, nrow=1, dimnames=list(NULL, newdataNames));
            }
            newdata <- as.data.frame(newdata);
        }
        else{
            dataOrders <- unlist(lapply(newdata,is.ordered));
            # If there is an ordered factor, remove the bloody ordering!
            if(any(dataOrders)){
                newdata[dataOrders] <- lapply(newdata[dataOrders],function(x) factor(x, levels=levels(x), ordered=FALSE));
            }
        }

        # Extract the formula and get rid of the response variable
        testFormula <- formula(object);
        testFormula[[2]] <- NULL;

        # Expand the data frame
        newdataExpanded <- model.frame(testFormula, newdata);
    }

    # Split the in-sample data into response and explanatory variables
    y <- object$data[,1];
    responseName <- colnames(object$data)[1];
    xreg <- object$data[,-1];
    xregNames <- colnames(object$data)[-1];
    newdataExpanded <- newdataExpanded[,xregNames];

    # Record the dimensions of the new data and the original one
    nRows <- nrow(newdataExpanded);
    nCols <- ncol(newdataExpanded);
    obsInSample <- nobs(object);

    ## Prepare stuff for calculating weights
    # Define, where the categorical variables are used
    factors <- sapply(xreg, is.factor);
    # Setup ranges of numerical variables
    ranges <- vector("numeric",nCols);
    rangesMatrix <- matrix(0,2,nCols);
    rangesMatrix[,!factors] <- apply(rbind(xreg,newdataExpanded)[,!factors,drop=FALSE], 2, range);
    ranges[] <- rangesMatrix[2,] - rangesMatrix[1,];

    ## Prepare the matrix of distances:
    # nrows - the nrows of the new data
    # ncols - the nrows of the original data
    distanceMatrix <- matrix(0,nRows,obsInSample);
    for(i in 1:nRows){
        for(j in 1:obsInSample){
            for(k in 1:nCols){
                distanceMatrix[i,j] <- distanceMatrix[i,j] + gower(newdataExpanded[i,k], xreg[j,k], ranges[k]);
            }
            # Produce weights based on the distances
            distanceMatrix[i,j] <- 1-distanceMatrix[i,j] / nCols;
        }
    }
    # Calculate weights based on the distance matrix
    weights <- distanceMatrix / apply(distanceMatrix,1,sum);

    # Produce point forecasts
    ourForecast <- weights %*% y;

    # if(interval=="prediction"){
    #     x <- sort(y);
    #
    #     levelLow
    #     levelUp
    #     w <- weights[1,]
    #     w[] <- w[order(y)]
    #
    #     v <- sum(0.95-w/2)
    #     x[floor(v)] + (v-floor(v)) * (x[floor(v)+1] - x[floor(v)])
    # }

    return(structure(list(model=object, mean=ourForecast, newdata=newdata,
                          weights=weights, level=level,
                          variances=NULL, newdataProvided=newdataProvided),
                     class="predict.wnn"));
}

#' @export
print.predict.wnn <- function(x, ...){
    ourMatrix <- as.matrix(x$mean);
    colnames(ourMatrix) <- "Mean";
    if(!is.null(x$lower)){
        ourMatrix <- cbind(ourMatrix, x$lower, x$upper);
        if(is.matrix(x$level)){
            level <- colMeans(x$level)[-1];
        }
        else{
            level <- x$level;
        }
        colnames(ourMatrix)[2:3] <- c(paste0("Lower ",round(level[1],3)*100,"%"),paste0("Upper ",round(level[2],3)*100,"%"));
    }
    print(ourMatrix);
}

#' @importFrom greybox graphmaker
#' @importFrom stats deltat fitted frequency start time ts
#' @export
plot.predict.wnn <- function(x, ...){
    yActuals <- actuals(x$model);
    yStart <- start(yActuals);
    yFrequency <- frequency(yActuals);
    yForecastStart <- time(yActuals)[length(yActuals)]+deltat(yActuals);

    if(!is.null(x$newdata)){
        yName <- all.vars(x$model$call$formula)[1];
        if(any(colnames(x$newdata)==yName)){
            yHoldout <- x$newdata[,yName];
            if(!any(is.na(yHoldout))){
                if(x$newdataProvided){
                    yActuals <- ts(c(yActuals,unlist(yHoldout)), start=yStart, frequency=yFrequency);
                }
                else{
                    yActuals <- ts(unlist(yHoldout), start=yForecastStart, frequency=yFrequency);
                }
            }
        }
    }

    # Change values of fitted and forecast, depending on whethere there was a newdata or not
    if(x$newdataProvided){
        yFitted <- ts(fitted(x$model), start=yStart, frequency=yFrequency);
        yForecast <- ts(x$mean, start=yForecastStart, frequency=yFrequency);
        vline <- TRUE;
    }
    else{
        yForecast <- ts(NA, start=yForecastStart, frequency=yFrequency);
        yFitted <- ts(x$mean, start=yStart, frequency=yFrequency);
        vline <- FALSE;
    }

    graphmakerCall <- list(...);
    graphmakerCall$actuals <- yActuals;
    graphmakerCall$forecast <- yForecast;
    graphmakerCall$fitted <- yFitted;
    graphmakerCall$vline <- vline;

    if(!is.null(x$lower)){
        if(x$newdataProvided){
            yLower <- ts(x$lower, start=yForecastStart, frequency=yFrequency);
            yUpper <- ts(x$upper, start=yForecastStart, frequency=yFrequency);
        }
        else{
            yLower <- ts(x$lower, start=yStart, frequency=yFrequency);
            yUpper <- ts(x$upper, start=yStart, frequency=yFrequency);
        }

        if(is.matrix(x$level)){
            level <- x$level[1];
        }
        else{
            level <- x$level;
        }
        graphmakerCall$level <- level;
        graphmakerCall$lower <- yLower;
        graphmakerCall$upper <- yUpper;

        if((any(is.infinite(yLower)) & any(is.infinite(yUpper))) | (any(is.na(yLower)) & any(is.na(yUpper)))){
            graphmakerCall$lower[is.infinite(yLower) | is.na(yLower)] <- 0;
            graphmakerCall$upper[is.infinite(yUpper) | is.na(yUpper)] <- 0;
        }
        else if(any(is.infinite(yLower)) | any(is.na(yLower))){
            graphmakerCall$lower[is.infinite(yLower) | is.na(yLower)] <- 0;
        }
        else if(any(is.infinite(yUpper)) | any(is.na(yUpper))){
            graphmakerCall$upper <- NA;
        }
    }

    do.call(graphmaker,graphmakerCall);
}

# @export
gower <- function(x, y, range=NULL){
    # The function compares two scalars, given the range and returns a value between 0 and 1
    if(is.factor(x) && is.factor(y)){
        return((x!=y)*1);
    }
    else{
        return(abs(x-y)/range);
    }
}

#' @importFrom stats sigma
#' @export
sigma.wnn <- function(object, ...){
    return(sqrt(sum((object$weights %*% residuals(object)-mean(residuals(object)))^2)));
}
