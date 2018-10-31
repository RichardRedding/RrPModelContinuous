wrapFunctionMacro <- function(labelArgPosition = 1, covArgPosition = integer(0), func = "f", collection = "Table[Double]", collectionDerefMethod = "apply"){
  paramArgDefinition <- paste(rep("Double", length(labelArgPosition)), collapse = ", ")
  argDefinition <- if(length(labelArgPosition) > 0){paste0('(', paste(paste0("param: (", paramArgDefinition, ")"), paste0('t: ', collection), sep = ", "), ')')} else {paste0('(t: ', collection, ')')}
  bodyDefinition <- character(sum(!is.na(labelArgPosition)) + sum(!is.na(covArgPosition)))
  bodyDefinition[labelArgPosition[!is.na(labelArgPosition)]] <- paste0("param._", seq_len(length(labelArgPosition)))[!is.na(labelArgPosition)]
  bodyDefinition[covArgPosition[!is.na(covArgPosition)]] <- paste0("t.", collectionDerefMethod, "(", seq_len(length(covArgPosition))-1, ")")[!is.na(covArgPosition)]
  bodyDefinition <- paste0("(", paste(bodyDefinition, collapse = ", "), ")")
  body <- paste0(func, bodyDefinition)
  if(length(func) > 1) {body <- paste0("Vector(", paste0(body, collapse = ", "), ")")}
  paste0(argDefinition, " => ", body)
}

as.NumericMatrix <- function(x){
  mat <- as.matrix(x)
  structure(as.numeric(mat), dim = dim(mat))
}

returnXtraForest <- function(covariates, staticTerm, observedTerm, numObservations,
                             numTrees = 30L, numSampledFeatures = ncol(covariates), 
                             greatestLeafSize = sqrt(nrow(covariates)), 
                             interval = rep(TRUE, ncol(covariates))){
  numParameters <- 2
  
  numCovariates <- ncol(covariates)
  
  estimator <- list(func = "mle", 
                    labelArgPosition = integer(0), 
                    covArgPosition = c(rep(NA, numCovariates), 1, 2, 3),
                    collection = "Table[Double]",
                    collectionDerefMethod = "col")
  
  scorer <- list(func = "maximumLogLikelihood", 
                 labelArgPosition = integer(0), 
                 covArgPosition = c(rep(NA, numCovariates), 1, 2, 3),
                 collection = "Table[Double]",
                 collectionDerefMethod = "col")
  
  predictor <- list(func = c("meanFutureNumObservations", "varianceFutureNumObservations"), 
                    labelArgPosition = seq_len(numParameters), 
                    covArgPosition = c(rep(NA, numCovariates), 3, 4),
                    collection = "IndexedSeq[Double]",
                    collectionDerefMethod = "apply")
  
  accumulator <- '(y1: Vector[Double], y2: Vector[Double], i: Int) => {val res1 = (y1(0) * i + y2(0)) / (i + 1); val res2 = ((y1(1) + y1(0) * y1(0)) * i + (y2(1) + y2(0) * y2(0))) / (i + 1) - res1 * res1; Vector(res1, res2)}'
  
  data <- as.NumericMatrix(cbind(covariates, staticTerm, observedTerm, numObservations))
  
  res <-
    s(data = data, numTrees = as.integer(numTrees), numSampledFeatures = as.integer(numSampledFeatures), 
      greatestLeafSize = as.integer(greatestLeafSize), interval = interval) ^ 
      paste('ExtraTrees(Table[Double](data, true)', 
            do.call(wrapFunctionMacro, estimator), 
            do.call(wrapFunctionMacro, scorer),
            do.call(wrapFunctionMacro, predictor),
            accumulator,
            'numTrees, numSampledFeatures, greatestLeafSize, interval)',
            sep = ", ")
  
  class(res) <- append(class(res), "XtraForest")
  
  res
}

predict.XtraForest <- function(forest, covariates, staticTerm, predictionTerm){
  data <- as.NumericMatrix(cbind(covariates, staticTerm, predictionTerm))
  structure(s(data = as.matrix(data), forest = forest) * 'Table(data, true).map(forest predict _).unpackColMajor.toArray', 
            dim = c(length(predictionTerm), 2),
            dimnames = list(list(), c("predictedMeanNumEvents", "predictedVarNumEvents")))
}