RFapply <- function(data, repeats, trees, mtry){
#data must have a first col named "Type" including the response and further cols containing the predictors
  require(caret)

  learning <- createDataPartition(data$Type, p = .8,
                          list = FALSE,
                          times = 1)

  myRFTrain <- data[ learning ,  ]
  myRFTest  <- data[-learning ,  ]

  testing<-as.integer(row.names(myRFTest))

  # Random Forest - Caret - Fit Model #

  rfControl <- trainControl(method = "repeatedcv",
                            number = 10, repeats = repeats,
                            classProbs = TRUE,
                            allowParallel = TRUE,
                            selectionFunction = "oneSE",
                            returnResamp = "final")

  rfGrid <- expand.grid(mtry = mtry)

  rfFit <- train(Type ~ ., data = myRFTrain,
                 method = "rf",
                 importance = TRUE, ntree=trees,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Kappa", maximize = TRUE)
  
  rfPred <- predict.train(rfFit, myRFTest[,-1], type = "raw")

  list(fit = rfFit,
       pred = predict.train(rfFit, myRFTest[,-1], type = "raw"),
       confusion = confusionMatrix(rfPred, myRFTest$Type),
       varImp = varImp(rfFit, scale = FALSE))
}
