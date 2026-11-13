## 12 July 2024
## Mia Venter

# RFES (Recursive Feature Discovery with Strikeout) -- a recursive feature selection procedure using Random Forest

# load libs 
library(tidyverse)
library(caret)
library(randomForest)
library(tictoc)

###########################################################################
#USER-INPUT
setwd("C:/Users/miatr/Other projects/Subtype_discovery/data/")

train.data <- readRDS("fd_vst_tcga526_353genes_scaled_train.rds")
test.data <- readRDS("fd_vst_tcga526_353genes_scaled_test.rds")
savename <- "TCGA526_"
strikeout <- 8
classPointer = 1
formula <- "cluster ~ ."
repeat_times <- 5
wd <- "27rfes_output/" # must end in /
mtry <- 1:40
ntree <- 300
cw <- "none" # else the name of the file
###########################################################################

set.seed(123)

# cp rfes.R into wd
x <- paste0("cp rfes.R ", wd)
system(x)

# create readme
readme <- "2025 11 19. TCGA526, using 353 GE genes , scaled. 300 trees"
x <- paste0("echo ", readme, " > " ,wd,"readme.txt")
system(x)

feature_discovery <- function(train.data, test.data, savename, formula, strikeout, repeat_index, ntree){
  
  # # split data into training and test sets
  # trainIndex <- createDataPartition(counts[,classPointer],
  #                                   p = train.test.split,
  #                                   list = FALSE,
  #                                   times = 1)
  # 
  # train.data <- counts[trainIndex, ]
  # test.data <- counts[-trainIndex, ]
  # 
  # # write training and testing data for repeatability
  # saveRDS(train.data,paste0(wd,savename,"_train.data.rds"))
  # saveRDS(test.data,paste0(wd,savename,"_test.data.rds"))
  
  
  # make names to avoid potential syntax errors
  names(train.data) <- make.names(names(train.data))
  names(test.data) <- make.names(names(test.data))
  
  # get y-variable
  yvar <- sapply(strsplit(formula, "\\s+"), "[", 1)
  
  # get y-var for when it's in a dataframe
  yvar.df <- paste0("test.data$",yvar)
  
  # get y-var for when it's in the testdata
  yvar.test <- paste0("test.data$",yvar)
  
  # reformulate the formula for later use 
  formula <- reformulate(".",yvar)
  
  
  # set s (a tracker for strikeout) and a performance metric value to compare against
  s = 0
  previous.metric = 0 #unreasonably small number so that the loop always runs at least a second time 
  
  # create dataframe to store results of feature discovery and an index to keep track of df rows
  df <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(df) <- c("run","num.genes.in.model","accuracy training",
                    "train_class_error","accuracy testing","mtry","nr_nodes","genes.in.model")
  index = 1
  
  # create dataframe to store accuracy (for accuracy vs num genes plot)
  x <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(x) <- c("Accuracy", "Num_genes")
  
  #while s < strikeout, perform recursive feature discovery
  while (s < strikeout){
    
    # set tuning params
    grid_rf <- expand.grid(mtry=mtry)
    train_ctrl <- trainControl(method="cv",
                               number=5, 
                               search = "grid" ,
    )
    
    # create directory to store imp_per_feature and random forest object
    dir.name <- paste0(wd,savename,"_Repeat_",repeat_index)
    if (!dir.exists(dir.name)){
      dir.create(dir.name)
    }
    
    if(cw == "none"){
      print("running RF")
      # run random forest
      print(dim(train.data))
      saveRDS(train.data,"train.data.rds")
      obj <- train(formula,
                   data = train.data,
                   method = "rf", # this will use the randomForest::randomForest function
                   metric = "Accuracy", # which metric should be optimized for 
                   trControl = train_ctrl, 
                   tuneGrid = grid_rf,
                   # options to be passed to randomForest
                   ntree = ntree,
                   keep.forest=T,
                   importance=TRUE,
                   savePredictions = T
      ) 
    }
    else{
      cw <- readRDS(cw)
      
      # run random forest with class weights
      obj <- train(formula,
                   data = train.data,
                   method = "rf", # this will use the randomForest::randomForest function
                   metric = "Accuracy", # which metric should be optimized for 
                   trControl = train_ctrl, 
                   tuneGrid = grid_rf,
                   # options to be passed to randomForest
                   ntree = ntree,
                   keep.forest=T,
                   importance=TRUE,
                   savePredictions = T,
                   classwt = cw
      ) 
      
    }
    
    #save importance per feature in a table
    print("dtermining feature imp")
    imp.per.feature <- varImp(obj, scale=T)
    imp.per.feature <- as.data.frame(imp.per.feature$importance)
    savename.imp.per.feature <- paste0(dir.name,"/",savename,"_importance_per_feature_Run",index,".tsv")
    write.table(imp.per.feature,savename.imp.per.feature,sep="\t",row.names = T)
    
    savename.rf.object <- paste0(dir.name,"/",savename,"RF_object_Run",index,".rds")
    saveRDS(obj,savename.rf.object)
    
    # get genes of interest as genes with importance in at least one class > 1st quartile 
    m <- apply(imp.per.feature,1,max)
    summary <- summary(m)[2]
    ind <- which(m > summary)
    goi <- rownames(imp.per.feature)[ind]
    goi <- append(goi,colnames(train.data)[classPointer])
    
    # get confusion matrix
    print("conf mat")
    
    conf.mat <- as.data.frame(obj$finalModel$confusion)
    
    # get prediction on testdata
    print("prediction on test data")
    
    p <- predict(obj, test.data)    
    test.conf.mat <- data.frame(predictions = p, true.values = test.data[,classPointer])
    test.conf.mat <- as.data.frame(table(test.conf.mat$predictions,test.conf.mat$true.values))
    test.conf.mat <- test.conf.mat %>% pivot_wider(names_from = "Var1",values_from = "Freq")
    test.conf.mat <- column_to_rownames(test.conf.mat,"Var2")
    
    # get accuracy metric for counts and testdata
    bestMtry <- obj$bestTune$mtry
    acc.train <- obj$results$Accuracy[obj$results$mtry == bestMtry]
    acc.test <-  diag(as.matrix(test.conf.mat))
    acc.test <- sum(acc.test)/length(rownames(test.data))
    
    #save all metrics in a dataframe
    df[index,1] <- index # run 
    df[index,2] <- length(colnames(train.data)) -1 # to account for y.var
    
    gene.list <- as.list(colnames(train.data)[-classPointer])
    df[index,3] <- acc.train #performance
    
    ce <- as.data.frame(obj$finalModel$confusion)
    ce$class.error <- round(ce$class.error,2)
    ce <- as.list(ce$class.error)
    df[index,4] <- paste(ce, collapse = ',')
    
    df[index,5] <- acc.test
    
    df[index,6] <- obj$finalModel$tuneValue
    df[index,7] <- obj$finalModel$forest$nrnodes
    df[index,8] <- paste(colnames(train.data)[-classPointer], collapse = ',') #list of genes used in model
    
    # an error occurs if no genes are significant
    
    # there's at least one significant gene
    if(length(goi) > 1){
      # reset data to only include the new genes of interest
      train.data <- as.data.frame(train.data[,colnames(train.data) %in% goi])
      test.data <- as.data.frame(test.data[,colnames(test.data) %in% goi])
      
      # store testing accuracy vs num.genes
      x[index,1] <- acc.test
      x[index,2] <- length(goi)
      
      # decide whether model performance improved
      decision = determine_improvement(acc.test,previous.metric)
      previous.metric = acc.test # replace old metric with new metric
      
      # if model did not improve, add 1 to the strikeout_tracker
      if (decision == F){
        s = s+1
      }
      index = index +1
      print(paste("s: ",s))
    }
    else{
      s=s+1
      print(paste("s: ",s)) 
      next()
    }
  }## end of loop
  
  #save 
  df.savename <- paste0(dir.name,"/",savename,".tsv")
  write.table(df,df.savename, sep = "\t",row.names = F)
  print("accuracy vs numgenes")
  
  # plot accuracy vs num.genes
  png(paste0(dir.name,"/",savename,"_accuracy_vs_numgenes.png"))
  print(ggplot(x, aes(x = Num_genes, y = Accuracy))+
          geom_line()+
          scale_x_reverse()+
          ylim(0,1)+
          labs(x = "Number of genes", y = "Accuracy")+
          theme_minimal()+
          theme(axis.title = element_text(size = 14, face="bold"),
                axis.text = element_text(size = 12, face = "bold", colour="black")))
  dev.off()
  
  return (df)
} ## end of function



determine_improvement <- function(metric, previous.metric){
  print("determining improvement")
  
  if (previous.metric < metric){ # model improves if previous metric is smaller
    return (T)
  }
  else{
    return (F)
  }
}



# repeat the feature discovery process repeat_times number of times to get the best run
# this is necessary because the RF algorithm has some randomness associated with it, so one might get better results
# in one run than the next
repeat_discovery <- function (repeat_times, train.data, test.data, savename, formula, strikeout, ntree) {
  
  # create dataframe to store the best run 
  
  best.run.df <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(best.run.df) <- c("run","num.genes.in.model","accuracy training",
                             "train_class_error","accuracy testing","mtry","nodesize","genes.in.model")
  i=1
  for (i in 1:repeat_times){
    df <- feature_discovery(train.data = train.data,test.data = test.data, savename = savename, formula = formula, 
                            strikeout = strikeout, repeat_index = i, ntree = ntree)
    
    # disregard the first row of the dataframe, since all input features may have the best performance, but 
    # is impractical for our purposes
    df <- df[-1,]
    
    # find all maximum performance values and get the index of the last maximum value 
    # last max value will be the run with fewest genes and best performance
    best.values <- which(df[,5] == max(df[,5]))
    best.values <- tail(best.values, n=1)
    
    
    
    # insert best run from this repeat into best.run.df
    best.run.df[i,] <- df[best.values,]
  } # end of for
  
  best.run.df <- add_column(best.run.df, Repeat = 1:repeat_times, .before = 1)
  
  return(best.run.df)
  
} # end of function

# time
tic("timer")

best.run.df <- repeat_discovery(train.data = train.data, savename = savename, formula = formula, test.data = test.data, strikeout = strikeout,
                                repeat_times = repeat_times, ntree = ntree)
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)
tic.clearlog()

dir.create(paste0(wd,savename,"_BestRun"))
best.run.savename <- paste0(wd,savename,"_BestRun/",savename,"_BestRun.tsv")

best.run.df[length(rownames(best.run.df))+1 ,1] <- paste(log.txt[[1]])

write.table(best.run.df,best.run.savename, sep = "\t",row.names = F)
