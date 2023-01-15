##################################################################
# Name: Brain Tumor Comparison
# Author: Junier Marrero Gutiérrez
# Descritpion: Here, we are grouping all functions that use in the
# project. Also, we explained specifications and places where were
# are used
##################################################################

########################################################
# Library

if(!require("pacman")){install.packages("pacman"); library("pacman")}

librerias <- c("devtools", "dplyr","openxlsx","spectacles", "ggplot2","DT",
               "mixOmics","MASS", "pROC", "caret", "factoextra","tidyverse",
               "doParallel")

p_load(char = librerias)





############################################################
### Processing Sample Number

cbind.fill <- function(...){
              nm <- list(...)
              nm <- lapply(nm, as.matrix)
              n <- max(sapply(nm, nrow))
              do.call(cbind, lapply(nm, function (x)
                rbind(x, matrix(,n-nrow(x), ncol(x)))))
}



############################################################
#### Sumary of Data

summary.iqr <- function(x){
  
  # x: data, onde a primeira coluna é ID e a segunda Label 
  
  aux1 <- as.data.frame(table(x$Label))
  
  aux4 <- NULL
  
  for (i in 3:ncol(x)) {
    
    min <- format(quantile(x[,i])[[1]], scientific = TRUE, digits = 2)
    Q1 <-  format(quantile(x[,i])[[2]], scientific = TRUE, digits = 2)
    median <- format(quantile(x[,i])[[3]], scientific = TRUE, digits = 2)
    Q2 <- format(quantile(x[,i])[[4]], scientific = TRUE, digits = 2)
    max <- format(quantile(x[,i])[[5]], scientific = TRUE, digits = 2)
    
    aux3 <- rbind(min,Q1, median,Q2, max)
    
    aux4 <- cbind(aux4, aux3)
    
  }
  
  aux4 <- data.frame(aux4)
  names(aux4) <- names(x)[3:ncol(x)]
  
  aux4 <- data.frame(IQR= rownames(aux4), aux4)
  rownames(aux4) <- NULL
  
  aux5 <- cbind.fill(aux1, aux4) 
  aux5 <- data.frame(aux5)
  names(aux5)[c(1,2)] <- c("Label","No")

return(aux5)
  
}


############################################################
#### Substitution of NA

substitution.na <- function(m) {
  
  x <- m
  
  c <- NULL  
  
  for (i in 3:ncol(x)){
    
    j <- x[,c(1,2,i)]
    
    if (length(which(is.na(j[,3])))!=0) {
      
      y <- j[-c(which(is.na(j[,3]))),]
      y <- min(y[,3]/5)
      j[,3][which(is.na(j[,3]))] <- y
      
    }else{
      
      j[,3]<- j[,3]
      
    }
    
    
    c <- cbind(c, j[,3])
    
  }  
  
  c <- data.frame(c)
  names(c) <- names(x)[3:ncol(x)]
  c <- data.frame(x[,c(1,2)],c)
  
  return(c)
  
}





############################################################
# Matrix de confusão de sPLS-DA 2x2
# To summarize the result the optimization threshold
Metric.Summary <- function(x, control, cases) {

  # x é o resultado do Best Prediction function
  # Label1 é a primiera classe dos dados
  # Label2 é a segunda classe dos dados

  Sens <- round(x$byClass["Sensitivity"], digits = 2)*100
  Sepc <- round(x$byClass["Specificity"], digits = 2)*100
  PPV  <- round(x$byClass["Pos Pred Value"], digits = 2)*100
  NPV  <- round(x$byClass["Neg Pred Value"], digits = 2)*100
  Acc  <- round(x$overall["Accuracy"], digits = 2)*100
  Acc.L<-round(x$overall["AccuracyLower"], digits = 2)
  Acc.U<-round(x$overall["AccuracyUpper"], digits = 2)

  mat <- cbind(Sens,Sepc, PPV,NPV, Acc, Acc.L,Acc.U)
  mat <- data.frame(mat)
  mat <- unite(data=mat,
               col=CI,
               sep="-",
               Acc.L,Acc.U)

  mat <- rbind(mat,mat)
  rownames(mat) <- NULL
  mat[1,] <- gsub("[0-9]","", mat[1,])
  mat[1,] <- gsub(".","", mat[1,])

  conf <- cbind(x$table)
  conf <- data.frame(conf)
  conf <- data.frame(Prediction = rownames(conf), 
                     conf, 
                     stringsAsFactors = FALSE)
 
  rownames(conf) <- NULL

  final <- data.frame (conf, mat)
  names(final)[4:8] <- paste0(names(final)[4:8], " (%)")
  names(final)[9] <- paste0(names(final)[9], " (95%)" )

  names(final)[1] <- ""
  final[1,1] <- cases
  names(final)[2] <- cases
  final[2,1] <- control
  names(final)[3] <- control

  lasso <- final

  #return (kable(lasso,align = rep("c", 5)))
  return(lasso)

}


############################################################
### Curvas ROC em função da intensidade dos espectros
### médios.

rocky.metabolites <- function(train.data,
                              test.data,
                              control, 
                              cases){
  
  # train.data: Label and metabolites in columns
  # test.data: Label and metabolites in columns
  # Control: representa os casos negativos, o que não se deseja detetar
  # Cases: representa os casos positivos, o que se deseja detetar
  # cutoff: o cutoof depende da distribuição dos valores
  # 0 e -1 podem ser valores recomendados

# train.data
  aux1 <- train.data[,c(1,2)] 

  ID <- names(aux1)[2]
  
  aux1 <- subset(aux1, Label== cases | Label == control)
  
  aux1[,1] <- ifelse(aux1[,1] == cases, "X0","X1")
  
  aux1[,1] <- as.character(aux1[,1])
  
  roc.ii <- roc(aux1[,1],
                aux1[,2],
                ci.thresholds=TRUE,
                levels=c("X1", "X0"))
  
  # plot.roc <-  plot.roc(roc.ii, 
  #                  print.auc=TRUE,
  #                  print.thres= TRUE, 
  #                  auc.polygon=TRUE,
  #                  plot = NULL)
  
  aux2.orig <- coords(roc.ii, 
                      best.method = "youden", 
                      x="best", 
                      transpose = TRUE)

  p1 <- ifelse(roc.ii$predictor <= aux2.orig[[1]],
                            "X0","X1") 
  p1 <- table(p1,aux1$Label)

  p2 <- ifelse(roc.ii$predictor >= aux2.orig[[1]],
               "X0","X1") 
  p2 <- table(p2,aux1$Label)
  

if (p1[1] > p2[1]){
  
  aux1$pred <- ifelse(roc.ii$predictor <= aux2.orig[[1]],
                      "X0","X1") 
}else{
  
  aux1$pred <- ifelse(roc.ii$predictor >= aux2.orig[[1]],
                      "X0","X1")
  
}  
 
  
f <- confusionMatrix(as.factor(aux1$pred),
                       as.factor(aux1$Label),
                       positive = 'X0')
  
f <- Metric.Summary(f, control, cases)


# Test Data

s <- test.data[,c(1,2)]
 
s <-  subset(s, Label== cases | Label == control)
  
s[,1] <- ifelse(s[,1] == cases, "X0","X1")
  
s$Label <- as.character(s$Label)
  
  
if (p1[1] > p2[1]){
  
  s$pred <- ifelse(s[,2] <= aux2.orig[[1]],
                      "X0","X1") 
}else{
  
  s$pred <- ifelse(s[,2] >= aux2.orig[[1]],
                      "X0","X1")
  
}  

g <- confusionMatrix(as.factor(s$pred),
                     as.factor(s$Label),
                     positive = "X0")

g <- Metric.Summary(g, control, cases)


r <- list(mz=ID, 
            ROC=roc.ii, 
            Cutoff= aux2.orig[[1]] ,
            #ROC.plot = plot.roc,
            MatrixConfusion.TrainData= f,
            MatrixConfusion.TestData= g)
  

  return(r)
    
}




