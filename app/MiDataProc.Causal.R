library(randomForest)
library(caret)
library(tidyverse)
library(ggplot2)
library(stringr)
library(edarf)
library(data.table)
library(dplyr)
library(reshape2)
library(aVirtualTwins)
library(grf)
library(MiVT)
library(rpart)
library(rpart.plot)

source("Source/MiDataProc.ML.Models.R")

double.sample.treatment.pred <- function(Feature, Response, Treatment, n.tree = 10000){
  s.time <- proc.time()
  cf.fit <- causal_forest(X = Feature, Y = Response, W = Treatment, num.trees = n.tree, 
                          sample.weights = NULL, clusters = NULL, equalize.cluster.weights = FALSE,
                          sample.fraction = 0.5, min.node.size = 5, 
                          honesty = TRUE, honesty.fraction = 0.5, honesty.prune.leaves = TRUE,
                          alpha = 0.05, imbalance.penalty = 0, stabilize.splits = TRUE, 
                          ci.group.size = 2, tune.parameters = "all",
                          tune.num.trees = n.tree/10, tune.num.reps = n.tree/10, tune.num.draws = n.tree/10, 
                          compute.oob.predictions = TRUE, num.threads = NULL, 
                          seed = 521)
  e.time <- proc.time()
  e.time - s.time
  
  pred <- predict(cf.fit, estimate.variance = TRUE)
  Treat.Effect <- pred$predictions 
  
  return(list(fit = cf.fit, Treat.Effect = pred$predictions))
}

propensity.treatment.pred <- function(Feature, Response, Treatment, Covariate, n.tree = 10000){
  
  reg.W.fit <- regression_forest(X = Feature, Y = Treatment, num.trees = n.tree, seed = 521)
  W.hat <- predict(reg.W.fit)$predictions
  
  reg.Y.fit <- regression_forest(X = Feature, Y = Response, num.trees = n.tree, seed = 521)
  Y.hat <- predict(reg.Y.fit)$predictions
  
  s.time <- proc.time()
  cf.fit <- causal_forest(X = Feature, Y = Response, W = Treatment, num.trees = n.tree, 
                          Y.hat = Y.hat, W.hat = W.hat, 
                          sample.weights = NULL, clusters = NULL, equalize.cluster.weights = FALSE,
                          sample.fraction = 0.5, min.node.size = 5, 
                          honesty = TRUE, honesty.fraction = 0.5, honesty.prune.leaves = TRUE,
                          alpha = 0.05, imbalance.penalty = 0, stabilize.splits = TRUE, 
                          ci.group.size = 2, tune.parameters = "all",
                          tune.num.trees = n.tree/10, tune.num.reps = n.tree/10, tune.num.draws = n.tree/10, 
                          compute.oob.predictions = TRUE, num.threads = NULL, 
                          seed = 521)
  e.time <- proc.time()
  e.time - s.time
  
  pred <- predict(cf.fit, estimate.variance = TRUE)
  Treat.Effect <- pred$predictions 
  
  return(list(fit = cf.fit, Treat.Effect = pred$predictions))
}


subgroup.id <- function(Treat.Effect, taxa.out, type, level.name){
  # Taxa <- scale(taxa.out[[type]][[level.name]])
  Taxa <- taxa.out[[type]][[level.name]]
  
  taxa.names <- character()
  for (i in 1:ncol(Taxa)) {
    taxon.name <- colnames(Taxa)[i]
    taxon.name <- unlist(strsplit(taxon.name, ";"))
    taxa.names[i] <- taxon.name[length(taxon.name)]
  }
  
  colnames(Taxa) <- paste(substr(level.name, 1, 1), 1:ncol(Taxa), sep = "")
  
  dat <- as.data.frame(cbind(Treat.Effect, Taxa))
  
  ### Regression Tree 
  
  control <- rpart.control(minsplit = 20, minbucket = round(20/3), cp = 1e-5, xval = nrow(dat)) # Do LOOCV cross-validation (xval = nrow(dat))
  dt.fit <- rpart(Treat.Effect ~ ., data = dat, method = "anova", control = control)
  
  if (identical(class(dt.fit$cptable[-1,]), c("matrix", "array"))) {
    cp <- dt.fit$cptable[-1,][which.min(dt.fit$cptable[-1,][,"xerror"]), "CP"]
  } else {
    cp <- dt.fit$cptable[-1,][[1]]
  }
  best.dt.fit <- prune(dt.fit, cp = cp)
  return(list(Taxa = Taxa, taxa.names = taxa.names, dat = dat, best.dt.fit = best.dt.fit, Treat.Effect = Treat.Effect))
}

subgroup.id.vis <- function(best.dt.fit){
  rpart.plot(best.dt.fit,type = 4, extra = 100, under = TRUE, fallen.leaves = TRUE, digits = 3, faclen = 3, cex = 0.85,
             tweak = 1.45, clip.right.labs = FALSE, box.palette = "Orange")
}

bort.func <- function(subgroup.id.result, level.name, n.tree = 10000){
  best.dt.fit <- subgroup.id.result$best.dt.fit
  taxa.names <- subgroup.id.result$taxa.names
  Treat.Effect <- subgroup.id.result$Treat.Effect
  ind <- grep(substr(level.name, 1, 1), as.character(best.dt.fit$frame$var))
  sel.taxa <- as.character(best.dt.fit$frame$var)[ind]
  taxa.num <- as.numeric(gsub(substr(level.name, 1, 1), "", as.character(best.dt.fit$frame$var)[ind]))
  Sel.Taxa <- t(as.data.frame(taxa.names[taxa.num]))
  rownames(Sel.Taxa) <- "Full names"
  colnames(Sel.Taxa) <- sel.taxa
  BoRT.out <- boot.test(Z.hat = Treat.Effect, best.dt.fit, n.boot = n.tree)
  BoRT.out <- round(BoRT.out, 3)
  colnames(BoRT.out) <- NULL
  out <- list(Sel.Taxa = Sel.Taxa, BoRT.out = BoRT.out)
  return(out)
}

bort.treatment.pred <- function(subgroup.id.result, level.name, n.tree = 10000){
  Taxa <- subgroup.id.result$Taxa
  Treat.Effect <- subgroup.id.result$Treat.Effect
  
  rf.cv <- rfcv(trainx = Taxa, trainy = Treat.Effect, cv.fold = 10, scale = "log", step = 0.8, 
                recursive = FALSE, ntree = n.tree/10)
  
  opt.mtry <- as.numeric(names(which.min(rf.cv$error.cv)))
  
  rf.fit <- randomForest(x = Taxa, y = Treat.Effect, mtry = opt.mtry, importance = TRUE, ntree = n.tree)
  return(rf.fit)
}

cf.imp.df <- function(fit, type){
  if(type == 0){
    imp <- randomForest::importance(fit)
  }
  else if(type == 1){
    imp <- randomForest::importance(fit, type = 1)
  }
  else if(type == 2){
    imp <- randomForest::importance(fit, type = 2)
  }
  if("%IncMSE" %in% colnames(imp)){
    ind <- which(colnames(imp) == "%IncMSE")
    colnames(imp)[ind] <- "IncMSE"
  }
  return(as.data.frame(imp))
}

cf.imp.plot <- function(fit, type, n = 30, subgroup.id.result, data.type, level.name){ # data, 

  imp.df1 <- cf.imp.df(fit, type = 1)
  imp.df2 <- cf.imp.df(fit, type = 2)
  # new.names <- subgroup.id.result$taxa.names
  # rownames(imp.df1) <- new.names
  # rownames(imp.df2) <- new.names
  
  ma <- get.mean.abundance(subgroup.id.result$Taxa, level.name)
  imp.df1 <- data.frame(imp.df1, ma)
  colnames(imp.df1)[2] <- "MeanAbundance"
  imp.df2 <- data.frame(imp.df2, ma)
  colnames(imp.df2)[2] <- "MeanAbundance"
  
  if(data.type == "clr"){
    data.type <- "CLR"
  }
  else if(data.type == "prop"){
    data.type <- "Proportion"
  }
  else if(data.type == "rare.count"){
    data.type <- "Rarefied Count"
  }
  else if(data.type == "arcsin"){
    data.type <- "Arcsine-Root"
  }
  
  names(imp.df1)[1] <- "IncMSE"
  imp.df1 <- imp.df1 %>%
    mutate(names = rownames(imp.df1)) %>%
    arrange(desc(IncMSE)) %>%
    top_n(n, IncMSE)
  
  imp.df1 <- imp.df1 %>%
    ggplot(aes(x=reorder(names, IncMSE), y=IncMSE)) +
    geom_segment( aes(x=reorder(names, IncMSE), xend=reorder(names, IncMSE), y=0, yend=IncMSE), color="grey") +
    # geom_point( color=ifelse(imp.df1$IncMSE > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
    geom_point(aes(color = MeanAbundance), size=5, alpha=0.8) +
    scale_color_gradient(low = "blue", high = "red", name = sprintf("Mean\nAbundance\n(%s)", data.type)) +
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 13)
    ) +
    xlab(element_blank()) +
    ylab("Decrease in Mean Squared Error")
  
  imp.df2 <- imp.df2 %>%
    mutate(names = rownames(imp.df2)) %>%
    arrange(desc(IncNodePurity)) %>%
    top_n(n, IncNodePurity)
  
  imp.df2 <- imp.df2 %>%
    ggplot(aes(x=reorder(names, IncNodePurity), y=IncNodePurity)) +
    geom_segment( aes(x=reorder(names, IncNodePurity), xend=reorder(names, IncNodePurity), y=0, yend=IncNodePurity), color="grey") +
    # geom_point( color=ifelse(imp.df2$IncNodePurity > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
    geom_point(aes(color = MeanAbundance), size=5, alpha=0.8) +
    scale_color_gradient(low = "blue", high = "red", name = sprintf("Mean\nAbundance\n(%s)", data.type)) +
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 13)
    ) +
    xlab(element_blank()) +
    ylab("Decrease in Node Impurity")
  
  if(type == 0){
    return(imp.df1 + imp.df2)
  }
  else if(type == 1){
    return(imp.df1)
  }
  else if(type == 2){
    return(imp.df2)
  }
}

cf.pdp.reg <- function(fit, X, n, data.type){
  rf.importance <- cf.imp.df(fit, type = 1) %>%
    mutate(taxon = rownames(cf.imp.df(fit, type = 1))) %>%
    arrange(-IncMSE) %>%
    rownames
  n <- min(length(rf.importance), n)
  feature <- rf.importance[1:n]
  result <- data.frame()
  
  if(data.type == "clr"){
    type = "CLR"
  }
  else if(data.type == "prop"){
    type = "Proportion"
  }
  else if(data.type == "rare.count"){
    type = "Rarefied Count"
  }
  else if(data.type == "arcsin"){
    type = "Arcsine-Root"
  }
  
  for(taxon.name in feature){
    val <- numeric()
    y_hat <- numeric()
    taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
    
    for(i in 1:length(taxon.val)){
      newX <- X
      newX[,taxon.name] <- rep(taxon.val[i], nrow(newX))
      y_pred <- predict(fit, newX)
      val <- c(val, taxon.val[i])
      y_hat <- c(y_hat, mean(y_pred))
    }
    pred.result <- data.frame(val, y_hat) %>%
      reshape2::melt(id.vars = "val", value.name = "Prediction")
    pred.result$title <- rep(taxon.name, nrow(pred.result))
    result <- rbind(result, pred.result)
  }
  result$title <- factor(result$title, levels = feature)
  
  p <- ggplot(result, aes(val, Prediction)) + 
    geom_line(size = 0.8) +
    theme_light() +
    xlab(type) + 
    ylab("Predicted Value") +
    theme(
      axis.text.x = element_text(size = 7.5),
      axis.text.y = element_text(size = 7.5),
      strip.text = element_text(size=12),
      panel.spacing.x = unit(1, "lines"),
      legend.title = element_blank()
    ) +
    facet_wrap(~ title, scales = "free_x", nrow = 5, dir = "v")
  p
}

cf.dt.used.var <- function(step2.result){
  new.name <- setdiff(step2.result$best.dt.fit$frame$var, "<leaf>")
  var.name.df <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
  
  ori.name <- c()
  for(name in new.name){
    ind <- which(name == var.name.df$sub)
    ori.name <- c(ori.name, var.name.df$ori[ind])
  }
  output <- data.frame(name = ori.name)
  rownames(output) <- new.name
  return(output)
}

rf.used.var <- function(step2.result, step4.result, n = 20){
  new.name <- cf.imp.df(step4.result, type = 1) %>% arrange(desc(IncMSE)) %>% head(n) %>% rownames
  var.name.df <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
  
  ori.name <- c()
  for(name in new.name){
    ind <- which(name == var.name.df$sub)
    ori.name <- c(ori.name, var.name.df$ori[ind])
  }
  output <- data.frame(name = ori.name)
  rownames(output) <- new.name
  return(output)
}
