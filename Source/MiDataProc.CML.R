double.sample.treatment.pred <- function(Feature, Response, Treatment, n.tree = 10000) {
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

propensity.treatment.pred <- function(Feature, Response, Treatment, Covariate, n.tree = 10000) {
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


subgroup.id <- function(Treat.Effect, taxa.out, type, level.name) {
  Taxa <- taxa.out[[type]][[level.name]]
  taxa.names <- character()
  for (i in 1:ncol(Taxa)) {
    taxon.name <- colnames(Taxa)[i]
    taxon.name <- unlist(strsplit(taxon.name, ";"))
    taxa.names[i] <- taxon.name[length(taxon.name)]
  }
  colnames(Taxa) <- paste(substr(level.name, 1, 1), 1:ncol(Taxa), sep = "")
  
  dat <- as.data.frame(cbind(Treat.Effect, Taxa))
  control <- rpart.control(minsplit = 20, minbucket = round(20/3), cp = 1e-5, xval = nrow(dat)) 
  dt.fit <- rpart(Treat.Effect ~ ., data = dat, method = "anova", control = control)
  
  if (identical(class(dt.fit$cptable[-1,]), c("matrix", "array"))) {
    cp <- dt.fit$cptable[-1,][which.min(dt.fit$cptable[-1,][,"xerror"]), "CP"]
  } else {
    cp <- dt.fit$cptable[-1,][[1]]
  }
  best.dt.fit <- prune(dt.fit, cp = cp)
  return(list(Taxa = Taxa, taxa.names = taxa.names, dat = dat, best.dt.fit = best.dt.fit, Treat.Effect = Treat.Effect))
}

subgroup.id.vis <- function(best.dt.fit) {
  rpart.plot(best.dt.fit,type = 4, extra = 100, under = TRUE, fallen.leaves = TRUE, digits = 3, faclen = 3, cex = 0.85,
             tweak = 1.45, clip.right.labs = FALSE, box.palette = "Orange")
}

bort.func <- function(subgroup.id.result, level.name, n.tree = 10000) {
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

bort.treatment.pred <- function(subgroup.id.result, level.name, n.tree = 10000) {
  Taxa <- subgroup.id.result$Taxa
  Treat.Effect <- subgroup.id.result$Treat.Effect
  rf.cv <- rfcv(trainx = Taxa, trainy = Treat.Effect, cv.fold = 10, scale = "log", step = 0.8, 
                recursive = FALSE, ntree = n.tree/10)
  opt.mtry <- as.numeric(names(which.min(rf.cv$error.cv)))
  rf.fit <- randomForest(x = Taxa, y = Treat.Effect, mtry = opt.mtry, importance = TRUE, ntree = n.tree)
  return(rf.fit)
}

cf.imp.df <- function(fit, type) {
  if(type == 0) {
    imp <- randomForest::importance(fit)
  } else if(type == 1) {
    imp <- randomForest::importance(fit, type = 1)
  } else if(type == 2) {
    imp <- randomForest::importance(fit, type = 2)
  }
  if("%IncMSE" %in% colnames(imp)) {
    ind <- which(colnames(imp) == "%IncMSE")
    colnames(imp)[ind] <- "IncMSE"
  }
  return(as.data.frame(imp))
}

cf.imp.plot <- function(fit, type, n = 30, subgroup.id.result, data.type, level.name) { 
  
  imp.df1 <- cf.imp.df(fit, type = 1)
  imp.df2 <- cf.imp.df(fit, type = 2)
  
  ma <- get.mean.abundance(subgroup.id.result$Taxa, level.name)
  imp.df1 <- data.frame(imp.df1, ma)
  colnames(imp.df1)[2] <- "MeanAbundance"
  imp.df2 <- data.frame(imp.df2, ma)
  colnames(imp.df2)[2] <- "MeanAbundance"
  
  if(data.type == "clr") {
    data.type <- "CLR"
  } else if(data.type == "prop") {
    data.type <- "Proportion"
  } else if(data.type == "rare.count") {
    data.type <- "Rarefied Count"
  } else if(data.type == "arcsin") {
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
  
  if(type == 0) {
    return(imp.df1 + imp.df2)
  } else if(type == 1) {
    return(imp.df1)
  } else if(type == 2) {
    return(imp.df2)
  }
}

cf.pdp.reg <- function(fit, X, n, data.type) {
  rf.importance <- cf.imp.df(fit, type = 1) %>%
    mutate(taxon = rownames(cf.imp.df(fit, type = 1))) %>%
    arrange(-IncMSE) %>%
    rownames
  n <- min(length(rf.importance), n)
  feature <- rf.importance[1:n]
  result <- data.frame()
  
  if(data.type == "clr") {
    type = "CLR"
  }
  else if(data.type == "prop") {
    type = "Proportion"
  } else if(data.type == "rare.count") {
    type = "Rarefied Count"
  } else if(data.type == "arcsin") {
    type = "Arcsine-Root"
  }
  
  for(taxon.name in feature) {
    val <- numeric()
    y_hat <- numeric()
    taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
    for(i in 1:length(taxon.val)) {
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

cf.dt.used.var <- function(step2.result) {
  new.name <- setdiff(step2.result$best.dt.fit$frame$var, "<leaf>")
  var.name.df <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
  ori.name <- c()
  for(name in new.name) {
    ind <- which(name == var.name.df$sub)
    ori.name <- c(ori.name, var.name.df$ori[ind])
  }
  output <- data.frame(name = ori.name)
  rownames(output) <- new.name
  return(output)
}

rf.used.var <- function(step2.result, step4.result, n = 20) {
  new.name <- cf.imp.df(step4.result, type = 1) %>% arrange(desc(IncMSE)) %>% head(n) %>% rownames
  var.name.df <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
  ori.name <- c()
  for(name in new.name) {
    ind <- which(name == var.name.df$sub)
    ori.name <- c(ori.name, var.name.df$ori[ind])
  }
  output <- data.frame(name = ori.name)
  rownames(output) <- new.name
  return(output)
}

colnames.to.ind <- function(data) { 
  origin <- list()
  new <- list()
  for(name in names(data)) {
    ind.dat <- data[[name]]
    col.names <- colnames(ind.dat)
    new.names <- character()
    first.letter <- str_to_upper(substr(name, 1, 1))
    for(i in 1:length(col.names)) {
      new.names <- c(new.names, sprintf("%s%i", first.letter, i))
    }
    origin[[name]] <- col.names
    new[[name]] <- new.names
  }
  return(list(origin = origin, new = new))
}

colnames.df <- function(colnames.list, name) {
  d <- data.frame(colnames.list$origin[[name]])
  rownames(d) <- colnames.list$new[[name]]
  colnames(d) <- str_to_title(name)
  return(d)
}

change.colnames <- function(data, new.names) {
  for(name in names(data)) {
    colnames(data[[name]]) <- new.names[[name]]
  }
  return(data)
}

"%notin%" <- Negate("%in%")

get.level.names <- function(include = TRUE) {
  if(!include) {
    return(c("phylum", "class", "order", "family", "genus"))
  }
  return(c("phylum", "class", "order", "family", "genus", "species"))
}

get.mean.abundance <- function(data, rank.name) {
  return(colMeans(data))
}

remove.na <- function(data, sam.dat, y.name, level.names) {
  ind1 <- which(is.na(sam.dat[[y.name]]))
  ind2 <- which(sam.dat[[y.name]] == -9.00 | sam.dat[[y.name]] == -99.00 | sam.dat[[y.name]] == -999.00)
  ind <- sort(c(ind1, ind2))
  if(length(ind) > 0) {
    sam.dat.na <- sam.dat[-ind,]
    for(name in level.names) {
      data[[name]] <- data[[name]][-ind,]
    }
  } else {
    sam.dat.na <- sam.dat
  }
  return(list(data = data, sam.dat.na = sam.dat.na))
}

category.names <- function(sam.dat, y.name) {
  return(names(table(sam.dat[[y.name]])))
}

str.check <- function(sam.dat, y.name) {
  out <- character()
  var <- sam.dat[[y.name]]
  len <- length(table(var))
  if(len == 2) {
    return("Binary")
  }
  else if(len >= 3 & len <= 8) {
    return("Multinomial")
  } else if(len > 8 & is.numeric(var)) {
    return("Continuous")
  } else {
    return("Neither")
  }
}

bmc.col.check <- function(sam.dat, type = c("Binary", "Multinomial", "Continuous")) {
  dtype.vector <- character()
  for(name in colnames(sam.dat)) {
    dtype.vector <- c(dtype.vector, str.check(sam.dat, name))
  }
  ind <- which(dtype.vector == type)
  return(colnames(sam.dat[,ind]))
}

col.str.check <- function(sam.dat, name) {
  dtype <- character()
  if(length(table(sam.dat[[name]])) == 1) {
    dtype <- "none"
  } else if((is.character(sam.dat[[name]])) | (is.factor(sam.dat[[name]]))) {
    if(length(unique(sam.dat[[name]])) == nrow(sam.dat)) {
      dtype <- "none"
    } else {
      dtype <- "factor"
    }
  } else if(is.numeric(sam.dat[[name]])) {
    if(length(table(sam.dat[[name]])) == 2) {
      dtype <- "factor"
    } else if(length(unique(sam.dat[[name]])) == nrow(sam.dat)) {
      dtype <- "none"
    } else {
      dtype = "numeric"
    }
  } else {
    dtype = "none"
  }
  return(dtype)
}

get.cov.col <- function(sam.dat) {
  dtype <- character()
  names <- colnames(sam.dat)
  for(name in names) {
    dtype <- c(dtype, col.str.check(sam.dat, name))
  }
  cov.col <- names[dtype != "none"]
  return(cov.col)
}

get.cat.levels <- function(sam.dat, y.name) {
  levels <- levels(as.factor(unlist(sam.dat[[y.name]])))
  if(length(levels) >= 2 & length(levels) <= 4) {
    return(levels)
  } else {
    stop(sprintf("%s is not categorical variable.", y.name))
  }
}

remove.symb <- function(X) {
  rep_str <- c("-" = "_", ";" = "_", ":" = "_", " " = "_", "\\(" = "", "\\)" = "", "\\]" = "", "\\[" = "", "\\*" = "_", "^" = "_", "&" = "_", "\\=" = "_", "\\." = "")
  colnames(X) <- str_replace_all(colnames(X), rep_str)
  colnames(X) <- substr(colnames(X), 2, nchar(colnames(X)))
  return(X)
}

cov.remove.na <- function(data, sam.dat, y.name, cov.name, level.names) {
  new.sam.dat <- sam.dat[,c(cov.name, y.name)]
  ind <- sort(as.vector(which(is.na(new.sam.dat), arr.ind = TRUE)[,1]))
  if(length(ind) > 0) {
    sam.dat.na <- new.sam.dat[-ind,]
    for(name in level.names) {
      data[[name]] <- data[[name]][-ind,]
    }
  } else {
    sam.dat.na <- new.sam.dat
  }
  return(list(data = data, sam.dat.na = sam.dat.na))
}

cov.linear.reg <- function(sam.dat, y.name) {
  f1 <<- as.formula(paste(y.name, "~", ".", sep=" "))
  fit <- lm(f1, data = data.frame(sam.dat))
  resid <- resid(fit)
  sam.dat.resid <- cbind(sam.dat, resid)
  return(sam.dat.resid)
}

cov.logistic.reg <- function(sam.dat, y.name) {
  f1 <<- as.formula(paste(y.name, "~", ".", sep=" "))
  fit <- glm(f1, data = data.frame(sam.dat), family = "binomial")
  resid <- residuals(fit, type = "pearson")
  sam.dat.resid <- cbind(sam.dat, resid)
  return(sam.dat.resid)
}

cov.mult.logistic.reg <- function(sam.dat, y.name) {
  f1 <<- as.formula(paste(y.name, "~", ".", sep=" "))
  fit <- vglm(f1, family = multinomial, data = data.frame(sam.dat))
  resid <- residuals(fit, type = "pearson")
  n <- length(category.names(sam.dat, y.name)) - 1
  resid.name <- character()
  for(i in 1:n) {
    resid.name <- c(resid.name, paste0("resid",i))
  }
  colnames(resid) <- resid.name
  sam.dat.resid <- cbind(sam.dat, resid)
  return(sam.dat.resid)
}
