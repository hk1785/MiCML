logistic.regression.no.cov <- function(taxa_out, include, treatment, treatment_levels, response){
  
  reg_list <- list() 
  
  for (j in 1:(5+include)){
    
    taxa <- taxa_out[[j]]
    
    coef.out <- as.data.frame(matrix(NA, ncol(taxa), 4))
    conf.out <- as.data.frame(matrix(NA, ncol(taxa), 2))
    
    for (i in 1:ncol(taxa)) {
      
      taxon <- as.numeric(scale(taxa[,i]))
      interaction <- taxon*treatment
      
      taxon.name <- names(taxa)[i]
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      fit <- glm(response ~ taxon + treatment + interaction, family = binomial())
      
      if (nrow(coef(summary(fit))) != 4) {
        coef.out[i,] <- rep(NA, 4)
        conf.out[i,] <- rep(NA, 2)
      } else {
        coef.out[i,] <- coef(summary(fit))[4,]
        conf.out[i,] <- confint.default(fit)[4,]
      }
      
      rownames(coef.out)[i] <- taxon.name
      rownames(conf.out)[i] <- taxon.name
      
    }
    
    colnames(coef.out) <- c("Estimate", "Std. Error", "Z", "P.value")
    colnames(conf.out) <- c("2.5 %", "97.5 %")
    
    coef.out$Q.value <- p.adjust(coef.out$P.value, method = "BH", n = length(coef.out$P.value))
    
    int.out <- cbind(coef.out, conf.out)
    int.out
    
    reg_list[[j]] <- int.out
  }
  names(reg_list) <- names(taxa_out)[1:(5+include)]
  
  return(reg_list)
}

logistic.regression.with.cov <- function(taxa_out, include, treatment, treatment_levels, covariate, response){
  
  reg_list <- list() 
  
  for (j in 1:(5+include)){
    
    taxa <- taxa_out[[j]]
    
    coef.out <- as.data.frame(matrix(NA, ncol(taxa), 4))
    conf.out <- as.data.frame(matrix(NA, ncol(taxa), 2))
    
    for (i in 1:ncol(taxa)) {
      
      taxon <- as.numeric(scale(taxa[,i]))
      interaction <- taxon*treatment
      
      taxon.name <- names(taxa)[i]
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      fit <<- glm(response ~ taxon + treatment + interaction + covariate, family = binomial())
      
      coef.out[i,] <- coef(summary(fit))[4,]
      conf.out[i,] <- confint.default(fit)[4,]
      
      rownames(coef.out)[i] <- taxon.name
      rownames(conf.out)[i] <- taxon.name
      
    }
    
    colnames(coef.out) <- c("Estimate", "Std. Error", "Z", "P.value")
    colnames(conf.out) <- c("2.5 %", "97.5 %")
    
    coef.out$Q.value <- p.adjust(coef.out$P.value, method = "BH", n = length(coef.out$P.value))
    
    int.out <- cbind(coef.out, conf.out)
    int.out
    
    reg_list[[j]] <- int.out
  }
  names(reg_list) <- names(taxa_out)[1:(5+include)]
  
  return(reg_list)
}

height_output <- function(test_result, num_taxa = "Q-value < 0.05", include){
  
  height <- list() 
  
  for (i in 1:(5+include)){
    result <- test_result[[i]]
    
    ind.sel <- which(result$Q.value < 0.05)
    sel.out <- result[ind.sel,]
    
    if (length(ind.sel) == 0){
      height[[i]] <- "250px"
    } else {
      height[[i]] <- paste0(ceiling(length(ind.sel)/3)*250, "px")
    }
  }
  return(height)
}

gaussian.regression.no.cov <- function(taxa_out, include, treatment, treatment_levels, response){
  
  reg_list <- list() 
  
  for (j in 1:(5+include)){
    
    taxa <- taxa_out[[j]]
    
    coef.out <- as.data.frame(matrix(NA, ncol(taxa), 4))
    conf.out <- as.data.frame(matrix(NA, ncol(taxa), 2))
    
    for (i in 1:ncol(taxa)) {
      
      taxon <- as.numeric(scale(taxa[,i]))
      interaction <- taxon*treatment
      
      taxon.name <- names(taxa)[i]
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      fit <- glm(response ~ taxon + treatment + interaction, family = gaussian())
      
      if (nrow(coef(summary(fit))) != 4) {
        coef.out[i,] <- rep(NA, 4)
        conf.out[i,] <- rep(NA, 2)
      } else {
        coef.out[i,] <- coef(summary(fit))[4,]
        conf.out[i,] <- confint.default(fit)[4,]
      }
      rownames(coef.out)[i] <- taxon.name
      rownames(conf.out)[i] <- taxon.name
      
    }
    colnames(coef.out) <- c("Estimate", "Std. Error", "Z", "P.value")
    colnames(conf.out) <- c("2.5 %", "97.5 %")
    
    coef.out$Q.value <- p.adjust(coef.out$P.value, method = "BH", n = length(coef.out$P.value))
    
    int.out <- cbind(coef.out, conf.out)
    int.out
    
    reg_list[[j]] <- int.out
  }
  names(reg_list) <- names(taxa_out)[1:(5+include)]
  
  return(reg_list)
}

gaussian.regression.with.cov <- function(taxa_out, include, treatment, treatment_levels, covariate, response){
  
  reg_list <- list() 
  
  for (j in 1:(5+include)){
    
    taxa <- taxa_out[[j]]
    
    coef.out <- as.data.frame(matrix(NA, ncol(taxa), 4))
    conf.out <- as.data.frame(matrix(NA, ncol(taxa), 2))
    
    for (i in 1:ncol(taxa)) {
      
      taxon <- as.numeric(scale(taxa[,i]))
      interaction <- taxon*treatment
      
      taxon.name <- names(taxa)[i]
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      fit <- glm(response ~ taxon + treatment + interaction + covariate, family = gaussian())
      
      if (nrow(coef(summary(fit))) != 4 + ncol(covariate)) {
        coef.out[i,] <- rep(NA, 4)
        conf.out[i,] <- rep(NA, 2)
      } else {
        coef.out[i,] <- coef(summary(fit))[4,]
        conf.out[i,] <- confint.default(fit)[4,]
      }
      
      rownames(coef.out)[i] <- taxon.name
      rownames(conf.out)[i] <- taxon.name
      
    }
    
    colnames(coef.out) <- c("Estimate", "Std. Error", "Z", "P.value")
    colnames(conf.out) <- c("2.5 %", "97.5 %")
    
    coef.out$Q.value <- p.adjust(coef.out$P.value, method = "BH", n = length(coef.out$P.value))
    
    int.out <- cbind(coef.out, conf.out)
    int.out
    
    reg_list[[j]] <- int.out
  }
  names(reg_list) <- names(taxa_out)[1:(5+include)]
  
  return(reg_list)
}

result_for_output <- function(taxa_out, test_result, treatment, treatment_levels, response, taxa_rank, num_taxa = "Q-value < 0.05", covariates = NULL, response.type = "binary"){
  
  taxa <- taxa_out[[taxa_rank]]
  result <- test_result[[taxa_rank]]
  
  ind.sel <- which(result$Q.value < 0.05)
  sel.out <- result[ind.sel,]
  
  list_result <- list() 
  
  print(ind.sel)
  
  for (i in ind.sel) {
    
    taxon <- as.numeric(scale(taxa[,i]))
    taxon.name <- names(taxa)[i]
    
    interaction <- taxon*treatment
    
    taxon.name <- unlist(strsplit(taxon.name, ";"))
    taxon.name <- taxon.name[length(taxon.name)]
    
    int.qval <- p.value.0.1(result$Q.value[i])
    
    if(response.type == "binary") {
      
      if (is.null(covariates)) {
        
        fit <- glm(response ~ taxon + treatment + interaction, family = binomial())
        
        pred.n <- 100
        
        pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
        pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
        pred.dat.0$response <- predict(fit, pred.dat.0, type = "response")
        
        pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
        pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
        pred.dat.1$response <- predict(fit, pred.dat.1, type = "response")
        
      } else {
        
        fit <- glm(response ~ taxon + treatment + interaction + covariates, family = binomial())
        
        pred.n <- 100
        
        pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
        pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
        
        pred.0 <- numeric()
        
        for (j in 1:pred.n) {
          p.dat_0 <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.0[j,]))), covariates))
          colnames(p.dat_0) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
          
          pred.0[j] <- mean(predict(fit, p.dat_0, type = "response"))
        }
        
        pred.dat.0$response <- pred.0
        
        pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
        pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
        pred.1 <- numeric()
        
        for (j in 1:pred.n) {
          p.dat_1 <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.1[j,]))), covariates))
          colnames(p.dat_1) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
          pred.1[j] <- mean(predict(fit, p.dat_1, type = "response"))
        }
        
        pred.dat.1$response <- pred.1
        
      }
      
    } else if (response.type == "continuous") {
      
      if (is.null(covariates)){
        
        fit <- glm(response ~ taxon + treatment + interaction, family = gaussian())
        
        pred.n <- 100
        
        pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
        pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
        pred.dat.0$response <- predict(fit, pred.dat.0, type = "response")
        
        pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
        pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
        pred.dat.1$response <- predict(fit, pred.dat.1, type = "response")
        
      } else {
        
        fit <- glm(response ~ taxon + treatment + interaction + covariates, family = gaussian())
        
        pred.n <- 100
        
        pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
        pred.dat.0$Interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
        pred.0 <- numeric()
        
        for (j in 1:pred.n) {
          
          p.dat <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.0[j,]))), covariates))
          colnames(p.dat) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
          pred.0[j] <- mean(predict(fit, p.dat, type = "response"))
          
        }
        
        pred.dat.0$response <- pred.0
        
        pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
        pred.dat.1$Interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
        pred.1 <- numeric()
        for (j in 1:pred.n) {
          p.dat <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.1[j,]))), covariates))
          colnames(p.dat) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
          pred.1[j] <- mean(predict(fit, p.dat, type = "response"))
        }
        
        pred.dat.1$response <- pred.1
      }
    }
    list_result[[i]] <- list(pred.dat.0, pred.dat.1)
  }
  
  return(list_result)
}

output.to.show <- function(taxa_out, test_result, treatment, treatment_levels, response, taxa_rank, num_taxa = "Q-value < 0.05", legend, covariates = NULL, response.type = "binary"){
  
  taxa <- taxa_out[[taxa_rank]]
  result <- test_result[[taxa_rank]]
  
  if (legend == "Top left") {
    leg.pos <- "topleft"
  } else if (legend == "Top right") {
    leg.pos <- "topright"
  } else if (legend == "Bottom left") {
    leg.pos <- "bottomleft"
  } else if (legend == "Bottom right") {
    leg.pos <- "bottomright"
  }
  
  ind.sel <- which(result$Q.value < 0.05)
  sel.out <- result[ind.sel,]
  
  if (length(ind.sel) == 0){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  } else {
    
    par(mfcol = c(ceiling(length(ind.sel)/3), 3))
    
    for (i in ind.sel) {
      
      taxon <- as.numeric(scale(taxa[,i]))
      taxon.name <- names(taxa)[i]
      
      interaction <- taxon*treatment
      
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      int.qval <- p.value.0.1(result$Q.value[i])
      
      if(response.type == "binary"){
        
        if (is.null(covariates)){
          
          fit <- glm(response ~ taxon + treatment + interaction, family = binomial())
          
          pred.n <- 100
          
          pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
          pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
          pred.dat.0$response <- predict(fit, pred.dat.0, type = "response")
          
          pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
          pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
          pred.dat.1$response <- predict(fit, pred.dat.1, type = "response")
          
          plot(taxon, response, ylim = c(0, 1), xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
          lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
          lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
          
          legend(leg.pos, title = "", legend = treatment_levels, 
                 fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
          title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
          
        } else {
          
          fit <- glm(response ~ taxon + treatment + interaction + covariates, family = binomial())
          
          pred.n <- 100
          
          pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
          pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
          
          pred.0 <- numeric()
          
          for (j in 1:pred.n) {
            p.dat_0 <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.0[j,]))), covariates))
            colnames(p.dat_0) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
            
            pred.0[j] <- mean(predict(fit, p.dat_0, type = "response"))
          }
          
          pred.dat.0$response <- pred.0
          
          pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
          pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
          pred.1 <- numeric()
          
          for (j in 1:pred.n) {
            p.dat_1 <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.1[j,]))), covariates))
            colnames(p.dat_1) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
            pred.1[j] <- mean(predict(fit, p.dat_1, type = "response"))
          }
          
          pred.dat.1$response <- pred.1
          
          plot(taxon, response, ylim = c(0, 1), xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
          
          lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
          lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
          
          legend(leg.pos, title = "", legend = treatment_levels, 
                 fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
          
          title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
        }
        
      } else if (response.type == "continuous") {
        
        if (is.null(covariates)){
          
          fit <- glm(response ~ taxon + treatment + interaction, family = gaussian())
          
          pred.n <- 100
          
          pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
          pred.dat.0$interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
          pred.dat.0$response <- predict(fit, pred.dat.0, type = "response")
          
          pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
          pred.dat.1$interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
          pred.dat.1$response <- predict(fit, pred.dat.1, type = "response")
          
          plot(taxon, response, xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
          
          lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
          lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
          legend(leg.pos, title = "", legend = treatment_levels, 
                 fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
          title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
          
        } else {
          
          fit <- glm(response ~ taxon + treatment + interaction + covariates, family = gaussian())
          
          pred.n <- 100
          
          pred.dat.0 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(0, pred.n))
          pred.dat.0$Interaction <- as.numeric(pred.dat.0$taxon * pred.dat.0$treatment)
          pred.0 <- numeric()
          
          for (j in 1:pred.n) {
            
            p.dat <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.0[j,]))), covariates))
            colnames(p.dat) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
            pred.0[j] <- mean(predict(fit, p.dat, type = "response"))
            
          }
          
          pred.dat.0$response <- pred.0
          
          pred.dat.1 <- data.frame(taxon = seq(min(taxon), max(taxon), len = pred.n), treatment = rep(1, pred.n))
          pred.dat.1$Interaction <- as.numeric(pred.dat.1$taxon * pred.dat.1$treatment)
          pred.1 <- numeric()
          for (j in 1:pred.n) {
            p.dat <- as.data.frame(cbind(t(replicate(length(response), as.numeric(pred.dat.1[j,]))), covariates))
            colnames(p.dat) <- c("taxon", "treatment", "interaction", colnames(covariates)) 
            pred.1[j] <- mean(predict(fit, p.dat, type = "response"))
          }
          
          pred.dat.1$response <- pred.1
          
          plot(taxon, response, xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
          
          lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
          lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
          legend(leg.pos, title = "", legend = treatment_levels, 
                 fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
          title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
          
        }
      }
    }
  }
}

output.to.show.2 <- function(taxa_out, test_result, result_list, treatment, treatment_levels, response, taxa_rank, num_taxa = "Q-value < 0.05", legend, covariates = NULL, response.type = "binary"){
  
  taxa <- taxa_out[[taxa_rank]]
  result <- test_result[[taxa_rank]]
  
  if (legend == "Top left") {
    leg.pos <- "topleft"
  } else if (legend == "Top right") {
    leg.pos <- "topright"
  } else if (legend == "Bottom left") {
    leg.pos <- "bottomleft"
  } else if (legend == "Bottom right") {
    leg.pos <- "bottomright"
  }
  
  ind.sel <- which(result$Q.value < 0.05)
  sel.out <- result[ind.sel,]
  
  if (length(ind.sel) == 0){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  } else {
    par(mfcol = c(ceiling(length(ind.sel)/3), 3))
    for (i in ind.sel) {
      taxon <- as.numeric(scale(taxa[,i]))
      taxon.name <- names(taxa)[i]
      
      interaction <- taxon*treatment
      
      taxon.name <- unlist(strsplit(taxon.name, ";"))
      taxon.name <- taxon.name[length(taxon.name)]
      
      int.qval <- p.value.0.1(result$Q.value[i])
      
      pred.dat.0 <- result_list[[i]][[1]]
      pred.dat.1 <- result_list[[i]][[2]]
      
      if(response.type == "binary") {
        plot(taxon, response, ylim = c(0, 1), xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
        
        lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
        lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
        
        legend(leg.pos, title = "", legend = treatment_levels, 
               fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
        
        title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
        
      } else if (response.type == "continuous") {
        plot(taxon, response, xlab = taxon.name, ylab = "Response", pch = 20, cex.lab = 1.2)
        
        lines(response ~ taxon, pred.dat.0, lwd = 2, col = "blue2")
        lines(response ~ taxon, pred.dat.1, lwd = 2, col = "orange2")
        legend(leg.pos, title = "", legend = treatment_levels, 
               fill = c("blue2", "orange2"), horiz = TRUE, bty = "n", cex = 1.2, pt.cex = 1.2, pt.lwd = 1.2)
        title(paste("Q-value: ", int.qval, sep = ""), adj = 1, cex.main = 1.3)
      }
    }
  }
}


cov.factorize.func <- function(sam.dat, cov_var){
  cov_mat <- c()
  i = 1
  
  for (cov in cov_var){
    covariate <- sam.dat[[cov]]
    
    if (is.character(covariate)){
      cov_mat <- cbind(cov_mat, model.matrix(~ as.factor(covariate))[,-1])
    } else {
      cov_mat <- cbind(cov_mat, as.numeric(covariate))
    }
    i = i+1
  }
  return(cov_mat)
}

