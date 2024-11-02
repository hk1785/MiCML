sel.var.binary <- function(sam_dat, not.sel.var){
  vars <- names(sam_dat)[names(sam_dat) != not.sel.var]
  vars_to_show <- vars[sapply(names(sam_dat), function(x){return(length(table(sam_dat[[x]])) == 2)})]
  return(vars_to_show)
}

select.covariates.func <- function(sam_dat, treat.var, response.var) {
  vars <- names(sam_dat)[!(names(sam_dat) %in% c(treat.var, response.var))]
  return(vars)
}

rename.bin.var <- function(vec, ori.cat, new.cat) {
  treat.vec <- replace(vec, vec == ori.cat[1], new.cat[1])
  treat.vec <- replace(treat.vec, treat.vec == ori.cat[2], new.cat[2])
  treat.var <- as.factor(treat.vec)
  return(treat.var)
}

bin.fisher.test <- function(Response, Treatment){
  count.tab <- table(Response, Treatment)
  fit <- fisher.test(count.tab)
  p.val <- round(fit$p.value, 3)
  return(p.val)
}

bin.chisq.test <- function(Response, Treatment){
  count.tab <- table(Response, Treatment)
  fit <- chisq.test(count.tab)
  p.val <- round(fit$p.value, 3)
  return(p.val)
}

con.wilcox.test <- function(Response, Treatment){
  fit <- wilcox.test(Response ~ Treatment)
  p.val <- round(fit$p.value, 3)
  return(p.val)
}

con.welch.test <- function(Response, Treatment){
  fit <- t.test(Response ~ Treatment)
  p.val <- round(fit$p.value, 3)
  return(p.val)
}

bin.des.plot <- function(Response, Treatment, legend, p_val, test){
  
  if (legend == "Top left") {
    panel <- "both"
    leg.pos <- "topleft"
  } else if (legend == "Top right") {
    panel <- "both"
    leg.pos <- "topright" 
  } else if (legend == "Top left on the left panel only") {
    panel <- "left"
    leg.pos <- "topleft"
  } else if (legend == "Top right on the left panel only"){
    panel <- "left"
    leg.pos <- "topright" 
  } else if (legend == "Top left on the right panel only"){
    panel <- "right"
    leg.pos <- "topleft"
  } else if (legend == "Top right on the right panel only"){
    panel <- "right"
    leg.pos <- "topright"
  }
    
  count.tab <- table(Response, Treatment)
  prop.tab <- prop.table(count.tab, 2)
  
  par(mfrow = c(1, 2))
  
  if (panel == "both") {
    
    barplot(count.tab, xlab = "Count", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
    legend(leg.pos, title = "Response", legend = substr(levels(Response), 1, 8), fill = c("lightyellow", "lightblue"), horiz = TRUE, bty = "n")
    barplot(prop.tab, xlab = "Proportion", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
    legend(leg.pos, title = "Response", legend = substr(levels(Response), 1, 8), fill = c("lightyellow", "lightblue"), horiz = TRUE, bty = "n")
  
    } else if (panel == "right") {
    
    barplot(count.tab, xlab = "Count", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
    barplot(prop.tab, xlab = "Proportion", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
    legend(leg.pos, title = "Response", legend = substr(levels(Response), 1, 8), fill = c("lightyellow", "lightblue"), horiz = TRUE, bty = "n")
  
    } else if (panel == "left") {
    
      barplot(count.tab, xlab = "Count", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
      legend(leg.pos, title = "Response", legend = substr(levels(Response), 1, 8), fill = c("lightyellow", "lightblue"), horiz = TRUE, bty = "n")
      barplot(prop.tab, xlab = "Proportion", col = c("lightyellow", "lightblue"), beside = TRUE, cex.lab = 1.2)
      
  }
  
  if (test == "Fisher's exact test (Default)") {
    title(paste("P-value (Fisher): ", p.value.0.1(p_val), sep = ""), cex.main = 1.2)
  } else if (test == "Pearson's Chi-squared test") {
    title(paste("P-value (Chi-sq): ", p.value.0.1(p_val), sep = ""), cex.main = 1.2)
  }
  
}

con.des.plot <- function(Response, Treatment, Treatment.cat, legend, p_val, test){
  
  if (legend == "Top left") {
    leg.pos <- "topleft"
  } else if (legend == "Top right") {
    leg.pos <- "topright" 
  } 
  
  par(mfrow = c(1, 2))
  boxplot(Response ~ Treatment, names = substr(levels(Treatment), 1, 8), col = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)), notch = FALSE)
  res.con <- Response[Treatment == Treatment.cat[1]]
  res.trt <- Response[Treatment == Treatment.cat[2]]
  
  res.con.ker.den.est <- density(res.con, kernel = "gaussian")
  res.trt.ker.den.est <- density(res.trt, kernel = "gaussian")
  
  plot(res.con.ker.den.est, xlab = "Response", ylab = "Density", 
       xlim = c(min(res.con.ker.den.est$x, res.trt.ker.den.est$x), max(res.con.ker.den.est$x, res.trt.ker.den.est$x)), 
       ylim = c(0, max(res.con.ker.den.est$y, res.trt.ker.den.est$y)*1.1),
       main = "", col = rgb(0,1,0,0.5), lwd = 2)
  
  lines(res.trt.ker.den.est, xlab = "Response", ylab = "Density", 
        xlim = c(min(res.con.ker.den.est$x, res.trt.ker.den.est$x), max(res.con.ker.den.est$x, res.trt.ker.den.est$x)), 
        ylim = c(0, max(res.con.ker.den.est$y, res.trt.ker.den.est$y)*1.1),
        main = "", col = rgb(1,0,0,0.5), lwd = 2)
  
  if (test == "Mann-Whitney test (Default)") {
    title(paste("P-value (Mann-Whitney): ", p.value.0.1(p_val), sep = ""), cex.main = 1.2)
  } else if (test == "Welch's t-test") {
    title(paste("P-value (Welch): ", p.value.0.1(p_val), sep = ""), cex.main = 1.2)
  }
  
  legend(leg.pos, title = "Treatment", legend = substr(levels(Treatment), 1, 8), fill = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)), horiz = TRUE, bty = "n")
  
}

