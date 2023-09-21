library(phyloseq)
library(proxy)
library(ecodist)
library(GUniFrac)
library(MiRKAT)
library(MiRKATMC)
library(vegan)

source("MiDataProc.Data.Input.R")
source("MiDataProc.Data.Upload.R")

# MiRKATMC ----------------------
Ds.Ks.func <- function(rare.biom, biom.after.qc, is.tree) {
  
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  
  if (is.tree == "withTree") {
    
    no.rare.tree <- phy_tree(biom.after.qc)
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
    u.unif <- unifs[, , "d_UW"]
    g.unif <- unifs[, , "d_0.5"]
    w.unif <- unifs[, , "d_1"]
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    u.unif.k <- D2K(u.unif)
    g.unif.k <- D2K(g.unif)
    w.unif.k <- D2K(w.unif)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
    rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
    rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
    
    return(list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
                Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k)))
    
  } 
  else if (is.tree == "withoutTree") {
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    
    return(list(Ds = list(Jaccard = jac, Bray.Curtis = bc), Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k)))
  }
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}

beta.PCoA.plot <- function(Ds, mirkatmc.pvs, y){
  par(mfrow = c(3, 2))
  for (i in 1:length(Ds)) {
    if (mirkatmc.pvs[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(mirkatmc.pvs[i]), sep="")
    }
    if (mirkatmc.pvs[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(mirkatmc.pvs[i]), sep="")
    }
    mod <- betadisper(as.dist(Ds[[i]]), y)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = 1:nlevels(y), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(y), fil = c(1:nlevels(y), cex=2.5, box.lty=0), 
         bty = "n", cex=1.5)
}