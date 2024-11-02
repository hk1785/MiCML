setSliderColor <- function(color, sliderId) {
  
  stopifnot(!is.null(color))
  stopifnot(is.character(color))
  stopifnot(is.numeric(sliderId))
  stopifnot(!is.null(sliderId))
  sliderId <- sliderId - 1
  
  sliderCol <- lapply(sliderId, FUN = function(i) {
    paste0(
      ".js-irs-", i, " .irs-single,",
      " .js-irs-", i, " .irs-from,",
      " .js-irs-", i, " .irs-to,",
      " .js-irs-", i, " .irs-bar-edge,",
      " .js-irs-", i,
      " .irs-bar{  border-color: transparent;background: ", color[i+1],
      "; border-top: 1px solid ", color[i+1],
      "; border-bottom: 1px solid ", color[i+1],
      ";}"
    )
  })
  custom_head <- tags$head(tags$style(HTML(as.character(sliderCol))))
  return(custom_head)
}

# Quality Control -----------

rem.tax.d <- c("", "metagenome", "gut metagenome", "mouse gut metagenome")
rem.tax.str.d <- c("uncultured", "incertae", "Incertae", "unidentified", "unclassified", "unknown")

tax.tab.clean <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- na.code
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  ind <- which(tax.tab.h[,2] != na.code)
  ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  ind <- which(tax.tab.h[,3] != na.code)
  ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  ind <- which(tax.tab.h[,4] != na.code)
  ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  ind <- which(tax.tab.h[,5] != na.code)
  ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  ind <- which(tax.tab.h[,6] != na.code)
  ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  ind <- which(tax.tab.h[,7] != na.code)
  ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

otu.tab.clean <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, tree.exists = TRUE) {
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- if(tree.exists) phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  if(tree.exists){
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  }
  else{
    biom <- merge_phyloseq(otu.tab, tax.tab, sam.dat)
  }
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}

biom.clean <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, tree.exists = TRUE) {
  
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  if(tree.exists) tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  if(tree.exists){
    ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  } else {
    ind.com.otu <- intersect(rownames(otu.tab), rownames(tax.tab))
  }
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- if(tree.exists) prune_taxa(ind.com.otu, tree)
  
  if(tree.exists) {
    if(!is.rooted(tree)) {
      tree <- phangorn::midpoint(tree)
    }
  }
  if (tax.tab.c) {
    tax.tab <- tax.tab.clean(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  if(tree.exists) {
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  } else {
    biom <- merge_phyloseq(otu.tab, tax.tab, sam.dat)
  }
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off, tree.exists = tree.exists)
  
  return(biom)  
}

num.tax.rank <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.cleaned <- tax.tab.clean(tax.tab, rem.tax, rem.tax.str, na.code = na.code)
  num.taxa <- c()
  for (i in 1:6) {
    taxa <- unique(tax.tab.cleaned[,i+1])
    uni.taxa <- sum(taxa == na.code)
    num.taxa[i] <- nrow(taxa) - uni.taxa
  }
  return(num.taxa)
}

lib.size.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  lib.size.sum <- c(mean(lib.size), quantile(lib.size))
  names(lib.size.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(lib.size = lib.size, lib.size.sum = lib.size.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

mean.prop.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  mean.prop.sum <- c(mean(mean.prop), quantile(mean.prop))
  names(mean.prop.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(mean.prop = mean.prop, mean.prop.sum = mean.prop.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

# Rarefying -----------

rarefy.func <- function(biom, cut.off, multi.rarefy = FALSE, tree.exists = TRUE) {
  if (!multi.rarefy | multi.rarefy == 1) {
    try(biom <- rarefy_even_depth(biom, cut.off, rngseed = 487), silent = TRUE)
  } else {
    otu.tab <- otu_table(biom)
    tax.tab <- tax_table(biom)
    if(tree.exists) tree <- phy_tree(biom)
    sam.dat <- sample_data(biom)
    
    otu.tab.list <- list()
    for (i in 1:multi.rarefy) {
      try(otu.tab.list[[i]] <- otu_table(rarefy_even_depth(biom, cut.off, rngseed = i), taxa_are_rows = TRUE), silent = TRUE)
    }
    
    sum.otu.tab <- otu.tab.list[[1]]
    for (i in 2:multi.rarefy) {
      sum.otu.tab <- sum.otu.tab + otu.tab.list[[i]]
    }
    otu.tab <- otu_table(round(sum.otu.tab/multi.rarefy), taxa_are_rows = TRUE)
    if(tree.exists) {
      biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat) 
    } else {
      biom <- merge_phyloseq(otu.tab, tax.tab, sam.dat) 
    }
  }
  
  return(biom)
}

is.mon.sin.rev.bin.con <- function(sam.dat) {
  
  n.var <- ncol(sam.dat)
  n.sam <- nrow(sam.dat)
  is.mon <- logical()
  is.rev <- logical()
  is.bin <- logical()
  is.con <- logical()
  is.sin <- logical()
  
  for (i in 1:n.var) {
    sam.var <- as.matrix(sam.dat[,i])
    if (length(table(sam.var)) == 1) {
      is.mon[i] <- TRUE
    }
    if (length(table(sam.var)) != 1) {
      is.mon[i] <- FALSE
    }
    if (length(table(sam.var)) == 2 & any(table(sam.var)==1)) {
      is.sin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2 | !any(table(sam.var)==1)) {
      is.sin[i] <- FALSE
    }
    if (length(table(sam.var)) == n.sam & sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.rev[i] <- TRUE
    }
    if (length(table(sam.var)) != n.sam | sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.rev[i] <- FALSE
    }
    if (length(table(sam.var)) == 2) {
      is.bin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2) {
      is.bin[i] <- FALSE
    }
    if (length(table(sam.var)) != 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- TRUE
    }
    if (length(table(sam.var)) == 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- FALSE
    }
    if (sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.con[i] <- FALSE
    }
    
  }
  return(list(is.mon = is.mon, is.sin = is.sin, is.rev = is.rev, is.bin = is.bin, is.con = is.con))
}

pri.func <- function(sam.dat, mon.sin.rev.bin.con) {
  colnames(sam.dat)[(mon.sin.rev.bin.con$is.bin | mon.sin.rev.bin.con$is.con) & !mon.sin.rev.bin.con$is.mon & !mon.sin.rev.bin.con$is.sin]
}

is.bin.con.pri <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind <- which(colnames(sam.dat) == sel.pri.var)
  if(length(ind) != 0){
    if (mon.sin.rev.bin.con$is.bin[ind]) {
      out <- "Binary"
    } else {
      out <- "Continuous"
    }
  }else {
    out = "Neither"
  }
  return(out)
}

# Diversity Calculation -----------

alpha.v1.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice, alpha.pd))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  
  return(alpha.div)
}

Ds.Ks.func <- function(rare.biom, biom.after.qc) {
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
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
  
  return(
    list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
         Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k))
  )
}

# Data Transformation -----------

tax.trans <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "NANANA") {
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    try(tax.clr <- compositions::clr(tax.imp.prop), silent = TRUE)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
}

taxa.names.rank <- function(taxa.out){ 
  taxon.names <- list()
  taxon.names <- lapply(taxa.out, function(x) str_split(names(x), ";"))
  dup.list <- list(NA,NA,NA,NA,NA,NA)
  
  ranks <- c("K_", "P_", "C_", "O_", "F_", "G_", "S_")
  
  taxon.names.rank <- list()
  for(rank in 1:6){
    taxon <- lapply(taxon.names[[rank]], function(x) str_sub(x,start = 3))
    taxon.names.rank[[rank]] <- sapply(taxon, tail, 1)
    
    if(length(taxon.names.rank[[rank]]) != length(unique(taxon.names.rank[[rank]]))){
      duplicated.taxons <- unique(taxon.names.rank[[rank]][duplicated(taxon.names.rank[[rank]])])
      
      for(i in 1:length(duplicated.taxons)){
        duplicated.taxon <- duplicated.taxons[i]
        ind.dup <- which(taxon.names.rank[[rank]] %in% duplicated.taxon)
        
        for(j in 1:length(ind.dup)){
          duplicated.taxon <- paste(duplicated.taxon,"*",collapse = "")
          taxon.names.rank[[rank]][ind.dup[j]] <- duplicated.taxon 
          dup.list[[rank]][j] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[j]]), sep = ""), collapse = " | "), sep = "")
        }
      }
    }
  }
  names(taxon.names.rank) <- names(taxa.out)
  return(list(names = taxon.names.rank, duplicates = dup.list))
}

taxa.bin.var.func <- function(sam.dat) {
  var.names <- colnames(sam.dat)
  return(var.names)
}

# Batch Effect Correction -----------

batch.correct.PCoA <- function(qc.biom, rare.biom, qc.biom.bat, rare.biom.bat, batch.var){
  par(mfrow = c(2, 4))
  
  qc.otu.tab <- otu_table(qc.biom)
  qc.sam.dat <- sample_data(qc.biom)
  
  ra.otu.tab <- otu_table(rare.biom)
  ra.sam.dat <- sample_data(rare.biom)
  
  batid <- as.factor(qc.sam.dat[[batch.var]])
  batid.first.char <- as.factor(substr(batid, 1, 1))
  
  cex <- 1.7
  cex.main <- 1.3
  cex.sub <- 1
  label.cex <- 1.3
  cex.lab <- 1.2
  
  ### Bray-Curtis
  
  batch.mod <- betadisper(bcdist(t(ra.otu.tab)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Bray-Curtis\n(Before)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Proportion
  
  batch.mod <- betadisper(compositions::dist(t(qc.otu.tab)/colSums(qc.otu.tab)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Proportion\n(Before)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Arcsine
  
  batch.mod <- betadisper(compositions::dist(asin(sqrt(t(qc.otu.tab)/colSums(qc.otu.tab)))), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Arcsine\n(Before)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Aitchison (CLR)
  
  batch.mod <- betadisper(compositions::dist(clr(t(qc.otu.tab) + 0.1)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Aitchison\n(Before)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  qc.otu.tab <- otu_table(qc.biom.bat)
  qc.sam.dat <- sample_data(qc.biom.bat)
  
  ra.otu.tab <- otu_table(rare.biom.bat)
  ra.sam.dat <- sample_data(rare.biom.bat)
  
  batid <- as.factor(qc.sam.dat[[batch.var]])
  batid.first.char <- as.factor(substr(batid, 1, 1))
  
  ### Bray-Curtis
  
  batch.mod <- betadisper(bcdist(t(ra.otu.tab)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Bray-Curtis\n(After)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Proportion
  
  batch.mod <- betadisper(compositions::dist(t(qc.otu.tab)/colSums(qc.otu.tab)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Proportion\n(After)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Arcsine
  
  batch.mod <- betadisper(compositions::dist(asin(sqrt(t(qc.otu.tab)/colSums(qc.otu.tab)))), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Arcsine\n(After)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
  
  ### Aitchison (CLR)
  
  batch.mod <- betadisper(compositions::dist(clr(t(qc.otu.tab) + 0.1)), batid.first.char)
  
  plot(batch.mod, ellipse = TRUE, hull = FALSE, xlab = "PC 1", ylab = "PC 2", main = "Aitchison\n(After)", sub = "",
       col = 1:nlevels(batid.first.char), mgp = c(2.5, 1, 0), pch = 1:nlevels(batid.first.char),
       cex = cex, cex.main = cex.main, cex.sub = cex.sub, label.cex = label.cex, cex.lab = cex.lab)
}

preprocess.tax.tab = function(tax.tab){
  trans.tax.tab <- matrix(NA, nrow(tax.tab), 7)
  tax.list <- strsplit(as.character(tax.tab$Taxon), ";")
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_0__", "", tax.list[[i]][grepl("D_0__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,1] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_1__", "", tax.list[[i]][grepl("D_1__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,2] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_2__", "", tax.list[[i]][grepl("D_2__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,3] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_3__", "", tax.list[[i]][grepl("D_3__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,4] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_4__", "", tax.list[[i]][grepl("D_4__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,5] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_5__", "", tax.list[[i]][grepl("D_5__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,6] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_6__", "", tax.list[[i]][grepl("D_6__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,7] <- taxon
    }
  }
  rownames(trans.tax.tab) <- tax.tab$Feature.ID
  colnames(trans.tax.tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  return(trans.tax.tab)
}

biom.check.samples <- function(otu.table, sample.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  sam.dat <- sample_data(sample.data)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  
  if (length(intersect(colnames(otu.tab), rownames(sam.dat))) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

biom.check.otu <- function(otu.table, tax.table) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  tax.tab <- tax_table(tax.table)
  
  
  if (length(intersect(rownames(otu.tab), rownames(tax.tab))) == 0) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

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
    
  } else if (is.tree == "withoutTree") {
    
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
  legend("center", title = NULL, legend = levels(y), fil = c(1:nlevels(y), cex = 2.5, box.lty = 0), 
         bty = "n", cex=1.5)
}

