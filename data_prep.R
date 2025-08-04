### load lobe data
dat <- read.csv("cephBrain.csv")

### standardize function
standardize <- function(x) {
  x <- scale(x)
  z <- as.numeric(x)
  attr(z,"scaled:center") <- attr(x,"scaled:center")
  attr(z,"scaled:scale") <- attr(x,"scaled:scale")
  return(z)
}

# columns to exclude from logging and standardizing
excl_cols <- c("phy.species", "species", "depth_cat", "habitat3", "benthic", "sociality", "sociality.3", "sociality.bin", "include", "source.1", "source.2", "source.3", "biolum")

# log and standardize continuous vars----
d <- dat
d[, !(names(d) %in% excl_cols)] <- lapply(d[, !(names(d) %in% excl_cols)], function(x) round(standardize(log(x)), digits=3))

### set factor variables
d$benthic <- factor(d$benthic)
d$sociality.bin <- factor(d$sociality.bin)
d$sociality.3 <- factor(d$sociality.3)
d$habitat3 <- factor(d$habitat3)
d$depth_cat <- factor(d$depth_cat)
d$biolum <- factor(d$biolum)

### load phylo data
cephtreeMCC <- ape::read.nexus("cephtreeMCC.tree") #consensus (maximum clade credibility) tree
cormat <- ape::vcv(cephtreeMCC, corr=TRUE) #correlation matrix for consensus phylogeny

### add uncerscores to align naming with phylogeny
d$phy.species <- gsub(" ", "_", d$phy.species)