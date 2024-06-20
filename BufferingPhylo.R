#Code to examine the phylogenetic signal in the sum stoch elasticities
#Rob Salguero-Gomez - rob.salguero@biology.ox.ac.uk

# clean up
rm(list=ls())

#Load up libraries
library(caper)
library(rotl)
library(phylosignal)
library(phytools)

#Read in data
d <- read.csv("buffering.csv")

#House-keeping
dim(d)
head(d)
colnames(d) <- c("Species","SumStochElas","MPMs")
d$Species <- gsub("_[0-9]","",d$Species)
d$Species <- gsub(" ","",d$Species)
d$Species <- gsub("_"," ",d$Species)
d$Species <- gsub("subsp. ","",d$Species)
dim(d)


#Get rid of duplicated populations per species, just taking in the first one into consideration
dim(d)
d <- d[-which(duplicated(d$Species)==T),]
dim(d)

#There's a "Gorilla beringei" and a ""Gorilla beringei beringei"
d <- d[-which(d$Species=="Gorilla beringei"),]
dim(d)

#Building phylogenetic tree
output_species_resolved <- tnrs_match_names(d$Species)
tree <- tol_induced_subtree(ott_ids = output_species_resolved$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = T)
tree$node.label <- NULL
length(tree$tip.label)
dim(d)

d$Species[which(!d$Species %in% tree$tip.label)]
#[1] "Macropus eugenii"            "Spermophilus armatus"       
#[3] "Spermophilus columbianus"    "Ursus americanus floridanus"

tree$tip.label[which(!tree$tip.label%in%d$Species)]
#[1] "Urocitellus columbianus"     "Urocitellus armatus"        
#[3] "Ursus americanus americanus" "Notamacropus eugenii"  

tree$tip.label[order(tree$tip.label)]
d$Species[order(d$Species)]

tree$tip.label[which(tree$tip.label=="Ursus americanus americanus")] <- "Ursus americanus floridanus"
tree$tip.label[which(tree$tip.label=="Urocitellus armatus")] <- "Spermophilus armatus"
tree$tip.label[which(tree$tip.label=="Urocitellus columbianus")] <- "Spermophilus columbianus"
tree$tip.label[which(tree$tip.label=="Notamacropus eugenii")] <- "Macropus eugenii"

tree <- compute.brlen(tree)
tree$edge.length
is.ultrametric(tree)

plot(tree, type='radial')
plot(tree)

#Phylogenetic signal
row.names(d) <- d$Species
d <- d[order(match(d$Species,tree$tip.label)),]
d$SumStochElas <- as.numeric(d$SumStochElas)
phylosig(tree, d$SumStochElas)
#Phylogenetic signal K : 0.225463 

#PGLS
d$logMPMs <- log(d$MPMs)
buffering <- comparative.data(tree, d, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pgls(SumStochElas ~ MPMs, buffering, lambda='ML')
summary(mod1)
#kappa  [Fix]  : 1.000
#lambda [ ML]  : 0.000
#lower bound : 0.000, p = 1    
#upper bound : 1.000, p = 0.068347
#95.0% CI   : (NA, NA)
#delta  [Fix]  : 1.000
#Estimate Std. Error t value Pr(>|t|)    
#MPMs        -0.0026128  0.0010336 -2.5279  0.01678 *  

mod2 <- pgls(SumStochElas ~ logMPMs, buffering, lambda='ML')
summary(mod2)










