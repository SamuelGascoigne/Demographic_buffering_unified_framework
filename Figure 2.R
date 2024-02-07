#################################################################
# Figure 2 Manuscript version from May 17th 2023
#################################################################

dev.off()
rm(list=ls())

set.seed(1)

require(popbio)
require(ggplot2)
require(scales)
require(demogR)
require(tidyverse)
require(popdemo)
require(magick)

#=======================================================================#
#			DATA LOADING							#
#=======================================================================#

# setwd("C:/Artigos e resumos publicados submetidos ideias/Submetidos/Demographic buffering - framework")

#load("COMADRE_v.2.0.1.RData")
load("COMADRE_v.3.0.0.RData")

#Check if is the currently version
comadre$version

#Mammals at Comadre
comadre$metadata %>%
  as_tibble()%>%
  filter(Class == "Mammalia") %>%
  distinct(SpeciesAccepted)

#=======================================================================#
#		AUTOMATIZATE and AUXILIARY FUNCTIONS				#
#=======================================================================#

#Load functions to perform analyses and automatizate the process

#file.edit("Script and data/Functions-Automatization.R")
source("Functions-Automatization (1).R")

# Main functions required:
# getMatA 		 - Return a list of MatA available per population
# there.is.na 	 - Verify missing value occurrences(Zero doens't count)
# check 		 - An easy way to check the population mean value
# stoch.sens_mean - Stochastic elasticity in respect of mean value
# stoch.sens_sig  - Stochastic elasticity in respect of variance
# array_to_matrix - Transform an array in a list of matrices
# reverselog_trans - Important step to plot negative values in log scale
# Bio_meaning     - Give a biological meaning for each matrix element
# Max_var_scale   - Correct Coefficient of variation by the maximum variance possible
# mysec - Automatize the extraction of the self-second derivatives

#=======================================================================#
#			DATA FILTERING							#
#=======================================================================#

#Giving IDs for each matrix
comadre$metadata$ID<-1:dim(comadre$metadata)[1]

#
Mammals.sub<-comadre$metadata %>%
  filter(MatrixCriteriaAge == "Yes")%>%as.tibble()%>%
  filter(MatrixFec == "Yes")%>%
  filter(Class == "Mammalia") %>%
  filter(MatrixComposite == "Individual") %>%
  #filter(MatrixTreatment == "Unmanipulated")%>%
  group_by(SpeciesAuthor,Order)%>%
  dplyr::summarise(n=n(),IDs=paste(ID, collapse=" "))%>%
  filter(n > 2)%>%print()

Mammals.all<-comadre$metadata %>%
  as.tibble()%>%
  filter(MatrixFec == "Yes")%>%
  filter(Class == "Mammalia") %>%
  filter(MatrixComposite %in% c("Individual","Pooled")) %>%
  #filter(MatrixTreatment == "Unmanipulated")%>%
  group_by(SpeciesAuthor,Order)%>%
  dplyr::summarise(n=n(),IDs=paste(ID, collapse=" "))%>%
  filter(n > 2)%>%print()


#-----------------------------------------------------------
#	NEW WAY TO GET MATRICES
#-----------------------------------------------------------
getAs.from.IDlist<-function(X){
  leng<-as.numeric(str_split(X, " ", n = Inf)[[1]])
  mxs<-NULL
  for(i in 1:length(leng)){
    mxs[[i]]<-comadre$mat[leng[i]][[1]]$matA}
  mxs<-mxs[unlist(there.is.na(mxs))]
  return(mxs)}
#-----------------------------------------------------------


#=================================================================================#
#		DEMOGRAPHIC BUFFERING-LABILITY CONTINUUM CONSTRUCTION
#=================================================================================#

DB.DL.all<-mxs<-NULL

for (i in 1:dim(Mammals.all)[1]){
  print(i)
  tryCatch({
    # mxs<-replicate(100,stoch.sens_mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    #I decide to standardize the usage of stoch.sens_mean to do not imply in bias in the variance
    #Rember that variance is determined by sampled size
    mxs<-replicate(50,stoch.sens_mean(sample(getAs.from.IDlist(Mammals.all$IDs[i]))[1:3]))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    DB.DL.all$E_Smean[i]<-mean(unlist(lapply(array_to_matrix(mxs),sum)))
    DB.DL.all$E_Smean.SD[i]<-sd(1-unlist(lapply(array_to_matrix(mxs),sum)))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    DB.DL.all$Lambda[i]<-lambda(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    DB.DL.all$Generation.time[i]<-generation.time(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    DB.DL.all$resilience.999[i]<- convt(mean(getAs.from.IDlist(Mammals.all$IDs[i])), accuracy=.01)[1]
    DB.DL.all$resilience.90[i]<- convt(mean(getAs.from.IDlist(Mammals.all$IDs[i])), accuracy=.1)[1]
    DB.DL.all$resilience.damping[i]<-damping.ratio(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  mxs<-NULL
}

DB.DL.all<-do.call(cbind.data.frame,DB.DL.all)
DB.DL.all<-cbind(as.data.frame(Mammals.all[,-4]),DB.DL.all)
DB.DL.all$E_Smean.SE<-DB.DL.all$E_Smean.SD/sqrt(DB.DL.all$n)

DB.DL.all$byVR<-DB.DL.all$SpeciesAuthor %in% Mammals.sub$SpeciesAuthor

DB.DL.all%>%
  filter(byVR=="TRUE")


#------------------------------------------------------------------
#			Some stats
#------------------------------------------------------------------
#Maximum value of mean value 
DB.DL.all[DB.DL.all$E_Smean==max(DB.DL.all$E_Smean,na.rm = TRUE),]
#Minimum value
DB.DL.all[DB.DL.all$E_Smean==min(DB.DL.all$E_Smean,na.rm = TRUE),]
#Number of orders
DB.DL.all%>%
  count(Order)

#populations
DB.DL.all%>%as.tibble()
names(DB.DL.all)

#===================================================================================
#			PAPER METADATA
#===================================================================================
metadata<-comadre$metadata[,c(1,2,3,6,13,14,15,16,18,25,34)]%>%
  filter(MatrixComposite == c("Individual","Pooled")) %>%
  filter(SpeciesAuthor %in% DB.DL.all$SpeciesAuthor)%>%
  mutate(NewSpeciesAccepted = str_replace_all(SpeciesAccepted, "_", " "))%>%
  mutate(NewSpeciesAccepted = str_replace_all(NewSpeciesAccepted, " subsp.*", ""))%>%
  group_by(SpeciesAuthor)%>%
  mutate(Matrices.N= n())%>% 
  distinct(SpeciesAuthor,.keep_all=TRUE)%>%
  as.data.frame()%>%print()
names(metadata)

#--------------------------------
etaler<-function(X){
  A<-str_split(X,";")[[1]]
  B<-ifelse(length(A) == 1, print(A[1]), 
            ifelse(length(A) == 2,
                   print(paste0(A[1]," &",A[2])),
                   print(paste0(A[1]," et al."))))}

#--------------------------------

metadata<-metadata%>%
  mutate(Authors = etaler(Authors))

metadata<-metadata%>%
  mutate(Authors = etaler(Authors))

metadata<-left_join(metadata,DB.DL.all,by="SpeciesAuthor")
colnames(metadata)
metadata<-metadata[,c(1,12,3,14,16,17,16,17,13,18,19,22,24)]
metadata[,6]<-metadata[,7]-1
colnames(metadata)<-c("SpeciesAuthorComadre","SpeciesName","CommonName","Order","E_Smu","E_Smu.SD","E_Ssig","E_Ssig.SD","# matrices","Lambda","Generation.time","Damping.ratio","ByAge")
#writeClipboard(knitr::kable(metadata))

head(metadata)

#=================================================================================#
#				SPECIES ANALIZED
#=================================================================================#

#Check that mean matrices below are singular
#	it means, they have at least one determinant equal to 0
#	So, I do not calculate second derivative for them and they should be removed from our analyes
#		FOR THIS REASON THE TOTAL AMOUNT OF POPULATIONS ARE 40 instead of 44
#		and the total amount of species are
DB.DL.all$SpeciesAuthor[c(3,10,19,33)]

comadre$metadata[,c(1,2,3)]%>%
  mutate(NewSpeciesAccepted = str_replace_all(SpeciesAccepted, "_", " "))%>%
  mutate(NewSpeciesAccepted = str_replace_all(NewSpeciesAccepted, " subsp.*", ""))%>%
  filter(SpeciesAuthor %in% DB.DL.all$SpeciesAuthor[-c(3,10,19,33)])%>%
  distinct(NewSpeciesAccepted)%>%print()

metadata%>%
  filter(SpeciesAuthorComadre %in% DB.DL.all$SpeciesAuthor[-c(3,10,19,33)])
head(metadata)
#=================================================================================#


#=================================================================================#
#		DEMOGRAPHIC BUFFERING-LABILITY STATS
#=================================================================================#

DB.DL.sub<-DB.DL.all%>%
  drop_na(.)%>%
  filter(byVR=="TRUE")%>%
  filter(Generation.time != -Inf)%>%print()

DB.DL.all%>%filter(byVR == "TRUE")
#=================================================================================#
#		PLOT FIGURE 2
#=================================================================================#

metadata <- metadata[!is.nan(metadata$E_Ssig), ]

metadata$E_SsigNEG <- NA
metadata$E_SsigNEG <- 1-metadata$E_Ssig

library(viridis)

# Define the color palette
order_colors <- viridis(length(unique(metadata$Order)))

gF2 <- ggplot(metadata, aes(x = E_SsigNEG, fill = as.factor(Order))) +
  geom_density(fill = "purple", alpha = 0.05) +
  geom_point(aes(y = 21), position = position_jitter(height = 10),
             shape = 21,
             size = 2.5,
             stroke = 1.5,
             alpha = 0.6) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_continuous(expand = c(0.002, 0.002)) +
  scale_fill_manual(values = order_colors) +
  xlab(expression(paste("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +"))) +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 18, color = "black")
  )
gF2



ordered_metadata <- metadata[order(metadata$E_SsigNEG, decreasing = TRUE), ]


# plot.db.dl <- ggplot(metadata, aes(x = E_SsigNEG, y = E_Smean - 1, color = Order)) +
#   geom_point(alpha = 0.7, size = 3) +
#   geom_errorbar(aes(ymin = (E_Smean - E_Smean.SE) - 1, ymax = (E_Smean + E_Smean.SE) - 1), alpha = 0.7) +
#   geom_errorbarh(aes(xmin = ((E_Smean - E_Smean.SE) - 1), xmax = ((E_Smean + E_Smean.SE) - 1)), alpha = 0.7) +
#   scale_color_brewer(palette = "Spectral") +
#   scale_x_continuous(
#     trans = log_trans(10),
#     breaks = c(0, 0.0001, 0.001, 0.001, 0.01, 0.1, 0.2, 0.4),
#     labels = c(1, 1.0001, 1.001, 1.001, 1.01, 1.1, 1.2, 1.4),
#     expand = c(0.1, 0.001)
#   ) +
#   scale_y_continuous(
#     trans = reverselog_trans(10),
#     breaks = c(0, 0.0001, 0.001, 0.001, 0.01, 0.1, 0.2, 0.4),
#     labels = c(0, "-0.0001", -0.001, -0.001, -0.01, -0.1, -0.2, -0.4),
#     expand = c(0.1, 0.001)
#   ) +
#   annotation_logticks() +
#   ylab(expression("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +")) +
#   xlab(expression("(More Buffered)            -       " %<-% "        "~Sigma~"E"^"s"^~mu~"        " %->%  "       +            (More Labile)")) +
#   theme_bw(base_family = "Times", base_size = 14) +
#   guides(
#     color = guide_legend(order = 1),
#     size = guide_legend(order = 3)
#   ) +
#   labs(
#     colour = "Order",
#     size = "# Matrices"
#   )
# 
# x11(9, 4)
# plot.db.dl

