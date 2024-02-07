#################################################################
# 		GABRIEL SANTOS' THESIS - CHAPTER 1
# To buffer or to be labile? A framework 1 to disentangle demographic 
#		patterns and evolutionary processes
#
#		wrote by Gabriel Santos 07 Jul 2022 
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
require(reshape2)

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
source("Functions-Automatization.R")

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

#=================================================================================#
#		DEMOGRAPHIC BUFFERING-LABILITY GRAPHIC
#=================================================================================#

DB.DL.all%>%filter(byVR == "TRUE")
# phylopics
# orca<- image_data("880129b5-b78b-40a9-88ad-55f7d1dc823f", size = "512")[[1]]
# gorilla<-image_data("d9af529d-e426-4c7a-922a-562d57a7872e", size = "512")[[1]]	#Source:http://phylopic.org/image/d9af529d-e426-4c7a-922a-562d57a7872e/
# Alce<-image_data("1a20a65d-1342-4833-a9dd-1611b9fb383c", size = "512")[[1]]	#Source:http://phylopic.org/image/1a20a65d-1342-4833-a9dd-1611b9fb383c/
# Vulpes<-image_data("d67d3bf6-3509-4ab6-819a-cd409985347e", size = "512")[[1]]	#Source:http://phylopic.org/image/d67d3bf6-3509-4ab6-819a-cd409985347e/
# Ursus.americanus<-image_data("5a5dafa2-6388-43b8-a15a-4fd21cd17594", size = "512")[[1]]	#Source:http://phylopic.org/image/5a5dafa2-6388-43b8-a15a-4fd21cd17594/
# Mustela<-image_data("592c3569-ccf3-4f11-920d-3765aa12343f", size = "512")[[1]]	#Source:http://phylopic.org/image/592c3569-ccf3-4f11-920d-3765aa12343f/
# Marmota<-image_data("eee50efb-40dc-47d0-b2cb-52a14a5e0e51", size = "512")[[1]]	#Source:http://phylopic.org/image/eee50efb-40dc-47d0-b2cb-52a14a5e0e51/
# Antechinus<-image_data("295cd9f7-eef2-441e-ba7e-40c772ca7611", size = "512")[[1]]	#Source:http://phylopic.org/image/295cd9f7-eef2-441e-ba7e-40c772ca7611/
# Spermophilus<-image_data("8de61ee7-89eb-49c0-85b7-bc25956544bc", size = "512")[[1]]	#Source:http://phylopic.org/image/8de61ee7-89eb-49c0-85b7-bc25956544bc/
# Cebus<-image_data("156b515d-f25c-4497-b15b-5afb832cc70c", size = "512")[[1]]	#Source:http://phylopic.org/image/156b515d-f25c-4497-b15b-5afb832cc70c/
# Macropus<-image_data("006f91fa-e49f-43f6-a02b-97c6d7b9178a", size = "512")[[1]]	#Source:http://phylopic.org/image/006f91fa-e49f-43f6-a02b-97c6d7b9178a/


plot.db.dl<-ggplot(DB.DL.all,aes(x=(E_Smean-1),y= E_Smean-1,color=Order))+
  #geom_point(aes(color=Order,shape=byVR,size=n),alpha=.7)+
  geom_point(aes(color=Order),alpha=.7,size=3)+
  geom_errorbar(aes(ymin =(E_Smean-E_Smean.SE)-1, ymax =(E_Smean+E_Smean.SE)-1),alpha=.7)+
  geom_errorbarh(aes(xmin = ((E_Smean-E_Smean.SE)-1), xmax = ((E_Smean+E_Smean.SE)-1)),alpha=.7)+
  scale_color_brewer(palette="Spectral")+
  scale_x_continuous(trans = log_trans(10),
                     breaks=c(0,0.0001,0.001,0.001,0.01,0.1,0.2,0.4),
                     labels=c(1,1.0001,1.001,1.001,1.01,1.1,1.2,1.4),expand = c(.1, .001))+
  scale_y_continuous(trans = reverselog_trans(10),
                     breaks=c(0,0.0001,0.001,0.001,0.01,0.1,0.2,0.4),
                     labels=c(0,"-0.0001",-0.001,-0.001,-0.01,-0.1,-0.2,-0.4),expand = c(.1, .001))+
  annotation_logticks() + 
  #ylab(bquote("Relatative effect of environemtal variation ("~sum(E[a[ij]]^S^""[sigma]~")")))+
  #xlab(bquote("Relatative effect of environemtal variation ("~sum(E[a[ij]]^S^""[mu]~")")))+
  ylab(expression("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +")) +
  xlab(expression("(More Buffered)            -       " %<-% "        "~Sigma~"E"^"s"^~mu~"        " %->%  "       +            (More Labile)")) +
  #scale_shape_manual(values=c(1,19),labels=c("Yes","No"))+
  theme_bw(base_family = "Times",base_size = 14)+
  #labs(Size= "#Matrices",color = "Order", Shape="")+
  guides(
    color = guide_legend(order = 1),
    colour = guide_legend(title="Order"),
    #   shape = guide_legend(order = 2),
    #	colour = guide_legend(title="#Vital rates"),
    size = guide_legend(order = 3),
    size= guide_legend(title="#Matrices"))+
  labs(
    #shape= "Vital rates",
    colour = "Order",
    size="# Matrices")

x11(9,4);plot.db.dl
#ggsave("DB-DL continuum.svg", width = 25, height = 15, units = "cm")
#dev.off()


#=================================================================================#
# SECOND DERIVATIVES
#=================================================================================#

Mammals.matrizes.all.SDs<-lapply(Mammals.all$IDs,getAs.from.IDlist)

outputs2<-outputs1<-list()

#------------------------------------------
# Updated Max_var_scale
#------------------------------------------

Max_var_scale2<-function(X){
  varmxcorrected<-(splitA(var2(X))$T*
                     (splitA(mean(X))$T*splitA(1-mean(X))$T))+
    splitA(var2(X))$F
  CV<-sqrt(varmxcorrected)/mean(X)*100
  return(list(A=varmxcorrected, CV=CV))
}


#------------------------------------------------------------------------------------
# Parte 1
#------------------------------------------------------------------------------------
for (i in 1:dim(Mammals.all)[1]){
  print(i)
  tryCatch({
    mean.mx<-mean(getAs.from.IDlist(Mammals.all$IDs[i]))
    outputs1[[i]]<-data.frame(
      rows = paste("a", 1:dim(mean.mx)[1], rep(1:dim(mean.mx)[1], each = dim(mean.mx)[1]), sep = ""),
      rows2 = paste("col.", 1:dim(mean.mx)[1]," / line.", rep(1:dim(mean.mx)[1], each = dim(mean.mx)[1]), sep = ""),
      age=paste(rep(1:dim(mean.mx)[1], each = dim(mean.mx)[1]), sep = ""),
      means=as.vector(mean.mx),
      Max.vars=as.vector(Max_var_scale2(getAs.from.IDlist(Mammals.all$IDs[i]))$CV),
      CV=as.vector(sqrt(var2(getAs.from.IDlist(Mammals.all$IDs[i])))/mean(getAs.from.IDlist(Mammals.all$IDs[i]))*100),
      SpeciesAuthor=Mammals.all$SpeciesAuthor[i],
      n=Mammals.all$n[i],
      IDs=Mammals.all$IDs[i])
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #mean.mx<-NULL
}
#------------------------------------------------------------------------------------
# Parte 2
#------------------------------------------------------------------------------------
for (i in 1:dim(Mammals.all)[1]){
  print(i)
  tryCatch({
    mean.mx<-mean(getAs.from.IDlist(Mammals.all$IDs[i]))
    outputs2[[i]]<-data.frame(
      rows = paste("a", 1:dim(mean.mx)[1], rep(1:dim(mean.mx)[1], each = dim(mean.mx)[1]), sep = ""),
      rows2 = paste("col.", 1:dim(mean.mx)[1]," / line.", rep(1:dim(mean.mx)[1], each = dim(mean.mx)[1]), sep = ""),
      Sensitivity = as.vector(sensitivity(mean.mx)),
      Elasticity = as.vector(elasticity( mean.mx)),
      Sec.derivative = as.vector(mysec(mean.mx)),
      Meanings = as.vector(Bio_meaning(mean.mx)),
      SpeciesAuthor = Mammals.all$SpeciesAuthor[i])
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  mean.mx<-NULL
}

SDs<-left_join(
  do.call(rbind,outputs2),do.call(rbind,outputs1),
  by=c("rows2","SpeciesAuthor"))

 # SDs<-SDs%>%
 #   filter(means>0)%>%
 #   drop_na(.)%>%
 #   filter(SpeciesAuthor %in% Mammals.sub$SpeciesAuthor)
 # 
 # SDs$Sec.derivative.bin<-SDs$Sec.derivative
 # SDs$Sec.derivative.bin[SDs$Sec.derivative < -1]=-1
 # SDs$Sec.derivative.bin[SDs$Sec.derivative > 1]=1
 # SDs$age<-as.numeric(as.character(SDs$age))
 # 
 # SDs$Sec.derivative.bin<-SDs$Sec.derivative
 # SDs$Sec.derivative.bin[SDs$Sec.derivative < -1]=-1
 # SDs$Sec.derivative.bin[SDs$Sec.derivative > 1]=1
 # SDs$age<-as.numeric(as.character(SDs$age))
###############################################
###############################################
###############################################

#Elasticities (Maja)

#par(mfrow=c(4,4))
 
#Blue Monkey

BM <- SDs %>%
  filter(SDs$SpeciesAuthor=="Cercopithecus_mitis")

elasticity_matrixBM <- matrix(BM$Elasticity, nrow = 9, byrow = TRUE)
elasticity_matrixBM
datagBM <- melt(elasticity_matrixBM)

hessian_matrixBM <- matrix(BM$Sec.derivative, nrow = 9, byrow = TRUE)
hessian_matrixBM
datagBM3 <- melt(hessian_matrixBM)

datagBM3$Elasticity <- NA
datagBM3$Elasticity <- datagBM$value
datagBM3 <- datagBM3 %>% rename(Sec.derivative = value)

datagBM3$absSecDer <- abs(datagBM3$Sec.derivative)
datagBM3 <- datagBM3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gBM3 <- ggplot(datagBM3, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) +
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.55)) +
  labs(x = "", y = "", title = "Blue monkey") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_reverse(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  geom_point(data = subset(datagBM3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(1, 2, 3, 4), limits = c(0, 4))

gBM3


CGS <- SDs %>%
  filter(SDs$SpeciesAuthor=="Spermophilus_columbianus_3")

elasticity_matrixCGS <- matrix(CGS$Elasticity, nrow = 9, byrow = TRUE)
elasticity_matrixCGS
datagCGS <- melt(elasticity_matrixCGS)

hessian_matrixCGS <- matrix(CGS$Sec.derivative, nrow = 9, byrow = TRUE)
hessian_matrixCGS
datagCGS3 <- melt(hessian_matrixCGS)

datagCGS3$Elasticity <- NA
datagCGS3$Elasticity <- datagCGS$value
datagCGS3 <- datagCGS3 %>% rename(Sec.derivative = value)

datagCGS3$absSecDer <- abs(datagCGS3$Sec.derivative)
datagCGS3 <- datagCGS3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gCGS3 <- ggplot(datagCGS3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.24)) +
  labs(x = "", y = "", title = "Columbian ground squirrel") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8,9), labels = c(1,2,3,4,5,6,7,8,9)) +
  geom_point(data = subset(datagCGS3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 10), breaks = c(0, 0.5, 1, 1.5), limits = c(0, 1.5))

gCGS3

#Eastern chimpanzee

EC <- SDs %>%
  filter(SDs$SpeciesAuthor=="Pan_troglodytes_subsp._schweinfurthii")

elasticity_matrixEC <- matrix(EC$Elasticity, nrow = 17, byrow = TRUE)
elasticity_matrixEC
datagEC <- melt(elasticity_matrixEC)

hessian_matrixEC <- matrix(EC$Sec.derivative, nrow = 17, byrow = TRUE)
hessian_matrixEC
datagEC3 <- melt(hessian_matrixEC)

datagEC3$Elasticity <- NA
datagEC3$Elasticity <- datagEC$value
datagEC3 <- datagEC3 %>% rename(Sec.derivative = value)

datagEC3$absSecDer <- abs(datagEC3$Sec.derivative)
datagEC3 <- datagEC3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gEC3 <- ggplot(datagEC3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Eastern chimpanzee") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)) +
  geom_point(data = subset(datagEC3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 6), breaks = c(0, 1, 2, 3), limits = c(0, 4.5))

gEC3

#Homo sapiens

HS <- SDs %>%
  filter(SDs$SpeciesAuthor=="Homo_sapiens_subsp._sapiens")

elasticity_matrixHS <- matrix(HS$Elasticity, nrow = 12, byrow = TRUE)
elasticity_matrixHS
datagHS <- melt(elasticity_matrixHS)

hessian_matrixHS <- matrix(HS$Sec.derivative, nrow = 12, byrow = TRUE)
hessian_matrixHS
datagHS3 <- melt(hessian_matrixHS)

datagHS3$Elasticity <- NA
datagHS3$Elasticity <- datagHS$value
datagHS3 <- datagHS3 %>% rename(Sec.derivative = value)

datagHS3$absSecDer <- abs(datagHS3$Sec.derivative)
datagHS3 <- datagHS3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gHS3 <- ggplot(datagHS3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.2)) +
  labs(x = "", y = "", title = "Human") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  geom_point(data = subset(datagHS3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 10), breaks = c(0.04, 0.06, 0.08, 0.1))

gHS3


#Killer whale

KW <- SDs %>%
  filter(SDs$SpeciesAuthor=="Orcinus_orca_2")

elasticity_matrixKW <- matrix(KW$Elasticity, nrow = 7, byrow = TRUE)
elasticity_matrixKW
datagKW <- melt(elasticity_matrixKW)

hessian_matrixKW <- matrix(KW$Sec.derivative, nrow = 7, byrow = TRUE)
hessian_matrixKW
datagKW3 <- melt(hessian_matrixKW)

datagKW3$Elasticity <- NA
datagKW3$Elasticity <- datagKW$value
datagKW3 <- datagKW3 %>% rename(Sec.derivative = value)

datagKW3$absSecDer <- abs(datagKW3$Sec.derivative)
datagKW3 <- datagKW3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gKW3 <- ggplot(datagKW3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Killer whale") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7), labels = c(1,2,3,4,5,6,7)) +
  geom_point(data = subset(datagKW3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 10), breaks = c(1, 2, 3, 4))

gKW3

#Moose

MS <- SDs %>%
  filter(SDs$SpeciesAuthor=="Alces_alces")

elasticity_matrixMS <- matrix(MS$Elasticity, nrow = 3, byrow = TRUE)
elasticity_matrixMS
datagMS <- melt(elasticity_matrixMS)

hessian_matrixMS <- matrix(MS$Sec.derivative, nrow = 3, byrow = TRUE)
hessian_matrixMS
datagMS3 <- melt(hessian_matrixMS)

datagMS3$Elasticity <- NA
datagMS3$Elasticity <- datagMS$value
datagMS3 <- datagMS3 %>% rename(Sec.derivative = value)

datagMS3$absSecDer <- abs(datagMS3$Sec.derivative)
datagMS3 <- datagMS3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gMS3 <- ggplot(datagMS3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0.14, 0.56)) +
  labs(x = "", y = "", title = "Moose") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3")) +
  scale_y_reverse(breaks = c(1,2,3), labels = c(1,2,3)) +
  geom_point(data = subset(datagMS3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 10), breaks = c(0.1, 0.2, 0.3, 0.4))

gMS3

#Mountain gorilla

MG <- SDs %>%
  filter(SDs$SpeciesAuthor=="Gorilla_beringei_subsp._beringei")

elasticity_matrixMG <- matrix(MG$Elasticity, nrow = 3, byrow = TRUE)
elasticity_matrixMG
datagMG <- melt(elasticity_matrixMG)

hessian_matrixMG <- matrix(MG$Sec.derivative, nrow = 3, byrow = TRUE)
hessian_matrixMG
datagMG3 <- melt(hessian_matrixMG)

datagMG3$Elasticity <- NA
datagMG3$Elasticity <- datagMG$value
datagMG3 <- datagMG3 %>% rename(Sec.derivative = value)

datagMG3$absSecDer <- abs(datagMG3$Sec.derivative)
datagMG3 <- datagMG3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gMG3 <- ggplot(datagMG3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.9)) +
  labs(x = "", y = "", title = "Mountain gorilla") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3")) +
  scale_y_reverse(breaks = c(1,2,3), labels = c(1,2,3)) +
  geom_point(data = subset(datagMG3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0, 0.2, 0.4, 0.6, 0.8))

gMG3

#Northern muriqui

NM <- SDs %>%
  filter(SDs$SpeciesAuthor=="Brachyteles_hypoxanthus_2")

elasticity_matrixNM <- matrix(NM$Elasticity, nrow = 3, byrow = TRUE)
elasticity_matrixNM
datagNM <- melt(elasticity_matrixNM)

hessian_matrixNM <- matrix(NM$Sec.derivative, nrow = 3, byrow = TRUE)
hessian_matrixNM
datagNM3 <- melt(hessian_matrixNM)

datagNM3$Elasticity <- NA
datagNM3$Elasticity <- datagNM$value
datagNM3 <- datagNM3 %>% rename(Sec.derivative = value)

datagNM3$absSecDer <- abs(datagNM3$Sec.derivative)
datagNM3 <- datagNM3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gNM3 <- ggplot(datagNM3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.8)) +
  labs(x = "", y = "", title = "Northern muriqui") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3")) +
  scale_y_reverse(breaks = c(1,2,3), labels = c(1,2,3)) +
  geom_point(data = subset(datagNM3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0, 0.2, 0.4, 0.6, 0.8))

gNM3

#Olive baboon

OB <- SDs %>%
  filter(SDs$SpeciesAuthor=="Papio_cynocephalus")

elasticity_matrixOB <- matrix(OB$Elasticity, nrow = 8, byrow = TRUE)
elasticity_matrixOB
datagOB <- melt(elasticity_matrixOB)

hessian_matrixOB <- matrix(OB$Sec.derivative, nrow = 8, byrow = TRUE)
hessian_matrixOB
datagOB3 <- melt(hessian_matrixOB)

datagOB3$Elasticity <- NA
datagOB3$Elasticity <- datagOB$value
datagOB3 <- datagOB3 %>% rename(Sec.derivative = value)

datagOB3$absSecDer <- abs(datagOB3$Sec.derivative)
datagOB3 <- datagOB3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gOB3 <- ggplot(datagOB3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Olive baboon") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8), labels = c(1,2,3,4,5,6,7,8)) +
  geom_point(data = subset(datagOB3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0, 0.2, 0.4, 0.6, 0.8))

gOB3

#Polar bear

PB <- SDs %>%
  filter(SDs$SpeciesAuthor=="Ursus_maritimus_2")

elasticity_matrixPB <- matrix(PB$Elasticity, nrow = 6, byrow = TRUE)
elasticity_matrixPB
datagPB <- melt(elasticity_matrixPB)

hessian_matrixPB <- matrix(PB$Sec.derivative, nrow = 6, byrow = TRUE)
hessian_matrixPB
datagPB3 <- melt(hessian_matrixPB)

datagPB3$Elasticity <- NA
datagPB3$Elasticity <- datagPB$value
datagPB3 <- datagPB3 %>% rename(Sec.derivative = value)

datagPB3$absSecDer <- abs(datagPB3$Sec.derivative)
datagPB3 <- datagPB3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gPB3 <- ggplot(datagPB3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.28)) +
  labs(x = "", y = "", title = "Polar bear") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6), labels = c(1,2,3,4,5,6)) +
  geom_point(data = subset(datagPB3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.1, 0.2, 0.3, 0.4, 0,5))

gPB3

#Rhesus macaque

RM <- SDs %>%
  filter(SDs$SpeciesAuthor=="Macaca_mulatta_3")

elasticity_matrixRM <- matrix(RM$Elasticity, nrow = 5, byrow = TRUE)
elasticity_matrixRM
datagRM <- melt(elasticity_matrixRM)

hessian_matrixRM <- matrix(RM$Sec.derivative, nrow = 5, byrow = TRUE)
hessian_matrixRM
datagRM3 <- melt(hessian_matrixRM)

datagRM3$Elasticity <- NA
datagRM3$Elasticity <- datagRM$value
datagRM3 <- datagRM3 %>% rename(Sec.derivative = value)

datagRM3$absSecDer <- abs(datagRM3$Sec.derivative)
datagRM3 <- datagRM3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))


gRM3 <- ggplot(datagRM3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Rhesus macaque") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5")) +
  scale_y_reverse(breaks = c(1,2,3,4,5), labels = c(1,2,3,4,5)) +
  geom_point(data = subset(datagRM3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.1, 0.2, 0.3, 0.4, 0,5))

gRM3

#Root vole

RV <- SDs %>%
  filter(SDs$SpeciesAuthor=="Microtus_oeconomus")

elasticity_matrixRV <- matrix(RV$Elasticity, nrow = 3, byrow = TRUE)
elasticity_matrixRV
datagRV <- melt(elasticity_matrixRV)

hessian_matrixRV <- matrix(RV$Sec.derivative, nrow = 3, byrow = TRUE)
hessian_matrixRV
datagRV3 <- melt(hessian_matrixRV)

datagRV3$Elasticity <- NA
datagRV3$Elasticity <- datagRV$value
datagRV3 <- datagRV3 %>% rename(Sec.derivative = value)

datagRV3$absSecDer <- abs(datagRV3$Sec.derivative)
datagRV3 <- datagRV3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gRV3 <- ggplot(datagRV3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 1)) +
  labs(x = "", y = "", title = "Root vole") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3")) +
  scale_y_reverse(breaks = c(1,2,3), labels = c(1,2,3)) +
  geom_point(data = subset(datagRV3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.1, 0.2, 0.3, 0.4, 0,5))

gRV3

#Soay sheep

SS <- SDs %>%
  filter(SDs$SpeciesAuthor=="Ovis_aries_2")

elasticity_matrixSS <- matrix(SS$Elasticity, nrow = 6, byrow = TRUE)
elasticity_matrixSS
datagSS <- melt(elasticity_matrixSS)

hessian_matrixSS <- matrix(SS$Sec.derivative, nrow = 6, byrow = TRUE)
hessian_matrixSS
datagSS3 <- melt(hessian_matrixSS)

datagSS3$Elasticity <- NA
datagSS3$Elasticity <- datagSS$value
datagSS3 <- datagSS3 %>% rename(Sec.derivative = value)

datagSS3$absSecDer <- abs(datagSS3$Sec.derivative)
datagSS3 <- datagSS3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gSS3 <- ggplot(datagSS3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Soay sheep") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6), labels = c(1,2,3,4,5,6)) +
  geom_point(data = subset(datagSS3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.1, 0.15, 0.2, 0.25, 0.3))

gSS3

#Tammar wallaby

TW <- SDs %>%
  filter(SDs$SpeciesAuthor=="Macropus_eugenii")

elasticity_matrixTW <- matrix(TW$Elasticity, nrow = 2, byrow = TRUE)
elasticity_matrixTW
datagTW <- melt(elasticity_matrixTW)

hessian_matrixTW <- matrix(TW$Sec.derivative, nrow = 2, byrow = TRUE)
hessian_matrixTW
datagTW3 <- melt(hessian_matrixTW)

datagTW3$Elasticity <- NA
datagTW3$Elasticity <- datagTW$value
datagTW3 <- datagTW3 %>% rename(Sec.derivative = value)

datagTW3$absSecDer <- abs(datagTW3$Sec.derivative)
datagTW3 <- datagTW3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gTW3 <- ggplot(datagTW3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.6)) +
  labs(x = "", y = "", title = "Tammar wallaby") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2")) +
  scale_y_reverse(breaks = c(1,2), labels = c(1,2)) +
  geom_point(data = subset(datagTW3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))

gTW3

#Verreaux's sifaka

VS <- SDs %>%
  filter(SDs$SpeciesAuthor=="Propithecus_verreauxi")

elasticity_matrixVS <- matrix(VS$Elasticity, nrow = 8, byrow = TRUE)
elasticity_matrixVS
datagVS <- melt(elasticity_matrixVS)

hessian_matrixVS <- matrix(VS$Sec.derivative, nrow = 8, byrow = TRUE)
hessian_matrixVS
datagVS3 <- melt(hessian_matrixVS)

datagVS3$Elasticity <- NA
datagVS3$Elasticity <- datagVS$value
datagVS3 <- datagVS3 %>% rename(Sec.derivative = value)

datagVS3$absSecDer <- abs(datagVS3$Sec.derivative)
datagVS3 <- datagVS3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gVS3 <- ggplot(datagVS3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.7)) +
  labs(x = "", y = "", title = "Verreaux's sifaka") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8")) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8), labels = c(1,2,3,4,5,6,7,8)) +
  geom_point(data = subset(datagVS3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.2, 0.4, 0.6, 0.8))

gVS3

#White faced capuchin monkey

CM <- SDs %>%
  filter(SDs$SpeciesAuthor=="Cebus_capucinus_2")

elasticity_matrixCM <- matrix(CM$Elasticity, nrow = 3, byrow = TRUE)
elasticity_matrixCM
datagCM <- melt(elasticity_matrixCM)

hessian_matrixCM <- matrix(CM$Sec.derivative, nrow = 3, byrow = TRUE)
hessian_matrixCM
datagCM3 <- melt(hessian_matrixCM)

datagCM3$Elasticity <- NA
datagCM3$Elasticity <- datagCM$value
datagCM3 <- datagCM3 %>% rename(Sec.derivative = value)

datagCM3$absSecDer <- abs(datagCM3$Sec.derivative)
datagCM3 <- datagCM3 %>% 
  mutate(new_var = ifelse(Sec.derivative == 0, NA, ifelse(Sec.derivative > 0, "P", "N")))

gCM3 <- ggplot(datagCM3, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill = ifelse(Elasticity > 0, Elasticity, NA))) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3, na.value = "transparent", limits = c(0, 0.8)) +
  labs(x = "", y = "", title = "White faced capuchin monkey") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 20, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, family = "Arial"),
    plot.title = element_text(size = 20, family = "Arial"),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_text(size = 20, family = "Arial")
  ) +
  scale_x_discrete(limits = c("1", "2", "3")) +
  scale_y_reverse(breaks = c(1,2,3), labels = c(1,2,3)) +
  geom_point(data = subset(datagCM3, Elasticity > 0), aes(x = Var1, y = Var2, size = absSecDer, color = as.factor(new_var))) +
  scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE) +
  scale_size_continuous(name = "|Self-second derivative|", range = c(1, 15), breaks = c(0.2, 0.4, 0.6, 0.8))

gCM3

max(datagBM3$Elasticity)
min(datagBM3$Elasticity)
max(datagBM3$Sec.derivative)
min(datagBM3$Sec.derivative)

max(datagCGS3$Elasticity)
min(datagCGS3$Elasticity)
max(datagCGS3$Sec.derivative)
min(datagCGS3$Sec.derivative)

max(datagEC3$Elasticity)
min(datagEC3$Elasticity)
max(datagEC3$Sec.derivative)
min(datagEC3$Sec.derivative)

max(datagHS3$Elasticity)
min(datagHS3$Elasticity)
max(datagHS3$Sec.derivative)
min(datagHS3$Sec.derivative)

max(datagKW3$Elasticity)
min(datagKW3$Elasticity)
max(datagKW3$Sec.derivative)
min(datagKW3$Sec.derivative)

max(datagMS3$Elasticity)
min(datagMS3$Elasticity)
max(datagMS3$Sec.derivative)
min(datagMS3$Sec.derivative)

max(datagMG3$Elasticity)
min(datagMG3$Elasticity)
max(datagMG3$Sec.derivative)
min(datagMG3$Sec.derivative)

max(datagNM3$Elasticity)
min(datagNM3$Elasticity)
max(datagNM3$Sec.derivative)
min(datagNM3$Sec.derivative)

max(datagOB3$Elasticity)
min(datagOB3$Elasticity)
max(datagOB3$Sec.derivative)
min(datagOB3$Sec.derivative)

max(datagPB3$Elasticity)
min(datagPB3$Elasticity)
max(datagPB3$Sec.derivative)
min(datagPB3$Sec.derivative)

max(datagRM3$Elasticity)
min(datagRM3$Elasticity)
max(datagRM3$Sec.derivative)
min(datagRM3$Sec.derivative)

max(datagRV3$Elasticity)
min(datagRV3$Elasticity)
max(datagRV3$Sec.derivative)
min(datagRV3$Sec.derivative)

max(datagSS3$Elasticity)
min(datagSS3$Elasticity)
max(datagSS3$Sec.derivative)
min(datagSS3$Sec.derivative)

max(datagTW3$Elasticity)
min(datagTW3$Elasticity)
max(datagTW3$Sec.derivative)
min(datagTW3$Sec.derivative)

max(datagVS3$Elasticity)
min(datagVS3$Elasticity)
max(datagVS3$Sec.derivative)
min(datagVS3$Sec.derivative)

max(datagCM3$Elasticity)
min(datagCM3$Elasticity)
max(datagCM3$Sec.derivative)
min(datagCM3$Sec.derivative)


# ordered_SDs <- arrange(SDs, Sec.derivative)
# 
# install.packages("gridExtra")
# library(gridExtra)
# library(cowplot)
# Fig3 <- plot_grid(gBM3,gCGS3, gEC3, gHS3, gKW3, gMS3, gMG3, gNM3, gOB3, gPB3,gRM3, gRV3, gSS3, gTW3, gVS3, gCM3, ncol = 4)
# Fig3


# library(ggplot2)
# library(reshape2)
# library(viridis)
# 
# # Melt the data
# melted_data <- melt(SDs.all, id.vars = c("Elasticity", "Max.vars", "Meanings", "Sec.derivative.bin"),
#                     variable.name = "Var1", value.name = "Var2")
# melted_data <- melted_data[!melted_data$Elasticity==0,]
# 
# # Plot the data as a raster plot
# g <- ggplot(melted_data, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill = Elasticity)) +
#   scale_fill_viridis(name = "Elasticity", begin = 0.3) +
#   labs(x = "", y = "", title = "") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(size = 15, angle = 0, vjust = 0.3),
#     axis.text.y = element_text(size = 15),
#     plot.title = element_text(size = 20)
#   ) +
#   scale_x_discrete(limits = c("1", "2", "3", "4")) +
#   scale_y_reverse()
# 
# # Plot the rest of the aesthetics
# g <- g +
#   geom_point(
#     data = SDs.all,
#     aes(x = Elasticity + 0.001, y = Max.vars, size = abs(Sec.derivative),
#         color = Sec.derivative.bin, shape = Meanings),
#     position = position_jitter(width = 0.1, height = 0.1),
#     alpha = 0.5
#   ) +
#   labs(size = "|Second derivative|", color = "Natural selection \n     linearity", shape = "Vital rate") +
#   scale_size(range = c(3, 6), breaks = c(0, 1, 2), labels = c("0", "1", "2+")) +
#   scale_shape_manual(values = c(17, 19), labels = c("Recruitment", "Stages transition")) +
#   scale_x_log10(breaks = c(0, 0.001, 0.01, 0.1, 1),
#                 labels = c(0, 0.001, 0.01, 0.1, 1),
#                 limits = c(0.0005, 1.1),
#                 expand = c(0.1, 0.1)) +
#   scale_y_sqrt(breaks = c(0, 0.005, 0.007, 0.01, 10, 50, 100, 200, 300),
#                labels = c(0, 0, "", "", 10, 50, 100, 200, 300),
#                expand = c(0.1, 0.1)) +
#   scale_color_manual(
#     values = c("skyblue", "grey", "tomato"),
#     labels = c("Stabilizing selection", "Non-informative", "Disruptive selection")
#   ) +
#   facet_wrap(. ~ CommonName + SpeciesAuthor, nrow = 5) +
#   labs(size = "|Second derivative|", color = "Natural selection \n     linearity", shape = "Vital rate") +
#   theme(
#     text = element_text(size = 14),
#     strip.text = element_text(size = 7),
#     axis.text.x = element_text(angle = 70, hjust = 1)
#   )
# 
# # Display the plot
# g





#----------------------------------------------------------------------
#		CHECKING AGE PATTERNS
#----------------------------------------------------------------------
#SDs$Sec.derivative.bin<-round(SDs$Sec.derivative,0)
SDs$Sec.derivative.bin<-SDs$Sec.derivative
SDs$Sec.derivative.bin[SDs$Sec.derivative < -1]=-1
SDs$Sec.derivative.bin[SDs$Sec.derivative > 1]=1
SDs$age<-as.numeric(as.character(SDs$age))

unique(SDs$SpeciesAuthor)

SDs.A<-ggplot(SDs,aes(y=Sec.derivative.bin,x=age))+
  geom_point()+
  ylim(-1,1)+
  xlab("Age")+
  ylab("Second derivative")+
  geom_hline(yintercept = 0)+
  facet_grid(Meanings~SpeciesAuthor,scale="free")

SDs.A

#----------------------------------------------------------------------
#		RESULT
#----------------------------------------------------------------------
SDs$Sec.derivative.bin<-round(SDs$Sec.derivative.bin,1)
SDs$Sec.derivative.bin[SDs$Sec.derivative.bin < 0]=-1
SDs$Sec.derivative.bin[SDs$Sec.derivative.bin > 0]=1
SDs$Sec.derivative.bin<-as.factor(SDs$Sec.derivative.bin)
SDs$age<-as.numeric(as.character(SDs$age))

#Remove duplicated data


SDs<-comadre$metadata%>%
  filter(SpeciesAuthor %in% unique(SDs$SpeciesAuthor))%>%
  select(c(SpeciesAuthor,SpeciesAccepted,CommonName))%>%
  unique()%>% 
  mutate(SpeciesAccepted = str_replace(SpeciesAccepted, "[_]"," "))%>%
  mutate(SpeciesAccepted = str_replace(SpeciesAccepted, "_subsp.*",""))%>%
  separate(SpeciesAccepted, c("A", "B"),sep=" ",remove = TRUE)%>%
  unite(NewSpeciesAccepted,c("A","B"),sep=" ",remove = TRUE)%>%
  left_join(SDs,.,by=c("SpeciesAuthor"))%>%
  filter(!(SpeciesAuthor %in% unique(SDs$SpeciesAuthor)[c(2,6,8)]))


plot.label<-unique(SDs$CommonName)
names(plot.label)<-unique(SDs$SpeciesAuthor)

plot.SDs<-SDs%>%
  filter( !SpeciesAuthor %in% 
            c("Callospermophilus_lateralis"))%>%
  ggplot(aes(x=Elasticity, y=Max.vars,label=Meanings,color=Sec.derivative.bin))+
  # geom_point(aes(x=Elasticity+0.001,y=Max.vars, size = abs(Sec.derivative),color=Sec.derivative.bin,shape=Meanings),alpha=.6)+
  geom_point(aes(x=Elasticity+0.001,y=Max.vars, size = abs(Sec.derivative),color=Sec.derivative.bin,shape=Meanings),
             position = position_jitter(width = 0.1, height = 0.1),alpha=.5)+
  #geom_abline(intercept = 250 , slope = -25)+
  #geom_smooth()+
  #geom_text(aes(x=elas_Det+0.001,y=CVcor+0.01,fill=Kind),size=3)+
  #scale_color_gradient2(midpoint=0, limits=c(-1.05,1.05),low="skyblue", mid="grey",
  #                 high="tomato", space ="Lab")+
  #labs(colour = " Second \nderivative")+
  ylab("Coefficient of Variation")+
  theme_bw()+
  scale_size(range=c(3,6),breaks=c(0,1,2),labels=c("0","1","2+"))+
  scale_shape_manual(values=c(17,19),labels=c("Recruitment","Survival"))+
  #scale_x_log10(breaks = c(0,0.001,.01,.1,1),labels = c(0,0.001,0.01,.1,1),limits=c(0.005,1.1),expand = c(.1, .1)) +
  scale_x_log10(breaks = c(0,0.001,.01,.1,1),labels = c(0,0.001,0.01,.1,1),limits=c(0.0005,1.1),expand = c(.1, .1)) +
  #scale_y_sqrt(breaks = c(0,0.01,10,50,100,200,300),labels = c(0,0,10,50,100,200,300),expand = c(.1, .1))+
  scale_y_sqrt(breaks = c(0,0.005,0.007,0.01,10,50,100,200,300),labels = c(0,0,"","",10,50,100,200,300),expand = c(.1, .1)) +
  scale_color_manual(values=c("skyblue", "grey", "tomato"),labels=c("Stabilizing selection","Non-informative","Disruptive selection"))+
  labs(size= "|Second derivative|",color = "Natural selection \n     linearity", shape="Vital rate")+
  facet_wrap(.~CommonName+" ")+
  theme( text=element_text(size=14),
         strip.text = element_text(size=9),
         axis.text.x = element_text(angle = 70, hjust = 1))

x11(20,15);plot.SDs
plot.SDs
#ggsave("Second derivatives MAIN.svg", width = 30, height = 20, units = "cm")
dev.off()


#---------------------------------------------------------------------
#	PLOTTING SECOND DERIVATIVES FOR ALL
#---------------------------------------------------------------------

SDs.all<-SDs%>%
  filter(means>0)%>%
  drop_na(.)

SDs.all$Sec.derivative.bin<-SDs.all$Sec.derivative
SDs.all$Sec.derivative.bin[SDs.all$Sec.derivative < -1]=-1
SDs.all$Sec.derivative.bin[SDs.all$Sec.derivative > 1]=1
SDs.all$age<-as.numeric(as.character(SDs.all$age))

SDs.all$Sec.derivative.bin<-round(SDs.all$Sec.derivative.bin,1)
SDs.all$Sec.derivative.bin[SDs.all$Sec.derivative.bin < 0]=-1
SDs.all$Sec.derivative.bin[SDs.all$Sec.derivative.bin > 0]=1
SDs.all$Sec.derivative.bin<-as.factor(SDs.all$Sec.derivative.bin)
SDs.all$age<-as.numeric(as.character(SDs.all$age))

SDs.all<-
  comadre$metadata%>%
  filter(SpeciesAuthor %in% unique(SDs.all$SpeciesAuthor))%>%
  select(c(SpeciesAuthor,SpeciesAccepted,CommonName))%>%
  unique()%>% 
  mutate(SpeciesAccepted = str_replace(SpeciesAccepted, "[_]"," "))%>%
  mutate(SpeciesAccepted = str_replace(SpeciesAccepted, "_subsp.*",""))%>%
  separate(SpeciesAccepted, c("A", "B"),sep=" ",remove = TRUE)%>%
  unite(NewSpeciesAccepted,c("A","B"),sep=" ",remove = TRUE)%>%
  left_join(SDs.all,.,by=c("SpeciesAuthor"))

SDs.all%>%
  ggplot(aes(x=Elasticity, y=Max.vars,label=Meanings,color=Sec.derivative.bin))+
  geom_point(aes(x=Elasticity+0.001,y=Max.vars, size = abs(Sec.derivative),color=Sec.derivative.bin,shape=Meanings),
             position = position_jitter(width = 0.1, height = 0.1),alpha=.5)+
  ylab("Coefficient of Variation")+
  theme_bw()+
  scale_size(range=c(3,6),breaks=c(0,1,2),labels=c("0","1","2+"))+
  scale_shape_manual(values=c(17,19),labels=c("Recruitment","Stages transition"))+
  scale_x_log10(breaks = c(0,0.001,.01,.1,1),labels = c(0,0.001,0.01,.1,1),limits=c(0.0005,1.1),expand = c(.1, .1)) +
  scale_y_sqrt(breaks = c(0,0.005,0.007,0.01,10,50,100,200,300),labels = c(0,0,"","",10,50,100,200,300),expand = c(.1, .1)) +
  scale_color_manual(values=c("skyblue", "grey", "tomato"),labels=c("Stabilizing selection","Non-informative","Disruptive selection"))+
  labs(size= "|Second derivative|",color = "Natural selection \n     linearity", shape="Vital rate")+
  facet_wrap(.~CommonName+SpeciesAuthor,nrow=5)+
  theme( text=element_text(size=14),
         strip.text = element_text(size=7),
         axis.text.x = element_text(angle = 70, hjust = 1))

#ggsave("Second derivatives Supplement.svg", width = 40, height = 45, units = "cm")
dev.off()

