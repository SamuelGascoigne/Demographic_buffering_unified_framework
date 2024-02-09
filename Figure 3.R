# Import necessary packages and helper functions.

source("setup.R")
source("functions.R")


# Import data.

load("COMADRE_v.3.0.0.RData")


#Give IDs for each matrix

comadre$metadata$ID<-1:dim(comadre$metadata)[1]


# Filter data

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


# Helper function to extract matrices. 

getAs.from.IDlist<-function(X){
  leng<-as.numeric(str_split(X, " ", n = Inf)[[1]])
  mxs<-NULL
  for(i in 1:length(leng)){
    mxs[[i]]<-comadre$mat[leng[i]][[1]]$matA}
  mxs<-mxs[unlist(there.is.na(mxs))]
  return(mxs)}



# Quantify measures of demographic buffering

DB.all<-mxs<-NULL

for (i in 1:dim(Mammals.all)[1]){
  print(i)
  tryCatch({
    mxs<-replicate(50,stoch.sens_mean(sample(getAs.from.IDlist(Mammals.all$IDs[i]))[1:3]))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    DB.all$E_Smean[i]<-mean(unlist(lapply(array_to_matrix(mxs),sum)))
    DB.all$E_Smean.SD[i]<-sd(1-unlist(lapply(array_to_matrix(mxs),sum)))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    DB.all$Lambda[i]<-lambda(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    DB.all$Generation.time[i]<-generation.time(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    DB.all$resilience.999[i]<- convt(mean(getAs.from.IDlist(Mammals.all$IDs[i])), accuracy=.01)[1]
    DB.all$resilience.90[i]<- convt(mean(getAs.from.IDlist(Mammals.all$IDs[i])), accuracy=.1)[1]
    DB.all$resilience.damping[i]<-damping.ratio(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  mxs<-NULL
}


# Coerce data into right format

DB.all<-do.call(cbind.data.frame,DB.all)
DB.all<-cbind(as.data.frame(Mammals.all[,-4]),DB.all)
DB.all$E_Smean.SE<-DB.all$E_Smean.SD/sqrt(DB.all$n)
DB.all$byVR<-DB.all$SpeciesAuthor %in% Mammals.sub$SpeciesAuthor



# Some summary statistics of interest

#Maximum value of mean value
DB.all[DB.all$E_Smean==max(DB.all$E_Smean,na.rm = TRUE),]

#Minimum value
DB.all[DB.all$E_Smean==min(DB.all$E_Smean,na.rm = TRUE),]

#Number of orders
DB.all%>%
  count(Order)

#populations
DB.all%>%as.tibble()


# Paper metadata

metadata<-comadre$metadata[,c(1,2,3,6,13,14,15,16,18,25,34)]%>%
  filter(MatrixComposite == c("Individual","Pooled")) %>%
  filter(SpeciesAuthor %in% DB.all$SpeciesAuthor)%>%
  mutate(NewSpeciesAccepted = str_replace_all(SpeciesAccepted, "_", " "))%>%
  mutate(NewSpeciesAccepted = str_replace_all(NewSpeciesAccepted, " subsp.*", ""))%>%
  group_by(SpeciesAuthor)%>%
  mutate(Matrices.N= n())%>% 
  distinct(SpeciesAuthor,.keep_all=TRUE)%>%
  as.data.frame()%>%print()
names(metadata)

# Function to convert strings of authors to et al. if more than 2 authors are present

etaler<-function(X){
  A<-str_split(X,";")[[1]]
  B<-ifelse(length(A) == 1, print(A[1]), 
            ifelse(length(A) == 2,
                   print(paste0(A[1]," &",A[2])),
                   print(paste0(A[1]," et al."))))}

metadata<-metadata%>%
  mutate(Authors = etaler(Authors))


# Merge datasets

metadata<-left_join(metadata,DB.all,by="SpeciesAuthor")
colnames(metadata)
metadata<-metadata[,c(1,12,3,14,16,17,16,17,13,18,19,22,24)]
metadata[,6]<-metadata[,7]-1
colnames(metadata)<-c("SpeciesAuthorComadre","SpeciesName","CommonName","Order","E_Smu","E_Smu.SD","E_Ssig","E_Ssig.SD","# matrices","Lambda","Generation.time","Damping.ratio","ByAge")


# Check dataset

head(metadata)


# Check that mean matrices below are singular - if not, we cannot calculate their second derivatives.
#	FOR THIS REASON THE TOTAL AMOUNT OF POPULATIONS ARE 40 instead of 44.

DB.all$SpeciesAuthor[c(3,10,19,33)]

comadre$metadata[,c(1,2,3)]%>%
  mutate(NewSpeciesAccepted = str_replace_all(SpeciesAccepted, "_", " "))%>%
  mutate(NewSpeciesAccepted = str_replace_all(NewSpeciesAccepted, " subsp.*", ""))%>%
  filter(SpeciesAuthor %in% DB.all$SpeciesAuthor[-c(3,10,19,33)])%>%
  distinct(NewSpeciesAccepted)%>%print()

metadata%>%
  filter(SpeciesAuthorComadre %in% DB.all$SpeciesAuthor[-c(3,10,19,33)])
head(metadata)


# Filter data

DB.sub<-DB.all%>%
  drop_na(.)%>%
  filter(byVR=="TRUE")%>%
  filter(Generation.time != -Inf)%>%print()

DB.all%>%filter(byVR == "TRUE")


# Graphics for figure 3 from phylopics

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


# Generate empty lists

outputs2<-outputs1<-list()


# Updated Max_var_scale function

Max_var_scale2<-function(X){
  varmxcorrected<-(splitA(var2(X))$T*
                     (splitA(mean(X))$T*splitA(1-mean(X))$T))+
    splitA(var2(X))$F
  CV<-sqrt(varmxcorrected)/mean(X)*100
  return(list(A=varmxcorrected, CV=CV))
}


# Calculate parameter matrices

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


# Calculate derivatives

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


# Build individual figures for all mammals
 
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

