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


# Plot figure 2

metadata <- metadata[!is.nan(metadata$E_Ssig), ]

metadata$E_SsigNEG <- NA
metadata$E_SsigNEG <- 1-metadata$E_Ssig


# Define the color palette

order_colors <- viridis(length(unique(metadata$Order)))

figure_2 <- ggplot(metadata, aes(x = E_SsigNEG, fill = as.factor(Order))) +
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

figure_2

