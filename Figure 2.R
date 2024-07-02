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
    
    mxs<-replicate(50, stoch.sens_sig(sample(getAs.from.IDlist(Mammals.all$IDs[i]))[1:3]))
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    
    DB.all$E_Smean[i]<-mean(unlist(lapply(array_to_matrix(mxs),sum)))
    DB.all$E_S_SD[i]<-sd(1-unlist(lapply(array_to_matrix(mxs),sum)))
    DB.all$E_S_SE[i]<-DB.all$E_S_SD[i]/sqrt(50)
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    
    DB.all$Lambda[i]<-lambda(mean(getAs.from.IDlist(Mammals.all$IDs[i])))
    DB.all$Stoch_lambda[i]<-stochastic_growth_rate_sim(getAs.from.IDlist(Mammals.all$IDs[i]), 
                                                       verbose = FALSE)
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  mxs<-NULL
  
}


# Coerce data into right format

DB.all<-do.call(cbind.data.frame,DB.all)
DB.all<-cbind(as.data.frame(Mammals.all[,-4]),DB.all)


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


# Filter and reorder metadata dataframe columns, then change column names

metadata<-metadata[,c(2, 3, 1, 4, 15, 19, 20, 16, 18)]
colnames(metadata)<-c("Species","Common_name","Species_COMADRE","Order","# matrices","Lambda", "Stochastic_lambda", "Stoch_elas_var", "Stoch_elas_var_SE")


# Check dataset

head(metadata)



# Plot figure 2

metadata <- metadata[!is.nan(metadata$Stoch_elas_var), ]



# Define the color palette

order_colors <- viridis(length(unique(metadata$Order)))

set.seed(1)

figure_2 <- ggplot(metadata, aes(x = Stoch_elas_var, fill = as.factor(Order))) +
  geom_density(fill = "purple", alpha = 0.05) +
  geom_errorbar(aes(xmin = Stoch_elas_var - Stoch_elas_var_SE, 
                    xmax = Stoch_elas_var + Stoch_elas_var_SE,
                    y = 15),
                position = position_jitter(height = 14, seed = 31),
                width = 0,
                alpha = 0.6) +
  geom_point(aes(y = 15, size = as.numeric(`# matrices`)), position = position_jitter(height = 14, seed = 31),
             shape = 21,
             stroke = 1.5,
             alpha = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 32)) +
  scale_x_continuous(expand = c(0.002, 0.002)) +
  scale_fill_manual(values = order_colors) +
  xlab(expression(paste("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +"))) +
  ylab("Density") +
  labs(fill = "Order", size = "# matrices") +
  theme_bw() +
  theme(
  axis.text = element_text(size = 24),
  axis.title = element_text(size = 30, face = "bold"),
  legend.title = element_text(size = 30, color = "black"),
  legend.text = element_text(size = 24, color = "black")
  )

figure_2 # 1,600x800

