#====================================================================
#			FIGURE 1 Ecol Lett 22nd May 2023 
#====================================================================


source("setup.R")


data(monkeyflower)
mim<- subset(monkeyflower, species == "cardinalis" &
               site == "Carlon" & year != "pooled", select = c(4:19))

## convert data frame to list of matrices using split
mim1<-split(mim, 2000:2002)
mim2<-lapply(mim1, matrix, nrow=4, byrow=TRUE)



Esig=rnorm(100, mean = -0.03, sd = 0.01)

#Respectivos valores da Elasticidade estoc?stica da m?dia relativa a cada Eso gerado
Emu= 1- Esig

#Raz?o entre os respectivos valores de elasticidade estoc?stica da vari?ncia e da m?dia
BCV<- Esig/Emu

# Organiza??o dos dados em uma ?nica tabela
data.teste<-data.frame(Esig,Emu,BCV)

#Primeiro gr?fico 
plot(Emu,Esig, xlab="Sum das Elasticidades da M?dia Su", ylab="Sum das Elasticidades da vari?ncia - So")
abline(v=.8, col="red")
text(.4,.5 ,
     "Entirematrixbuffered")

E.min<-min(Emu)
E.max<-max(Emu)

#Gr?fico figura 6
g1<-ggplot(data.teste, aes(x=Emu, y=Esig, color=Emu)) + 
  geom_point(size=5,aes(color=Emu),alpha=.6,size=4) +
  scale_colour_gradient(guide=guide_colourbar(reverse = TRUE),
                        low = "blue", high = "red",
                        breaks=c(E.min,E.max),labels=c("   More  \nbuffered"," More \n labile"))+
  labs(colour = "")+
  scale_alpha(range=c(0.3))+
  ylab(expression("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +")) +
  xlab(expression("(More Buffered)            -       " %<-% "        "~Sigma~"E"^"s"^~mu~"        " %->%  "       +            (More Labile)")) +
  #geom_hline(aes(yintercept = 0),linetype="dashed",size=1.3) +
  theme_bw()+
  theme(axis.text=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 18, colour = "black")
  )
g1

set.seed(1234)

g1 <- ggplot(data.teste, aes(x = Esig)) +
  geom_density(fill = muted("purple"), alpha = 0.05) +
  geom_point(aes(fill = Esig, y = 21), position = position_jitter(height = 10),
             shape = 21,
             size = 2.5,
             stroke = 1.5,
             alpha = 0.6) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_fill_viridis_c(begin = 0.3) +
  xlab(expression("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +")) +
  ylab("Density") +
  #lab(fill = expression(~Sigma~"E"^"s"^~sigma~)) +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.title= element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18, colour = "black")
  )
g1

# the correct one:
g1 <- ggplot(data.teste, aes(x = Esig)) +
  geom_density(fill = "purple", alpha = 0.05) +
  geom_point(aes(fill = Esig, y = 21), position = position_jitter(height = 10),
             shape = 21,
             size = 2.5,
             stroke = 1.5,
             alpha = 0.6) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_fill_viridis_c(begin = 0.3) +
  xlab(expression(paste("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +"))) +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 18, color = "black")
  )
g1

# g1 <- ggplot(data.teste, aes(x = Esig)) +
#   geom_density(fill = "purple", alpha = 0.05) +
#   geom_point(aes(fill = Esig, y = 21), position = position_jitter(height = 10),
#              shape = 21,
#              size = 2.5,
#              stroke = 1.5,
#              alpha = 0.6) +
#   scale_y_continuous(expand = c(0, 1), breaks = "all") +
#   scale_x_continuous(expand = c(0.001, 0.001), breaks = "all") +
#   scale_fill_viridis_c(begin = 0.3) +
#   xlab(expression(paste("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +"))) +
#   ylab("Density") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 18, face = "bold"),
#     legend.title = element_text(size = 18, color = "black"),
#     legend.text = element_text(size = 18, color = "black"),
#     axis.ticks = element_line(color = "black")
#   )
# 
# #g1 <- ggplot(data.teste, aes(x = Esig)) +  # Create a ggplot object with "Esig" variable on the x-axis
#   geom_density(fill = "blue", alpha = 0.05) +  # Add a density plot with blue fill and 0.05 alpha
#   geom_jitter(aes(color = Esig), height = 0.1, width = 0.1) +  # Add jittered points with color mapped to "Esig"
#   scale_color_viridis_c(begin = 0.3) +  # Use "viridis" color palette for the points
#   xlab(expression("-    " %<-% "  "~Sigma~"E"^"s"^~sigma~"  " %->%  "  +")) +  # Add x-axis label with mathematical expression
#   ylab("Density") +  # Add y-axis label
#   theme_bw() +  # Use a black and white theme
#   theme(axis.text=element_blank(),  # Remove axis text
#         axis.title=element_text(size=18,face="bold"),  # Customize axis titles
#         legend.title= element_text(size = 18, colour = "black"),  # Customize legend title
#         legend.text = element_text(size = 18, colour = "black")  # Customize legend text
#   )
# 
#   
# #g1

#===============================================================
#	Second derivative Calculator
#===============================================================
teste<-function (A){
  q<- A != 0
  size<- dim(A)
  qq<- matrix(q, ncol = 1)
  D <- NULL
  for (j in 1:size[2]) {
    for (i in 1:size[1]) {
      d2 <- secder(A, i, j)
      D <- cbind(D, matrix(d2, ncol = 1) * qq)}}
  D}
#===============================================================
#-----------------------------------------------------------------------------------
#					REPLACING
#		Add labels for each matrix's element
#-----------------------------------------------------------------------------------
replacing<-function(A){
  A<-A
  A[upper.tri(A)]<-"T"
  A[1,]<-"R"
  A[lower.tri(A,diag=T)]<-"T"
  return(A)
}
#===============================================================teste<-function (A){

z<-diag(teste(mean(mim2)))

zlog<-log10(z)

zlog(is.nan[zlog])<-0

zlog[is.nan(zlog)]<-0
zlog[is.nan(zlog)]<-0
zlog[zlog == -Inf]<-0

#image(matrix(zlog,4))

zmx<-matrix(zlog,4)

mx<-melt(zmx)
names(mx) <- c("x", "y", "z")
elas<-elasticity(mean(mim2))
mx.CV<-sqrt(var2(mim2))/mean(mim2)

data.teste<-data.frame(
  as.vector(zlog),
  as.vector(elas),
  as.vector(mx.CV),
  rotulo<-as.vector(replacing(mean(mim2))))

names(data.teste)<-c("SD","elas","CV","rotulo")
zmin<-min(zlog)
zmax<-max(zlog)

#data.teste2<-data.teste
#data.teste2$SD<-data.teste2$SD+1
#data.teste2$SD[data.teste2$SD<= -1]<- -1
#data.teste2$SD[data.teste2$SD+1 >= 1]<- 1

g3<-ggplot(data.teste, aes(x = CV, y = elas,color=SD)) +
  geom_point(size=7,alpha=.3) +
  geom_text(size=5,aes(label=rotulo),color="black")+		#hjust=1.3,vjust=1.3
  #scale_color_gradient2(midpoint=0, limits=c(min(data.teste2$SD),max(data.teste2$SD)),low="skyblue", mid="grey",
  #high="tomato", space ="Lab",
  scale_color_viridis_c(begin = 0.3,
                        breaks=c(min(data.teste$SD),0,max(data.teste$SD)),labels=c("-   \n More \n buffered","0","Less \n buffered \n+  "))+
  labs(colour = " Second \nderivative \n \n")+
  scale_x_log10()+scale_y_log10()+
  xlab(expression("-   " %<-% "  Importance to"  ~ lambda ~" " %->%"  +"))+
  ylab(expression("-   " %<-% "  CV  "  %->%"  +"))+
  geom_smooth(method = "lm", se = FALSE,color="black")+
  theme_bw()+
  theme(axis.text=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.title= element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18, colour = "black")
  )
g3

elas <- as.matrix(elas)
datag2 <- melt(elas)



g2 <- ggplot(datag2, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill = value)) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3) +
  labs(x = "Column", y = "Row", title = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.3),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_x_discrete(limits = c("1", "2", "3", "4")) +
  scale_y_reverse()

g2


g3 <- ggplot(datag2, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3) +  # Change the name of the fill legend
  labs(x="Column", y="Row", title="") +
  theme_bw() + theme(axis.text.x=element_text(size=15, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=15),
                     plot.title=element_text(size=20)) +
  scale_x_discrete(limits=c("1","2","3","4")) +
  scale_y_reverse() +
  geom_point(aes(size = data.teste$SD)) +
  scale_size_continuous(name = "Second-order derivative", range = c(15, 1))  # Change the name of the size legend

g3


