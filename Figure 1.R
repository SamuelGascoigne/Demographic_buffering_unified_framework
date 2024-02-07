# Import necessary packages and helper functions.

source("setup.R")
source("functions.R")


# Generate a vector of randomly drawn stochastic elasticity values and coerce into a dataframe.

Esig <- rnorm(100, mean = -0.03, sd = 0.01)
Esig_df<-data.frame(Esig)


# Set seed for figure 1.

set.seed(1234)


# Build figure 1a.

figure_1a <- ggplot(Esig_df, aes(x = Esig)) +
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
figure_1a



# Now it is time to build figure 1b and 1c.

# Subset the monkeyflower data from the popbio package.
# These will be used later in the script.

mim<- subset(monkeyflower, species == "cardinalis" &
               site == "Carlon" & year != "pooled", select = c(4:19))

## convert data frame to list of matrices using split

mim1 <- split(mim, 2000:2002)
mim2 <- lapply(mim1, matrix, nrow=4, byrow=TRUE)


# Calculate elasticity values.

elas <- elasticity(mean(mim2))


# Calculate an scale second derivative values.

z <- diag(secder_calculator(mean(mim2)))
zlog <- log10(z)
zlog(is.nan[zlog]) <- 0
zlog[is.nan(zlog)] <- 0
zlog[is.nan(zlog)] <- 0
zlog[zlog == -Inf] <- 0


# Story second derivative values in a dataframe.

zlog_df<-data.frame(as.vector(zlog))
names(zlog_df)<-c("SD")


# Store elasticity values

elas <- as.matrix(elas)
elas_df <- melt(elas)


# Build figures 1b and 1c.

figure_1b <- ggplot(elas_df, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill = value)) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3) +
  labs(x = "Column", y = "Row", title = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.3),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_x_discrete(limits = c("1", "2", "3", "4")) +
  scale_y_reverse()

figure_1b


figure_1c <- ggplot(elas_df, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_viridis_c(name = "Elasticity", begin = 0.3) +
  labs(x="Column", y="Row", title="") +
  theme_bw() + theme(axis.text.x=element_text(size=15, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=15),
                     plot.title=element_text(size=20)) +
  scale_x_discrete(limits=c("1","2","3","4")) +
  scale_y_reverse() +
  geom_point(aes(size = zlog_df$SD)) +
  scale_size_continuous(name = "Second-order derivative", range = c(15, 1))

figure_1c


