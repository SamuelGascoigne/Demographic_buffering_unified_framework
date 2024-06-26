# Return the matrix complete cases

there.is.na<-function(your_list){
  there<-NULL
  for( i in 1:length(your_list)){
    there[[i]]<-all(complete.cases(your_list[[i]]))
  }
  return(as.list(there))
}


#This function returns a list with all MatA matrices with non-zero values

getMatA<-function(sp){
	#sp = "SpeciesAuthor" values
MEANS<-Matrices2<-NULL
sp<-as.character(sp)
 id<-grep(paste(paste("^",sp,"$", sep="")),comadre$metadata$SpeciesAuthor)
	#Take all matrices in a hit and than reduce to matA only
 Matrices<-comadre$mat[(id)]
for(i in 1:length(Matrices)){
Matrices2[[i]]<-Matrices[[i]][[1]]
}
return(Matrices2[unlist(there.is.na(Matrices2))])
}


# Function to extract key properties of matrices

check<-function(sp){
  
  resultado<-NULL
  resultado$population<-as.character(sp)
  resultado$Numero.matrices=length(getMatA(sp))
  resultado$IDs.Comadre<-grep(paste(paste("^",sp,"$", sep="")),comadre$metadata$SpeciesAuthor)
  resultado$matrix.media=mean(getMatA(sp))
  
  return(resultado)
}


# Calculate stochastic elasticities with respect to means

stoch.sens_mean<-function (A, tlimit = 100) 
{
    if (!is.list(A) && !is.matrix(A[[1]])) 
        stop("A should be a list of matrices")
    k <- ncol(A[[1]])
    A <- sample(A, tlimit, replace = TRUE)
    tlimit <- length(A)
    wvec <- rep(1/k, k)
    w <- cbind(wvec)
    r <- rep(0, tlimit)
    for (i in 1:tlimit) {
        a <- A[[i]]
        wvec <- a %*% wvec
        r[i] <- sum(wvec)
        wvec <- wvec/r[i]
        w <- cbind(w, wvec)
    }
    vvec <- rep(1/k, k)
    v <- cbind(vvec)
    for (i in rev(1:tlimit)) {
        a <- A[[i]]
        vvec <- vvec %*% a
        v <- cbind(t(vvec), v)
    }
    sensmat <- matrix(0, nrow = k, ncol = k)
    elasmat <- matrix(0, nrow = k, ncol = k)
    elasmean <- matrix(0, nrow = k, ncol = k)
elassig <- matrix(0, nrow = k, ncol = k)
for (i in 1:tlimit) {
      sensmat <- sensmat + ((v[, i + 1] %*% t(w[, i]))/as.numeric(r[i] * 
            t(v[, i + 1]) %*% w[, i + 1]))
        a <- A[[i]]
     elasmat <- elasmat + ((v[, i + 1] %*% t(w[, i]) * a)/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
elasmean <- elasmean + ((v[, i + 1] %*% t(w[, i]) * mean(A))/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
elassig <- elassig + ((v[, i + 1] %*% t(w[, i]) * (a-mean(A)))/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
    }
    sensmat <- sensmat/tlimit
    elasmat <- elasmat/tlimit
    elasmean <- elasmean/tlimit
    elassig <- elassig/tlimit
    out <- elasmean
    out
}

# The following function to calculate the stochastic population growth rate by simulation
# was taken from the popbio package.

stochastic_growth_rate_sim <- function(matrices, prob = NULL, maxt = 50000, verbose = FALSE){
  if (is.list(matrices)) {
    matrices <- matrix(unlist(matrices), ncol = length(matrices))
  }
  s <- sqrt(dim(matrices)[1]) ## number of stage classes
  n <- dim(matrices)[2] ## number of matrixes
  # default equal probabilities
  if (is.null(prob)) {
    prob <- rep(1 / n, n)
  }

  ##  Simulation
  r <- numeric(maxt)
  n0 <- rep(1, times = s)
  for (t in 1:maxt)
  {
    if (verbose) {
      if (t == 1 || t %% 10000 == 0) {
        message("Calculating stochastic growth at time ", t)
      }
    }
    col <- sample(1:n, 1, prob = prob)
    A <- matrix(matrices[, col], nrow = s)
    n0 <- A %*% n0
    N <- sum(n0)
    r[t] <- log(N)
    n0 <- n0 / N
  }
  loglsim <- mean(r)
  dse <- 1.96 * sqrt(var(r) / maxt)
  CI <- c(loglsim - dse, loglsim + dse)

  ## output...
  stoch <- exp(loglsim)
  
  return(stoch)
  
}


# Calculate stochastic elasticities with respect to variances

stoch.sens_sig<-function (A, tlimit = 100) 
{
    if (!is.list(A) && !is.matrix(A[[1]])) 
        stop("A should be a list of matrices")
    k <- ncol(A[[1]])
    A <- sample(A, tlimit, replace = TRUE)
    tlimit <- length(A)
    wvec <- rep(1/k, k)
    w <- cbind(wvec)
    r <- rep(0, tlimit)
    for (i in 1:tlimit) {
        a <- A[[i]]
        wvec <- a %*% wvec
        r[i] <- sum(wvec)
        wvec <- wvec/r[i]
        w <- cbind(w, wvec)
    }
    vvec <- rep(1/k, k)
    v <- cbind(vvec)
    for (i in rev(1:tlimit)) {
        a <- A[[i]]
        vvec <- vvec %*% a
        v <- cbind(t(vvec), v)
    }
    sensmat <- matrix(0, nrow = k, ncol = k)
    elasmat <- matrix(0, nrow = k, ncol = k)
    elasmean <- matrix(0, nrow = k, ncol = k)
elassig <- matrix(0, nrow = k, ncol = k)
for (i in 1:tlimit) {
        sensmat <- sensmat + ((v[, i + 1] %*% t(w[, i]))/as.numeric(r[i] * 
            t(v[, i + 1]) %*% w[, i + 1]))
        a <- A[[i]]
        elasmat <- elasmat + ((v[, i + 1] %*% t(w[, i]) * a)/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
elasmean <- elasmean + ((v[, i + 1] %*% t(w[, i]) * mean(A))/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
elassig <- elassig + ((v[, i + 1] %*% t(w[, i]) * (a-mean(A)))/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
    }
    sensmat <- sensmat/tlimit
    elasmat <- elasmat/tlimit
    elasmean <- elasmean/tlimit
    elassig <- elassig/tlimit
    out <- elassig
    out
}


# Transform an array into a list of matrices

array_to_matrix<-function(A){
  
  lapply(seq(dim(A)[3]), function(x) A[ , , x])
  
  }


# A function to log transform negative values

 reverselog_trans<-function(base = exp(1)) {
   
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
    
    }

# Function to scale CV by max CV
 
Max_var_scale<-function(X){
  
  varmxcorrected<-var2(X)*(mean(X)*(1-mean(X)))
  CV<-sqrt(varmxcorrected)/mean(X)*100
  return(list(A=varmxcorrected, CV=CV))
  
}




# Define biological meaning of matrix elements

Bio_meaning<-function(A){
A<-A
A[upper.tri(A)]<-"T"
A[1,]<-"R"
A[lower.tri(A,diag=T)]<-"T"
return(A)
}


# Automate the extraction of self second derivative values

mysec <- function(A){
B<-numeric(dim(A)[[1]])%o%numeric(dim(A)[[1]])
for(i in 1:dim(A)[1]){
	for(j in 1:dim(A)[1]){
B[i,j]<-secder(A,i,j)[i,j]
	}
   }
return(B)
}


# Calculate all second derivatives associated with an MPM

secder_calculator<-function (A){
  q<- A != 0
  size<- dim(A)
  qq<- matrix(q, ncol = 1)
  D <- NULL
  for (j in 1:size[2]) {
    for (i in 1:size[1]) {
      d2 <- secder(A, i, j)
      D <- cbind(D, matrix(d2, ncol = 1) * qq)}}
  D}


# Replace names of demographic rates

replacing<-function(A){
  A<-A
  A[upper.tri(A)]<-"T"
  A[1,]<-"R"
  A[lower.tri(A,diag=T)]<-"T"
  return(A)
}


