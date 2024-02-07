#################################################################
#		CHAPTER 1 - FUNCTIONS/AUTOMATIZATION
#	TWENTY YEARS OF DEMOGRAPHIC BUFFERING - A NEW FRAMEWORK
#		wrote by Gabriel Santos 04 Oct 2019
#################################################################

#=======================================================================#
#			there.is.na Function
#=======================================================================#
#Return the matrix complete cases

there.is.na<-function(your_list){
there<-NULL
for( i in 1:length(your_list)){
there[[i]]<-all(complete.cases(your_list[[i]]))
}
return(as.list(there))
}

#USAGE:
#Create some data - Putting a example here!


#=======================================================================#
#			getMatA
#=======================================================================#

#This function return a list with all MatA with non-zero values
#------------------------------------------------------------------------------
#Matrices2<-Matrices<-NULL

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


#USAGE
#getMatA(unique(Mammals.all$SpeciesAuthor)[159])

#====================================================================================
#			CHECK MATRICES ONE BY ONE
#=====================================================================================

check<-function(sp){
resultado<-NULL
resultado$population<-as.character(sp)
resultado$Numero.matrices=length(getMatA(sp))
resultado$IDs.Comadre<-grep(paste(paste("^",sp,"$", sep="")),comadre$metadata$SpeciesAuthor)
resultado$matrix.media=mean(getMatA(sp))
return(resultado)
}

#USAGE:
#check(Mammals.list[[1]])

#Pode ajudar tamb?m
#comadre$mat[(764:789)]




#====================================================================================
#	Stochastic elasticity from the MEAN value for each matrix value
#=====================================================================================
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

#====================================================================================
#	Stochastic elasticity from the VARIANCE value for each matrix value
#====================================================================================
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
#====================================================================================
#====================================================================================
#					Array to matrix function
#====================================================================================
	#Transform an array in a list of matrices

array_to_matrix<-function(A){
lapply(seq(dim(A)[3]), function(x) A[ , , x])}

#====================================================================================
#====================================================================================
#					REVERSELOG_TRANS
#====================================================================================

 reverselog_trans<-function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))}

#====================================================================================
#		CORRECT COEFFICIENT OF VARIATION BY MAXIMUM VARIANCE
#====================================================================================
Max_var_scale<-function(X){
varmxcorrected<-var2(X)*(mean(X)*(1-mean(X)))
CV<-sqrt(varmxcorrected)/mean(X)*100
return(list(A=varmxcorrected, CV=CV))
}

#CORRECTING ONLY TRANSITIONS ELEMENTS - DEPRECATED
#Max_var_scale<-function(X){
#varmxcorrected<-splitA(var2(X))$T*(mean(X)*(1-mean(X))) +
# splitA(var2(X))$F
#CV<-sqrt(varmxcorrected)/mean(X)*100
#return(list(A=varmxcorrected, CV=CV))
#}



#====================================================================================
#		DEFINE BIOLOGICAL MEANING OF MATRIX ELEMENTS
#====================================================================================

Bio_meaning<-function(A){
A<-A
A[upper.tri(A)]<-"T"
A[1,]<-"R"
A[lower.tri(A,diag=T)]<-"T"
return(A)
}



#====================================================================================
#		AUTOMATIZE THE EXTRACTION OF THE SELF-SECOND DERIVATIVE
#====================================================================================
mysec <- function(A){
B<-numeric(dim(A)[[1]])%o%numeric(dim(A)[[1]])
for(i in 1:dim(A)[1]){
	for(j in 1:dim(A)[1]){
B[i,j]<-secder(A,i,j)[i,j]
	}
   }
return(B)
}





