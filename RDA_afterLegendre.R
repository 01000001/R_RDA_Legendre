#######################################
#RDA from scratch after Legendre (1998)
#
#RDA steps on the covariance matrix:
#
#    1, Compute multivariate linear regression of the centered data matrix on the
#    standradized matrix of explinatory variables. This gives the matrix of fitted values.
#
#    2, Compute PCA on the matrix of fitted values
#
#    3, Compute the two types of site scores
#
#    4, Output the results



RDA <- function(Y,X) {
                             
###########################  
#1. Preperation of the data
###########################  
  
  Y.mat <- as.matrix(Y)
  Yc <- scale(Y.mat, scale = FALSE)
  
  X.mat <- as.matrix(X)
  Xcr <- scale(X.mat)

##########################################  
#2. Compute multivariate linear regression  
##########################################

#Compute B, the matrix of regression coefficients of all response variables Y on the
#regressors X - after eq. 11.4; B = [X' X]^-1 [X' Y]

#?solve
#solve(A) = Inverse of A where A is a square matrix.
  B <- solve(t(Xcr) %*% Xcr) %*% (t(Xcr) %*% Yc)

#Matrix of fitted values - eq. 11.5; Yhat = X * b
  Yhat <- Xcr %*% B

# Dimensions
  n <- nrow(Y)
  p <- ncol(Y)
  m <- ncol(X)

########################
#3. PCA on fitted values
########################

#Covariance matrix - eq. 11.7 S = [1/(n-1)] Yhat' Yhat

  S <- (1/(n-1)) * t(Yhat) %*% Yhat
  #S <- cov(Yhat)

# Eigenvalue decomposition of covariance matrix of Yhat
  eigenS <- eigen(S)

#Number of canonical axes
  ka <- length(which(eigenS$values > 0.0000000001))

#Eigenvalues of canonical axes
  ev <- eigenS$values[1:ka]

#Total variance of the centered matrix Yc
  trace <- sum(diag(cov(Yc)))

#Orthonormal eigenvectors (contributions of response variables)
#"Matrix u, of size (p * p), contains only min[p,m, n-1] eigenvectors with nonzero
#eigenvalues
#The canonical coefficients in the normalized matrix U give the contributions of the
#variables of to the canonical axes."

  U <- eigenS$vectors[,1:ka]
  row.names(U) <- colnames(Y)

#"The ordination of objects in the space of the response variables Y can be
#obtained directly from the centred matrix Y, using the standard equation for principal
#components (eq. 9.4) and matrix U of the eigenvectors"

#"The ordination vectors (columns of F) defined in eq. 11.12 are called the vectors of
#"site scores". They have variances that are close, but not equal to the corresponding
#eigenvalues."

#"Site" scores (vegan's 'wa' scores, scaling 1 - eq. 11.12) / Matrix of principal components
  Fmatrix <- Yc %*% U
  row.names(F) <- row.names(Y)

#"Site" constrains /"fitted site scores" (vegan's 'lc' scores, scaling 1 - eq. 11.13)
  Z <- Yhat %*% U
  row.names(Z) <- row.names(Y)

#Canonical coefficients - eq. 11.14
  CC <- B %*% U
  row.names(CC) <- colnames(X)

#Unadjusted R2
  R2 <- sum(ev/trace)

########################
#4. Output
########################


  result <- list(B,       #"regression_coefficients"
                 Yhat,    #"fitted_values"
                 U,       #"Matric of eigenvectors"
                 CC,      #"Canonical Coefficients"
                 trace,   #"Total variance"
                 R2,      #"R2"
                 ev,      #"Eigenvalues of canonical axes"
                 Z,       #"Site constraints"
                 Fmatrix, #"Site scores"
                 S       #"Covariance_matrix_of_fitted_values"
                 )

  names(result) <- c("regression_coefficients", 
                     "fitted_values", 
                     "Matric_of_eigenvectors", 
                     "Canonical_Coefficients",
                     "Total_variance",
                     "R2",
                     "Eigenvalues_of_canonical_axes",
                     "Site_constraints",
                     "Site_scores",
                     "Covariance_matrix_of_fitted_values")
  result
  
}


########################
#Playaround data
########################
C = matrix( c(1,2,2,1,3,-1,1,1,1,1,1,1), nrow = 3, ncol = 4, byrow = TRUE)

D = matrix( c(1,2,2,1,3,-1), nrow = 3, ncol = 2, byrow = TRUE)

##################################################################
#Data from Table 10.5 for chapter 11.1.2 Numerical examples of RDA
##################################################################

book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]


RDA_output <- RDA(Y,X)
#summary(RDA_output)
RDA_output$Eigenvalues_of_canonical_axes

#################
RDA_output
RDA_output$Eigenvalues_of_canonical_axes

library(vegan)
vegan_rda <- rda(Y,X)
vegan_rda
