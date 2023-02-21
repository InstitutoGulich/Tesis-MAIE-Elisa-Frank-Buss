"
This function calculates the geodetic distance (GD) between two matrices (equation obtained from page 2 of Ratha et al. 2019,
also on page 3 of Mandal et al 2020)

Input data:
- two K matrixes

Author: Elisa Frank Buss
Encoding: UTF-8
Creation date: Mar 2022
Last modification: May 2022

For the bibliography cited in this file see .../Theses-MAIE-Elisa-Frank-Buss/Texts/PDF/Tesis_Maestria_Elisa_Frank_Buss.pdf

"

GDFunction <- function(K1, K2){
    
    # Trace operator calculation
    traceFunction <- function(M) {
        # Obtain matrix dimensions
        n <- dim(M)[1]
        
        tr <- 0
        
        # Go through the diagonal elements of the matrix and sum them to tr
        for (k in 1:n) {
            l <- M[k,k]
            tr <- tr + l
        }
        # Save the result
        return(tr[[1]])
    }
    
    # Numerator of the GD equation: is the trace of the matrix resulting from the product between the transpose of K1 and K2
    num <- traceFunction(t(K1)%*%K2)
    
    # One of the denominator terms of the GD equation, is the square root of the trace of the matrix resulting from the product
    # of the transpose of K1 and the K1 matrix
    den1 <- sqrt(traceFunction(t(K1)%*%K1))
    
    # The another term in the denominator of the GD equation, is the square root of the trace of the matrix resulting from
    # the product of the transpose of K2 and the K2 matrix
    den2 <- sqrt(traceFunction(t(K2)%*%K2))
    
    # Calculation of the geodetic distance, 2/pi normalises the range to [0, 1]
    GD <- (2/pi)*acos(num/(den1*den2))
    
    # Save the result
    return(GD)
}

