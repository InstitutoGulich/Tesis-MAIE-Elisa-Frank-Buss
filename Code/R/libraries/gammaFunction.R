"
This function calculates the copol ratio, γ = <|SHH|^2>/<|SVV|^2>, according to 
the model of Antropov et al (2011), for each pixel of the input matrix, 
applying multilooking, with x,y = 2.

Input data:
_the polarizations HH, HV, VV readed with the myread function
_the values for x and y, in this case 2.

Author: Elisa Frank Buss
Encoding: UTF-8
Creation date: Feb 2022
Last modification: May 2022

For the bibliography cited in this file see .../Theses-MAIE-Elisa-Frank-Buss/Texts/PDF/Tesis_Maestria_Elisa_Frank_Buss.pdf

"

gammaFunction <- function(SHH, SVV, x, y) {
    
    # Obtain matrix dimensions
    dim.input <- dim(SHH)
    
    # Determine the dimensions of the output matrix
    dim.output <- c(floor(dim.input[1]/x), floor(dim.input[2]/y))
    
    # Save gamma (γ) values in an array.
    gammaArray <- array(rep(0, prod(dim.output)),
                        dim = c(dim.output, 1))
    
    xy <- x*y
    
    # The following rules indicate not to analyse the last row if the number of rows is odd.
    # The same for the columns.
    # WARNING, this part is designed with the assumption that x,y = 2.
    indexx <- 0
    for(br in seq(from=1, to=(dim.input[1]-1), by=x)) {
        indexx <- indexx + 1
        indexy <- 0
        for(bc in seq(from=1, to=(dim.input[2]-1), by=y)) {
            indexy <- indexy + 1
            
            # start indexes
            I <- 0 
            J <- 0
            BR <- 0
            BC <- 0
            
            numerador <- 0
            denominador <- 0
            
            # Multilooking
            for(i in br:(br+(x-1))) { 
                for(j in bc:(bc+(y-1))) {
                    
                    # update indexes
                    I <- I + i
                    J <- J + j
                    BR <- BR + br
                    BC <- BC + bc
                    
                    numerador <- numerador + (Re(SHH[i,j]))^2
                    denominador <- denominador + (Re(SVV[i,j]))^2
                }
            }
            
            # Store gamma values
            gammaArray[indexx, indexy,] <- (numerador/xy)/(denominador/xy) # same as "(numerador / denominador) / xy
                                                                           
        }
    }
    # Save the result
    return(gammaArray)
}
