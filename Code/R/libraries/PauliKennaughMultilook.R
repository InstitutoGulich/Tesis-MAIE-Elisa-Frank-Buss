"
This function obtains the coherence matrix (T) for the Pauli decomposition 
and the Kennaugh matrix (K). At the same time, it applies a multilooking of x*y. 
The result is an average of 4 pixels for each element of the T matrix and for each 
element of the K matrix.
NoData are not taken into account in the analysis. When a NoData is present in 
at least one of the pixels, the result is NoData.

Input data:
_the polarizations HH, HV, VV readed with the myread function
_the values for x and y, in this case 2.

Author: Alejandro Frery and Elisa Frank Buss
Encoding: UTF-8
Creation date: May 2021
Last modification: Sep 2022

"

PauliKennaughMultilook <- function(SHH, SHV, SVV, x, y) {
  # Constant
  sqrt2 <- sqrt(2)
  xy <- x*y
  
  # Obtain matrix dimensions
  dim.input <- dim(SHH)
  
  # Determine the dimensions of the output matrix
  dim.output <- c(floor(dim.input[1]/x), floor(dim.input[2]/y))
  
  # Create an array to store the T matrix
  Pauli <- array(rep(complex(real=0, imaginary=0), prod(dim.output)*9),
                  dim=c(dim.output, 9))
  
  # Create an array to store the K matrix
  Kennaugh <- array(rep(0, prod(dim.output)*16), 
                    dim = c(dim.output, 16))
  
  # The following rules indicate not to analyse the last row if the number of rows is odd.
  # The same for the columns.
  # WARNING, this part is designed with the assumption that x,y = 2.
  indexx <- 0
    for(br in seq(from=1, to=(dim.input[1]-1), by=x)) {
    indexx <- indexx + 1
    
    indexy <- 0
    for(bc in seq(from=1, to=(dim.input[2]-1), by=y)) {
      indexy <- indexy + 1 
      
      TP <- matrix(rep(complex(real=0, imaginary=0), 9),
                   nrow=3, ncol=3)
      
      # start indexes
      I <- 0 
      J <- 0
      BR <- 0
      BC <- 0
      
      # Multilooking
      for(i in br:(br+(x-1))) { 
          for(j in bc:(bc+(y-1))) {
              # update indexes
              
              I <- I + i
              J <- J + j
              BR <- BR + br
              BC <- BC + bc
              
              # kp vector generation
              P <- c(SHH[i,j]+SVV[i,j], SHH[i,j]-SVV[i,j], 2*SHV[i,j])
              
              # The T matrix is obtained by multiplying kp by its conjugate transpose, resulting in a 3x3 matrix.
              # In turn, the 4 T matrices of the 4 pixels involved in the multilooking are added together.
              TP <- TP + (P %*% t(Conj(P)))
              }
          }
      
          # T matrix (9 elements) is obtained. It divides by 4 (xy) completing the multilooking result.
          # It is necessary to transpose because R organises the result in columns.
          # It is possible to get the Pauli decomposition from this result (because of that, the name Pauli).
          Pauli[indexx, indexy,] <- sqrt2 * as.vector(t(TP)) / xy
          
          # K matrix (16 elements) is obtained. It divides by 4 (xy) completing the multilooking result.
          # Also, it is necessary to transpose because R organises the result in columns.
          Kennaugh[indexx, indexy,] <- c(Re((t(TP)[1])+(t(TP)[5])+(t(TP)[9]))/2,     #1 element
                                             Re((t(TP)[2])),                         #2       
                                             Re((t(TP)[3])),                         #3           
                                             Im((t(TP)[6])),                         #4           
                                             Re((t(TP)[2])),                         #5   
                                             Re((t(TP)[1])+(t(TP)[5])-(t(TP)[9]))/2, #6 
                                             Re((t(TP)[6])),                         #7
                                             Im((t(TP)[3])),                         #8 
                                             Re((t(TP)[3])),                         #9
                                             Re((t(TP)[6])),                         #10
                                             Re((t(TP)[1])-(t(TP)[5])+(t(TP)[9]))/2, #11
                                             (-Im((t(TP)[2]))),                      #12
                                             Im((t(TP)[6])),                         #13
                                             Im((t(TP)[3])),                         #14
                                             (-Im((t(TP)[2]))),                      #15
                                             Re(-(t(TP)[1])+(t(TP)[5])+(t(TP)[9]))/2 #16 element
                                             ) /xy
    }
    }
     
  # Save the result
  return(list(Pauli=Pauli, Kennaugh=Kennaugh))
}
