"
This function calculates the K matrix (as a vector) of the volume scattering model proposed by Antropov et al. 2011. 
The parameters for a volume scattering model are γ and ρ. According to Antropov et al 2011, γ = <|SHH|^2>/<|SVV|^2> 
(copol ratio for each pixel), while ρ = 1/3 (implicit in the matrix equation).

Input data:
- γ (from gammaFunction)

Author: Elisa Frank Buss
Encoding: UTF-8
Creation date: Feb 2022
Last modification: May 2022

For the bibliography cited in this file see .../Theses-MAIE-Elisa-Frank-Buss/Texts/PDF/Tesis_Maestria_Elisa_Frank_Buss.pdf

"

vectorVFunction <- function(g) {
    
    # Obtain array dimensions
    dim.input <- dim(g)
    
    # Convert array to matrix
    g <- matrix(g, dim.input[1], dim.input[2])
    
    # Obtain matrix dimensions
    dim.input2 <- dim(g)
    
    # Create a matrix with zeros (4x4) for the matrix Kv 
    Kv <- array(rep(0, prod(dim.input2)*16), 
                dim = c(dim.input2, 16))
    
    
    # Obtain Kv
    indexx <- 0 
    for(i in 1:(dim.input2[1])) {  
        indexx <- indexx + 1                     
        indexy <- 0 
        for(j in 1:(dim.input2[2])) {
            indexy <- indexy + 1 
            
            # Generate Kv (equation obtained from page 2 of Ratha et al 2019)
            vV <- c(3/2*(1+g[i,j])-(sqrt(g[i,j])/3),
                    g[i,j]-1,
                    0,
                    0,
                    g[i,j]-1,
                    1/2*(1+g[i,j])+(sqrt(g[i,j])/3),
                    0,
                    0,
                    0,
                    0,
                    1/2*(1+g[i,j])+(sqrt(g[i,j])/3),
                    0,
                    0,
                    0,
                    0,
                    1/2*(1+g[i,j])-sqrt(g[i,j]))
            mul <- 1/((3*(1+g[i,j]))/4 - sqrt(g[i,j])/6)
            
            # Store Kv matrix
            Kv[indexx, indexy,] <- mul * vV
        }
    }
    
    # Save the result
    return(Kv)
}
