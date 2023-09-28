"
This function calculates the Generalised Radar Vegetation Index (GRVI) and applies a prior a multilooking 
(x*y) of the radar bands.

Input data:
- the polarizations HH, HV, VV readed with the myread function
- the values for x and y, in this case 2.

Author: Elisa Frank Buss
Encoding: UTF-8
Creation date: Mar 2022
Last modification: May 2022

For the bibliography cited in this file see .../Theses-MAIE-Elisa-Frank-Buss/Texts/PDF/Tesis_Maestria_Elisa_Frank_Buss.pdf

"

GRVIFunction <- function(SHH, SHV, SVV, x, y) {
    
    ### First step: generation of the observed Kennaugh matrix with multilooking (x*y)
    Resultado <- PauliKennaughMultilook(SHH, SHV, SVV, x, y)
    Kennaugh <- Resultado$Kennaugh
    
    # Obtain array dimensions
    dim.input <- dim(Kennaugh)
    dim.input2 <- c(dim.input[1], dim.input[2])
    
    ### Second step: formulation of the elementary target vectors to generate their K matrixes
    vt <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1)
    vc <- c(5/8, 3/8, 0, 0, 3/8, 5/8, 0, 0, 0, 0, 1/2, 0, 0, 0, 0, -1/2)
    vd <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1)
    vnd <- c(5/8, 3/8, 0, 0, 3/8, 5/8, 0, 0, 0, 0, -1/2, 0, 0, 0, 0, 1/2)
    
    # K matrixes of the elementary targets
    Kt <- matrix(vt, nrow = 4, ncol = 4, byrow = TRUE)
    Kc <- matrix(vc, nrow = 4, ncol = 4, byrow = TRUE) 
    Kd <- matrix(vd, nrow = 4, ncol = 4, byrow = TRUE)
    Knd <- matrix(vnd, nrow = 4, ncol = 4, byrow = TRUE)
    
    ### Third step: calculation of γ for the volume scattering model according to Antropov et al 2011,
    # and from γ the Kv matrix
    
    # Gamma (γ) is calculated for each pixel, also doing a multilooking x*y
    g <- gammaFunction(SHH, SVV, x, y)
    
    # Calculate the K matrix of the volume scattering model (per pixel)
    Kv <- vectorVFunction(g)
    
    ### Fourth step: calculation of the GRVI
    
    # Create an array to store the GRVI
    GRVI <- array(rep(0, prod(dim.input2)), dim = c(dim.input2, 1))
    
    # Calculate GRVI
    indexx <- 0 
    for(i in 1:dim.input2[1]) {  
        indexx <- indexx + 1
        indexy <- 0 
            for(j in 1:dim.input2[2]) {
                indexy <- indexy + 1 
                
                # The observed K matrix per pixel is generated
                Kennaugh0 <- matrix(0, 4, 4)
                Kennaugh0[1,1] <- Kennaugh[i, j, 1]
                Kennaugh0[1,2] <- Kennaugh[i, j, 2]
                Kennaugh0[1,3] <- Kennaugh[i, j, 3]
                Kennaugh0[1,4] <- Kennaugh[i, j, 4]
                Kennaugh0[2,1] <- Kennaugh[i, j, 5]
                Kennaugh0[2,2] <- Kennaugh[i, j, 6]
                Kennaugh0[2,3] <- Kennaugh[i, j, 7]
                Kennaugh0[2,4] <- Kennaugh[i, j, 8]
                Kennaugh0[3,1] <- Kennaugh[i, j, 9]
                Kennaugh0[3,2] <- Kennaugh[i, j, 10]
                Kennaugh0[3,3] <- Kennaugh[i, j, 11]
                Kennaugh0[3,4] <- Kennaugh[i, j, 12]
                Kennaugh0[4,1] <- Kennaugh[i, j, 13]
                Kennaugh0[4,2] <- Kennaugh[i, j, 14]
                Kennaugh0[4,3] <- Kennaugh[i, j, 15]
                Kennaugh0[4,4] <- Kennaugh[i, j, 16]
                K <- Kennaugh0
                
                # The geodesic distances between the observed K matrix and the elementary target matrixes are calculated
                GDs <- c(GDt <- GDFunction(K, Kt),
                         GDc <- GDFunction(K, Kc),
                         GDd <- GDFunction(K, Kd),
                         GDnd <- GDFunction(K, Knd))
                
                # The values of p and q, which are required for the calculation of β in GRVI, are determined
                p <- min(GDs)
                q <- max(GDs)
                
                # The K matrix of the volume scattering model is generated per pixel
                Kv0 <- matrix(0, 4, 4)
                Kv0[1,1] <- Kv[i, j, 1]
                Kv0[1,2] <- Kv[i, j, 2]
                Kv0[1,3] <- Kv[i, j, 3]
                Kv0[1,4] <- Kv[i, j, 4]
                Kv0[2,1] <- Kv[i, j, 5]
                Kv0[2,2] <- Kv[i, j, 6]
                Kv0[2,3] <- Kv[i, j, 7]
                Kv0[2,4] <- Kv[i, j, 8]
                Kv0[3,1] <- Kv[i, j, 9]
                Kv0[3,2] <- Kv[i, j, 10]
                Kv0[3,3] <- Kv[i, j, 11]
                Kv0[3,4] <- Kv[i, j, 12]
                Kv0[4,1] <- Kv[i, j, 13]
                Kv0[4,2] <- Kv[i, j, 14]
                Kv0[4,3] <- Kv[i, j, 15]
                Kv0[4,4] <- Kv[i, j, 16]
                KV <- Kv0
                
                # The geodesic distance between the observed K matrix and the K matrix of the volume scattering model (Kv)
                # needed to determine both β and fv in the GRVI is calculated.
                GDv <- GDFunction(K, KV)
                
                # Calculate β and fv
                fv <- 1 - GDv
                B <- (p/q)^(2*GDv)
                
                # Store GRVI
                GRVI[indexx, indexy,] <- B * fv
            }
        }
    # Save the result (this function returns the GRVI as a multilooking (x*y) array)
    return(GRVI)
}
