"
This script tests the GRVIFunction with an clip of a SAOCOM image. The image is not one of the images used in the thesis, 
because SAOCOM data cannot be made publicly available. Instead, a cropping of an image authorised by the institution that 
manages the SAOCOM data was used. Find information about the SAOCOM image in the TXT file accompanying the image files.

Input data:
-real and imaginary parts of the SAOCOM image polarizations

Author: Elisa Frank Buss
Encoding: UTF-8
Creation date: May 2022
Last modification: Jul 2022

"
### Calculate GRVI

# Directory
wd <- setwd("D:/Documents/GitHub/Tesis-MAIE-Elisa-Frank-Buss/")

# Required functions
source("Code/R/libraries/imagematrix.R")
source("Code/R/libraries/myread.ENVI.R")
source("Code/R/libraries/PauliKennaughMultilook.R")
source("Code/R/libraries/gammaFunction.R")
source("Code/R/libraries/vectorVFunction.R")
source("Code/R/libraries/GDFunction.R")
source("Code/R/libraries/GRVIFunction.R")

# Path
path <- "Data/SAOCOM/subset_0_of_S1A_OPER_SAR_EOSSP__CORE_L1A_OLVF_20200703T145917_TC.data"

# Read the real and imaginary bands
R_HH <- myread.ENVI(paste(path,"/i_HH.img", sep=""), paste(path,"/i_HH.hdr", sep=""))
I_HH <- myread.ENVI(paste(path,"/q_HH.img", sep=""), paste(path,"/q_HH.hdr", sep=""))
R_HV <- myread.ENVI(paste(path,"/i_HV.img", sep=""), paste(path,"/i_HV.hdr", sep=""))
I_HV <- myread.ENVI(paste(path,"/q_HV.img", sep=""), paste(path,"/q_HV.hdr", sep=""))
R_VV <- myread.ENVI(paste(path,"/i_VV.img", sep=""), paste(path,"/i_VV.hdr", sep=""))
I_VV <- myread.ENVI(paste(path,"/q_VV.img", sep=""), paste(path,"/q_VV.hdr", sep=""))

# Obtain matrix dimensions
dimensions <- dim(R_HH)

# Determine the dimensions of the output matrix
dimensions_2 <- floor(dimensions/2)

# Join the real and imaginary bands together to form the matrices SHH, SHV, SVV
SHH <- matrix(complex(real=R_HH, imaginary=I_HH), nrow=dim(R_HH)[1])
SHV <- matrix(complex(real=R_HV, imaginary=I_HV), nrow=dim(R_HV)[1])
SVV <- matrix(complex(real=R_VV, imaginary=I_VV), nrow=dim(R_VV)[1])

rm(I_HH,I_VV,I_HV,R_HH,R_VV,R_HV)

# Calculate GRVI
GRVI <- GRVIFunction(SHH, SHV, SVV, 2, 2)

# Plot the GRVI result
plot(
    imagematrix(
        equalize(GRVI[1:floor(dimensions_2[1]),1:floor(dimensions_2[2]),1])
    )
)


### Export GRVI

#Import library
library(raster)

# Open a raster (.tif) of the clipping to obtain its extent and CRS
cuadricula <- raster(paste0(wd, "/Data/SAOCOM/TIF/subset_0_of_S1A_OPER_SAR_EOSSP__CORE_L1A_OLVF_20200703T145917_TC_HH.tif"))

# Obtain CRS
crs <- proj4string(cuadricula)

# Obtain extent
coor1 = cuadricula@extent@xmin
coor2 = cuadricula@extent@xmax
coor3 = cuadricula@extent@ymin
coor4 = cuadricula@extent@ymax

rm(cuadricula)

# Convert GRVI from array to matrix
GRVI_m <- GRVI[,,1]

# Convert GRVI from matrix to raster
GRVI_raster <- raster(GRVI_m)

# Apply a 5x5 focal filter
GRVI_5x5 <- focal(GRVI_raster, w = matrix(1, 5, 5), mean)

# Assign the raster extent
extent(GRVI_5x5) <- c(coor1, coor2, coor3, coor4)

# Assign the CRS of the raster
projection(GRVI_5x5) <- crs

#Export GRVI as.tif
writeRaster(GRVI_5x5, 'Data/SAOCOM/GRVI_5x5.tif', overwrite = TRUE)
