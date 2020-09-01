# Make sure to first install the soundecology and data.table packages from CRAN repository

# Load the required packages and libraries
p <- c("parallel", "pracma", "oce", "ineq", "gsw", "seewave", "soundecology", "data.table", "openxlsx", "tuneR", "vegan")
lapply(p, require, character.only = TRUE);

# Set the working directory to where your WAV data is stored and retrieve the wave files.  Replace \ with / in pathways.
oldwd <- getwd()
setwd("insert WAV directory pathway")
files <- dir(pattern = "wav$");

# Delete the existing csv file from the target folder
toDelete <- list.files("C:/SALVE/DAP Winter Remaining WAVs/DAP Winter Remaining WAVs", pattern = ".xlsx")
file.remove(toDelete);

# Prepare the function to calculate eight alpha indices
indices <- function(x){
    x <- readWave(x)
    return(
             c(
             sh(meanspec(x, plot=FALSE)),
             acoustic_evenness(x)$aei_left,
             ACI(x),
             acoustic_diversity(x)$adi_left,
             bioacoustic_index(x)$left_area,
             H(x),
             nrow(fpeaks(meanspec(x, plot=TRUE))),
             ndsi(x, fft_w = 512, anthro_min = 1000, anthro_max = 2000, bio_min = 2000, bio_max = 8000)$ndsi_left))};

# Prepare the first object where the results will be written in
n <- length(files)
num <- rep(NA, n)
res <- data.frame(Hf=num, AEI=num, ACI=num, ADI=num, BIO=num, H=num, NP=num, NDSI=num, row.names=files);

# Run a for loop with the prepared function to calculate the eight indices for all wav files in the directory
for(i in 1:n) try(res[i,] <- indices(files[i]));

# Calculate and store the remaining three alpha indices AR, M and Ht in a second object
ARI <- AR(getwd(), datatype="files");

# Join the three alpha indices from the second object to the first object containing the rest of the indices
res$M <- ARI$M
res$Ht <- ARI$Ht
res$AR <- ARI$AR;

setDT(res, keep.rownames = "Filename");

# Export all eleven alpha indices to an excel file
write.xlsx(res, file="insert WAV directory pathway/WAM_Alpha.xlsx", row.names=TRUE);


