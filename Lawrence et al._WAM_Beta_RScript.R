# Make sure to first install the soundecology, data.table, FactoMineR, openxlsx packages from CRAN repository

# Load the required packages and libraries
p <- c("parallel", "pracma", "oce", "ineq", "gsw", "seewave", "soundecology", "data.table", "tuneR", "vegan", "openxlsx", "FactoMineR")
lapply(p, require, character.only = TRUE);

oldwd <- getwd()
setwd("Insert WAV data directory")
files <- dir(pattern = "wav$")
n <- length(files);

# Get the mean spectrum of all wav files
wl <- 512
mspectra <- matrix(NA, nrow=wl/2, ncol=n)
for(i in 1:n) mspectra[,i] <- meanspec(readWave(files[i]),
                                        wl=wl, plot=FALSE)[,2]
comb <- combn(1:n, 2);


# Create a matrix for DCF for all combinations of wav files
dcf <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dcf[comb[2,i], comb[1,i]] <- diffcumspec(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            plot=FALSE)};

# Put filenames as row and column headings
colnames(dcf) <- rownames(dcf) <- unlist(strsplit(files, split=".wav"))

# Fill the upper triangle of the matrix
dcf[upper.tri(dcf)] <- t(dcf)[upper.tri(dcf)]

# Assign zeroes to diagonals of the matrix
diag(dcf) <- 0

# Convert the matrix to a data frame
mydataframedcf <- as.data.frame(dcf)
setDT(mydataframedcf, keep.rownames = "Filename")



dks <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dks[comb[2,i], comb[1,i]] <- ks.dist(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            f=44100, plot=FALSE)$D};
colnames(dks) <- rownames(dks) <- unlist(strsplit(files, split=".wav"))

dks[upper.tri(dks)] <- t(dks)[upper.tri(dks)]

diag(dks) <- 0

mydataframedks <- as.data.frame(dks)
setDT(mydataframedks, keep.rownames = "Filename")



dis <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dis[comb[2,i], comb[1,i]] <- itakura.dist(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            )$D};
colnames(dis) <- rownames(dis) <- unlist(strsplit(files, split=".wav"))

dis[upper.tri(dis)] <- t(dis)[upper.tri(dis)]

diag(dis) <- 0

mydataframedis <- as.data.frame(dis)
setDT(mydataframedis, keep.rownames = "Filename")



dkl <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dkl[comb[2,i], comb[1,i]] <- kl.dist(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            )$D};
colnames(dkl) <- rownames(dkl) <- unlist(strsplit(files, split=".wav"))

dkl[upper.tri(dkl)] <- t(dkl)[upper.tri(dkl)]

diag(dkl) <- 0

mydataframedkl <- as.data.frame(dkl)
setDT(mydataframedkl, keep.rownames = "Filename")



mi <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     mi[comb[2,i], comb[1,i]] <- 1-symba(round(mspectra[,comb[1,i]],2),
                                            round(mspectra[,comb[2,i]],2),
                                            plot=FALSE)$I};
colnames(mi) <- rownames(mi) <- unlist(strsplit(files, split=".wav"))

mi[upper.tri(mi)] <- t(mi)[upper.tri(mi)]

diag(mi) <- 0

mydataframemi <- as.data.frame(mi)
setDT(mydataframemi, keep.rownames = "Filename")



dls <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dls[comb[2,i], comb[1,i]] <- logspec.dist(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            scale=TRUE)};
colnames(dls) <- rownames(dls) <- unlist(strsplit(files, split=".wav"))

dls[upper.tri(dls)] <- t(dls)[upper.tri(dls)]

diag(dls) <- 0

mydataframedls <- as.data.frame(dls)
setDT(mydataframedls, keep.rownames = "Filename")



s <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     s[comb[2,i], comb[1,i]] <- 100-simspec(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            plot=FALSE)};
colnames(s) <- rownames(s) <- unlist(strsplit(files, split=".wav"))

s[upper.tri(s)] <- t(s)[upper.tri(s)]

diag(s) <- 0

mydataframes <- as.data.frame(s)
setDT(mydataframes, keep.rownames = "Filename");



r <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     r[comb[2,i], comb[1,i]] <- sqrt(1-cor(mspectra[,comb[1,i]]/sum(mspectra[,comb[1,i]]),
                                            mspectra[,comb[2,i]]/sum(mspectra[,comb[2,i]])))};
colnames(r) <- rownames(r) <- unlist(strsplit(files, split=".wav"))

r[upper.tri(r)] <- t(r)[upper.tri(r)]

diag(r) <- 0

mydataframer <- as.data.frame(r)
setDT(mydataframer, keep.rownames = "Filename")



df <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     df[comb[2,i], comb[1,i]] <- diffspec(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            plot=FALSE)};
colnames(df) <- rownames(df) <- unlist(strsplit(files, split=".wav"))

df[upper.tri(df)] <- t(df)[upper.tri(df)]

diag(df) <- 0



dt <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     dt[comb[2,i], comb[1,i]] <- diffenv(mspectra[,comb[1,i]],
                                            mspectra[,comb[2,i]],
                                            f=44100)};
colnames(dt) <- rownames(dt) <- unlist(strsplit(files, split=".wav"))

dt[upper.tri(dt)] <- t(dt)[upper.tri(dt)]

diag(dt) <- 0


# Calculate the matrix for D by multiplying Df and Dt matrices
d <- df*dt;


mydataframedf <- as.data.frame(df)
setDT(mydataframedf, keep.rownames = "Filename")
mydataframedt <- as.data.frame(dt)
setDT(mydataframedt, keep.rownames = "Filename")
mydataframed <- as.data.frame(d)
setDT(mydataframed, keep.rownames = "Filename");



# Increase the hanning window to calculate RV
wl <- 1024
mspectra <- matrix(NA, nrow=wl/2, ncol=n)
for(i in 1:n) mspectra[,i] <- meanspec(readWave(files[i]),
                                        wl=wl, plot=FALSE)[,2]
rv <- matrix(NA, nrow=n, ncol=n)
for(i in 1:ncol(comb)){
     rv[comb[2,i], comb[1,i]] <- 1-coeffRV(spectro(mspectra[,comb[1,i]], f=44100, plot=FALSE)$amp,
                                           spectro(mspectra[,comb[2,i]], f=44100, plot=FALSE)$amp)$rv};
colnames(rv) <- rownames(rv) <- unlist(strsplit(files, split=".wav"))

rv[upper.tri(rv)] <- t(rv)[upper.tri(rv)]

diag(rv) <- 0

mydataframerv <- as.data.frame(rv)
setDT(mydataframerv, keep.rownames = "Filename")



# Create an excel table with different index matrices in different sheets
list_of_datasets <- list("DCF" = mydataframedcf, "DKS" = mydataframedks, "DIS" = mydataframedis, "DKL" = mydataframedkl, 
"MI" = mydataframemi, "DLS" = mydataframedls, "S" = mydataframes, "R" = mydataframer, "DF" = mydataframedf, 
"DT" = mydataframedt, "D" = mydataframed, "RV" = mydataframerv);

write.xlsx(list_of_datasets, file = "Insert WAV data pathway/WAM_Beta.csv")
