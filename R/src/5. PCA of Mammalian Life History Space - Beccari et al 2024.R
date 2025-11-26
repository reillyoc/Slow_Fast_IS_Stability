# The below code is directly from Beccari et al 2024., code to create the PC1 used to represent the mammalian slow-fast continuum

##### Beccari et al 2024 Code #####

# In the folowing script we will compute global mammalian life history space and we will check for the reliability of the so built space. 

# To built mammalian life history space: 
# 1. We performed a series of linear regressions between single log-transformed life history traits and log-transformed species' adult body mass.
# 2. We extracted body mass–corrected residuals for each species from each regression. 
# 3. We constructed the mammalian life history space by performing a principal component analysis (PCA). 
#
# We estimated the reliability of our methodological choices by comparing this life history space (i.e., the scores of species in the selected number of principal components) with three other spaces: 
# 1. a space built considering only complete set of records (i.e., non-imputed data; N = 1,293); 
# 2. a space built considering the imputed dataset without body mass correction (N = 3,438); 
# 3. a space built considering the same body mass–corrected residuals used to build the life history space (N = 3,438) but accounting for the phylogenetic history of species
# see Supplementary material 2 for further details. 

library(tidyverse)
library(ape)
library(phytools)


## Loading imputed set of Traits contained in the dataset mammalTraitsImputed
MammaUse <- read.table( "../R/Data/Beccari et al 2024/mammalTraitsImputed.txt")
#ls   = litter size
#ly   = litter per year
#bm   = body mass
#long = longevity
#gest = gestation length
#wea  = weaning length 
#fmat = age at female maturity
#General_system = classification in environmental realms: 
# SA = semi-aquatic; F = Aerial; W = Aquatic; T = Terrestrial 

# Let's remove the Phylogenetic Eigenvectors that are no more needed
MammTraitReadyScaled <- MammaUse[, c("ly","ls", "bm", "long", "gest", "wea", "fmat", "General_system")]



## Function to extract angles among PCA's eigenvectors
angle <- function(x, y){
  norm_x <- sqrt(sum(x^2))#sd of eigenvector
  norm_y <- sqrt(sum(y^2))
  cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8)#dot product/magnitude of vectors
  #round introduced to avoid numerical problems with the acos function
  angle <- acos(cosXY) * 180 / pi
  return(angle)
}

#################################################################################################################################
#################################################################################################################################
########### Creating mammalian life history space ##########
#################################################################################################################################
#################################################################################################################################
#### 1. Account for the effect of body mass on life history traits  __________________________________________________________#####

# We performed a series of linear regressions between single log-transformed life history traits and log-transformed species' adult body mass. 
# We extracted body mass–corrected residuals for each species from each regression.

trait <- c("ly","ls", "long", "gest", "wea", "fmat")
Resid.BM_TraitsLM <- MammTraitReadyScaled # Object to store the residuals 
Resid.BM_TraitsLM[, trait] <- rep(NA, times=nrow(Resid.BM_TraitsLM))
BM_model <- list()
for(t in trait){
  BM_model_Aux <- lm(MammTraitReadyScaled[,t] ~ MammTraitReadyScaled$bm) 
  Resid.BM_TraitsLM[,t] <- BM_model_Aux$residuals
  BM_model[[t]] <- BM_model_Aux
}
str(Resid.BM_TraitsLM)



#### 2. Creating the Global Life History space _______________________________________________________________________________####

# Let's construct the mammalian life history space by performing a principal component analysis (PCA) on the complete set of mammalian species.
# For that, we will use body mas corrected residuals
colnames(Resid.BM_TraitsLM) 

RawPCA_BMCorr <- princomp(scale(Resid.BM_TraitsLM[,c("ly","ls", "long", "gest", "wea", "fmat")]))
### SAVING PC Scores  ________________________________________________
PCA_BM_Corr <- Resid.BM_TraitsLM
PCA_BM_Corr[,"PC1"] <- RawPCA_BMCorr$scores[,1]
PCA_BM_Corr[,"PC2"] <- RawPCA_BMCorr$scores[,2]
PCA_BM_Corr$General_system <- as.factor(PCA_BM_Corr$General_system)
## These scores are the one that will be use across all set of analyses. 

### Let's extract the Loadings ___________________
loadingsPCA_Use <- RawPCA_BMCorr$loadings[, 1:6]
for(i in 1:ncol(loadingsPCA_Use)){
  loadingsPCA_Use[, i] <- loadingsPCA_Use[, i] * (sd(RawPCA_BMCorr$scores[,i]))
}
loadingsPCA_Use


### Let's extract the angles between life history traits  ________________________________________________
angMatPCA_Use <- matrix(NA, nrow = nrow(loadingsPCA_Use), ncol = ncol(loadingsPCA_Use), 
                        dimnames = list(rownames(loadingsPCA_Use), rownames(loadingsPCA_Use)))

for(i in 1:nrow(loadingsPCA_Use)){
  for(j in 1:nrow(loadingsPCA_Use)){
    angMatPCA_Use [i, j] <- angle(loadingsPCA_Use[i, 1:2], loadingsPCA_Use[j, 1:2])  ## we select the number of dimensions retained = 2 PCs
  }
}
round(angMatPCA_Use , 2)

###### Plotting Code #####
summary(RawPCA_BMCorr)

# Prepare data
pca_data <- data.frame(PC1 = PCA_BM_Corr$PC1,
                       PC2 = PCA_BM_Corr$PC2,
                       General_system = PCA_BM_Corr$General_system)

loadings <- data.frame(Variable = rownames(loadingsPCA_Use),
                       PC1 = loadingsPCA_Use[, 1],
                       PC2 = loadingsPCA_Use[, 2])

# Create the PCA plot
beccari_pca <- ggplot(data = pca_data, aes(x = PC1, y = PC2)) +
  # Scatter points for species
  geom_point(size = 2.5, alpha = 0.65, aes(color = PC1)) +
  scale_color_gradient(low = "#3CA373", high = "#EB8F00") +
  # Arrows for loadings
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4),
               arrow = arrow(length = unit(0.3, "cm")), size = 0.8, color = "black") +
  # Labels for loadings
  #geom_text(data = loadings, aes(x = PC1 * 5, y = PC2 * 5, label = Variable),
          #  size = 6, color = "black", fontface = "bold") +
  # Aesthetic enhancements
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid", size = 2),
          axis.ticks = element_line(linetype = "blank"),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "none",
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y =element_text(size = 16), 
          axis.title.x = element_text(size = 16), 
          text = element_text(family = "Arial")) +
  labs(
    x = paste0("PC1 (", "60.54% variance explained)"),
    y = paste0("PC2 (", "14.53% variance explained)")
  )

# Display the plot
print(beccari_pca)

#ggsave("../R/Figures/Figure 1 - PCA - Beccari et al 2024.jpeg", plot = beccari_pca, width = 8, height = 6, dpi = 300)

