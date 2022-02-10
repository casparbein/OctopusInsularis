## TESS3R analysis
## load libraries 
library(tess3r) 
library(LEA) 
library(maps)
library(raster) 
library(marmap) 
library(rworldmap)
library(ggplot2) 
library("ggspatial")

## get and set wd
setwd("~/octopus/tess3r")

## In order to use the Tess3 clustering algorithm, input vcf files have to converted into the
## lfmm format also used by the LEA-package. This is done by the function vcf2lfmm from the LEA package.
## vcf files filtered for singletons, biallelic sites and one SNP per locus:
## 40MD: 69inds_40MD_filt.vcf, name in this script: perm (for permissive)
## 20MD: 64inds_20MD_filt.vcf, name in this script: stri (for stringent)

vcf2lfmm("~/clean_vcfs/69inds_40MD.recode.vcf", 
         "69_40MD.lfmm")

vcf2lfmm("~/clean_vcfs/64inds_20MD.recode.vcf", 
         "64_20MD.lfmm")

## After conversion, the lfmm-files can be read into R by the read.lfmm function, 
## also from the LEA package.

## load the stri-dataset
stri_genotype <- read.lfmm("~/tess3r/64inds_20MD.recode.lfmm")
## take a look at part of the matrix
stri_genotype[1:64,1:100]

## load the perm_dataset
perm_genotype <- read.lfmm("~/tess3r/69inds_40MD.recode.lfmm")
## take a look at part of the matrix
perm_genotype[1:69,1:100]

## Tess3R cannot interpret 9 (for missing value in the Matrix), so I changed 9 to NA.
stri_genotype[stri_genotype==9] <- NA
perm_genotype[perm_genotype==9] <- NA

## check number of columns of genotype matrix
length(colnames(stri_genotype))
length(colnames(perm_genotype))

## Tess3R needs a coordinates file to set geographic "priors". 
## load coordinates file
## check if CEARA GPS data is correct??

stri_coordinates <- read.table("~/R scripts/coordinates_64.txt", 
                          header = TRUE)
View(stri_coordinates)

perm_coordinates <- read.table("~/R scripts/coordinates_69.txt", 
                          header = TRUE)
View(perm_coordinates)


## plot the coordinates in a mock map to see if they are correct
plot(perm_coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

## create a tess3 object for perm:
## "keep all" ensures that all runs, not only the best one, are kept

perm_tess3_nofil_11.obj <- tess3(X = perm_genotype, coord = as.matrix(perm_coordinates), 
                               K = 1:11, 
                               method = "projected.ls", 
                               ploidy = 2, 
                               openMP.core.num = 4,
                               max.iteration = 1000,
                               rep = 20,
                               keep ="all"
)

## and for stri
stri_tess3_nofil_11.obj <- tess3(X = stri_genotype, coord = as.matrix(stri_coordinates), 
                        K = 1:11, 
                        method = "projected.ls", 
                        ploidy = 2, 
                        openMP.core.num = 4,
                        max.iteration = 1000,
                        rep = 20,
                        keep ="all"
)


## plot Cross validation scores of the tess3 objects 
## (indication of which K might best fit the data).
cv_64 <- plot(stri_tess3_nofil_11.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main = "Dataset: 64inds_20MD")

cv_69 <- plot(perm_tess3_nofil_11.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main = "Dataset: 69inds_40MD"
     )


## create q-matrices:
## q-matrices are stored in tess3 objects. To work on a particular q-matrix outside of tess3, 
## one has to first define the q_matrices for each K separately. 

## q.matrices for K=1 - K=8, perm-dataset

perm_q.matrix_1 <- qmatrix(perm_tess3_nofil_11.obj, K = 1)
perm_q.matrix_2 <- qmatrix(perm_tess3_nofil_11.obj, K = 2)
perm_q.matrix_3 <- qmatrix(perm_tess3_nofil_11.obj, K = 3)
perm_q.matrix_4 <- qmatrix(perm_tess3_nofil_11.obj, K = 4)
perm_q.matrix_5 <- qmatrix(perm_tess3_nofil_11.obj, K = 5)
perm_q.matrix_6 <- qmatrix(perm_tess3_nofil_11.obj, K = 6)
perm_q.matrix_7 <- qmatrix(perm_tess3_nofil_11.obj, K = 7)
perm_q.matrix_8 <- qmatrix(perm_tess3_nofil_11.obj, K = 8)
perm_q.matrix_9 <- qmatrix(perm_tess3_nofil_11.obj, K = 9)
perm_q.matrix_10 <- qmatrix(perm_tess3_nofil_11.obj, K = 10)
perm_q.matrix_11 <- qmatrix(perm_tess3_nofil_11.obj, K = 11)

## same for stri-dataset.

stri_q.matrix_1 <- qmatrix(stri_tess3_nofil_11.obj, K = 1)
stri_q.matrix_2 <- qmatrix(stri_tess3_nofil_11.obj, K = 2)
stri_q.matrix_3 <- qmatrix(stri_tess3_nofil_11.obj, K = 3)
stri_q.matrix_4 <- qmatrix(stri_tess3_nofil_11.obj, K = 4)
stri_q.matrix_5 <- qmatrix(stri_tess3_nofil_11.obj, K = 5)
stri_q.matrix_6 <- qmatrix(stri_tess3_nofil_11.obj, K = 6)
stri_q.matrix_7 <- qmatrix(stri_tess3_nofil_11.obj, K = 7)
stri_q.matrix_8 <- qmatrix(stri_tess3_nofil_11.obj, K = 8)
stri_q.matrix_9 <- qmatrix(stri_tess3_nofil_11.obj, K = 9)
stri_q.matrix_10 <- qmatrix(stri_tess3_nofil_11.obj, K = 10)
stri_q.matrix_11 <- qmatrix(stri_tess3_nofil_11.obj, K = 11)

## STRUCTURE-like barplot for the Q-matrix (ancestry coefficients), plotted in tess3R.
## Problems: Colors change for each K, no consistency in clusters.
## Final plotting will be done with the pophelper package.
## example plot for K = 4, perm_dataset
barplot(stri_q.matrix_4, border = "black", space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix", sort.by.Q = F) -> bp

## One feature of Tess3R is geographic extrapolation, so basically a map, where different color shades
## indicate assignment probability to a certain cluster, if an individual had been sampled from that area.
## Usually, tess3R only displays these shades on land (elevation above sea level >=0), but if I invert 
## the elevation pattern so that areas below 0 are interpreted as above 0 and viceversa, 
## it also works with oceans and other water bodies. 
## The method we used has been proposed by the package developer on github.

## get the map from NOAA (it is free access), filename of file I created:
## 'marmap_coord_-85;-22;10;15_res_10.csv'

map_bathy <- marmap::getNOAA.bathy(lon1 = -85, lon2 = 10, lat1 = -22,
                                   lat2 = 15, res=10,
                                   keep=T)

## invert land/ocean surface elevation values

map_bathy1 = - map_bathy

## It is possible to exclude certain depths and areas manually so they are not taken into account 
## by tess3 R. Since there are no Octopus insularis in the Pacific, I excluded the Pacific by assigning
## elevation values of -1 to latitudes/longitudes pairs that fall into this region.

## exclude Pacific
map_bathy1[1:100,1:180] <- -1 

## depth control: exclude areas that are deeper than 4000 meters 
## (in the end, we did not use this feature)
for (i in 1:570) {
        for (j in 1:220) {
                if (map_bathy1[i,j] > 4000 | map_bathy1[i,j] < 0) {
                        map_bathy1[i,j] <- -1  
                }
        }
}

summary(map_bathy1)

## convert the bathymetry map to a raster

asc_raster <- marmap::as.raster(map_bathy1)

## Write raster to file for later use

raster::writeRaster(asc_raster, "~/octopus/raster files/brazil_test_10.asc", 
                    overwrite=TRUE)

## Test plot of the perm-dataset with adapted raster and K = 4

plot(perm_q.matrix_7, perm_coordinates, method="map.max", cex=1,
     raster.filename = "~/octopus/raster files/brazil_test_10.asc", 
     interpol=FieldsKrigModel(9.5),
     main=paste0("Ancestry coefficients"), resolution=c(300,300),
     xlab="longitude", ylab="Latitude")

## I assembled the final plot with ggplot, and the ggtess3Q function from the tess3R package.
## first, I called a costum color palette. Also, I need a world map as a backdrop for the plot,
## which can be created using the rworldmap package. 

## Get a map
testmap <- getMap(resolution = "high", projection=NA)

## create a color palette and shades of the color palette with CreatePalette
my.colors <- c("#1D72F5","#DF0101","#77CE61","#FF9326","#A945FF",
                "#0089B2","#FDF060","#FFA6B2")
my.palette <- CreatePalette(my.colors, 15)


## test plot for K = 4, perm-dataset
pl <- ggtess3Q(perm_q.matrix_4, perm_coordinates, 
               raster.filename = "~/octopus/raster files/brazil_test_10.asc",
               col.palette = my.palette)
pl4 <- pl + geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
        xlim(-80, -10) + 
        ylim(-21, 12) + 
        coord_equal()  +
        geom_point(data = as.data.frame(perm_coordinates), aes(x = V1, y = V2), size = 3, 
                   color = "black") + 
        xlab("Longitute") +
        ylab("Latitude") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
        scale_colour_continuous(my.palette) +
        annotation_scale(location = "bl", width_hint = 0.15) 

pl4

## Save q-matrices as text files and align the clusters across files,
## so that in the plots one set colour always indicates the same cluster.
## For the perm-dataset

## get all q-matrices as data frames for handling
perm_tess1  <- as.data.frame(perm_q.matrix_1)
perm_tess2  <- as.data.frame(perm_q.matrix_2)
perm_tess3  <- as.data.frame(perm_q.matrix_3)
perm_tess4  <- as.data.frame(perm_q.matrix_4)
perm_tess5  <- as.data.frame(perm_q.matrix_5)
perm_tess6  <- as.data.frame(perm_q.matrix_6)
perm_tess7  <- as.data.frame(perm_q.matrix_7)
perm_tess8  <- as.data.frame(perm_q.matrix_8)
perm_tess9  <- as.data.frame(perm_q.matrix_9)
perm_tess10  <- as.data.frame(perm_q.matrix_10)
perm_tess11 <- as.data.frame(perm_q.matrix_11)


## function to plot a q-matrix
q_mat <- function(dat.frame, labs = F){
        ## assign names
        ## remove filtered out individuals: SPS3, FN5
        dat.frame$inds <- s4_frame$name[c(-40,-57)]
        
        ## reorder individuals based on geographic origin:
        ## OIC - CE - RN - AR - FN - SPS - ASC - STH - AL - BA - TM
        dat.frame <- dat.frame[c(42:44, 26:35, 45:53,  5:14, 36:41,54:56, 15:16, 57, 1:4, 17:25, 58:69),]
        dat.frame$inds <- factor(dat.frame$inds, levels = as.character(dat.frame$inds))
        
        ## melt dataframe
        cols_name <- colnames(dat.frame)
        dat.frame_longer <- pivot_longer(dat.frame, cols = cols_name[1]:cols_name[ncol(dat.frame)-1] , 
                                         names_to = "cluster", 
                                         values_to = "ancestry")
        
        ## clusters as factors
        dat.frame_longer$cluster <- factor(dat.frame_longer$cluster, levels= unique(dat.frame_longer$cluster))
        
        ## plot
        if (labs==T) { ## individual labels on x-axis
        plot <- ggplot(data = dat.frame_longer, aes(inds, ancestry)) + 
                geom_col(aes(fill = cluster), color=NA) + 
                theme_classic() + 
                theme(axis.text.x = element_text(angle = 90, 
                                                 vjust = 0.5, 
                                                 hjust=1, size = 5),
                      axis.title.x=element_blank(),
                      legend.position="none" #,
                      # axis.line=element_blank()
                      ) +
                scale_fill_discrete(type=c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                                           "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                                           "#4a805e", "#cdc0b0"))
        } else { ## no labels
        plot <- ggplot(data = dat.frame_longer, aes(inds, ancestry)) + 
                geom_col(aes(fill = cluster)) + 
                theme_classic() + 
                theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x=element_blank(),
                      axis.line=element_blank(),
                      legend.position="none") + 
                scale_fill_discrete(type=c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                                           "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                                           "#4a805e", "#cdc0b0"))
                
        }
        
        return(plot)
        
        
        
        
}


## helper function for a data frame in long/melted format
q_frame <- function(dat.frame){
        ## assign names
        s4_data <- read.delim("~/octopus/final stats clean/s4_joint_estimate.txt",sep="")
        s4_frame$name[c(-40,-57)]
        dat.frame$inds <- s4_frame$name[c(-40,-57)]
        dat.frame <- dat.frame[c(42:44, 26:35, 45:53,  5:14, 36:41,54:56, 15:16, 57, 1:4, 17:25, 58:69),]
        dat.frame$inds <- factor(dat.frame$inds, levels = as.character(dat.frame$inds))
        
        ## melt dataframe
        cols_name <- colnames(dat.frame)
        dat.frame_longer <- pivot_longer(dat.frame, cols = cols_name[1]:cols_name[ncol(dat.frame)-1] , 
                                         names_to = "cluster", 
                                         values_to = "ancestry")
        
        return(dat.frame_longer)
        
        
        
        
}

## plot preliminary q-matrices
qmat_tess1<- q_mat(perm_tess1)
qmat_tess2<- q_mat(perm_tess2)
qmat_tess3<- q_mat(perm_tess3)
qmat_tess4<- q_mat(perm_tess4)
qmat_tess5<- q_mat(perm_tess5)
qmat_tess6<- q_mat(perm_tess6)
qmat_tess7<- q_mat(perm_tess7)
qmat_tess8<- q_mat(perm_tess8)
qmat_tess9<- q_mat(perm_tess9)
qmat_tess10 <- q_mat(perm_tess10)
qmat_tess11 <- q_mat(perm_tess11)

qmat_tess1
qmat_tess2
qmat_tess3
qmat_tess4
qmat_tess5
qmat_tess6
qmat_tess7
qmat_tess8
qmat_tess9
qmat_tess10
qmat_tess11

## create a list with all q-matrices to loop over
all_tess <- list(perm_tess1, perm_tess2, perm_tess3, perm_tess4,
                 perm_tess5, perm_tess6, perm_tess7, perm_tess8,
                 perm_tess9, perm_tess10, perm_tess11)
length(all_tess)

## a helper list for the loop
helper <- list()

for (i in 1:11) {
        temp <- all_tess[[i]] ## grab the i-th q-matrix
        temp_cols <- ncol(temp) ## get the number of clusters in this q-matrix
        checklist <- c(which(temp[26,] > 0.5), ## cluster 1
                       which(temp[58,] > 0.5), ## cluster 2
                       which(temp[56,] > 0.5), ## cluster 3
                       which(temp[1,] > 0.5), ## cluster 4
                       which(temp[43,]>0.5), ## cluster 5 
                       which(temp[15,]>0.5), ##cluster 6
                       which(temp[5,]>0.5), ## cluster 7 
                       which(temp[44,]>0.5), ## cluster 8
                       which(temp[17,]>0.5), ## cluster 9
                       which(temp[45,]>0.5), ## cluster 10
                       which(temp[36,]>0.5)) ## cluster 11   
        temp_index <- checklist[1:i] ## find the i-th item in checklist, 
        ##which corresponds to the i-th cluster
        helper[[i]] <- as.data.frame(temp)[,temp_index] ## the i-th entry of helper is 
        ## the reordered i-th input
        if (i>1){
        colnames(helper[[i]]) <- c(seq(1,i,1)) ## column names are numeric and will be ordered correctly in 
        ## overview figures 
        }
        if (i<10){
        filename <- paste("tess69_new0",i, ".txt", sep = "") ## ordering in the folder needs a 0 
        ## in front of one digit numbers
        } else {
        filename <- paste("tess69_new",i, ".txt", sep = "")
        }
        write.table(helper[[i]], ## write output file
                    file = paste("~/octopus/tess3r/tess_manual_69_40/cluster_11/", filename, sep = ""), 
                    append = FALSE, sep = "\t", dec = ".",
                    row.names = F, col.names = F)
}

## ordered q-matrices are stored in the "helper" list
## example:
helper[[11]]
q_mat(helper[[8]], labs = T)


## plot all Q-matrices

##69
library(cowplot)
ggdraw() + 
        draw_plot(q_mat(helper[[2]], labs = F), x = 0, y = 0.91, height = 0.09, width =0.9) +
        draw_plot(q_mat(helper[[3]], labs = F), x = 0, y = 0.82, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[4]], labs = F), x = 0, y = 0.73, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[5]], labs = F), x = 0, y = 0.64, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[6]], labs = F), x = 0, y = 0.55, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[7]], labs = F), x = 0, y = 0.46, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[8]], labs = F), x = 0, y = 0.37, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[9]], labs = F), x = 0, y = 0.28, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[10]], labs = F), x = 0, y = 0.19, height = 0.09, width =0.9) + 
        draw_plot(q_mat(helper[[11]], labs = T), x = 0, y = 0.05, height = 0.14, width =0.9) +
        draw_plot_label(
                c("K= 2", "K= 3", "K= 4", "K= 5", "K= 6", "K= 7", "K= 8", "K= 9", "K=10", "K=11"),
                c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)-0.01,
                c(0.99, 0.90, 0.81, 0.72, 0.63, 0.54, 0.45, 0.36, 0.27, 0.18)-0.03,
                size = 15
        )


## same for the stri-dataset

stri_tess1  <- as.data.frame(stri_q.matrix_1)
stri_tess2  <- as.data.frame(stri_q.matrix_2)
stri_tess3  <- as.data.frame(stri_q.matrix_3)
stri_tess4  <- as.data.frame(stri_q.matrix_4)
stri_tess5  <- as.data.frame(stri_q.matrix_5)
stri_tess6  <- as.data.frame(stri_q.matrix_6)
stri_tess7  <- as.data.frame(stri_q.matrix_7)
stri_tess8  <- as.data.frame(stri_q.matrix_8)
stri_tess9  <- as.data.frame(stri_q.matrix_9)
stri_tess10  <- as.data.frame(stri_q.matrix_10)
stri_tess11  <- as.data.frame(stri_q.matrix_11)


## design q_mat function for 64 inds

q_mat_64 <- function(dat.frame, labs = F){
        ## assign names
        s4_data <- read.delim("~/octopus/final stats clean/s4_joint_estimate.txt",sep="")
        
        ## remove individuals filtered out (OIC1, STH2, RN2A, RN13, BA18)
        dat.frame$inds <- s4_frame$name[c(-22, -40, -43, -49, -52, -57, -59 )]
        
        ## reorder individuals based on geographic origin
        ## ## OIC - CE - RN - AR - FN - SPS - ASC - AL - BA - TM
        dat.frame <- dat.frame[c(41:42, 25:34, 43:49, 5:14, 35:40, 50:52, 15:16, 1:4, 17:24, 53:64),]
        dat.frame$inds <- factor(dat.frame$inds, levels = as.character(dat.frame$inds))
        
        ## melt dataframe
        cols_name <- colnames(dat.frame)
        dat.frame_longer <- pivot_longer(dat.frame, cols = cols_name[1]:cols_name[ncol(dat.frame)-1] , 
                                         names_to = "cluster", 
                                         values_to = "ancestry")
        
        ## clusters as factors
        dat.frame_longer$cluster <- factor(dat.frame_longer$cluster, levels= unique(dat.frame_longer$cluster))
        
        ## plot
        
        if (labs == F) {
                plot <- ggplot(data = dat.frame_longer, aes(inds, ancestry)) + 
                        geom_col(aes(fill = cluster)) + 
                        theme_classic() + 
                        theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x=element_blank(),
                              axis.line=element_blank(),
                              legend.position="none") + 
                        scale_fill_discrete(type=c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                                                   "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                                                   "#4a805e", "#cdc0b0"))
                
        } else {
                plot <- ggplot(data = dat.frame_longer, aes(inds, ancestry)) + 
                        geom_col(aes(fill = cluster)) + 
                        theme_classic() + theme(#text = element_text(size=15),
                                        axis.text.x = element_text(angle = 90, 
                                                                   vjust = 0.5, 
                                                                   hjust=1, size = 5),
                                        #axis.line=element_blank(),
                                        legend.position="none"
                                        ) +
                        scale_fill_discrete(type=c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                                                   "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                                                   "#4a805e", "#cdc0b0"))
                
        }

        return(plot)

        
}


## Order Q-matrices for consistency 

all_tess_64 <- list(stri_tess1, stri_tess2, stri_tess3, stri_tess4,
                 stri_tess5, stri_tess6, stri_tess7, stri_tess8,
                 stri_tess9, stri_tess10, stri_tess11)
length(all_tess_64)

helper_64 <- list()

## for some reason, cluster 9 disappears at K = 10 

for (i in 1:11) {
        if (i < 10 | i > 10) {
        temp <- all_tess_64[[i]]
        temp_cols <- ncol(temp)
        checklist <- c(which(temp[26,]>0.5),  ##1
                       which(temp[53,]>0.5),  ##2
                       which(temp[50,]>0.5),  ##3
                       which(temp[1,]>0.5),   ##4
                       which(temp[41,]>0.4),  ##5
                       which(temp[15,]>0.5),  ##6
                       which(temp[5,]>0.5),   ##7
                       which(temp[42,]>0.5),  ##8
                       which(temp[17,]>0.5),  ##9
                       which(temp[35,]>0.5),  ##10
                       which(temp[43,]>0.5)   ##11
                       )
        } else {
                temp <- all_tess_64[[i]]
                temp_cols <- ncol(temp)
                checklist <- c(which(temp[26,]>0.5),  ##1
                               which(temp[53,]>0.5),  ##2
                               which(temp[50,]>0.5),  ##3
                               which(temp[1,]>0.5),   ##4
                               which(temp[41,]>0.4),  ##5
                               which(temp[15,]>0.5),  ##6
                               which(temp[5,]>0.5),   ##7
                               which(temp[42,]>0.5),  ##8
                               which(temp[35,]>0.5),  ##10
                               which(temp[43,]>0.5),   ##11
                               which(temp[17,]>0.5))  ##9
                
        }
        temp_index <- checklist[1:i]
        helper_64[[i]] <- as.data.frame(temp)[,temp_index]
        if (i>1){
                colnames(helper_64[[i]]) <- c(seq(1,i,1))
        }
        if (i<10){
                filename <- paste("tess64_new0",i, ".txt", sep = "")
        } else {
                filename <- paste("tess64_new",i, ".txt", sep = "")
        }
        write.table(helper_64[[i]],
                    file = paste("~/octopus/tess3r/tess_manual_64_20/cluster_11/", filename, sep = ""), 
                    append = FALSE, sep = "\t", dec = ".",
                    row.names = F, col.names = F)
}



## plot all Q-matrices for 64 inds
ggdraw() + 
        draw_plot(q_mat_64(helper_64[[2]], labs = F), x = 0, y = 0.91, height = 0.09, width =0.9) +
        draw_plot(q_mat_64(helper_64[[3]], labs = F), x = 0, y = 0.82, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[4]], labs = F), x = 0, y = 0.73, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[5]], labs = F), x = 0, y = 0.64, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[6]], labs = F), x = 0, y = 0.55, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[7]], labs = F), x = 0, y = 0.46, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[8]], labs = F), x = 0, y = 0.37, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[9]], labs = F), x = 0, y = 0.28, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[10]], labs = F), x = 0, y = 0.19, height = 0.09, width =0.9) + 
        draw_plot(q_mat_64(helper_64[[11]], labs = T), x = 0, y = 0.05, height = 0.14, width =0.9) +
        draw_plot_label(
                c("K= 2", "K= 3", "K= 4", "K= 5", "K= 6", "K= 7", "K= 8", "K= 9", "K=10", "K=11"),
                c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)-0.01,
                c(0.99, 0.90, 0.81, 0.72, 0.63, 0.54, 0.45, 0.36, 0.27, 0.18)-0.03,
                size = 15
        )

## After realigning and writing the q-matrices to txt.-files, now I want to plot them next to each other.
## I use the pophelper package to draw a stacked barplot of K=1 to K=11. 
## Individually designing plots seems to be a bit more versatile.

library(pophelper)
library(stringr)

## read in all q-matrix files from a particular folder
## perm-dataset

setwd("~/octopus/tess3r/tess_manual_69_40/cluster_11")
str_list <- list.files(path=getwd(), full.names=T)
str_list

## create a tess3R box with all q-matrices stored in it
prem_tess3_q_manual <- readQ(files = str_list,
                        filetype="basic")

## create a plot object for this box
bars <- plotQ(prem_tess3_q_manual,
              returnplot=T,
              imgoutput = "join",
              exportplot=F,
              basesize=11,
              barbordercolour="grey",
              barbordersize = 0.05,
              showindlab=F,
              useindlab=T,
              panelspacer = 0.3,
              splab = c("K = 1", "K = 2", "K = 3", "K = 4", "K = 5", "K = 6", "K = 7", "K = 8",
                        "K = 9", "K = 10", "K = 11"),
              clustercol = c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF",
                             "#ae5d5d","#FDF060","#FFA6B2","#992525",
                              "#4a805e", "#cdc0b0")) ## colors will be reused in plotting maps

## show final plot
bars$plot[[1]]

## same for the stri-dataset

setwd("~/octopus/tess3r/tess_manual_64_20/")
str_list <- list.files(path=getwd(), full.names=T)
str_list

stri_tess3_q_manual <- readQ(files = str_list,
                             filetype="basic")

bars_stri <- plotQ(stri_tess3_q_manual,
              returnplot=T,
              imgoutput = "join",
              exportplot=F,
              basesize=11,
              barbordercolour="black",
              barbordersize = 0.05,
              showindlab=F,
              useindlab=T,
              panelspacer = 0.3,
              splab = c("K = 1", "K = 2", "K = 3", "K = 4", "K = 5", "K = 6", "K = 7", "K = 8"),
              clustercol = c("#1D72F5","#DF0101","#77CE61", "#FF9326", "#ae5d5d","#A945FF",
                             "#FFA6B2","#FDF060")) ## colors will be reused in plotting maps

bars_stri$plot[[1]]

## Interpolation 
## The last part of this analysis is plotting the maps where the probability of an individual to belong
## to a certain cluster is indicated by the shade of color at a certain area.


## I wrote two functions sort_colors_69/64 for this, taking a q-matrix from the helper list as input 
## and outputting an interpolated ggplot map.

sort_color_69 <- function(matrix, point_size = 0.5, ## point size on the map
                          labs = T ## indicate if "Latitude and Longitude" will be plotted
                          ) {
                testmap <- getMap(resolution = "high", projection=NA) ## map as backdrop 
                temp <- as.data.frame(matrix) ## q-matrix to dataframe
                temp_matrix <- matrix
                my.colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                               "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                               "#4a805e", "#cdc0b0")
                my.palette <- CreatePalette(my.colors, 15)
                ## create the ggtess object for plotting
                pl_temp <- ggtess3Q(temp_matrix, perm_coordinates, 
                                            raster.filename = "~/octopus/raster files/brazil_test_10.asc",
                                            col.palette = my.palette,
                                    interpolation.model = FieldsKrigModel(9.5))
                if (labs == T) {
                pl_final <- pl_temp + 
                        geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
                        xlim(-80, -0) +
                        ylim(-22, 12) +
                        coord_equal()  +
                        geom_point(data = as.data.frame(perm_coordinates), 
                                   aes(x = V1, y = V2), size = point_size,
                                   color = "black") +
                        xlab("Longitute") +
                        ylab("Latitude") +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=0.5)
                              )
                } else {
                pl_final <- pl_temp + 
                        geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
                        xlim(-80, -0) +
                        ylim(-22, 12) +
                        coord_equal()  +
                        geom_point(data = as.data.frame(perm_coordinates), 
                                        aes(x = V1, y = V2), size = point_size,
                                        color = "black") +
                        theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                axis.title.x=element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y=element_blank()
                                )
                        
                }
                
                return(pl_final)
}


## for 64 individuals
sort_color_64 <- function(matrix, point_size = 0.5, ## point size on the map
                          labs = T ## indicate if "Latitude and Longitude" will be plotted
                          ) {
        testmap <- getMap(resolution = "high", projection=NA) ## map as backdrop 
        temp <- as.data.frame(matrix) ## q-matrix to dataframe
        temp_matrix <- matrix
       
        my.colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326",
                       "#ae5d5d","#A945FF","#FDF060","#FFA6B2","#992525",
                       "#4a805e", "#cdc0b0")
        my.palette <- CreatePalette(my.colors, 15)
        
        ## create the ggtess object for plotting
        pl_temp <- ggtess3Q(temp_matrix, stri_coordinates, 
                            raster.filename = "~/octopus/raster files/brazil_test_10.asc",
                            interpolation.model = FieldsKrigModel(8.0),
                            col.palette = my.palette)
        if (labs == T){
        pl_final <- pl_temp + 
                geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
                xlim(-80, -0) +
                ylim(-22, 12) +
                coord_equal()  +
                geom_point(data = as.data.frame(stri_coordinates), 
                           aes(x = V1, y = V2), size = 3,
                           color = "black") +
                xlab("Longitute") +
                ylab("Latitude") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=0.5))
        } else {
        pl_final <- pl_temp + 
                geom_path(data = testmap, aes(x = long, y = lat, group = group)) +
                xlim(-80, -0) +
                ylim(-22, 12) +
                coord_equal()  +
                geom_point(data = as.data.frame(perm_coordinates), 
                                aes(x = V1, y = V2), size = point_size,
                                color = "black") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y=element_blank()
                )
                
                
        }
        
        return(pl_final)
}

## test it
sort_color_69(helper[[6]])

## combine maps and q-matrices for 69 inds

map2 <- sort_color_69(helper[[2]], labs = F)
map3 <- sort_color_69(helper[[3]], labs = F)
map4 <- sort_color_69(helper[[4]], labs = F)
map5 <- sort_color_69(helper[[5]], labs = F)
map6 <- sort_color_69(helper[[6]], labs = F)
map7 <- sort_color_69(helper[[7]], labs = F)
map8 <- sort_color_69(helper[[8]], labs = F)
map9 <- sort_color_69(helper[[9]], labs = F)
map10 <- sort_color_69(helper[[10]], labs = F)
map11 <- sort_color_69(helper[[11]], labs = F)

ggdraw() +
        draw_plot(map2, x = 0, y = 0.9, height = 0.1 ,width =0.5) +
        draw_plot(map3, x = 0, y = 0.8, height = 0.1 ,width =0.5) + 
        draw_plot(map4, x = 0, y = 0.7, height = 0.1 ,width =0.5) + 
        draw_plot(map5, x = 0, y = 0.6, height = 0.1 ,width =0.5) + 
        draw_plot(map6, x = 0, y = 0.5, height = 0.1 ,width =0.5) + 
        draw_plot(map7, x = 0, y = 0.4, height = 0.1 ,width =0.5) + 
        draw_plot(map8, x = 0, y = 0.3, height = 0.1 ,width =0.5) + 
        draw_plot(map9, x = 0, y = 0.2, height = 0.1 ,width =0.5) + 
        draw_plot(map10, x = 0, y = 0.1, height = 0.1, width =0.5) + 
        draw_plot(map11, x = 0, y = 0.0, height = 0.1, width =0.5) +
        draw_plot(q_mat(helper[[2]], labs = F), x = 0.4, y = 0.91, height = 0.09, width =0.5) +
        draw_plot(q_mat(helper[[3]], labs = F), x = 0.4, y = 0.81, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[4]], labs = F), x = 0.4, y = 0.71, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[5]], labs = F), x = 0.4, y = 0.61, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[6]], labs = F), x = 0.4, y = 0.51, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[7]], labs = F), x = 0.4, y = 0.41, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[8]], labs = F), x = 0.4, y = 0.31, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[9]], labs = F), x = 0.4, y = 0.21, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[10]], labs = F), x = 0.4, y = 0.11, height = 0.09, width =0.5) + 
        draw_plot(q_mat(helper[[11]], labs = T), x = 0.4, y = 0.0, height = 0.11, width =0.5) +
        draw_plot_label(
                c("K= 2", "K= 3", "K= 4", "K= 5", "K= 6", "K= 7", "K= 8", "K= 9", "K=10", "K=11"),
                c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)-0.01,
                c(0.99, 0.89, 0.79, 0.69, 0.59, 0.49, 0.39, 0.29, 0.19, 0.09)-0.02,
                size = 15
        )


## combine maps and q-matrices for 64 inds

map2_64 <- sort_color_64(helper_64[[2]], labs = F)
map3_64 <- sort_color_64(helper_64[[3]], labs = F)
map4_64 <- sort_color_64(helper_64[[4]], labs = F)
map5_64 <- sort_color_64(helper_64[[5]], labs = F)
map6_64 <- sort_color_64(helper_64[[6]], labs = F)
map7_64 <- sort_color_64(helper_64[[7]], labs = F)
map8_64 <- sort_color_64(helper_64[[8]], labs = F)
map9_64 <- sort_color_64(helper_64[[9]], labs = F)
map10_64 <- sort_color_64(helper_64[[10]], labs = F)
map11_64 <- sort_color_64(helper_64[[11]], labs = F)


ggdraw() +
        draw_plot(map2_64, x = 0, y = 0.9, height = 0.1 ,width =0.5) +
        draw_plot(map3_64, x = 0, y = 0.8, height = 0.1 ,width =0.5) + 
        draw_plot(map4_64, x = 0, y = 0.7, height = 0.1 ,width =0.5) + 
        draw_plot(map5_64, x = 0, y = 0.6, height = 0.1 ,width =0.5) + 
        draw_plot(map6_64, x = 0, y = 0.5, height = 0.1 ,width =0.5) + 
        draw_plot(map7_64, x = 0, y = 0.4, height = 0.1 ,width =0.5) + 
        draw_plot(map8_64, x = 0, y = 0.3, height = 0.1 ,width =0.5) + 
        draw_plot(map9_64, x = 0, y = 0.2, height = 0.1 ,width =0.5) + 
        draw_plot(map10_64, x = 0, y = 0.1, height = 0.1, width =0.5) + 
        draw_plot(map11_64, x = 0, y = 0.0, height = 0.1, width =0.5) +
        draw_plot(q_mat_64(helper_64[[2]], labs = F), x = 0.4, y = 0.91, height = 0.09, width =0.5) +
        draw_plot(q_mat_64(helper_64[[3]], labs = F), x = 0.4, y = 0.81, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[4]], labs = F), x = 0.4, y = 0.71, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[5]], labs = F), x = 0.4, y = 0.61, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[6]], labs = F), x = 0.4, y = 0.51, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[7]], labs = F), x = 0.4, y = 0.41, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[8]], labs = F), x = 0.4, y = 0.31, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[9]], labs = F), x = 0.4, y = 0.21, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[10]], labs = F), x = 0.4, y = 0.11, height = 0.09, width =0.5) + 
        draw_plot(q_mat_64(helper_64[[11]], labs = T), x = 0.4, y = 0.0, height = 0.11, width =0.5) +
        draw_plot_label(
                c("K= 2", "K= 3", "K= 4", "K= 5", "K= 6", "K= 7", "K= 8", "K= 9", "K=10", "K=11"),
                c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)-0.01,
                c(0.99, 0.89, 0.79, 0.69, 0.59, 0.49, 0.39, 0.29, 0.19, 0.09)-0.02,
                size = 15
        )
