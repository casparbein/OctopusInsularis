## load libraries
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
setwd("~/LMU/octopus/IQtree/")
getwd()

## load ML-tree for 69inds_40MD dataset
miss40_69 <- read.tree("69inds_40MD.newick")
plot(miss40_69, type = "phylogram", edge.width = 2)

## root the tree, outgroup: ALBA, TM
## To choose these clusters as outgroup is arbitrary, but a rooted tree is more visually appealing
miss40_69$tip.label[10:53] ## other outgroup
outgroup <- c("AL1" ,  "AL2"  , "BA10" , "BA17" , "BA12",  "AL7",  
              "BA16" , "BA18" , "BA3" ,  "BA14" ,
             "BA2"  , "BA5" ,  "TM1"  , "TM10" , "TMV5" , "TM2" ,  "TM11" , "TMV4"  ,"TM14" , 
             "TMV6" , "TM3" , "TM21"  ,"TM23" , "TM15"  ,"AL5")

rooted_tree <- root(miss40_69, outgroup= outgroup, resolve.root=TRUE)

## plot the rooted tree
plot(rooted_tree, type = "phylogram", edge.width = 2,
     show.node.label = F,
     align.tip.label = T, font =2, cex = 0.8)

## with ggtree

## group clades according to clusters inferred with PCA,Tess3R,fineRADstructure

## To be able to plot the three with bootstrap values as circles and clades colored 
## according to the clusters the animals were assigned to, 
## one has to get the node numbers of respective clades in the tree,
## which is possible by plotting them and taking notes (the old fashioned way with pen and paper),
## using the ggtree commands geom_text2(aes(subset=!isTip, label=node)) and 
## geom_nodelab()

## daylight (unrooted):

miss40_69_nr <- groupClade(miss40_69, c( 82, ##CEN
                                            121, ##ATL
                                            123, ##OIC
                                            116, ##SPS
                                            126 ##TMV
))

ggtree(miss40_69_nr, layout = "daylight",
       size = 1) + #aes(color=group)
  #geom_tiplab() +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-1) +
  geom_tippoint(aes(color=group), size=3, alpha=.75) + 
  scale_color_manual(name = "cluster", labels = c("S-Coastal","N-Coastal", "N-Atlantic","N-Car", "N-SPS", "S-Oceanic"), 
                     values = c("#FF9326","#1D72F5", "#A945FF","#ae5d5d","#77CE61",  
                               "#DF0101")) + 
  geom_point2(fill = "black", color = "white",  size = 3, shape = 21, aes(
  subset = (node %in% c(123, 116, 117, 126, 133, 134, 135, 129, 130, 132,
                        121, 104, 101, 112, 83, 82, 94, 89, 87))), ## nodes with bootstrap support 
  ## higher than 90
  show.legend = FALSE) +
  # geom_cladelabel(node=82, label="northern clusters",  offset.text  = 0.005, color="#DF0101", angle = 270,
  #                 barsize = 1, fontsize = 5, offset = 0.01) + #, align = TRUE +
  # geom_strip("BA12", "TMV6", label="southern clusters", align = TRUE, offset.text  = 0.005, color="#1D72F5", angle = 270,
  #            barsize = 1, hjust = 0.5, fontsize = 5) +
  # #theme(plot.margin = unit(c(10,10,10,10), "mm")) 
  #geom_treescale() +
  ggtitle("Unrooted Maximum Likelihood tree for the 69inds_40MD dataset")
  

## rooted tree
miss40_69_r <- groupClade(rooted_tree, c( 82, ##CEN
                                          121, ##ATL
                                       123, ##OIC
                                       116, ##SPS
                                       126 ##TMV
                                       ))

ggtree(miss40_69_r,#aes(color=group), 
       size = 1) + 
  geom_rootedge(rootedge = 0.025) +
  geom_tippoint(aes(color=group), size=3, alpha=.75) + 
  #geom_text2(aes(subset=!isTip, label=node), hjust=-1) + ## get node names for assignments of clades
  #geom_tiplab() + 
   scale_color_manual(name = "cluster", labels = c("ALBA","CEN", "ATL","OIC", "SPS", "TMV"), 
                      values = c("#FF9326","#1D72F5","#A945FF","#ae5d5d","#77CE61",  
                                   "#DF0101")) + 
  #guides(color = guide_legend(override.aes = list(size = 5))) 
  #geom_hilight(node=128, fill="blue", alpha=.04) +
  #geom_hilight(node = 75, fill = "blue", alpha = .04)+
  #geom_hilight(node = 76, fill = "blue", alpha = .04) +
  #geom_nodelab(hjust = 1, fontface = 2, size = 4)+ ## for bootstrap support
  geom_point2(fill = "black", color = "white",  size = 3, shape = 21, aes(
             subset = (node %in% c(123, 116, 117, 126, 133, 134, 135, 129, 130, 132,
                                  121, 104, 101, 112, 83, 82, 94, 89, 87))), ## nodes with bootstrap support 
             ## higher than 90
             show.legend = FALSE) +
  geom_cladelabel(node=82, label="North", align=TRUE,  offset.text  = 0.005, color="#1D72F5", angle = 270,
                  barsize = 1, fontsize = 5, offset = 0.01) +
  geom_strip("BA3", "BA12", label="South", align = TRUE, offset.text  = 0.005, color="#DF0101", angle = 270,
             barsize = 1, hjust = 0.5, fontsize = 5) +
  # #theme(plot.margin = unit(c(10,10,10,10), "mm")) 
  geom_treescale() +
  ggtitle("Rooted Maximum Likelihood tree for the 69inds_40MD dataset")

## colors for a proper legend

library(RColorBrewer)
myColors <- c("#FF9326","#DF0101","#A945FF","#77CE61", "#ae5d5d" ,
"#1D72F5") 
names(myColors) <- c("CEN", "SPS", "ATL", "OIC", "TMV", "ALBA")
colScale <- scale_fill_manual(name = "locs",values = myColors)

## for the 64inds_20MD dataset

miss20_64 <- read.tree("~/LMU/octopus/IQtree/64inds_20MD.newick")
plot(miss20_64, type = "phylogram", edge.width = 2)

## root the tree, outgroup: ALBA, TM
outgroup <- c("AL1" ,  "AL2"  , "BA10" , "BA17" , "BA12",  "AL7",
              "BA16" , "BA3" ,  "BA14" ,
              "BA2"  , "BA5" , "AL5","TM1"  , "TM10" , "TMV5" , "TM2" ,  "TM11" , "TMV4"  ,"TM14" ,
              "TMV6" , "TM3" , "TM21"  ,"TM23" , "TM15" )

# outgroup <- c("AR1"  , "CE21" , "CE3"  , "RN12" , "RN3",  
#               "AR5"  , "FN15" , "FN7"  , "CE10" , "OIC2" , "CE16" , "AR10" , "RN13" , "FN03" , "AR11" ,
#               "FN02" , "RN4"  , "AR12" , "FN8"  , "CE12" , "CE17" , "ARO28", "CE19" , "CE20" , "AR23" ,
#               "CE11" , "RN11" , "AR4"  , "CE1"  , "FN01" , "RN2"  , "AR24" , "RN10" , "RN14" , "AR25" ,
#               "ASC1" , "ASC2" , "STH2" , "SPS10", "SPS8" , "SPS2" , "RN2A" , "OIC1" , "OIC3" , "OIC4")

rooted_tree <- root(miss20_64, outgroup, resolve.root=TRUE)

## plot the rooted tree
plot(rooted_tree, type = "phylogram", edge.width = 2,
     show.node.label = F,
     align.tip.label = T, font =2, cex = 0.8)

## with ggtree

## unrooted
miss20_64_nr <- groupClade(miss20_64, c( 74, ##CEN
                                           98, ##ATL
                                           112, ##OIC
                                           99, ##SPS
                                           116 ##TMV
))

ggtree(miss20_64_nr, layout = "daylight", #aes(color=group)
       size = 1) +
  
  geom_tippoint(aes(color=group), size=3, alpha=.75) + 
  #geom_rootedge(rootedge = 0.025) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-1) + ## get node names for assignments of clades
  #geom_tiplab() + 
  # scale_color_manual(name = "cluster", labels = c("ALBA","CEN", "ATL","OIC", "SPS", "TMV"),
  #                    values = c("#FF9326","#DF0101","#A945FF","#ae5d5d","#77CE61",
  #                               "#1D72F5")) +
  scale_color_manual(name = "cluster", labels = c("S-Coastal","N-Coastal", "N-Atlantic","N-Car", "N-SPS", "S-Oceanic"), 
                     values = c("#FF9326","#1D72F5", "#A945FF","#ae5d5d","#77CE61",  
                                "#DF0101")) + 
  #guides(color = guide_legend(override.aes = list(size = 5)))
  #geom_hilight(node=128, fill="blue", alpha=.04) +
  #geom_hilight(node = 75, fill = "blue", alpha = .04)+
  #geom_hilight(node = 76, fill = "blue", alpha = .04) +
  #geom_nodelab(hjust = 1, fontface = 2, size = 4)+ ## for bootstrap support
  geom_point2(fill = "black", color = "white",  size = 3, shape = 21, aes(
    subset = (node %in% c(97,99,100,
                          98,112,100,
                          75,116,92,
                          80,123,126,
                          124,125,74))), ## nodes with bootstrap support
    ## higher than 90
    show.legend = FALSE) +
  # geom_cladelabel(node=74, label="northern clusters", align=TRUE,  offset.text  = 0.005, color="#DF0101", angle = 270,
  #                 barsize = 1, fontsize = 5, offset = 0.01) +
  # geom_strip("BA3", "BA14", label="southern clusters", align = TRUE, offset.text  = 0.005, color="#1D72F5", angle = 270,
  #            barsize = 1, hjust = 0.5, fontsize = 5) +
  #theme(plot.margin = unit(c(10,10,10,10), "mm"))
  #geom_treescale() +
  ggtitle("Unrooted Maximum Likelihood tree for the 64inds_20MD dataset")


## rooted

miss20_64_r <- groupClade(rooted_tree, c( 74, ##CEN
                                          98, ##ATL
                                          112, ##OIC
                                          99, ##SPS
                                          116 ##TMV
))

ggtree(miss20_64_r,#aes(color=group), 
       size = 1) + 
  geom_rootedge(rootedge = 0.025) +
  geom_tippoint(aes(color=group), size=3, alpha=.75) + 
  scale_color_manual(name = "cluster", labels = c("S-Coastal","N-Coastal", "N-Atlantic","N-Car", "N-SPS", "S-Oceanic"), 
                     values = c("#FF9326","#1D72F5", "#A945FF","#ae5d5d","#77CE61",  
                                "#DF0101")) +
  
  #geom_text2(aes(subset=!isTip, label=node), hjust=-1) + ## get node names for assignments of clades
  #geom_tiplab() + 
  # scale_color_manual(name = "cluster", labels = c("ALBA","CEN", "ATL","OIC", "SPS", "TMV"),
  #                    values = c("#FF9326","#DF0101","#A945FF","#ae5d5d","#77CE61",
  #                               "#1D72F5")) +
  #guides(color = guide_legend(override.aes = list(size = 5)))
  #geom_hilight(node=128, fill="blue", alpha=.04) +
  #geom_hilight(node = 75, fill = "blue", alpha = .04)+
  #geom_hilight(node = 76, fill = "blue", alpha = .04) +
  #geom_nodelab(hjust = 1, fontface = 2, size = 4)+ ## for bootstrap support
  geom_point2(fill = "black", color = "white",  size = 3, shape = 21, aes(
    subset = (node %in% c(97,99,100,
                          98,112,100,
                          75,116,92,
                          80,123,126,
                          124,125,74))), ## nodes with bootstrap support
    ## higher than 90
    show.legend = FALSE) +
  geom_cladelabel(node=74, label="northern clusters", align=TRUE,  offset.text  = 0.005, color="#1D72F5", angle = 270,
                  barsize = 1, fontsize = 5, offset = 0.01) +
  geom_strip("BA3", "BA14", label="southern clusters", align = TRUE, offset.text  = 0.005, color= "#DF0101", angle = 270,
             barsize = 1, hjust = 0.5, fontsize = 5) +
  #theme(plot.margin = unit(c(10,10,10,10), "mm"))
  geom_treescale() +
  ggtitle("Rooted Maximum Likelihood tree for the 64inds_20MD dataset")


## Tree for the tetrad comparison, design slightly different (this was used in the manuscript)
miss40_69_r <- groupClade(rooted_tree, c( 82, ##CEN
                                          121, ##ATL
                                          123, ##OIC
                                          116, ##SPS
                                          126 ##TMV
))

?ggtree()

ggtree(miss40_69_r,#aes(color=group), 
       size = 0.5) + 
  geom_rootedge(rootedge = 0.025) +
  geom_tippoint(aes(color=group), size=3, alpha=1) + 
  #geom_text2(aes(subset=!isTip, label=node), hjust=-1) + ## get node names for assignments of clades
  #geom_tiplab() + 
  scale_color_manual(name = "cluster", labels = c("ALBA","CEN", "ATL","OIC", "SPS", "TMV"), 
                     values = c("#FF9326","#1D72F5","#A945FF","#ae5d5d","#77CE61",  
                                "#DF0101"), guide = "none") + 
  #guides(color = guide_legend(override.aes = list(size = 5))) 
  #geom_hilight(node=128, fill="blue", alpha=.04) +
  #geom_hilight(node = 75, fill = "blue", alpha = .04)+
  #geom_hilight(node = 76, fill = "blue", alpha = .04) +
  #geom_nodelab(hjust = 1, fontface = 2, size = 4)+ ## for bootstrap support
  geom_point2(fill = "black", color = "white",  size = 3, shape = 21, aes(
    subset = (node %in% c(123, 116, 117, 126, 133, 134, 135, 129, 130, 132,
                          121, 104, 101, 112, 83, 82, 94, 89, 87))), ## nodes with bootstrap support 
    ## higher than 90
    show.legend = FALSE) +
  #geom_cladelabel(node=82, label="North", align=TRUE,  offset.text  = 0.005, color="#1D72F5", angle = 270,
  #                barsize = 1, fontsize = 5, offset = 0.01) +
  #geom_strip("BA3", "BA12", label="South", align = TRUE, offset.text  = 0.005, color="#DF0101", angle = 270,
  #           barsize = 1, hjust = 0.5, fontsize = 5) +
  # #theme(plot.margin = unit(c(10,10,10,10), "mm")) 
  geom_treescale() #+
  #ggtitle("Rooted Maximum Likelihood tree for the 69inds_40MD dataset")
