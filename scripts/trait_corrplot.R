library(tidyverse)
library(corrplot)
library(ape)
library(ggtree)
library(ggtreeExtra)

# corrplot requires a matrix of continuous values
# to add the significance stars, I need a matrix formatted the exact same way but with significance values

################
# wrangle data #
################

# read in phenotype data
xipho_phen <- read_csv("20230824_phy_test_info_revised_PIC_samples.csv")
xipho_phen <- xipho_phen %>%
  mutate(swd_index = swd_length/std_length)

# make subset for males
xipho_phen_m <- xipho_phen %>%
  filter(sex == "M") %>%
  # distinguish between sneakers and courters
  mutate(species = case_when(str_detect(fish_ID, "CRTR") == T ~ paste(species,"_large",sep = ""),
                             str_detect(fish_ID, "SNKR") == T ~ paste(species,"_small",sep = ""),
                             str_detect(fish_ID, "SNKR") == F & str_detect(fish_ID, "CRTR") == F ~ species)) %>%
  group_by(species) %>%
  # group by and average by species
  summarise_at(vars(std_length,swd_length,vertical_bars,upper_sword_edge_width,upper_sword_edge_length,
                    lower_sword_edge_width,lower_sword_edge_length,dorsal_fin_length,dorsal_fin_height,
                    body_depth,peduncle_depth,caudal_fin_length,caudal_fin_height,peduncle_edge,lower_edge_2_peduncle,swd_index),mean)
  

# read in tree data
swdtree <- read.nexus("tobir3concat_pruned_v4_revised_August2022_branchlengths.tre")
swdtree <- root(swdtree, outgroup = "Psjonesii")
species <- swdtree$tip.label

#remove duplicate ids/missing species
species <- subset(species, species != "Xbirref" &
                   species != "Priapella")
pruned_swdtree <- drop.tip(swdtree,swdtree$tip.label[-match(species, swdtree$tip.label)])

# match up phenotype and tree data:
# the tree species and the phenotype species need to match exactly, no missing or extra species


rownames(xipho_phen_m) <- xipho_phen_m$species
# this seems to undo the rowname setting
xipho_phen_m <- xipho_phen_m[match(pruned_swdtree$tip.label,rownames(xipho_phen_m)),]

#####################
#  helper functions #
#####################

# test whether trait needs to be log transformed, and if so, transform it
pic_diagnostics <- function(){
  
}


# create a matrix of PIC values
pic_matrix<- function(df, tree){
  
  # create matrix without the species info (just numeric data)
  df_dbl <- df %>%
    select(is.double)
  
  # initialize n x n matrix where n = number of variables in phenotype data
  mat <- matrix(nrow=length(colnames(df_dbl)),
                ncol = length(colnames(df_dbl)),
                dimnames = list(colnames(df_dbl)))
  colnames(mat) <- colnames(df_dbl)
  
  # iterate through every combination of variables and run PIC
  for (i in colnames(mat)){
    for (j in colnames(mat)){
      
      vec_df_i <- with(df, setNames(df[[i]], df[["species"]]))
      vec_df_j <- with(df, setNames(df[[j]], df[["species"]]))
      pic_i <- pic(vec_df_i, tree)
      pic_j <- pic(vec_df_j, tree)
      
      mat[i,j] <- cor(pic_i,pic_j)
      
    }
  }
  return(mat)
}

#############
# plot data #
#############

# create correlation matrix (no PIC)
xpm_corr <- cor(xipho_phen_m %>%
                  select(is.double))
xpm_corr_res <- cor.mtest(xipho_phen_m %>%
                            select(is.double),conf.level=0.95)
#make corrplot
corrplot(xpm_corr,
         p.mat = xpm_corr_res$p,
         sig.level=c(0.001,0.01,0.05),
         pch.cex=0.9,
         insig = 'label_sig',
         pch.col = 'magenta',
         method = 'color',
         diag = FALSE,
         type = 'upper',
         col=COL2("BrBG"),
         tl.col = "black")

# create PIC matrix
xipho_phen_m_pic <- pic_matrix(xipho_phen_m,pruned_swdtree)

xipho_phen_m_pic_res <- cor.mtest(xipho_phen_m_pic, conf.level=0.95)

# make corrplot for pic
corrplot(xipho_phen_m_pic,
         p.mat = xipho_phen_m_pic_res$p,
         sig.level=c(0.001,0.01,0.05),
         pch.cex=0.9,
         insig = 'label_sig',
         pch.col = 'magenta',
         method = 'color',
         diag = FALSE,
         type = 'upper',
         col=COL2("BrBG"),
         tl.col = "black")