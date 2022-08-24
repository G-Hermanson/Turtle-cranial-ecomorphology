#SCRIPT 1 - FULL LANDMARK DATASET
set.seed(123)

#load packages
library(AICcmodavg)
library(ape)
library(class)
library(geiger)
library(geomorph)
library(landvR)
library(lattice)
library(mda)
library(Morpho)
library(nlme)
library(nnet)
library(phylolm)
library(phytools)
library(rr2)

#Clear workspace
rm( list = ls() )

#Set working directory
#This working directory should contain the following folder/files:
# - "Landmarks_full_dataset" folder
# - "Landmarks_partial_dataset" folder
# - "Pereira_tree_pruned_full.tre" tree file
# - "Pereira_tree_pruned_partial.tre" tree file
# - "Evers_tree_pruned.tre" tree file
# - "Sterli_tree_pruned.tre" tree file
# - Scripts 1-3 to run analyses and "Custom R functions.R" script

#Set working directory
setwd ( "Working directory" )

#load custom functions

source ( 'Custom R functions.txt' )

#ANALYSES

# LOAD 3D COORDINATES AND SLIDERS INFO; SET 3D ARRAY FOR DOWNSTREAM STEPS

#Load 3D coordinates from folder containing the data
dir.temp <- paste0(getwd(),'/Landmarks_full_dataset/')
land.temp <- paste0(dir.temp,list.files(dir.temp,pattern = '.txt'))
land.temp <- lapply ( land.temp , read.table)
land.rownames <- read.csv(paste0(dir.temp,'rownames.csv'), sep=';')[,2]

landmarks.full <- array(unlist(land.temp), c(dim(land.temp[[1]]), length(land.temp))) 
names.temp <- gsub('.txt','',list.files(dir.temp,pattern = '.txt'))
dimnames(landmarks.full)[[1]] <- land.rownames
dimnames(landmarks.full)[[2]] <- c('x','y','z')
dimnames(landmarks.full)[[3]] <- names.temp

landmarks.full

sliders.full <- read.csv(paste0(dir.temp,'sliders.csv'),sep = ' ')
sliders.full <- as.matrix(sliders.full)

#plot specimen
i = 1
plot3d(landmarks.full[,,i], aspect = 'iso', box='n', size=5)
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( landmarks.full[ lines.plot.temp[j,] , ,i ] , 
                                                 lwd = 2 , col = 'black' ) }
title3d(main=dimnames(landmarks.full)[[3]][i],line = 1)


# LOAD TREES

#Molecular-based tree of Pereira et al. (2017) pruned to 'full landmark dataset' taxa sample
extant_tree_full <- read.tree ( 'Pereira_tree_pruned_full.tre' )
#Molecular-based tree of Pereira et al. (2017) pruned to 'partial landmark dataset' taxa sample
extant_tree_partial <- read.tree ( 'Pereira_tree_pruned_partial.tre' )
#Time-scaled composite topology based on Evers et al. (2019)
Evers_tree <- read.tree ( 'Evers_tree_pruned.tre' )
#Time-scaled composite topology based on Sterli et al. (2018)
Sterli_tree <- read.tree ( 'Sterli_tree_pruned.tre' )

# LOAD SPECIMENS INFO

turtle_data <- read.csv ( "SupportingInformationS1.csv" , 
                          header = T , row.names = 1 , sep = ",")

full_taxa <- rownames (turtle_data) [ turtle_data$full_dataset == 1 ] 

### DATA ORDINATION ANALYSES ###

#GPA

GPA.full <- gpagen(landmarks.full , curves = sliders.full , 
                   surfaces = as.matrix(399:nrow(landmarks.full)),  ProcD = F )

#get corrected centroid sizes and fix some names

GPA.full$Csize[ GPA.full$Csize > 3000 ] <- GPA.full$Csize[ GPA.full$Csize > 3000 ] / 1000

## change remaining
GPA.full$Csize["Chelonoidis_sp_SMF67582"] <- GPA.full$Csize["Chelonoidis_sp_SMF67582"] * 2 ##Chelonoidis sp
GPA.full$Csize["Cycloderma_frenatum_NHMUK84241"] <- GPA.full$Csize["Cycloderma_frenatum_NHMUK84241"] / 2 ##Cycloderma
GPA.full$Csize["Heosemys_grandis_unnumbered"] <- GPA.full$Csize["Heosemys_grandis_unnumbered"] * 100 ##Heosemys

size.full <- GPA.full$Csize/10
names( size.full ) <- unlist( lapply( lapply( strsplit( names( size.full ) , "_" ) , function(X){X[c(1,2)]} ) , paste , collapse = "_" ) )
names( size.full )[ names( size.full ) == "Chelonoidis_sp" ]  <- "Chelonoidis_nigra"
names( size.full )[ names( size.full ) == "Gopherus_agassizi" ] <- "Gopherus_agassizii"
names( size.full )[ names( size.full ) == "Deirochelys_reticularis" ] <- "Deirochelys_reticularia"
names( size.full )[ names( size.full ) == "Kinosternon_suburum" ] <- "Kinosternon_subrubrum"

#PCA
PCA.results.full <- plotTangentSpace ( GPA.full$coords , warpgrids = F ) #old version of geomorph PCA function
rownames(PCA.results.full$pc.scores) <- sort(full_taxa)

#### ECOMORPHOLOGICAL ANALYSES ####

#create 'hardness index'

# temporary data frame containing only food items
turtle_data.temp_food <- turtle_data[full_taxa,-c(1:5,19:ncol(turtle_data))]

# define food hardness categories based on Vanhooydonck et al. 2007 (See main text)
food_hardness <- list('soft' = c('Flowers','Terrestrial_leaves','Aquatic_leaves','Fungi','Jellyfish','Worms'),
                      'interm' = c('Stems','Vertebrates','Aquatic_insects','Terrestrial_arthropods'),
                      'hard' = c('Seeds_fruits','Mollusks','Crustaceans'))

w_matrix <- cbind (turtle_data.temp_food[,food_hardness[['soft']]]*0,
                   turtle_data.temp_food[,food_hardness[['interm']]]*0.5,
                   turtle_data.temp_food[,food_hardness[['hard']]]*1)

hardness_index <- round(apply(w_matrix,1,sum) / apply(turtle_data.temp_food,1,sum),2)


#create 'evasiveness index'

# define food evasiveness categories based on Vanhooydonck et al. 2007 (See main text)
food_evasiveness <- list('sedentary' = c('Seeds_fruits','Flowers','Stems','Terrestrial_leaves','Aquatic_leaves','Fungi','Jellyfish','Worms'),
                         'interm' = c('Terrestrial_arthropods','Mollusks'),
                         'evasive' = c('Vertebrates','Aquatic_insects','Crustaceans'))

w_matrix <- cbind (turtle_data.temp_food[,food_evasiveness[['sedentary']]]*0,
                   turtle_data.temp_food[,food_evasiveness[['interm']]]*0.5,
                   turtle_data.temp_food[,food_evasiveness[['evasive']]]*1)
evasiveness_index <- round(apply(w_matrix,1,sum) / apply(turtle_data.temp_food,1,sum),2)

### D-PGLS models ###

# create geomorph data frame

gdf.full <- geomorph.data.frame( "shape" = GPA.full$coords,
                                 "skull_size" = log10(size.full)[extant_tree_full$tip.label] ,
                                 "feeds_on_water" = turtle_data[extant_tree_full$tip.label,"Feed_on_water"] , 
                                 "feeds_on_land" = turtle_data[extant_tree_full$tip.label,'Feed_on_land'] , 
                                 "suction" = turtle_data[extant_tree_full$tip.label,'Suction_feeding'],
                                 "durophagous" = turtle_data[extant_tree_full$tip.label,'Mostly_hard_food..durophagy.'],
                                 "plant" = turtle_data[extant_tree_full$tip.label,'Mostly_vegetable_matter..herbivory.'],
                                 "meat" = turtle_data[extant_tree_full$tip.label,'Mostly_animal_matter..carnivory.'] ,
                                 "marine" = turtle_data[extant_tree_full$tip.label,'Marine'] ,
                                 "swim" = turtle_data[extant_tree_full$tip.label,'Open_swimmer'],
                                 "neck_retraction" = turtle_data[extant_tree_full$tip.label,'Neck_retraction'],
                                 "pleurodira" = turtle_data[extant_tree_full$tip.label,'Side_neck'],
                                 "hardness" = hardness_index[extant_tree_full$tip.label],
                                 "evasiveness" = evasiveness_index[extant_tree_full$tip.label],
                                 "phy" = extant_tree_full )

dimnames(gdf.full$shape)[[3]] <- sort(full_taxa)
gdf.full$shape <- gdf.full$shape[,,extant_tree_full$tip.label]

# create vector of models

right.sides.full <- c("skull_size",
                      "skull_size + feeds_on_water","skull_size + feeds_on_land","skull_size + suction","skull_size+durophagous",
                      "skull_size+plant","skull_size+meat","skull_size+swim","skull_size+neck_retraction","skull_size+pleurodira",
                      "skull_size+hardness","skull_size+evasiveness",
                      "skull_size+suction+durophagous+plant+meat+neck_retraction+hardness+evasiveness",
                      "skull_size+suction+durophagous+neck_retraction+hardness+evasiveness",
                      "skull_size+feeds_on_water+suction+durophagous+neck_retraction+hardness+evasiveness",
                      "skull_size+feeds_on_water+suction+durophagous+swim+neck_retraction+hardness+evasiveness")

models.full <- paste( "shape ~" , right.sides.full )
models.full <- lapply( models.full , as.formula )

# Run D-PGLS 

procD.fit.full <- list()

for( i in 1:length(models.full) ) {
  procD.fit.full[[ i ]] <- procD.pgls( models.full[[ i ]] , data = gdf.full , phy = phy , SS.type = "II" )
  names(procD.fit.full)[[i]] <- as.character(models.full[i])
  
}

#inspect results
lapply( procD.fit.full , summary )

#get corrected p-values using FDR
p.corr <- lapply( lapply ( lapply ( procD.fit.full , get.procD.p) , 
                   p.adjust, method="fdr") , round , 4)


# retrieve R2 values from models and order them

#custom function to get R2 from ProcD models
models.R2.full <- sort(unlist(lapply(procD.fit.full,get.R2)),decreasing = T)

# save and export results      
turtle.aov <- do.call( what = rbind ,
                       arg = lapply( procD.fit.full[order(models.R2.full,decreasing = T)] , 
                                     function( X ) { X$aov.table } ) )

write.csv(turtle.aov , "turtle_aov_full.csv" )

# Define best model (according to R2 scores, and with all predictors significant)

model.temp <- names (models.R2.full)[3] #change number in brackets until conditions above are met
procD.fit.full[[model.temp]]$aov.table  #check

bestmodel.full <- procD.fit.full[[model.temp]] ## best model for full dataset

# Get regression scores
reg.full <- procD.scores ( bestmodel.full , plot = F )


### 2B-PLS analyses between emarginations ###

## Get rows with emargination landmarks
n.emarg.temp <-  grep('Temporal', x = land.rownames)
n.emarg.cheek <-  grep('Cheek', x = land.rownames)

# subset temporal emargination landmarks
turtle_coords.full_temporal <- gdf.full$shape[n.emarg.temp,,]

## subset cheek emargination landmarks
turtle_coords.full_cheek <- gdf.full$shape[n.emarg.cheek,,]

# Run 2B-PLS

PLS.emarg <- phylo.integration(turtle_coords.full_cheek , turtle_coords.full_temporal , phy= gdf.full$phy ,
                               iter=999)

# inspect results
PLS.emarg

### pGLS models ###

# Create 'skull' dataset

skull_dataset <- data.frame ( skull_size = gdf.full$skull_size , 
                              carapace_size = log10(turtle_data[extant_tree_full$tip.label,'Body_size']) ,
                              feeds_on_water = gdf.full$feeds_on_water , 
                              feeds_on_land = gdf.full$feeds_on_land , 
                              suction = gdf.full$suction , durophagous = gdf.full$durophagous,
                              plant = gdf.full$plant , meat = gdf.full$meat , 
                              marine = gdf.full$marine , neck = gdf.full$neck_retraction,
                              pleurodira = gdf.full$pleurodira, swim = gdf.full$swim ,
                              hardness = gdf.full$hardness, evasiveness = gdf.full$evasiveness)

# create vector of models

skull_models <- c('carapace_size' , 'feeds_on_water', 'feeds_on_land',
                  'suction' , 'durophagous' , 'plant', 'meat',
                  'marine' , 'neck', 'swim', 'hardness','evasiveness',
                  'carapace_size + marine' , 'carapace_size + swim' , 'carapace_size + neck',
                  'carapace_size + marine + swim + neck' ,
                  'carapace_size + neck + durophagous', 'carapace_size + neck + marine' , 'carapace_size + neck + swim', 
                  'carapace_size + neck + hardness')

skull_models <- paste('skull_size ~' , skull_models)
skull_models <- lapply(skull_models, as.formula)

# Run pGLS

pgls_models <- list()

for ( i in 1:length(skull_models)){
  
  pgls_models[[i]] <- gls ( skull_models[[i]] , data = skull_dataset , 
                            correlation = corPagel(1,extant_tree_full),
                            method = 'ML' )
  names(pgls_models)[[i]] <- as.character(skull_models[i])
}

# inspect results
lapply ( pgls_models , summary )

#corrected p-values
pgls.p.corr <- lapply( lapply ( lapply ( pgls_models, get.pgls.p), p.adjust,
                        method="fdr") , round , 4)

#save and export pgls results

pgls.results <- list()
for ( i in 1:length(pgls_models)){
  
  pgls.results[[i]] <- round(data.frame(summary(pgls_models[[i]])$tTable,
                                  "corr p.value"=pgls.p.corr[[i]]) ,4)
}

to_export <- do.call(rbind, pgls.results)

write.csv(to_export , file=paste0('PGLS_skullsize_results_',Sys.Date(),'.txt') , sep=' ')  

# save and export AIC comparison table

AICc_table <- cbind(aictab(pgls_models, sort = F) , R2 = unlist(lapply(pgls_models, R2.pred)),
                    lambda = unlist(lapply(pgls_models,function(x) x$modelStruct)))
rownames(AICc_table) <- AICc_table$Modnames
AICc_table$Modnames <- NULL
AICc_table <- round(AICc_table , 3)
#order table so best model (i.e. lowest AICc) is on top
AICc_table <- AICc_table[order(AICc_table$Delta_AICc,decreasing = F),]
write.csv(AICc_table , file=paste0('PGLS_skullsize_AICc_',Sys.Date(),'.txt') , sep=' ')  

# test skull size vs. neck/carapace ratio (as proxy for neck length)

taxa_neck_info <- rownames(turtle_data)[!is.na(turtle_data$Neck_ratio)]
taxa_to_use <-  full_taxa %in% taxa_neck_info

skull_dataset_reduced <- data.frame(skull_dataset[taxa_neck_info,] ,
                                    neck_ratio = na.omit(turtle_data[taxa_neck_info,"Neck_ratio"]) )

skull_dataset_reduced <- skull_dataset_reduced[ complete.cases(skull_dataset_reduced),]

extant_tree_full_neck <- keep.tip(extant_tree_full ,rownames(skull_dataset_reduced))

skull_dataset_reduced <- skull_dataset_reduced[extant_tree_full_neck$tip.label,]

skull_models_reduced <- c('carapace_size' , 'feeds_on_water', 'feeds_on_land',
                  'suction' , 'durophagous' , 'plant', 'meat',
                  'marine' , 'neck', 'swim', 'hardness','evasiveness',
                  'neck_ratio',
                  'carapace_size + marine' , 'carapace_size + swim' , 
                  'carapace_size + neck', 'carapace_size + neck_ratio',
                  'carapace_size + marine + swim + neck' ,
                  'carapace_size + neck + swim',
                  'carapace_size + neck + durophagous', 
                  'carapace_size + neck + marine' , 
                  'carapace_size + neck + hardness',
                  'carapace_size + swim + neck_ratio',
                  'carapace_size + marine + neck_ratio')

skull_models_reduced <- paste('skull_size ~' , skull_models_reduced)
skull_models_reduced <- lapply(skull_models_reduced, as.formula)

# pgls_reduced test #

pgls_models_reduced <- list()

for ( i in 1:length(skull_models_reduced)){
  
  pgls_models_reduced[[i]] <- phylolm(skull_models_reduced[[i]], 
                                      data = skull_dataset_reduced ,
                                      phy = extant_tree_full_neck, 
                                      model="lambda")
  names(pgls_models_reduced)[[i]] <- as.character(skull_models_reduced[i])
}


# inspect results
lapply ( pgls_models_reduced , summary )

#get corrected p-values
pgls_reduced.p.corr <- lapply( lapply ( lapply ( pgls_models_reduced, get.phylolm.p), p.adjust,
                                method="fdr") , round , 4)

#save and export pgls results

pgls_reduced.results <- list()
for ( i in 1:length(pgls_models_reduced)){
  
  pgls_reduced.results[[i]] <- round(data.frame(summary(pgls_models_reduced[[i]])$coefficients,
                                        "corr p.value"=pgls_reduced.p.corr[[i]]) ,4)
}

to_export <- do.call(rbind, pgls_reduced.results)

write.csv(to_export , file=paste0('PGLS_skullsize_reduced_results_',Sys.Date(),'.txt') , sep=' ')  

# save and export AIC comparison table

AICc_table_reduced <- cbind(AIC=unlist(lapply(pgls_models_reduced, AIC) ), 
                            R2 = unlist(lapply(pgls_models_reduced, R2.pred, phy=extant_tree_full_neck)),
                    lambda = unlist(lapply(pgls_models_reduced,function(x) x$optpar)))
AICc_table_reduced <- data.frame( round(AICc_table_reduced , 3))

AICc_table_reduced$deltaAIC <- AICc_table_reduced$AIC - min(AICc_table_reduced$AIC)

#export
write.csv(AICc_table_reduced , file=paste0('PGLS_skullsize_reduced_AICc_',Sys.Date(),'.txt') , sep=' ')  


#D-PGLS with "neck ratio" info

coords.neck <- GPA.full$coords
dimnames(coords.neck)[[3]] <- sort(full_taxa)
coords.neck <- coords.neck[,,extant_tree_full_neck$tip.label]

gdf.neck <- geomorph.data.frame( "shape" = coords.neck,
                                 "skull_size" = log10(size.full)[extant_tree_full_neck$tip.label] ,
                                 "feeds_on_water" = turtle_data[extant_tree_full_neck$tip.label,"Feed_on_water"] , 
                                 "feeds_on_land" = turtle_data[extant_tree_full_neck$tip.label,'Feed_on_land'] , 
                                 "suction" = turtle_data[extant_tree_full_neck$tip.label,'Suction_feeding'],
                                 "durophagous" = turtle_data[extant_tree_full_neck$tip.label,'Mostly_hard_food..durophagy.'],
                                 "plant" = turtle_data[extant_tree_full_neck$tip.label,'Mostly_vegetable_matter..herbivory.'],
                                 "meat" = turtle_data[extant_tree_full_neck$tip.label,'Mostly_animal_matter..carnivory.'] ,
                                 "marine" = turtle_data[extant_tree_full_neck$tip.label,'Marine'] ,
                                 "swim" = turtle_data[extant_tree_full_neck$tip.label,'Open_swimmer'],
                                 "neck_retraction" = turtle_data[extant_tree_full_neck$tip.label,'Neck_retraction'],
                                 "pleurodira" = turtle_data[extant_tree_full_neck$tip.label,'Side_neck'],
                                 "hardness" = hardness_index[extant_tree_full_neck$tip.label],
                                 "evasiveness" = evasiveness_index[extant_tree_full_neck$tip.label],
                                 "neck_ratio" = turtle_data[extant_tree_full_neck$tip.label,"Neck_ratio"],
                                 "phy" = extant_tree_full_neck )

# create vector of models

right.sides.neck <- c("skull_size",
                      "skull_size + feeds_on_water","skull_size + feeds_on_land","skull_size + suction","skull_size+durophagous",
                      "skull_size+plant","skull_size+meat","skull_size+swim","skull_size+neck_retraction","skull_size+pleurodira",
                      "skull_size+hardness","skull_size+evasiveness","skull_size+neck_ratio",
                      "skull_size+suction+durophagous+plant+meat+neck_retraction+hardness+evasiveness+neck_ratio",
                      "skull_size+suction+durophagous+neck_retraction+hardness+evasiveness+neck_ratio",
                      "skull_size+feeds_on_water+suction+durophagous+neck_retraction+hardness+evasiveness+neck_ratio",
                      "skull_size+feeds_on_water+suction+durophagous+swim+neck_retraction+hardness+evasiveness+neck_ratio")

models.neck <- paste( "shape ~" , right.sides.neck )
models.neck <- lapply( models.neck , as.formula )

# Run D-PGLS 

procD.fit.neck <- list()

for( i in 1:length(models.neck) ) {
  procD.fit.neck[[ i ]] <- procD.pgls( models.neck[[ i ]] , data = gdf.neck , phy = phy , SS.type = "II" )
  names(procD.fit.neck)[[i]] <- as.character(models.neck[i])
  
}

#inspect results
lapply( procD.fit.neck , summary )

#get corrected p-values using FDR
p.corr.neck <- lapply( lapply ( lapply ( procD.fit.neck , get.procD.p) , 
                           p.adjust, method="fdr") , round , 4)


# retrieve R2 values from models and order them

#custom function to get R2 from ProcD models
models.R2.neck <- sort(unlist(lapply(procD.fit.neck,get.R2)),decreasing = T)

model.temp.neck <- names (models.R2.neck)[2] #change number in brackets until conditions above are met
procD.fit.neck[[model.temp.neck]]$aov.table  #check

bestmodel.neck <- procD.fit.neck[[model.temp.neck]] ## best model for full dataset using "neck ratio" info

# save and export results      
turtle.aov.neck <- do.call( what = rbind ,
                       arg = lapply( procD.fit.neck[order(models.R2.neck,decreasing = T)] , 
                                     function( X ) { X$aov.table } ) )

write.csv(turtle.aov.neck , "turtle_aov_full_neck.csv" )

