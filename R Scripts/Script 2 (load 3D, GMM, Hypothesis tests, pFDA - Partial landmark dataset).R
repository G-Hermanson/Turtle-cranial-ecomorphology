#SCRIPT 2 - PARTIAL LANDMARK DATASET

# LOAD 3D COORDINATES AND SLIDERS INFO; SET 3D ARRAY FOR DOWNSTREAM STEPS

dir.temp <- paste0(getwd(),'/Landmarks_partial_dataset/')
land.temp <- paste0(dir.temp,list.files(dir.temp,pattern = '.txt'))
land.temp <- lapply ( land.temp , read.table, col.names=c('x','y','z'))
land.rownames <- read.csv(paste0(dir.temp,'rownames.csv'), sep=';')[,2]

landmarks.partial <- array(unlist(land.temp), c(dim(land.temp[[1]]), length(land.temp))) 
names.temp <- gsub('.txt','',list.files(dir.temp,pattern = '.txt'))
dimnames(landmarks.partial)[[1]] <- land.rownames
dimnames(landmarks.partial)[[2]] <- c('x','y','z')
dimnames(landmarks.partial)[[3]] <- names.temp

sliders.partial <- read.csv(paste0(dir.temp,'sliders.csv'),sep = ' ')
sliders.partial <- as.matrix(sliders.partial)

#plot specimen
i = 1
plot3d(landmarks.partial[,,i], aspect = 'iso', box='n', size=5)
lines.plot.temp <- rbind( sliders.partial[ , 1:2 ] , sliders.partial[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( landmarks.partial[ lines.plot.temp[j,] , ,i ] , 
                                                 lwd = 2 , col = 'black' ) }
title3d(main=dimnames(landmarks.partial)[[3]][i],line = 1)

partial_taxa <- rownames (turtle_data) [ turtle_data$partial_dataset == 1 ]

### DATA ORDINATION ANALYSES ###

#GPA

GPA.partial <- gpagen(landmarks.partial , curves = sliders.partial , 
                      surfaces = as.matrix(140:nrow(landmarks.partial)),  ProcD = F )


#get corrected centroid sizes and fix some names

GPA.partial$Csize[ GPA.partial$Csize > 3000 ] <- GPA.partial$Csize[ GPA.partial$Csize > 3000 ] / 1000

## change remaining
GPA.partial$Csize["Chelonoidis_sp_SMF67582"] <- GPA.partial$Csize["Chelonoidis_sp_SMF67582"] * 2 ##Chelonoidis sp
GPA.partial$Csize["Cycloderma_frenatum_NHMUK84241"] <- GPA.partial$Csize["Cycloderma_frenatum_NHMUK84241"] / 2 ##Cycloderma
GPA.partial$Csize["Heosemys_grandis_unnumbered"] <- GPA.partial$Csize["Heosemys_grandis_unnumbered"] / 10 ##Heosemys
GPA.partial$Csize["Argillochelys_antiqua_NHMUKR38955"] <- GPA.partial$Csize["Argillochelys_antiqua_NHMUKR38955"] * 10 ##Argillochelys
GPA.partial$Csize["Puppigerus_camperi_IRSNBR0076"] <- GPA.partial$Csize["Puppigerus_camperi_IRSNBR0076"] * 2 ##Puppigerus
GPA.partial$Csize["Rhinochelys_pulchriceps_CAMSMB55783"] <- GPA.partial$Csize["Rhinochelys_pulchriceps_CAMSMB55783"] * 10 ##Rhinochelys

size.partial <- GPA.partial$Csize/10
names( size.partial ) <- unlist( lapply( lapply( strsplit( names( size.partial ) , "_" ) , function(X){X[c(1,2)]} ) , paste , collapse = "_" ) )
names( size.partial )[ names ( size.partial ) == "Chelonoidis_sp" ]  <- "Chelonoidis_nigra"
names( size.partial )[ names ( size.partial ) == "Gopherus_agassizi" ] <- "Gopherus_agassizii"
names( size.partial )[ names ( size.partial ) == "Deirochelys_reticularis" ] <- "Deirochelys_reticularia"
names( size.partial )[ names ( size.partial ) == "Kinosternon_suburum" ] <- "Kinosternon_subrubrum"
names( size.partial )[ names ( size.partial ) == "Geoclemys_hamiltoni" ] <- "Geoclemys_hamiltonii"
names( size.partial )[ names ( size.partial ) == "Hieremys_annandalei" ] <- "Hieremys_annandalii"
names( size.partial )[ names ( size.partial ) == "Annemys_IVPPV18106" ] <- "Annemys_sp"

#PCA

PCA.results.partial <- plotTangentSpace ( GPA.partial$coords , warpgrids = F ) #old version of geomorph PCA function
rownames(PCA.results.partial$pc.scores) <- sort(partial_taxa)

#### ECOMORPHOLOGICAL ANALYSES ####

#create 'hardness index'

# temporary data frame containing only food items
turtle_data.temp_food <- turtle_data[partial_taxa,-c(1:5,19:ncol(turtle_data))]

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
GPA.partial$coords <- GPA.partial$coords[,,sort(dimnames(GPA.partial$coords)[[3]])]
dimnames(GPA.partial$coords)[[3]] <- sort(partial_taxa)

gdf.partial <- geomorph.data.frame( "shape" = GPA.partial$coords[,,extant_tree_partial$tip.label],
                                    "skull_size" = log10(size.partial)[extant_tree_partial$tip.label] ,
                                    "feeds_on_water" = turtle_data[extant_tree_partial$tip.label,"Feed_on_water"] , 
                                    "feeds_on_land" = turtle_data[extant_tree_partial$tip.label,'Feed_on_land'] , 
                                    "suction" = turtle_data[extant_tree_partial$tip.label,'Suction_feeding'],
                                    "durophagous" = turtle_data[extant_tree_partial$tip.label,'Mostly_hard_food..durophagy.'],
                                    "plant" = turtle_data[extant_tree_partial$tip.label,'Mostly_vegetable_matter..herbivory.'],
                                    "meat" = turtle_data[extant_tree_partial$tip.label,'Mostly_animal_matter..carnivory.'] ,
                                    "marine" = turtle_data[extant_tree_partial$tip.label,'Marine'] ,
                                    "swim" = turtle_data[extant_tree_partial$tip.label,'Open_swimmer'],
                                    "neck_retraction" = turtle_data[extant_tree_partial$tip.label,'Neck_retraction'],
                                    "pleurodira" = turtle_data[extant_tree_partial$tip.label,'Side_neck'],
                                    "hardness" = hardness_index[extant_tree_partial$tip.label],
                                    "evasiveness" = evasiveness_index[extant_tree_partial$tip.label],
                                    "phy" = extant_tree_partial )

# create vector of models

right.sides.partial <- c("skull_size","feeds_on_water","feeds_on_land","suction","durophagous","plant","meat","swim",
                         "neck_retraction","pleurodira","hardness","evasiveness",
                         "skull_size + feeds_on_water","skull_size + feeds_on_land","skull_size + suction",
                         "skull_size+durophagous","skull_size+plant","skull_size+meat","skull_size+swim",
                         "skull_size+neck_retraction","skull_size+pleurodira","skull_size+hardness","skull_size+evasiveness",
                         "skull_size+suction+durophagous+plant+meat+neck_retraction+hardness+evasiveness",
                         "skull_size+suction+durophagous+neck_retraction+evasiveness",
                         "skull_size+suction+neck_retraction+hardness+evasiveness",
                         "skull_size+suction+durophagous+neck_retraction+evasiveness+feeds_on_water",
                         "skull_size+suction+durophagous+neck_retraction+evasiveness+swim")

models.partial <- paste( "shape ~" , right.sides.partial )
models.partial <- lapply( models.partial , as.formula )

# Run D-PGLS 

procD.fit.partial <- list()

for( i in 1:length(models.partial) ) {
  procD.fit.partial[[ i ]] <- procD.pgls( models.partial[[ i ]] , data = gdf.partial , phy = phy , SS.type = "II" )
  names(procD.fit.partial)[[i]] <- as.character(models.partial[i])
  
}

#inspect results
lapply( procD.fit.partial , summary )

# retrieve R2 values from models and order them

models.R2.partial <- sort(unlist(lapply(procD.fit.partial,get.R2)),decreasing = T)

# save and export results      
turtle.aov <- do.call( what = rbind ,
                       arg = lapply( procD.fit.partial[order(models.R2.partial,decreasing = T)] , 
                                     function( X ) { X$aov.table } ) )

write.csv(turtle.aov , "turtle_aov_full.csv" )

# Define best model (according to R2 scores, and with all predictors significant)

model.temp <- names ( models.R2.partial)[4] #change number in brackets until conditions above are met
procD.fit.partial[[model.temp]]$aov.table  

bestmodel.partial <- procD.fit.partial[[model.temp]] ## best model for partial dataset

# Get regression scores
reg.partial <- procD.scores ( bestmodel.partial , plot = F )

#Define fossils to get their regressions scores

fossils <- rownames ( turtle_data ) [  turtle_data$Type == 'Fossil' ]
fossils

reg.fossil <- predict.procD.scores ( bestmodel.partial , 
                                     newdata= GPA.partial$coords[,,fossils] )

# Bind extant and fossil scores into a single data frame
# Add predictors from the best D-PGLS model to this same data frame

scores.all <- data.frame ( rbind ( reg.partial , reg.fossil ) ,  
                           suct = c(bestmodel.partial$data$suction , rep(NA,length(fossils))),
                           durop = c(bestmodel.partial$data$durophagous , rep(NA,length(fossils))),
                           neck = c(bestmodel.partial$data$neck_retraction , rep(NA,length(fossils))),
                           ev = c(bestmodel.partial$data$evasiveness , rep(NA,length(fossils))),
                           type= c( rep ('extant',Ntip(extant_tree_partial)) , rep('fossil',length(fossils))))

scores.all$sizes <- log10(size.partial)[rownames(scores.all)]

#Predict evasiveness index for fossils

#stepwise PGLS with 'phylostep' function of 'phylolm'
pgls.temp <-  phylostep(ev~skull_size+suction+durophagous+neck_retraction+evasiveness,
                        data=scores.all[extant_tree_partial$tip.label,],phy = extant_tree_partial,model='lambda') 

#pGLS with 'gls' function of 'nlme', using predictors of best model according to phylostep results

pgls_ev <- gls(ev~skull_size+durophagous+neck_retraction,
               data=scores.all[extant_tree_partial$tip.label,],
               correlation = corPagel(1,extant_tree_partial), method='ML')

#predict index for extinct taxa
ev_preds <- predict(pgls_ev,newdata = data.frame(skull_size=scores.all[fossils,'skull_size'],
                                                 durophagous=scores.all[fossils,'durophagous'],
                                                 neck_retraction=scores.all[fossils,'neck_retraction']))

names(ev_preds) <- fossils
ev_preds <- setNames(as.numeric(ev_preds),fossils)

#Optional step

#Outputs from analyses that require a phylogenetic framework (e.g. pPCA, pGLS) still "need to be
#analyzed using phylogenetic methods" (Revell 2009, p. 3259) in downstream procedures. This following (optional) step verifies that
#residual calculation using 'procD.pgls' is the same as in 'phyl.resid' function from Revell (2009), and justifies our use
#of phylogenetic discriminant analyses on regression scores provided by 'procD'.pgls tests.

resid.temp <- phyl.resid (tree = extant_tree_partial, x = bestmodel.partial$X[,-1], Y = bestmodel.partial$Y, method = 'BM')
resid.pls <- two.b.pls (bestmodel.partial$pgls.residuals , resid.temp$resid)
#resid.pls #r-PLS=1
#plot(resid.pls,bty='n')
#text(resid.pls$XScores[,1],resid.pls$YScores[,1],
#     rownames(resid.pls$XScores),cex=.4,
#     pos=ifelse(resid.pls$XScores[,1]<(-0.05),4,2))


## Phylogenetic Flexible Discriminant Analyses ##
# Use pFDA to predict presence/absence of binary ecological and functional traits from the best D-PGLS model

set.seed(123)
start <- Sys.time() # this is optional; just to check how long this analyis took
reps = 100 #number of iterations

#From this point on, pFDA will run using composite toplogy based on Evers et al. (2019); 
#To run it using the topology from Sterli et al. (2018), add the '#' prior to line 233, and remove the '#' prior to line 234;

full.tree <- Evers_tree
#full.tree <- Sterli_tree


#create a vector of numbers representing the row indices of extant turtles that DON'T have the trait (i.e. 'trait' = 0)
row.no <- row(scores.all)[scores.all[['suct']]==0 & scores.all$type=='extant',1]
#create a vector of numbers representing the row indices of extant turtles that HAVE the trait (i.e. 'trait' = 1)
row.yes <- row(scores.all)[scores.all[['suct']]==1 & scores.all$type=='extant',1]
#create a vector of numbers representing the row indices of extant turtles
row.extant <- row(scores.all)[scores.all$type=='extant',1]
#create a vector of numbers representing the row indices of fossil turtles
row.fossils <- row(scores.all)[scores.all$type=='fossil',1]


#Lists to store matrices of 100 columns (number of iterations). The number of rows correspond to how many
#extant turtles HAVE a given trait. This is important to set an equal prior probability of 
#having/not having it. These lists will all have the same length as 'row.extant', because we're also
#predicting the probability of extant turtles to have it as a training step, to evaluate how accurate the pFDA
#can classify them
ind.no <- list() 
ind.yes <- list()

#List to store the 100 pFDA runs for the ith extant turtle
fda.temp <- list()
#List to store the results of every extant turtle. So the list is supposed to have length=i (number of extant species), and every ith element has 100 pFDA runs (i.e. number of 'reps')
fda.suct <- list()

#List to store the predictions for fossils (and the ith extant turtle) based on 
#the pFDA results of fda.temp[[i]]
fda.pred.temp <- list()
#List of length=i. For every ith extant turtle, you have their predicted probability of having a trait, along
#with the predicted probability of the fossils to have it as well.
fda.pred.suct <- list()

# Run pFDA

for ( i in 1:length(row.extant)){    
  for ( j in 1:reps){  
    
    ind.yes[[i]]<- replicate(reps, 
                             c(i , sample(row.yes[which(row.yes!=i)],size=length(row.yes)-1,replace=F  )))
    ind.no[[i]]<- replicate(reps, 
                            c(i , sample(row.no[which(row.no!=i)],size=length(row.yes)-1,replace=F  ))) 
    
    #From the 'scores.all' object, pick data from X turtles that don't use suction,
    # and X that do, and combine them into this new temporary object, 'score.temp'.
    score.temp <- rbind(scores.all[ sort(ind.no[[i]][-1,j]) ,] , scores.all[ sort(ind.yes[[i]][-1,j]) ,])
    #Then add to this new 'score.temp' object the data from the extinct ones and the ith extant taxon,
    #that will be predicted too as if it were a fossil.
    score.temp <- rbind ( score.temp[,1:5] , as.matrix(scores.all[c(row.extant[i],row.fossils),1:5]))
    #Now create the discrete binary variable that will be predicted, adding 'unknown' for the fossils and
    #the ith extant taxon
    g.temp <- c(scores.all[ sort(ind.no[[i]][-1,j]) ,'suct'] , scores.all[ sort(ind.yes[[i]][-1,j]) ,'suct'] ,
                rep('unkown',length(row.fossils)+1))
    names(g.temp) <- rownames(score.temp)
    
    #Set a temporary tree containing only the taxa from the 'score.temp' object.
    full.tree.temp <- keep.tip(full.tree , rownames(score.temp))
    #This is just to leave the row names of 'score.temp' in the same order as the tip.labels
    score.temp <- score.temp[full.tree.temp$tip.label,]
    #And this to leave the discrete variable object 'g.temp' in the same order as the tip.labels
    g.temp <- g.temp[full.tree.temp$tip.label]
    #Define row number of those that will be predicted
    n.temp <- row(score.temp)[g.temp=='unkown',1]
    #Define all taxa names from your temporary tree
    taxaA.temp <- full.tree.temp$tip.label
    #Get the lambda value necessary for the pFDA
    lambda.temp <- optLambda.pred(score.temp[,1:5] , g.temp , taxaA.temp, full.tree.temp, n.temp)
    lambda.temp <- lambda.temp$optlambda[1,1]
    
    #Run jth iteration of pFDA for the ith taxon+fossils
    fda.temp[[j]] <- phylo.fda.pred(score.temp[,1:5],g.temp,taxtaxA = taxaA.temp,full.tree.temp,
                                    testlistn = n.temp,val=lambda.temp)
    #Store results
    fda.suct[[i]] <- fda.temp
    
    #Predict posterior probabilities for the ith taxon+fossils
    fda.pred.temp[[j]] <- predict(fda.temp[[j]] , newdata = fda.temp[[j]]$DATAtest,
                                  type = 'posterior')
    rownames(fda.pred.temp[[j]]) <- full.tree.temp$tip.label[n.temp]
    
    #Store results
    fda.pred.suct[[i]] <- fda.pred.temp
    
    #This was created just to keep track of the ongoing progress of the analysis (Hello anxiety :D)
    setTxtProgressBar(txtProgressBar(0,reps,style = 3),j)
    
    cat('\nRunning replicates for taxon',i,'/',length(row.extant))  
    
  }
}

# Get posterior probabilities for extant turtles
post.probs.suct <- matrix ( , ncol=reps, nrow=length(row.extant),
                            dimnames=list(rownames(scores.all)[1:76]))

for ( i in 1:length(fda.pred.suct)){
  
  name.temp <- rownames(scores.all)[i]
  
  post.probs.suct[i,] <-  unlist(lapply ( fda.pred.suct[[i]] , function(x) x[name.temp,2]))
}

#View(post.probs.suct)

mean.post.probs.suct <- apply ( post.probs.suct , 1 , mean) #calculate mean value

#set class (0/1) based on 0.66 cutoff of Chapelle et al. 2020 (See main text)
class.temp.suct <- ifelse ( mean.post.probs.suct>0.66,1,0)

#Check how accurate the pFDA predicted the classes of extant turtles
mean ( class.temp.suct == scores.all[extant_tree_partial$tip.label,'suct']) #check % of correct predictions

## Get posterior probabilities for fossil turtles
post.probs.suct.fos <- list()

for ( i in 1:length(fda.pred.suct)){
  post.probs.suct.fos[[i]] <-  do.call(rbind,lapply ( fda.pred.suct[[i]] , function(x) x[fossils,2]))
}

post.probs.suct.fos <- do.call(rbind,post.probs.suct.fos)

#View(post.probs.suct.fos)

mean.post.probs.suct.fos <- apply ( post.probs.suct.fos , 2 , mean ) #calculate mean
mean.post.probs.suct.fos


## If desired, create PDFs to help visualising some results ##

## POSTERIOR PROBABILITY FOR EXTANT TURTLES ##
pdf ( 'posterior_probs_suction_extant.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:nrow(post.probs.suct)){
  hist ( sort(as.numeric(post.probs.suct[i,])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', rownames(post.probs.suct)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,reps),
         font.main = 4)
  
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.suct[i] , col='red', lty=2, lwd=1)
}
dev.off()


## POSTERIOR PROBABILITY FOR FOSSIL TURTLES ##
pdf ( 'posterior_probs_suction_fossils.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:ncol(post.probs.suct.fos)){
  hist ( sort(as.numeric(post.probs.suct.fos[,i])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', colnames(post.probs.suct.fos)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,nrow(post.probs.suct.fos)),
         font.main = 4)
  
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.suct.fos[i] , col='red', lty=2, lwd=1)
}
dev.off()  


# pFDA to predict 'durophagy' and 'neck retraction' follow the same rationale of previous code. So,
# there are fewer comments in the following lines

### DUROPHAGOUS ###

set.seed(123)

ind.no <- list()
ind.yes <- list()
fda.durop <- list()
fda.temp <- list()
fda.pred.durop <- list()
fda.pred.temp <- list()

row.no <- row(scores.all)[scores.all[['durop']]==0 & scores.all$type=='extant',1]
row.yes <- row(scores.all)[scores.all[['durop']]==1 & scores.all$type=='extant',1]
row.extant <- row(scores.all)[scores.all$type=='extant',1]
row.fossils <- row(scores.all)[scores.all$type=='fossil',1]

for ( i in 1:length(row.extant)){    
  for ( j in 1:reps){  
    
    ind.yes[[i]]<- replicate(reps, 
                             c(i , sample(row.yes[which(row.yes!=i)],size=length(row.yes)-1,replace=F  )))
    ind.no[[i]]<- replicate(reps, 
                            c(i , sample(row.no[which(row.no!=i)],size=length(row.yes)-1,replace=F  ))) #}
    
    score.temp <- rbind(scores.all[ sort(ind.no[[i]][-1,j]) ,] , scores.all[ sort(ind.yes[[i]][-1,j]) ,])
    score.temp <- rbind ( score.temp[,1:5] , as.matrix(scores.all[c(row.extant[i],row.fossils),1:5]))
    g.temp <- c(scores.all[ sort(ind.no[[i]][-1,j]) ,'durop'] , scores.all[ sort(ind.yes[[i]][-1,j]) ,'durop'] ,
                rep('unk',length(row.fossils)+1))
    names(g.temp) <- rownames(score.temp)
    #score.temp$g <- g.temp
    
    full.tree.temp <- keep.tip(full.tree , rownames(score.temp))
    score.temp <- score.temp[full.tree.temp$tip.label,]
    g.temp <- g.temp[full.tree.temp$tip.label]
    n.temp <- row(score.temp)[g.temp=='unk',1]
    taxaA.temp <- full.tree.temp$tip.label
    lambda.temp <- optLambda.pred(score.temp[,1:5] , g.temp , taxaA.temp, full.tree.temp, n.temp)
    lambda.temp <- lambda.temp$optlambda[1,1]
    
    fda.temp[[j]] <- phylo.fda.pred(score.temp[,1:5],g.temp,taxtaxA = taxaA.temp,full.tree.temp,
                                    testlistn = n.temp,val=lambda.temp)
    fda.durop[[i]] <- fda.temp
    
    fda.pred.temp[[j]] <- predict(fda.temp[[j]] , newdata = fda.temp[[j]]$DATAtest,
                                  type = 'posterior')
    rownames(fda.pred.temp[[j]]) <- full.tree.temp$tip.label[n.temp]
    
    fda.pred.durop[[i]] <- fda.pred.temp
    
    setTxtProgressBar(txtProgressBar(0,reps,style = 3),j)
    
    cat('\nRunning replicates for taxon',i,'/',length(row.extant))  
    
  }
}

# get posterior probs for extant turtles
post.probs.durop <- matrix ( , ncol=reps, nrow=length(row.extant),
                             dimnames=list(rownames(scores.all)[1:76]))

for ( i in 1:length(fda.pred.durop)){
  
  name.temp <- rownames(scores.all)[i]
  
  post.probs.durop[i,] <-  unlist(lapply ( fda.pred.durop[[i]] , function(x) x[name.temp,2]))
}

#View(post.probs.durop)

mean.post.probs.durop <- apply ( post.probs.durop , 1 , mean) #calculate mean value

#set class (0/1) based on 0.66 cutoff
class.temp.durop <- ifelse ( mean.post.probs.durop>0.66,1,0) 

mean ( class.temp.durop == scores.all[extant_tree_partial$tip.label,'durop']) #check % of correct predictions

## Get probs for fossil turtles
post.probs.durop.fos <- list()

for ( i in 1:length(fda.pred.durop)){
  post.probs.durop.fos[[i]] <-  do.call(rbind,lapply ( fda.pred.durop[[i]] , function(x) x[fossils,2]))
}

post.probs.durop.fos <- do.call(rbind,post.probs.durop.fos)

#View(post.probs.durop.fos)

mean.post.probs.durop.fos <- apply ( post.probs.durop.fos , 2 , mean ) #calculate mean
mean.post.probs.durop.fos

## PDF ##
## POSTERIOR PROB FOR EXTANT TURTLES ##
pdf ( 'posterior_probs_durophagous_extant.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:nrow(post.probs.durop)){
  hist ( sort(as.numeric(post.probs.durop[i,])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', rownames(post.probs.durop)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,reps),
         font.main = 4)
  
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.durop[i] , col='red', lty=2, lwd=1)
}
dev.off()


## POSTERIOR PROB FOR FOSSIL TURTLES ##
pdf ( 'posterior_probs_durophagous_fossils.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:ncol(post.probs.durop.fos)){
  hist ( sort(as.numeric(post.probs.durop.fos[,i])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', colnames(post.probs.durop.fos)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,nrow(post.probs.durop.fos)),
         font.main = 4)
  
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.durop.fos[i] , col='red', lty=2, lwd=1)
}
dev.off()  



### NECK RETRACTION ###

set.seed(123)

ind.no <- list()
ind.yes <- list()
fda.neck <- list()
fda.temp <- list()
fda.pred.neck <- list()
fda.pred.temp <- list()

row.no <- row(scores.all)[scores.all[['neck']]==0 & scores.all$type=='extant',1]
row.yes <- row(scores.all)[scores.all[['neck']]==1 & scores.all$type=='extant',1]
row.extant <- row(scores.all)[scores.all$type=='extant',1]
row.fossils <- row(scores.all)[scores.all$type=='fossil',1]

for ( i in 1:length(row.extant)){    
  for ( j in 1:reps){  
    
    ind.yes[[i]]<- replicate(reps, 
                             c(i , sample(row.yes[which(row.yes!=i)],size=length(row.no)-1,replace=F  )))
    ind.no[[i]]<- replicate(reps, 
                            c(i , sample(row.no[which(row.no!=i)],size=length(row.no)-1,replace=F  ))) #}
    
    score.temp <- rbind(scores.all[ sort(ind.no[[i]][-1,j]) ,] , scores.all[ sort(ind.yes[[i]][-1,j]) ,])
    score.temp <- rbind ( score.temp[,1:5] , as.matrix(scores.all[c(row.extant[i],row.fossils),1:5]))
    g.temp <- c(scores.all[ sort(ind.no[[i]][-1,j]) ,'neck'] , scores.all[ sort(ind.yes[[i]][-1,j]) ,'neck'] ,
                rep('unk',length(row.fossils)+1))
    names(g.temp) <- rownames(score.temp)
    #score.temp$g <- g.temp
    
    full.tree.temp <- keep.tip(full.tree , rownames(score.temp))
    score.temp <- score.temp[full.tree.temp$tip.label,]
    g.temp <- g.temp[full.tree.temp$tip.label]
    n.temp <- row(score.temp)[g.temp=='unk',1]
    taxaA.temp <- full.tree.temp$tip.label
    lambda.temp <- optLambda.pred(score.temp[,1:5] , g.temp , taxaA.temp, full.tree.temp, n.temp)
    lambda.temp <- lambda.temp$optlambda[1,1]
    
    fda.temp[[j]] <- phylo.fda.pred(score.temp[,1:5],g.temp,taxtaxA = taxaA.temp,full.tree.temp,
                                    testlistn = n.temp,val=lambda.temp)
    fda.neck[[i]] <- fda.temp
    
    fda.pred.temp[[j]] <- predict(fda.temp[[j]] , newdata = fda.temp[[j]]$DATAtest,
                                  type = 'posterior')
    rownames(fda.pred.temp[[j]]) <- full.tree.temp$tip.label[n.temp]
    
    fda.pred.neck[[i]] <- fda.pred.temp
    
    setTxtProgressBar(txtProgressBar(0,reps,style = 3),j)
    
    cat('\nRunning replicates for taxon',i,'/',length(row.extant))  
    
  }
}

# get posterior probs for extant turtles
post.probs.neck <- matrix ( , ncol=reps, nrow=length(row.extant),
                            dimnames=list(rownames(scores.all)[1:76]))

for ( i in 1:length(fda.pred.neck)){
  
  name.temp <- rownames(scores.all)[i]
  
  post.probs.neck[i,] <-  unlist(lapply ( fda.pred.neck[[i]] , function(x) x[name.temp,2]))
}

#View(post.probs.neck)

mean.post.probs.neck <- apply ( post.probs.neck , 1 , mean) #calculate mean value

#set class (0/1) based on 0.66 cutoff
class.temp.neck <- ifelse ( mean.post.probs.neck>0.66,1,0) 

mean ( class.temp.neck == scores.all[extant_tree_partial$tip.label,'neck']) #check % of correct predictions

## Get probs for fossil turtles
post.probs.neck.fos <- list()

for ( i in 1:length(fda.pred.neck)){
  post.probs.neck.fos[[i]] <-  do.call(rbind,lapply ( fda.pred.neck[[i]] , function(x) x[fossils,2]))
}

post.probs.neck.fos <- do.call(rbind,post.probs.neck.fos)

#View(post.probs.neck.fos)

mean.post.probs.neck.fos <- apply ( post.probs.neck.fos , 2 , mean ) #calculate mean
mean.post.probs.neck.fos

## PDF ##
## POSTERIOR PROB FOR EXTANT TURTLES ##
pdf ( 'posterior_probs_neck_extant.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:nrow(post.probs.neck)){
  hist ( sort(as.numeric(post.probs.neck[i,])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', rownames(post.probs.neck)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,reps),
         font.main = 4)
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.neck[i] , col='red', lty=2, lwd=1)
}
dev.off()


## POSTERIOR PROB FOR FOSSIL TURTLES ##
pdf ( 'posterior_probs_neck_fossils.pdf' , height = 5, width = 5, useDingbats = F)
par(mfrow=c(2,2), cex.main=0.7, cex.lab=.6, cex.axis=0.5,
    mar=c(5,4.5,4,2))

for ( i in 1:ncol(post.probs.neck.fos)){
  hist ( sort(as.numeric(post.probs.neck.fos[,i])) ,
         col=rgb(0.8,0.8,0.8,0.25),
         main=gsub('_',' ', colnames(post.probs.neck.fos)[i]),
         xlab = 'Predicted\n posterior probability', xlim = c(0,1), 
         ylab = 'Frequency\n (# of replicates)', ylim=c(0,nrow(post.probs.neck.fos)),
         font.main = 4)
  abline ( v = c(0.33,0.66) , col='black', lty=2, lwd=.8)
  abline ( v = mean.post.probs.neck.fos[i] , col='red', lty=2, lwd=1)
}
dev.off()  


### PDF ###
### MEAN PP FOR EACH PREDICTOR ###


mean.pp.all <- list ( mean.post.probs.suct.fos ,
                      mean.post.probs.durop.fos ,
                      mean.post.probs.neck.fos)

names(mean.pp.all) <- colnames(scores.all)[2:4]

## colour gradients
fossil.colour.suction <- colorRamp(c('grey90','darkblue'))
fossil.colour.durop <- colorRamp(c('grey90','tomato'))
fossil.colour.neck <- colorRamp(c('#80bdcf','grey90'))
fossil.colour.ev <- colorRamp(c('grey90','grey90','#882255'))

fossil.colour.funcs <- list ( suction = fossil.colour.suction,
                              durop = fossil.colour.durop ,
                              neck = fossil.colour.neck ,
                              ev = fossil.colour.ev)

pdf ( 'fossil_mean_probs_barplot_scores.pdf' , height = 5.5 , width = 6.5 , useDingbats = F)

layout ( matrix ( c (1,1,2,2,
                     3,3,4,4),2,byrow = T))
par(mar=c(5,8,4,1.8))

for ( i in 1:length(mean.pp.all)){
  barplot(mean.pp.all[[i]], 
          col=rgb(fossil.colour.funcs[[i]](mean.pp.all[[i]])/255),
          names.arg = gsub('_',' ', fossils), las=1,
          cex.axis = .65, cex.names = .65 , cex.lab = .7,
          xlab = 'Mean PP (pFDA)', ylab = '',
          xlim = c(0,1) , main = colnames(scores.all)[2:4][i] , cex.main = 1, horiz = T )
  abline(v  = c(0.33,0.66) , lty=2 , lwd=.8, col='black')
}

barplot(as.numeric(ev_preds), 
        col=rgb(fossil.colour.ev(ev_preds)/255),
        names.arg = gsub('_',' ', fossils), las=1,
        cex.axis = .65, cex.names = .65 , cex.lab = .7,
        xlab = 'Evasiveness index (predicted)', ylab = '',
        xlim = c(0,1) , main = 'evasiveness index' , cex.main = 1, horiz = T )
dev.off()


## SAVE AND EXPORT RESULTS ##

#If using Evers et al. (2019) topology, use following lines:
#For extant
extant_PP_Evers <- data.frame (suction = mean.post.probs.suct,
                               durophagy = mean.post.probs.durop,
                               neck_retraction = mean.post.probs.neck)

#For fossils
pfda.class.Evers <- matrix(round(unlist(mean.pp.all),2),ncol=3,nrow=length(fossils),
                           dimnames=list(fossils,
                                         colnames(scores.all)[2:4]))

write.table(pfda.class, 'pFDA_table_Evers.txt', sep = '\t')

save(fda.suct,fda.durop,fda.neck, file='fda_runs_Evers.RData')
remove(fda.suct,fda.durop,fda.neck)
save(fda.pred.suct,fda.pred.durop,fda.pred.neck, file='fda_predictions_Evers.RData')
remove(fda.pred.suct,fda.pred.durop,fda.pred.neck)

#If using Sterli et al. (2018) topology, use following lines:
#For extant
#extant_PP_Sterli <- data.frame (suction = mean.post.probs.suct,
#                               durophagy = mean.post.probs.durop,
#                              neck_retraction = mean.post.probs.neck)

#For fossils
#pfda.class.Sterli <- matrix(round(unlist(mean.pp.all),2),ncol=3,nrow=length(fossils),
#                           dimnames=list(fossils,
#                                        colnames(scores.all)[2:4]))

#write.table(pfda.class, 'pFDA_table_Sterli.txt', sep = '\t')

#save(fda.suct,fda.durop,fda.neck, file='fda_runs_Sterli.RData')
#remove(fda.suct,fda.durop,fda.neck)
#save(fda.pred.suct,fda.pred.durop,fda.pred.neck, file='fda_predictions_Sterli.RData')
#remove(fda.pred.suct,fda.pred.durop,fda.pred.neck)


end <- Sys.time()
end - start # check elapsed time ;)

