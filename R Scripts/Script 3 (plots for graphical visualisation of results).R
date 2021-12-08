#SCRIPT 3 - VARIOUS PLOTS
#This scripts contains code to generate images present in both the main text and the Supplementary Information text
#All plots were exported as vectors and later modified using Adobe Illustrator to add extra graphical elements (e.g. silhouettes, 3D skull models)

library(geomorph)
library(phytools)

####1 SCATTERPLOTS####

cols.clade <- c("#AA4499","#332288","#88ccee",
                "#44aa99","#999933","#cc6677",
                'grey88',"#117733","#ddcc77")
cols.clade2 <- cols.clade[-7]

clades <- sort(unique(turtle_data$Clade))
clades2 <- clades[-7]


#'FULL DATASET' PLOTS

size.temp <- exp(log10(size.full))

#Principal Component Analysis plot (PCA w/phylomorphospace)

X <- "PC1"
Y <- "PC2"
XLAB <- paste( X , " (" , round( PCA.results.full $ pc.summary $ importance[ 2 , X ] * 100 , 1 ) , "%)" , sep = "" )
YLAB <- paste( Y , " (" , round( PCA.results.full $ pc.summary $ importance[ 2 , Y ] * 100 , 1 ) , "%)" , sep = "" )

#par( cex.axis = .75, cex.lab=0.95)

phylomorphospace(extant_tree_full,PCA.results.full$pc.scores[extant_tree_full$tip.label,c(X,Y)], label="off",
                 node.size=c(NA,NA),lwd=.25,colors=rgb(0.6,0.6,0.6,0.5),
                 xlab = XLAB , ylab = YLAB  ,bty='l')

points(PCA.results.full$pc.scores[extant_tree_full$tip.label,c(X,Y)],
       pch=21,
       bg= cols.clade2 [as.numeric(as.factor(turtle_data[extant_tree_full$tip.label,'Clade']))] ,
       cex=size.temp[extant_tree_full$tip.label] / (min(size.temp)) )


#To plot individual clades (as shown in Supp Info)

par(mfrow=c(2,4), cex.axis = .75, cex.lab=0.95)

for ( i in 1:length(cols.clade2)){
  
  X <- "PC1"
  Y <- "PC2"
  XLAB <- paste( X , " (" , round( PCA.results.full $ pc.summary $ importance[ 2 , X ] * 100 , 1 ) , "%)" , sep = "" )
  YLAB <- paste( Y , " (" , round( PCA.results.full $ pc.summary $ importance[ 2 , Y ] * 100 , 1 ) , "%)" , sep = "" )
  
  phylomorphospace(extant_tree_full,PCA.results.full$pc.scores[extant_tree_full$tip.label,c(X,Y)], label="off",
                   node.size=c(NA,NA),lwd=.25,control=list(col.edge='grey92'),
                   xlab = XLAB , ylab = YLAB  ,bty='l')
  
  points(PCA.results.full$pc.scores[extant_tree_full$tip.label,c(X,Y)],
         pch=21, col=rgb(0.3,0.3,0.3,0.5),
         bg= rgb(0.9,0.9,0.9,0.7),lwd=.25,
         cex=size.temp / (min(size.temp)) )
  
  clade.temp <- rownames(PCA.results.full$pc.scores) %in% rownames(turtle_data)[turtle_data$Clade==clades2[i]]
  
  
  points(PCA.results.full$pc.scores[clade.temp,c(X,Y)],
         pch=21,
         bg= cols.clade2[i],
         cex=size.temp[clade.temp] / (min(size.temp)) )
  
  
  text(PCA.results.full$pc.scores[clade.temp,c(X,Y)],
       cex=.45)
  
  mtext(LETTERS[i], side = 3,at = -0.2,adj = 1,padj = -1.4,
        cex=1)
  
}


# PGLS + 2B-PLS plot

par(bty='l', mfrow=c(1,2), mar=c(5,4,6,1), asp=1, cex.lab=.7, cex.axis=.65 )

plot ( y = skull_dataset$skull_size , x = skull_dataset$carapace_size , bty='l', pch=21,
       col=ifelse(!bestmodel.full$data$neck_retraction,'black',rgb(0.5,0.5,0.5,0.3)),
       bg=ifelse(!bestmodel.full$data$neck_retraction,'#80bdcf',rgb(0.8,0.8,0.8,0.2)),
       cex=size.temp / min(size.temp) , 
       ylab = 'skull size (log10)' , 
       xlab = 'carapace size (log10)')

text(y = skull_dataset$skull_size , x = skull_dataset$carapace_size, cex=.4)

##

plot ( PLS.emarg , bty='l' , cex=1.5)

points ( PLS.emarg$XScores[,1],PLS.emarg$YScores[,1], bty='l', pch=21,
         col=ifelse(!bestmodel.full$data$neck_retraction,'black',rgb(0.5,0.5,0.5,0.3)),
         bg=ifelse(!bestmodel.full$data$neck_retraction,'#80bdcf',rgb(0.9,0.9,0.9,0.5)),
         cex=1.5)

text ( PLS.emarg$XScores[,1],PLS.emarg$YScores[,1] , cex=.4)


dev.off()

# Regression scores plot

par(mfcol=c(2,2), bty='l', cex.axis=0.65, cex.lab=0.7)

size.temp <- size.temp[extant_tree_full$tip.label]

#1 aquatic feeding x suction
plot(reg.full[,c(3,2)], pch=NA,
     ylab = paste0('aquatic feeding scores'),
     xlab = paste0('suction scores'))
points(reg.full[,3], reg.full[,2],
       pch=ifelse(bestmodel.full$data$feeds_on_water==0,22,21),
       col=rgb(0.5,0.5,0.5,0.3),
       bg=rgb(0.8,0.8,0.8,0.2),
       cex= size.temp/min(size.temp)  ) #plot all points
points(reg.full[bestmodel.full$data$feeds_on_water==0,3], reg.full[bestmodel.full$data$feeds_on_water==0,2],
       pch=22,bg='beige',
       cex=size.temp[bestmodel.full$data$feeds_on_water==0]/min(size.temp)) #points for non-aquatic feeders
points(reg.full[bestmodel.full$data$durophagous==1,3], reg.full[bestmodel.full$data$durophagous==1,2],
       pch=21,bg='tomato',
       cex=size.temp[bestmodel.full$data$durophagous==1]/min(size.temp)) #points for durophages
points(reg.full[bestmodel.full$data$suction==1,3], reg.full[bestmodel.full$data$suction==1,2],
       pch=21,bg='darkblue',
       cex=size.temp[bestmodel.full$data$durophagous==1]/min(size.temp)) #points for suction-feeders
points(reg.full['Cycloderma_frenatum',3], reg.full['Cycloderma_frenatum',2],
       pch=21,bg=colorRampPalette(c('darkblue','tomato'))(3)[2],
       cex=size.temp['Cycloderma_frenatum']/min(size.temp)) #points for Cycloderma

text(reg.full[,c(3,2)] , cex=.45)


legend('bottomright', legend=c('suction','durophage','non-aquatic feeder'), 
       pch=c(21,21,22),pt.bg=c('darkblue','tomato','beige'),
       bty='n',cex=.6, pt.cex=.8)


#2 neck retraction vs suction
plot(reg.full[,c(5,3)], 
     pch=ifelse(bestmodel.full$data$feeds_on_water==0,22,21), 
     col=ifelse(bestmodel.full$data$suction,'black',rgb(0.5,0.5,0.5,0.3)),
     bg=ifelse(bestmodel.full$data$suction,'darkblue',rgb(0.8,0.8,0.8,0.2)),
     cex=size.temp/min(size.temp),
     ylab = paste0('suction scores'),
     xlab = paste0('neck retraction scores'))
points(reg.full[,5][bestmodel.full$data$neck_retraction==0] , 
       reg.full[,3][bestmodel.full$data$neck_retraction==0] ,
       pch=21, bg='#80bdcf',
       cex=size.temp[bestmodel.full$data$neck_retraction==0]/min(size.temp))
points(reg.full['Dermochelys_coriacea',5], reg.full['Dermochelys_coriacea',3],
       pch=21, bg=colorRampPalette(c('darkblue' , '#80bdcf'))(3)[2],
       cex=size.temp['Dermochelys_coriacea']/min(size.temp))

text(reg.full[,c(5,3)] , cex=.45)

legend('bottomright', legend=c('suction','no neck retraction','non-aquatic feeder'), 
       pch=c(21,21,22),pt.bg=c('darkblue','#80bdcf','grey95'),
       bty='n',cex=.6, pt.cex=.8)

#3 hardness x durophagy (highlight durophages)
plot(reg.full[,c(6,4)], pch=NA,
     ylab = paste0('durophagy scores'),
     xlab = paste0('hardness scores'), 
     cex=size.temp/min(size.temp))
points(reg.full[bestmodel.full$data$feeds_on_water==0,c(6,4)],
       pch=22,
       cex=size.temp[bestmodel.full$data$feeds_on_water==0]/min(size.temp),
       col=rgb(0.5,0.5,0.5,0.3),
       bg=rgb(0.8,0.8,0.8,0.2))
points(reg.full[bestmodel.full$data$durophagous==0 & bestmodel.full$data$feeds_on_water==1,c(6,4)],
       pch=21,
       cex=size.temp[bestmodel.full$data$durophagous==0 & bestmodel.full$data$feeds_on_water==1]/min(size.temp),
       col=rgb(0.5,0.5,0.5,0.3),
       bg=rgb(0.8,0.8,0.8,0.2))
points(reg.full[bestmodel.full$data$durophagous==1,c(6,4)],
       pch=21,
       col='black',
       cex=size.temp[bestmodel.full$data$durophagous==1]/min(size.temp),
       bg='tomato')

text(reg.full[,c(6,4)] , cex=.45)

legend('bottomright', legend=c('non-aquatic feeder',
                               'durophage'), pch=c(22,21), pt.bg=c('grey95','tomato'),
       bty='n',cex=.6, pt.cex=.8)

#4 evasiveness x suction (highlight suction- and non-aquatic feeders)
plot(reg.full[,c(7,3)], pch=NA,
     xlab = paste0('evasiveness scores'),
     ylab = paste0('suction scores'))
points(reg.full[bestmodel.full$data$suction==0,c(7,3)],
       pch=ifelse(bestmodel.full$data$feeds_on_water==0,22,21),
       cex=size.temp[bestmodel.full$data$suction==0] / min(size.temp),
       col=rgb(0.5,0.5,0.5,0.3),
       bg=rgb(0.8,0.8,0.8,0.2))
points(reg.full[bestmodel.full$data$suction==1,c(7,3)],
       pch=21,
       cex=size.temp[bestmodel.full$data$suction==1] / min(size.temp),
       col='black',
       bg='darkblue')

text(reg.full[,c(7,3)] , cex=.45)

legend('bottomright', legend=c('suction',
                               'non-aquatic feeder'), pch=c(21,22), pt.bg=c('darkblue','grey95'),
       bty='n',cex=.6, pt.cex=.8)

dev.off()

###

#'PARTIAL DATASET' PLOTS

size.temp <- exp(log10(size.partial)) [Evers_tree$tip.label]

#Principal Component Analysis plot (PCA w/phylomorphospace)

X <- "PC1"
Y <- "PC2"
XLAB <- paste( X , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , X ] * 100 , 1 ) , "%)" , sep = "" )
YLAB <- paste( Y , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , Y ] * 100 , 1 ) , "%)" , sep = "" )

par( cex.axis = .75, cex.lab=0.95)

phylomorphospace(Evers_tree,
                 PCA.results.partial$pc.scores[Evers_tree$tip.label,c(X,Y)], label="off",
                 node.size=c(NA,NA),lwd=.25,colors=rgb(0.6,0.6,0.6,0.5),
                 xlab = XLAB , ylab = YLAB  ,bty='l')

points(PCA.results.partial$pc.scores[Evers_tree$tip.label,c(X,Y)],
       pch=ifelse(turtle_data[Evers_tree$tip.label,"Type"]=="Extant",21,22),
       bg= cols.clade [as.numeric(as.factor(turtle_data[Evers_tree$tip.label,'Clade']))] ,
       cex=size.temp / min(size.temp) )


#To plot individual clades (as shown in Supp Info)

#Show only living turtles

tree.temp <- drop.tip(Evers_tree,fossils)

par(mfrow=c(3,3), cex.axis = .75, cex.lab=0.95)

for ( i in 1:length(cols.clade2)){
  
  X <- "PC1"
  Y <- "PC2"
  XLAB <- paste( X , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , X ] * 100 , 1 ) , "%)" , sep = "" )
  YLAB <- paste( Y , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , Y ] * 100 , 1 ) , "%)" , sep = "" )
  
  phylomorphospace(tree.temp,PCA.results.partial$pc.scores[tree.temp$tip.label,c(X,Y)], label="off",
                   node.size=c(NA,NA),lwd=.25,control=list(col.edge='grey92'),
                   xlab = XLAB , ylab = YLAB  ,bty='l',
                   xlim=range(PCA.results.partial$pc.scores[,X]),
                   ylim=range(PCA.results.partial$pc.scores[,Y]))
  
  points(PCA.results.partial$pc.scores[tree.temp$tip.label,c(X,Y)],
         pch=21, col=rgb(0.3,0.3,0.3,0.5),
         bg= rgb(0.9,0.9,0.9,0.7),lwd=.25,
         cex=size.temp[tree.temp$tip.label] / min(size.temp) )
  
  clade.temp <- rownames(PCA.results.partial$pc.scores) %in% rownames(turtle_data)[turtle_data$Clade==clades2[i] & turtle_data$Type=='Extant']
  
  
  points(PCA.results.partial$pc.scores[clade.temp,c(X,Y)],
         pch=21,
         bg= cols.clade2[i],
         cex=size.temp[rownames(PCA.results.partial$pc.scores)][clade.temp] / (min(size.temp)) )
  
  
  text(PCA.results.partial$pc.scores[clade.temp,c(X,Y)],
       cex=.45)
  
  mtext(LETTERS[i], side = 3,at = -0.2,adj = 1,padj = -1.4,
        cex=1)
  
}


#Show living and extinct turltes

par(mfrow=c(3,3), cex.axis = .75, cex.lab=0.95)

for ( i in 1:length(cols.clade)){
  
  X <- "PC1"
  Y <- "PC2"
  XLAB <- paste( X , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , X ] * 100 , 1 ) , "%)" , sep = "" )
  YLAB <- paste( Y , " (" , round( PCA.results.partial $ pc.summary $ importance[ 2 , Y ] * 100 , 1 ) , "%)" , sep = "" )
  
  phylomorphospace(Evers_tree,PCA.results.partial$pc.scores[ Evers_tree$tip.label,c(X,Y)], label="off",
                   node.size=c(NA,NA),lwd=.25,control=list(col.edge='grey92'),
                   xlab = XLAB , ylab = YLAB  ,bty='l',
                   xlim=range(PCA.results.partial$pc.scores[,X]),
                   ylim=range(PCA.results.partial$pc.scores[,Y]))
  
  points(PCA.results.partial$pc.scores[Evers_tree$tip.label,c(X,Y)],
         pch=ifelse(turtle_data[Evers_tree$tip.label,'Type']=='Extant',21,22), 
         col=rgb(0.3,0.3,0.3,0.5),
         bg= rgb(0.9,0.9,0.9,0.7),lwd=.25,
         cex=size.temp[Evers_tree$tip.label] / (min(size.temp)) )
  
  clade.temp <- rownames(PCA.results.partial$pc.scores) %in% rownames(turtle_data)[turtle_data$Clade==clades[i]]
  
  points(PCA.results.partial$pc.scores[clade.temp,c(X,Y)],
         pch=ifelse(rownames(PCA.results.partial$pc.scores[clade.temp,c(X,Y)]) %in% fossils, 22,21),
         bg= cols.clade[i],
         cex=size.temp[rownames(PCA.results.partial$pc.scores)][clade.temp] / (min(size.temp)) )
  
  
  text(PCA.results.partial$pc.scores[clade.temp,c(X,Y)],
       cex=.45)
  
  mtext(LETTERS[i], side = 3,at = -0.2,adj = 1,padj = -1.4,
        cex=1)
  
}


#Regression scores plotted against predicted probabilities (pFDA) and evasiveness index

par(mfrow=c(2,2))

#SUCTION

suction.preds <- c(mean.post.probs.suct, mean.post.probs.suct.fos)

plot ( y=suction.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'suction'], bty="l" , 
       pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',21,22) ,
       col= rgb(0.5,0.5,0.5,0.3), 
       bg = rgb(fossil.colour.funcs$suction(suction.tips[Evers_tree$tip.label])/255, alpha=.85),
       cex = ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',1.25,1.75) ,
       ylab='probability of suction-feeding', xlab='suction scores', cex.axis=.65, lwd=.5)

text ( y=suction.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'suction'],
       cex=.45)


#DUROPHAGY

durop.preds <- c(mean.post.probs.durop, mean.post.probs.durop.fos)

plot ( y=durop.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'durophagous'] , bty="l" , 
       pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',21,22) ,
       col= rgb(0.5,0.5,0.5,0.3), 
       bg = rgb(fossil.colour.funcs$durop(durop.tips[Evers_tree$tip.label])/255, alpha=.85),
       cex = ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',1.25,1.75),
       ylab='probability of durophagy', xlab='durophagy scores', cex.axis=.65, lwd=.5)

text ( y=durop.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'durophagous'],
       cex=.45)


#EVASIVENESS

plot ( y=ev.tips[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'evasiveness'] , bty="l" , 
       pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',21,22) ,
       col= rgb(0.5,0.5,0.5,0.3), 
       bg = rgb(fossil.colour.funcs$ev(ev.tips[Evers_tree$tip.label])/255, alpha=.85),
       cex = ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',1.25,1.75),
       ylab='food evasiveness index', xlab='evasiveness scores', cex.axis=.65, lwd=.5)

text ( y=ev.tips[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'evasiveness'],
       cex=.45)


#NECK RETRACTION

neck.preds <- c(mean.post.probs.neck, mean.post.probs.neck.fos)

plot ( y=neck.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'neck_retraction'] , bty="l" , 
       pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',21,22) ,
       col= rgb(0.5,0.5,0.5,0.3), 
       bg = rgb(fossil.colour.funcs$neck(neck.tips[Evers_tree$tip.label])/255, alpha=.85),
       cex = ifelse(scores.all[Evers_tree$tip.label,'type']=='extant',1.25,1.75),
       ylab='probability of neck retraction', xlab='neck retraction scores', cex.axis=.65, lwd=.5)

text ( y=neck.preds[Evers_tree$tip.label],
       x=scores.all[Evers_tree$tip.label,'neck_retraction'],
       cex=.45)

dev.off()

# To plot correlation between PPs calculated using either topology (as in Supplementary Figure 13)

pfda.results <- list(pfda.class.Evers , pfda.class.Sterli)

par(mfrow=c(1,3))

for ( i in 1:3){
  plot(pfda.results[[1]][,i],pfda.results[[2]][,i],cex.axis=.65,cex.lab=.75,
       bty='l', bg=rgb(0.8,0.8,0.8,0.4),pch=21,cex=1.5,lwd=.6,
       xlab=paste0(colnames(pfda.results[[1]])[i],'_Evers'),
       ylab=paste0(colnames(pfda.results[[2]])[i],'_Sterli'))
  abline(v=0.66,h=0.66,col=rgb(1,0,0,0.5),lwd=.6,lty=2)
  text(pfda.results[[1]][,i],pfda.results[[2]][,i],
       cex=.45)
  pfda.cor <- round(cor(do.call(cbind,lapply(pfda.results,function(X)X[,i]))),3)[1,2]*100
  legend('topleft',legend=paste0('cor = ',pfda.cor,'%'),
         bty='n',cex=.8)
}
dev.off()

####

####2 SHAPE DEFORMATION PLOTS ####

#PCA - maximum and minimum PC1-2 scores shapes superimposed on the consensus shape

# Switch between the 'GPA' objects (i.e. 'GPA.full' or 'GPA.partial') and 'sliders'
#objects (i.e. 'sliders.full' or 'sliders.partial' to visualise landmark configurations of either dataset

open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 2, nc = 2, byrow = TRUE, sharedMouse = TRUE)
choose <- "PC1min"
plot3d( GPA.full$consensus , col = "grey80" , size = 4 , box = "n",aspect = 'iso')	
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( GPA.full$consensus[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'grey80' ) }
plot3d( PCA.results.full $ pc.shapes[[ choose ]] , col = "black" , size = 5 , box = "n", add = T)
aspect3d( "iso" )
#title3d( main = choose, line=1 , cex = 1)
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( PCA.results.full $ pc.shapes[[ choose ]][ lines.plot.temp[j,], ] , 
                                                 lwd = 2 , col = 'black' ) }

next3d()
choose <- "PC1max"
plot3d( GPA.full$consensus , col = "grey80" , size = 4 , box = "n",aspect = 'iso')	
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( GPA.full$consensus[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'grey80' ) }
plot3d( PCA.results.full $ pc.shapes[[ choose ]] , col = "black" , size = 5 , box = "n", add = T)
aspect3d( "iso" )
#title3d( main = choose, line=1 , cex = 1)
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( PCA.results.full $ pc.shapes[[ choose ]][ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'black' ) }	

next3d()
choose <- "PC2min"
plot3d( GPA.full$consensus , col = "grey80" , size = 4 , box = "n",aspect = 'iso')	
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( GPA.full$consensus[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'grey80' ) }
plot3d( PCA.results.full $ pc.shapes[[ choose ]] , col = "black" , size = 5 , box = "n", add = T)
aspect3d( "iso" )
#title3d( main = choose, line=1 , cex = 1)
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( PCA.results.full $ pc.shapes[[ choose ]][ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'black' ) }	

next3d()
choose <- "PC2max"
plot3d( GPA.full$consensus , col = "grey80" , size = 4 , box = "n",aspect = 'iso')	
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( GPA.full$consensus[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'grey80' ) }
plot3d( PCA.results.full $ pc.shapes[[ choose ]] , col = "black" , size = 5 , box = "n", add = T)
aspect3d( "iso" )
#title3d( main = choose, line=1 , cex = 1)
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( PCA.results.full $ pc.shapes[[ choose ]][ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'black' ) }

#D-PGLS : plot shape changes associated with individual predictors in the best D-PGLS model

# Best D-PGLS model using the 'full dataset' (i.e. object 'bestmodel.full')

#Get coefficients from 'procD.pgls' object
coefficients.temp <- bestmodel.full$pgls.coefficients 

#Get 1st and 3rd quartiles of continuous variables to avoid potential outliers (minimum and maximum values)
skull.size.quart <-  quantile(bestmodel.full$data$skull_size,prob=c(.25,.75))
hard.quart <-  quantile(bestmodel.full$data$hardness,prob=c(.25,.75))
ev.quart <-  quantile(bestmodel.full$data$evasiveness,prob=c(.25,.75))

# With this next step we define which variables to show the shape deformation for.
# Create two lists. The elements of the list correspond to the variables in the model + the Intercept.
# Let's say we want to see the deformation of a binary variable. In the first list we set the variable as absent,
# and in the second as present. The other variables in the model are kept constant in both lists. In the case of 
# e.g. allometry, set 'skull_size' values as the minimum/maximum of the 'skull.size.quart' object created.

#Allometry variation
variables1 <- list( Intercept = 1  , skull_size = min(skull.size.quart) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,   neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness) , evasiveness = mean(bestmodel.full$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = max(skull.size.quart) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0  ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness) ,evasiveness = mean(bestmodel.full$data$evasiveness) )

# Now to show the deformation for the ecological/functional variables, set the 'skull_size' value as the mean of your vector of sizes 
# contained within the model object

#Aquatic feeding (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 0 , suction = 0  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
#Suction (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 1  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
#Durophagy (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 1 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
#Neck retraction (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,   neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 1 ,   neck_retraction = 0, hardness = mean(bestmodel.full$data$hardness),evasiveness = mean(bestmodel.full$data$evasiveness) )
#Hardness index variation
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,   neck_retraction = 1, hardness = min(hard.quart) ,evasiveness = mean(bestmodel.full$data$evasiveness))
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,   neck_retraction = 1, hardness = max(hard.quart) ,evasiveness = mean(bestmodel.full$data$evasiveness))
#Evasiveness index variation
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,   neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness) ,evasiveness = min(ev.quart))
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.full$data$skull_size) ,  feeds_on_water = 1 , suction = 0  , durophagous = 0 ,  neck_retraction = 1, hardness = mean(bestmodel.full$data$hardness) ,evasiveness = max(ev.quart))

# Change the names of these lists to match the names in your 'coefficients.temp' object

names( variables1 ) <- rownames ( coefficients.temp )
names( variables2 ) <- rownames ( coefficients.temp )

# Create lists to which you will save the values calculated for each landmark in each shape deformation	
coefficient.list.temp1 <- list()
coefficient.list.temp2 <- list()

# Calculate the response values in two steps. First: multiplying the coefficients of your best model by the values in your lists
# 'variables1' and 'variables2' ...

# (step 1)
for( row.temp in 1:nrow( coefficients.temp ) ) {
  coefficient.list.temp1[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables1[[ rownames( coefficients.temp )[ row.temp ] ]]
  coefficient.list.temp2[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables2[[ rownames( coefficients.temp )[ row.temp ] ]]
}

# (step 2) 
# ... and second: by summing the values you obtained, in a way that you have such formula: 
# e.g. Y (shape) = Intercept + eco_trait1*coefficient_eco_trait1 + eco_trait2*coefficient_eco_trait2
# and putting these results in a 3-column matrix (X,Y,Z)

shape1 <- matrix( apply( do.call( rbind , coefficient.list.temp1 ) , 2 , sum ) , ncol = 3 , byrow = T ) 
shape2 <- matrix( apply( do.call( rbind , coefficient.list.temp2 ) , 2 , sum ) , ncol = 3 , byrow = T )

## Create colour scheme to visualise distance between landmarks of the two shape configurations

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

point.distance.scale <- colorRamp(c("lightgrey" ,'red')) # define the colours you want

point.distances <- matrix(nrow=dim(shape2)[1],ncol=dim(shape2)[2])
for ( i in 1:nrow(shape2) ){
  for ( j in 1:ncol(shape2)){
    
    # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
    point.distances[i,j] <- Edist(shape1[i,j] , shape2[i,j]) 
  }
}	

#point.distances

# normalise the distances so they range from 0 to 1
point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

# and then you're able to create a colorRamp that goes from 0 (grey) to 1 (red)
point.colours <- point.distance.scale(point.distances.norm)

# 3D plot of shape deformation

open3d()
plot3d(shape2, size = 9, col=rgb(point.colours,maxColorValue = 255) , box = "n" , aspect="iso")
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape2[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = rgb(point.colours,maxColorValue = 255) ) }

points3d(shape1, size = 7, col='lightgrey' )
lines.plot.temp <- rbind( sliders.full[ , 1:2 ] , sliders.full[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape1[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'lightgrey' ) }


# Best D-PGLS model using the 'partial dataset' (i.e. object 'bestmodel.partial')
# This follows the same rationale as with the 'full dataset' results visualisation, so the following lines don't have as many comments

#Get coefficients from 'procD.pgls' object
coefficients.temp <- bestmodel.partial$pgls.coefficients 

#Get 1st and 3rd quartiles of continuous variables to avoid potential outliers (minimum and maximum values)
skull.size.quart <-  quantile(bestmodel.partial$data$skull_size,prob=c(.25,.75))
ev.quart <-  quantile(bestmodel.partial$data$evasiveness,prob=c(.25,.75))

# Inspect individual variables

#Allometry variation
variables1 <- list( Intercept = 1  , skull_size = min(skull.size.quart) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = max(skull.size.quart) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
#Suction (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 1 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
#Durophagy (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 1  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
#Neck retraction (presence/absence)
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 0 , evasiveness = mean(bestmodel.partial$data$evasiveness) )
#Evasiveness index variation
variables1 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = min(ev.quart) )
variables2 <- list( Intercept = 1  , skull_size = mean(bestmodel.partial$data$skull_size) , suction = 0 , durophagous = 0  , neck_retraction = 1 , evasiveness = max(ev.quart) )

# Change the names of these lists to match the names in your 'coefficients.temp' object

names( variables1 ) <- rownames ( coefficients.temp )
names( variables2 ) <- rownames ( coefficients.temp )

# Create lists to save the values calculated for each landmark in each shape deformation

coefficient.list.temp1 <- list()
coefficient.list.temp2 <- list()

# Calculate the response values

for( row.temp in 1:nrow( coefficients.temp ) ) {
  coefficient.list.temp1[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables1[[ rownames( coefficients.temp )[ row.temp ] ]]
  coefficient.list.temp2[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables2[[ rownames( coefficients.temp )[ row.temp ] ]]
}

shape1 <- matrix( apply( do.call( rbind , coefficient.list.temp1 ) , 2 , sum ) , ncol = 3 , byrow = T ) 
shape2 <- matrix( apply( do.call( rbind , coefficient.list.temp2 ) , 2 , sum ) , ncol = 3 , byrow = T )

# Visualise landmark configurations for each condition

point.distance.scale <- colorRamp(c("lightgrey" ,'red')) # define the colours you want

point.distances <- matrix(nrow=dim(shape2)[1],ncol=dim(shape2)[2])
for ( i in 1:nrow(shape2) ){
  for ( j in 1:ncol(shape2)){
    
    point.distances[i,j] <- Edist(shape1[i,j] , shape2[i,j]) 
  }
}	

#point.distances

point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

point.colours <- point.distance.scale(point.distances.norm)

# 3D plot

open3d()
plot3d(shape2, size = 9, col=rgb(point.colours,maxColorValue = 255) , box = "n" , aspect="iso")
lines.plot.temp <- rbind( sliders.partial[ , 1:2 ] , sliders.partial[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape2[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = rgb(point.colours,maxColorValue = 255) ) }

points3d(shape1, size = 7, col='lightgrey' )
lines.plot.temp <- rbind( sliders.partial[ , 1:2 ] , sliders.partial[ , 2:3 ] )
for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape1[ lines.plot.temp[j,] , ] , 
                                                 lwd = 2 , col = 'lightgrey' ) }


####3 TREE PLOTS (mapping traits onto the topology of Evers et al. 2019)####

# Suction #

suction.tips <- setNames ( c(scores.all$suct[1:Ntip(extant_tree_partial)],
                             mean.post.probs.suct.fos),rownames(scores.all))

suction.tips <- round(suction.tips,3)
suction.tips <- suction.tips [Evers_tree$tip.label]

plotBranchbyTrait(Evers_tree,suction.tips,mode='tips', edge.width=2, legend = F,
                  palette = colorRampPalette(c('grey90','darkblue')), show.tip.label = F)

cols.temp <- colorRamp(c('grey90','darkblue'))

tiplabels(tip=1:Ntip(Evers_tree), pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='fossil',15,16) , 
          col=rgb(cols.temp(suction.tips[Evers_tree$tip.label])/255),
          cex= 1.25  )
anc.temp <- fastAnc(Evers_tree,x=suction.tips)

nodelabels(node=1:Evers_tree$Nnode+Ntip(Evers_tree), pch=16,
           col = rgb(cols.temp(anc.temp)/255) )

tiplabels(tip=1:Ntip(Evers_tree), text=1:Ntip(Evers_tree),frame = 'none',offset = 7,
          cex= .45  )


# Durophagy #

durop.tips <- setNames ( c(scores.all$durop[1:Ntip(extant_tree_partial)],
                           mean.post.probs.durop.fos),rownames(scores.all))

durop.tips <- round(durop.tips,3)
durop.tips <- durop.tips [Evers_tree$tip.label]

plotBranchbyTrait(Evers_tree, durop.tips, mode='tips', edge.width=2, legend = F,
                  palette = colorRampPalette(c('grey90','tomato')), show.tip.label = F,
                  direction='leftwards')

cols.temp <- colorRamp(c('grey90','tomato'))

tiplabels(tip=1:Ntip(Evers_tree), pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='fossil',15,16) , 
          col=rgb(cols.temp(durop.tips[Evers_tree$tip.label])/255),
          cex= 1.25  )
anc.temp <- fastAnc(Evers_tree,x=durop.tips)

nodelabels(node=1:Evers_tree$Nnode+Ntip(Evers_tree), pch=16,
           col = rgb(cols.temp(anc.temp)/255) )

#tiplabels(tip=1:Ntip(Evers_tree), text=1:Ntip(Evers_tree),frame = 'none',offset = 4.5,
#         cex= .45  )


# Evasiveness #

ev.tips <- setNames ( c(scores.all$ev[1:Ntip(extant_tree_partial)],
                        ev_preds),rownames(scores.all))

ev.tips <- round(ev.tips,3)
ev.tips <- ev.tips [Evers_tree$tip.label]

plotBranchbyTrait(Evers_tree,ev.tips,mode='tips', edge.width=2, legend = F,
                  palette = colorRampPalette(c('grey90','grey90','#882255')), show.tip.label = F)

cols.temp <- colorRamp(c('grey90','grey90','#882255'))

tiplabels(tip=1:Ntip(Evers_tree), pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='fossil',15,16) , 
          col=rgb(cols.temp(ev.tips[Evers_tree$tip.label])/255),
          cex= 1.25  )
anc.temp <- fastAnc(Evers_tree,x=ev.tips)

nodelabels(node=1:Evers_tree$Nnode+Ntip(Evers_tree), pch=16,
           col = rgb(cols.temp(anc.temp)/255) )
tiplabels(tip=1:Ntip(Evers_tree), text=1:Ntip(Evers_tree),frame = 'none',offset = 7,
          cex= .45  )


# Neck retraction #

neck.tips <- setNames ( c(scores.all$neck[1:Ntip(extant_tree_partial)],
                          mean.post.probs.neck.fos),rownames(scores.all))

neck.tips <- round(neck.tips,3)
neck.tips <- neck.tips [Evers_tree$tip.label]

plotBranchbyTrait(Evers_tree,neck.tips,mode='tips', edge.width=2, legend = F,
                  palette = colorRampPalette(c('#80bdcf' , 'grey90')), show.tip.label = F,
                  direction='leftwards')

cols.temp <- colorRamp(c('#80bdcf' , 'grey90'))

tiplabels(tip=1:Ntip(Evers_tree), pch=ifelse(scores.all[Evers_tree$tip.label,'type']=='fossil',15,16) , 
          col=rgb(cols.temp(neck.tips[Evers_tree$tip.label])/255),
          cex= 1.25  )
anc.temp <- fastAnc(Evers_tree,x=neck.tips)

nodelabels(node=1:Evers_tree$Nnode+Ntip(Evers_tree), pch=16,
           col = rgb(cols.temp(anc.temp)/255) )

#tiplabels(tip=1:Ntip(Evers_tree), text=1:Ntip(Evers_tree),frame = 'none',offset = 4.5,
#         cex= .45  )

dev.off()
