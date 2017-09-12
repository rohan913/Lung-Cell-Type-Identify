##Scoring attempt R##
rm(list=ls()) #clean environment
library(dplyr)
setwd("~/Desktop/Test of Scoring/")

###Training Dataset (more can be added)
markers=read.csv('~/Desktop/MRGCCAMarkers/CCAMRGmarkers.csv')


last=last(markers$cluster)
clusters=c(0:last)

# # indexes of clusters  
# for (i in clusters){ print(i)}
# for (i in clusters){
#   first=(first(grep(paste("\\b",i,"\\b",sep = ""),markers$cluster)))
#   last=(last(grep(paste("\\b",i,"\\b",sep=""),markers$cluster)))
  
  
#}

# #this works alone to get you all the subset files lol-idt I need to writ them all though 
#last=last(markers$cluster)
# for (i in c(1:last)){
#   tmp=subset(markers,cluster==i)
#   filename=paste('cluster',i, "markers.csv")
#   write.csv(tmp,filename)
# }  

###1 split long file into shorter tables not needed anymore
# grep_SC14=grep("SC14_",c9@cell.names)
# SC14 <- SubsetData(CCAMRG, cells.use=c9@cell.names[first(grep_SC14):last(grep_SC14)])
# grep_SC15=grep("SC15_",c9@cell.names)
# SC15 <- SubsetData(CCAMRG, cells.use=c9@cell.names[first(grep_SC15):last(grep_SC15)])
#   
# #2 develop scoring system thinking list all cell types and simultaneously keep track of scores for each one 
# ###Fingerprints for each one: 
# AT2genes=c("Sftpa1","Sftpb","Cxcl15","Scd1","Slc34a2","Sftpd","Sfta2","S100g","Chil1","Lamp3")
# ClassMon=c('Plac8','S100a4','Ly6c2',"Ms4a6c","Ifitm3","Ccr2","F13a1","Ifi27l2a","Tgfbi","Ccl9","Fn1","Ms4a4c","Smpdl3a","Lgals3","Ifitm6","Lst1","Alox5ap","Pld4","Clec4a1","Ms4a6b","Vcan","Ifi204","Tyrobp","Emilin2","Ms4a6d","Prdx5","Gsr","Gpx1","Plbd1","Fcer1g","Clec4a3","Zeb2","Cybb","Vim","Pirb","Igsf6","S100a6","Fam96a","Ifitm2","Sirpb1c","Crip1","Ifi30","F10","Tpd52","Ly86","Emb","AB124611","Gm9733","Vsir","Arpc1b")

##easier way:

positions<-c()
for (i in clusters){
  first=(first(grep(paste("\\b",i,"\\b",sep = ""),markers$cluster)))
  positions<-c(positions,first)
}
celltypes=c('AT2','ClassMon','AMs','BCell','Neutrophils','Tcell','CD8Tcell','NonClassMon','NK Cell','IMs','Club/Cill','Airway Epith','Lymph endo','DC1','nonlymph endo','cellcycle','naiveT','PDC','Fibroblasts','DC2','Mesothelial Cell','AT1','lymph progen','smooth musc','mast cell')

genesss=c()
for (i in positions){
  top50=(as.character(markers$gene[i:(i+50)]))
  genesss=c(genesss,top50)
}


for (i in 1:length(celltypes)){assign(as.character(celltypes[i]),genesss[(1+(50*(i-1))):(50*i)])}

### now trying to set the scoring areas into 10 entries for x5 ie top 50 genes

###markers2 is the nonref data you already have the top 50 in each celltype to compare
setwd('~/Desktop/1502/')
markers2=read.csv("~/Desktop/1502/1502markers.csv")
last=last(markers2$cluster)
clusters=c(0:last)
score=c()
for (i in clusters) {
  tmp=subset(markers2,cluster==i)
  tmp1=as.character(tmp$gene[1:10])
  tmp2=as.character(tmp$gene[11:20])
  tmp3=as.character(tmp$gene[21:30])
  tmp4=as.character(tmp$gene[31:40])
  tmp5=as.character(tmp$gene[41:50])
  
  for (m in 1:length(celltypes)){
    assign(paste0("score",celltypes[m]),0)
    checking=as.array(get(celltypes[m]))
      for (n in checking){
        if (is.element(n,tmp1) == TRUE){assign(paste0('score',celltypes[m]),get(paste0('score',celltypes[m]))+100)}
        if (is.element(n,tmp2) == TRUE){assign(paste0('score',celltypes[m]),get(paste0('score',celltypes[m]))+75)}
        if (is.element(n,tmp3) == TRUE){assign(paste0('score',celltypes[m]),get(paste0('score',celltypes[m]))+50)}
        if (is.element(n,tmp4) == TRUE){assign(paste0('score',celltypes[m]),get(paste0('score',celltypes[m]))+25)}
        if (is.element(n,tmp5) == TRUE){assign(paste0('score',celltypes[m]),get(paste0('score',celltypes[m]))+10)}
      }
    score=c(score,(paste(i,",",celltypes[m],paste(',',get(paste0('score',celltypes[m]))))))
  }
}

write.table(score, file = "allscores1502.csv", quote = FALSE, sep = " ", 
            eol = "\n", row.names = FALSE, col.names = FALSE)
#write colnames
scoretable=read.csv('allscores1502.csv',header=FALSE)
colnames(scoretable)<-c("cluster","celltype","score")
write.table(scoretable, file = "allscores1502.csv",sep=',',row.names = FALSE)

#final score table for further analysis
scoretable=read.csv('allscores1502.csv')


for (i in clusters){
  first=(first(grep(paste("\\b",i,"\\b",sep = ""),scoretable$cluster)))
  last=(last(grep(paste("\\b",i,"\\b",sep=""),scoretable$cluster)))
  subsetted=(scoretable[first:last,])
  #ordered=tail(order(subsetted$score))
  #for (i in ordered){print(subsetted[i,])}
  ordered=tail(order(subsetted$score),n=3)
  for (i in ordered){print(subsetted[i,])}
}
