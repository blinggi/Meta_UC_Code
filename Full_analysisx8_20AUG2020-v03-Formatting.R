# Full Analysis
# Meta Deg manuscript ###
### Versions
# 13APR2020, v2 to remove Smillie parts; reorganize to show meta with singles
# 1MAY2020, got list of datasets and sample from KULeuven. a lot of overlap in GSE16879, GSE59071, GSE73661
# 15May2020, fix symbols to remove any ///
# can only use 1 set because controls used completely overlap (12 in GSE73661 contain 11 from GES59071, 6 in GSE16879 
# are in the 11 and 12 from other). Cannot use same samples as it messes up variance estimate of controls (GY message)
# will use only GSE73661
# 17JUNE2020- x10, remove GSE42911 because ages <18
# 10July2020 drop GSE22619- some uninflamed
# 11July202- remove ileal sample from 97012
# 4Aug2020-add sample size plot
# 13Aug2020-remove GSE97012- ibd and non-ibd were processed separatly acc to methods- probable reason for so many deg
# 20AUG2020-v1 -update to remove GSE9452 from enrichment top pathways because is 0 and chosing by alphebetics, also clean
# v2- fix sc to >= from > in volcano plot analysis
#v3-fix formatting for publication

rm(list=ls())
require(limma);require(Biobase);require(BiocGenerics);require(ggplot2)
wd= setwd('C:\\Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Projects/')
source('C:/Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/Code_General/DEGfunctions_v01.R')
#load previous run

# Get DEG datasets, output from DEGout codes#####
#gse16=readRDS('GSE16879/Output/GSE16879_coloncontrol_DEG.rds')#redid 23mar2020 to include all probes, (before had removed na symbol, but had not done same in other datasets)
#gse22= readRDS('GSE22619/Output/GSE22619_DEG.RDS')
gse38=readRDS('GSE38713/Output/GSE38713_DEG.rds')
#gse59=readRDS('GSE59071/Output/GSE59071_DEG.RDS')
gse73=readRDS('GSE73661/Output/GSE73661_DEG.rds')
gse87=readRDS('GSE87466/Output/GSE87466_DEG.RDS')
gse9452=readRDS('GSE9452/Output/GSE9452_DEG.RDS')
#gse97=readRDS('GSE97012/Output/GSE97012_DEG_noIleum.rds')
gse53=readRDS('GSE53306/Output/GSE53306_DEGrds')
#gse48=readRDS('GSE48634/Output/GSE48634_DEGrds')#uninflamed
gse47=readRDS('GSE47908/Output/GSE47908_DEGrds')
#gse42=readRDS('GSE42911/Output/GSE42911_DEGrds')
#gse36=readRDS('GSE36807/Output/GSE36807_DEGrds') # is uninflamed acc to paper
gse133= readRDS('GSE13367/Output/GSE13367_DEGrds')
gse114=readRDS('GSE114527/Output/GSE114527_DEGrds')# , is agilent, cant get in in , tried a lot, # is illumina looks ok as geo download 23mar2020

 
par(mfrow=c(1,1))
gse=list(   'GSE13367'=gse133,# 11 datasets
            'GSE9452'= gse9452,
       
           # 'GSE22619'=gse22,
            'GSE53306'=gse53,
            'GSE38713'=gse38,
            'GSE47908'=gse47, 
           # 'GSE42911'=gse42,
             
           # 'GSE97012'=gse97,
            'GSE73661'=gse73,
            'GSE114527'=gse114,
            'GSE87466'=gse87 )
   
gsett=sapply(gse, function(x) (x$toptable) )# 
sapply(gse,function(x) colnames(x$toptable))
for (j in 1:length(gsett)){
   colnames(gsett[[j]])[grep('symbol',tolower(colnames(gsett[[j]])))]='SYMBOL'}

# Info on dataset #####
pdata= print(sapply(gse,function(x) summary(pData(x$input))))#used either Disease or INDC, use this to create summary table
fdata= print(sapply(gse,function(x) summary(fData(x$input))))#
adata= print(sapply(gse,function(x) annotation(x$input)))#


# Volcano plots  ####
#& volcano plots were produced with the R package 'EnhancedVolcano'
require('EnhancedVolcano')
runit=F  # to suppress printing

if (runit){
for (i in 1:length(gse)){

res1=gse[[i]]$toptable
symix=grep('symbol',tolower(colnames(res1)))
tiff(paste("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/",names(gse[i]),".tiff",sep=''), units="in", width=8, height=6, res=1600)
pl=EnhancedVolcano(res1,
                lab = res1[,symix],
                x = 'logFC',
                y = 'P.Value',
               ,caption='',subtitle = '',transcriptLabSize = 2,
                pCutoff = .05,FCcutoff = log2(1.5), drawConnectors = F,legendPosition = 'right',legend=c('NS','Log2FC','P','Log2FC and P'))
print(pl)
dev.off()

}
   
}

# Enrichment on individ, using toptable ####
  #& Gene set enrichment analysis was performed with the Reactome Pathways database (reactome.org) using the the R packages
  #... ReactomePA and clusterProfiler using genes that have abs(log2(FC) > 1.5 and adjusted pvalue (benjamini) < .05
  # options************************************************
lfc=log2(1.5)
apv= .05
  #  *******************************************************

   enrupTt=  vector(mode = "list", length = length(gse))
   enrdownTt=vector(mode = "list", length = length(gse))
   #k=1
   require(dplyr)
   for (k in 1:length(gse)){
     res1     =gse[[k]]$toptable
     symix=    grep('symbol',tolower(colnames(res1)))
     apvix=    grep('adj',   tolower(colnames(res1)))
     lfcix=    grep('logfc' ,tolower(colnames(res1)))

     ixfcup=   res1[,lfcix] >   lfc
     ixfcdown= res1[,lfcix] <  -lfc
     ixp= res1[,apvix]      <   apv
    
    if ( sum(ixfcup & ixp)    >0 ){
     enrupTt[[k]]= enrichuniv( res1[ ixfcup & ixp, symix],univ=res1[,symix]  )
     require(ggplot2)
     enrdownTt[[k]]=enrichuniv(res1[ ixfcdown & ixp, symix],univ=res1[,symix]   )
     tp=  enrdownTt[[k]][1:10,]
       }}

   #plot together
#  Compile up and down enriched pathways for individual into single table
   #for Up
   require(plyr)
   a= sapply(enrupTt, length)
   ix=which(a!=0)# some with 0, need to skip
   enrupTt2       =  enrupTt[ix]#
   names(enrupTt2)=names(gse)[ix]#
   allenrTtup     = Reduce(function(x,y) merge(x,y,all.x=T,all.y=T,by='Description'),enrupTt2)
   allenrupTtsmall= allenrTtup[,c(1,grep('p.adjust',tolower(colnames(allenrTtup))))]
   colnames(allenrupTtsmall)[2:ncol(allenrupTtsmall)]= names(gse)[ix]#
   # feed above into plotting after adding meta results

   #down enriched 
   a=sapply(enrdownTt, length)
   ix=which(a!=0)# some with 0, need to skip
   enrTt2=enrdownTt[ix]#
   names(enrTt2)=names(gse)[ix]#
   allenrTt = Reduce(function(x,y) merge(x,y,all.x=T,all.y=T,by='Description'),enrTt2)
   allenr_down_Ttsmall=allenrTt[,c(1,grep('p.adjust',tolower(colnames(allenrTt))))]
   colnames(allenr_down_Ttsmall)[2:ncol(allenr_down_Ttsmall)]= names(gse)[ix]#*** 
   # feed above into plotting after adding meta results
   
   
# Run Meta Deg  ######
#check colnames
sapply(gse,function(x) colnames(x$toptable))
for (j in 1:length(gsett)){
colnames(gsett[[j]])[grep('symbol',tolower(colnames(gsett[[j]])))]='SYMBOL'
}
sapply(gsett,function(x) colnames(x))

# ^ 15may change symbols with /// to _ to allow plotting output
# also is a '.' SW-CL.36 17may

test='CH /// CH'
test2=gsub('[ /// ]','_',test)

for (i in 1:length(gse)){
#print(gsub(' /// ','_' ,gsett[[i]]$SYMBOL))
   gsett[[i]]$SYMBOL= gsub(' /// ','_' ,gsett[[i]]$SYMBOL)
  # print(gsub('[.]','_' ,gsett[[i]]$SYMBOL))
   gsett[[i]]$SYMBOL= gsub('[.]','_' ,gsett[[i]]$SYMBOL)

  # print(gsub('[.]','_' ,gsett[[i]]$SYMBOL))
   gsett[[i]]$SYMBOL= gsub('[.]','_' ,gsett[[i]]$SYMBOL)
   gsett[[i]]$SYMBOL= gsub('[/]','_' ,gsett[[i]]$SYMBOL)
}

# also is a '.' SW-CL.36
#TEMP=data.frame(gsett$GSE42911$SYMBOL)
# # & DEG metanalysis was perfomred using the R. Package 'MetavolcanoR' using the paramaters cvar=T, collaps=T (yes, is mispelled in code) , and metathr =0.1
 pcriteria="P.Value"
 foldchangecol='logFC'
 genenamecol1="SYMBOL"
 geneidcol=NULL
 collaps=T
 names(gsett)
 require(MetaVolcanoR)
 require(metafor)
 #comment out if loading previous run
# meta_degs_rem <- rem_mv(diffexp=gsett,
#                        pcriteria=pcriteria,
#                         foldchangecol=foldchangecol,
#                         genenamecol=genenamecol1,
#                         geneidcol=NULL,
#                         collaps=T,
#                         llcol='CI.L',
#                         rlcol='CI.R',
#                         vcol = NULL,
#                         cvar=T,
#                         metathr=0.1,
#                         jobname="MetaVolcano",
#                         outputfolder=".",
#                         draw='HTML',
#                         ncores=1)

#saveRDS(meta_degs_rem,'C:\\Users/bryan.linggi/Box Sync/PMED_Data(bryan.linggi@robartsinc.com 2)/MetaProjects/MetaDEG_v01/Output/mremx8_aug2020.rds')
#  # # # #load previous run
meta_degs_rem=readRDS('../MetaProjects/MetaDEG_v01/Output/mremx8_aug2020.rds')
#result tables
mr=meta_degs_rem@metaresult
#write.csv(mr,'../MetaProjects/MetaDEG_v01/Output/mr.csv')
#& the adjusted p value for meta results was calculated using the r p.adjust function with the parameter ('BH')
mrinput=  meta_degs_rem@input
mr$adjP=  p.adjust(mr$randomP,'BH')# add adjusted pvalue to table

# Graph DEG from indiv and Meta

# Meta analysisis and diagnositics ####
sc=.75*length(gse) # for Sign consistency variable, 8.25 for 11
pthresh=.05#*** 

mrup=   mr[mr$signcon     >=  sc & mr$adjP < pthresh & mr$randomSummary >  log2(1.5),];nrow(mrup)
mrdown= mr[mr$signcon     <= -sc & mr$adjP < pthresh & mr$randomSummary < -log2(1.5),];nrow(mrdown)
write.csv(mr,'../MetaProjects/MetaDEG_v01/Output/mremresultx8.csv',row.names = F)
summary(mrup)
summary(mrdown)
dput(colnames(mrup))
#export top 10 up and down DEG
write.csv(mrup[1:10,c("SYMBOL", "signcon", "randomSummary", "randomCi.lb", "randomCi.ub", 
                      "randomP", "adjP")],'../MetaProjects/MetaDEG_v01/Output/mrup.csv',row.names = F)
write.csv(mrdown[1:10,c("SYMBOL", "signcon", "randomSummary", "randomCi.lb", "randomCi.ub", 
                        "randomP", "adjP")],'../MetaProjects/MetaDEG_v01/Output/mrdown.csv',row.names = F)

# get counts
length(mr$SYMBOL[mr$error==F])
length(mrup$SYMBOL[mrup$signcon==length(gse)])
length(mrdown$SYMBOL[mrdown$signcon==-length(gse)])

#*********************************************************************************
# Stats on individual by tt####
# input is gse (list of gse lists)
# output is table of tc deg
#& Differential gene expression was determine using limma with abs(log2(FC) > 1.5 and adjust pvalue (benjamini) < .05
lfc= log2(1.5)  #log2(1) is 0
apv= .05
#
DEGtable=data.frame(matrix(nrow=length(gsett)))
unique=T #collapse probesets to gene level
for (i in 1:length(gsett))   {
   In=       gsett[[i]]
   apval=    In[which(In$adj.P.Val   < apv),]
   up=       apval[which(apval$logFC > lfc ),]
   down=     apval[which(apval$logFC < -lfc ),]
   
   rownames(DEGtable)[i]=names(gsett[i])
   colnames(DEGtable)[1]='Up'
   if (!unique){
   ifelse( nrow(up)   >0,DEGtable$Up  [i] <-   nrow(up),   DEGtable$Up[i]   <- 0)
   ifelse( nrow(down) >0,DEGtable$Down[i] <- nrow(down), DEGtable$Down[i]   <- 0)}
   
   if (unique){
      ifelse( nrow(up)   >0,DEGtable$Up  [i] <-   length(unique(up$SYMBOL)),   DEGtable$Up[i]   <- 0)
      ifelse( nrow(down) >0,DEGtable$Down[i] <- length(unique(down$SYMBOL)), DEGtable$Down[i]   <- 0)}
}

summary(DEGtable$Up)
summary(DEGtable$Down)


#DEGtable$Sum=DEGtable$Up+DEGtable$Down
# R Differential gene expression analysis was performed using R limma to compare UC patients to normal control patients, using the same settings to allow comparison acrross dattasets (lfc=log2(1.5), BH apv<.05.  there was a wide range of 

print(DEGtable) ## 
#add meta results
DEGtable_wMeta= rbind(DEGtable,c(nrow(mrup),nrow(mrdown)))
# compare deg ind to meta
median(DEGtable$Up)/nrow(mrup)
median(DEGtable$Down)/nrow(mrdown)

rownames(DEGtable_wMeta)[length(gse)+1]='Meta-x8'
DEGtable_wMeta$Dataset= rownames(DEGtable_wMeta)
require(reshape2)

melt_deg=melt(DEGtable_wMeta)

melt_deg$Dataset = factor(melt_deg$Dataset,levels =DEGtable_wMeta$Dataset)
colnames(melt_deg)=c('Dataset','Direction','Number')

tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/DEG.tiff", units="in", width=7, height=7, res=1600)

ggplot(melt_deg, aes(y=Number, x=Dataset,fill=Direction))+geom_bar(stat='identity')+theme_minimal()+scale_fill_grey() +
   theme( axis.ticks = element_line(size = .3),panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(
      color="black",size=10, angle=90,vjust = .5, hjust=1))+ylab('Number of DEG')+scale_y_continuous(expand = c(0,0)) 
dev.off()

## add sample size vs deg
DEGtable$Total = DEGtable$Up+DEGtable$Down
DEGtable$SampleSize=c(18,13,24,28,54,79,21,108)#manually entered

tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/DEGcor-total.tiff", units="in", width=7, height=7, res=1600)

ggplot(DEGtable, aes(x=SampleSize,y=Total))+geom_point()+theme_minimal()+ xlab('Sample size')+ theme( axis.ticks = element_line(size = .3),panel.border = element_blank(),
                                                                                panel.grid.major = element_blank(),
                                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(
                                                                                   color="black",size=10, angle=90,vjust = .5, hjust=1))+ylab('Number of DEG')
dev.off()
tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/DEGcor-down.tiff", units="in", width=7, height=7, res=1600)

ggplot(DEGtable, aes(x=SampleSize,y=Down))+geom_point()+theme_minimal()+ xlab('Sample size')+ theme( axis.ticks = element_line(size = .3),panel.border = element_blank(),
                                                                                                      panel.grid.major = element_blank(),
                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(
                                                                                                         color="black",size=10, angle=90,vjust = .5, hjust=1))+ylab('Number of DEG')
dev.off()
tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/DEGcor-up.tiff", units="in", width=7, height=7, res=1600)

ggplot(DEGtable, aes(x=SampleSize,y=Up))+geom_point()+theme_minimal()+ xlab('Sample size')+ theme( axis.ticks = element_line(size = .3),panel.border = element_blank(),
                                                                                                      panel.grid.major = element_blank(),
                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(
                                                                                                         color="black",size=10, angle=90,vjust = .5, hjust=1))+ylab('Number of DEG')
dev.off()

cor(DEGtable$SampleSize,DEGtable$Total,method = 'spearman')
cor(DEGtable$SampleSize,DEGtable$Up,method = 'spearman')
cor(DEGtable$SampleSize,DEGtable$Down,method = 'spearman')

# Volcano plot meta results ####
#use other volcano plotter
require(EnhancedVolcano)
# # all genes!!
# EnhancedVolcano(mr,
#                 lab = mr$SYMBOL,
#                 x = 'randomSummary',
#                 y = 'adjP',# is the column used 
#                 xlim = c(-5, 8), title = 'Meta11',subtitle = '',transcriptLabSize = 2,
#                 pCutoff = .05,FCcutoff = log2(1.5), drawConnectors = F)

#use other volcano plotter
require(EnhancedVolcano)
#genes with signcon filter
mrsc=mr[abs(mr$signcon) >= sc,] # removed = version 3* change back to >=, like in paper
   tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/volc2.tiff", units="in", width=8, height=5, res=1600)
EnhancedVolcano(mrsc,
                lab = mrsc$SYMBOL,
                x = 'randomSummary',
                y = 'adjP',
                xlim = c(-4.2, 5),ylim=c(0,25), title = '',caption='',subtitle = '',transcriptLabSize = 2,
                pCutoff = .05,FCcutoff = log2(1.5), drawConnectors = F,legendPosition = 'right',legend=c('NS','Log2FC','P','Log2FC and P'))
dev.off()

#no y axis limit
EnhancedVolcano(mrsc,
                lab = mrsc$SYMBOL,
                x = 'randomSummary',
                y = 'adjP',
                xlim = c(-4.2, 5), title = '',subtitle = '',transcriptLabSize = 2,
                pCutoff = .05,FCcutoff = log2(1.5), drawConnectors = F,legendPosition = 'right',legend=c('NS','Log2FC','P','Log2FC and P'))
length(mr$SYMBOL[mr$error=='FALSE'])

#count numbers in each group on volcano plot
fcandp=print(sum(abs(mrsc$randomSummary) >log2(1.5) & mrsc$adjP<.05)) # log and p
fconly=print(sum(abs(mrsc$randomSummary) >log2(1.5) & !mrsc$adjP<.05))# logfc only
ponly=print(sum((mrsc$adjP<.05) & !abs(mrsc$randomSummary) >log2(1.5)))# p only
none=print(sum( abs(mrsc$randomSummary) < log2(1.5) & !mrsc$adjP<.05))# neither

(total=(fcandp+fconly+ponly+none))

#diagnostics
# 
# hist( mr$signcon       ,breaks = 25)
# hist((mr$randomSummary),breaks = 100)
# hist(-log10(mr$adjP)   ,breaks=20000)
# #combo=merge(mrinput,mr)
# 
# 
head(mrup,10)
head(mrdown,10)

# Forest plots ####
#plot.new();meta_degs_rem@MetaVolcano
pcriteria="P.Value"
foldchangecol='logFC'
genenamecol1="SYMBOL"
getwd()
##up
source('../MetaProjects/MetaDEG_v01/Code/draw_forest_mod_v2.R') # mod allows manual pick of axis length


for (i in 1:2){
   require(MetaVolcanoR);require(dplyr)
   tiff(paste("..//MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/forest_",mrup$SYMBOL[i],'.tiff',sep=""), units="in", width=7, height=5, res=1600)
   
   pl=draw_forest_mod(remres=meta_degs_rem,
                  gene=mrup$SYMBOL[i],
                  genecol=genenamecol1,
                  foldchangecol=foldchangecol,
                  llcol="CI.L",
                  rlcol="CI.R",
                  jobname="MetaVolcano",
                  outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
                  draw="PDF",xlim=c(-2.5,12))
   
  
   #output stats
   print(pl)
   dev.off()
   print(mrup$SYMBOL[i])
   print( max(mrinput[mrinput$SYMBOL %in% mrup$SYMBOL[i], grep('logFC',colnames(mrinput))],na.rm=T))
   print( min(mrinput[mrinput$SYMBOL %in% mrup$SYMBOL[i], grep('logFC',colnames(mrinput))]))
   
   
}
dev.off()
draw_forest
##down
for (i in 1:2){#*** Numbers%@#@$
   tiff(paste("..//MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/forest_",mrdown$SYMBOL[i],'.tiff',sep=""), units="in", width=7, height=5, res=1600)
   
   pl=draw_forest_mod(remres=meta_degs_rem,
                  gene=mrdown$SYMBOL[i],
                  genecol=genenamecol1,
                  foldchangecol=foldchangecol,
                  llcol="CI.L",
                  rlcol="CI.R",
                  jobname="MetaVolcano",
                  outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
                  draw="PDF",xlim=c(-12, 2.5))
   print(mrdown$SYMBOL[i])
  print( min(mrinput[mrinput$SYMBOL %in% mrdown$SYMBOL[i], grep('logFC',colnames(mrinput))],na.rm=T))

  print( max(mrinput[mrinput$SYMBOL %in% mrdown$SYMBOL[i], grep('logFC',colnames(mrinput))],na.rm=T))
   print(pl)
   dev.off()
   # print(mrdown[i,])
   # print(mrinput[mrinput$SYMBOL %in% mrdown$SYMBOL[i],])
}
require(MetaVolcanoR)

# Selected forests ####
genes=c('OSM','TREM1','IL13RA2','IL1RN')
#genes=c('CMTM2','C5AR1','FGF2','GK','HGF','IL1RN','LILRA2','NAMPT','PAPPA','SNCA','SOD2','STEAP4',"ZBED3")
for (ea in genes){
   tiff(paste("..//MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/forest_",ea,'.tiff',sep=""), units="in", width=7, height=5, res=1600)
   
pl=draw_forest_mod(remres=meta_degs_rem,
               gene=ea,
               genecol=genenamecol1,
               foldchangecol=foldchangecol,
               llcol="CI.L",
               rlcol="CI.R",
               jobname="MetaVolcano",
               outputfolder="../MetaProjects/MetaDEG_v01/Output/Forests",
               draw="PDF",xlim=c(-2.5,7))
print(pl)
dev.off()
print(ea)
print( (mrup[mrup$SYMBOL %in% ea, grep('randomSummary',colnames(mrup))]))
print( (mrup[mrup$SYMBOL %in% ea, grep('randomCi.lb',colnames(mrup))]))
print( (mrup[mrup$SYMBOL %in% ea, grep('randomCi.ub',colnames(mrup))]))


print( min(mrinput[mrinput$SYMBOL %in% ea, grep('logFC',colnames(mrinput))],na.rm=T))
print( max(mrinput[mrinput$SYMBOL %in% ea, grep('logFC',colnames(mrinput))],na.rm=T))
print('----------------------------------------------------------------------')
}
# ranges of selected from input- replace 'sel' for other genes
sel='IL1RN'
mrinput[mrinput$SYMBOL %in% sel, grep('logFC',colnames(mrinput))]
min(mrinput[mrinput$SYMBOL %in% sel, grep('logFC',colnames(mrinput))],na.rm=T)
max(mrinput[mrinput$SYMBOL %in% sel, grep('logFC',colnames(mrinput))],na.rm=T)
  

# Enrichment analysis meta #####
require(vctrs)
require(clusterProfiler)
require(ReactomePA)

menrUp=enrichuniv(mrup$SYMBOL,univ = mr$SYMBOL)# mrup and mrdown defined above 
menrDown=enrichuniv(mrdown$SYMBOL, univ=mr$SYMBOL)
#& tables were created with the top 10 (lowest p.adjust) for enriched reactome pathways


#plot Enr with Individual#####
#Compile up and down enriched pathways for individual into single table
#                                for Up 
require(plyr)
#from indiv output 
head(enrupTt2)
enr_ind_meta= merge(allenrupTtsmall,menrUp[,c('Description','p.adjust')],all.y=T,all.x=T,by='Description')
colnames(enr_ind_meta)[ncol(enr_ind_meta)]='Meta-x8'

#select top 3 from each dataset (unless there are 0 enriched)
topDs=list()
k=3
for (k in 2:(length(gse)+1)){
x= enr_ind_meta[,c(1,k)]
y= x %>% arrange(x[,2]) 
topDs=c(topDs,y$Description[1:3])
}

c=enr_ind_meta[enr_ind_meta$Description %in% topDs,]
require(reshape2)
mc=melt(c)

mc$Description=factor(mc$Description)
levels(mc$Description)
mc$Desc_Group= 'Other'  
#                                check for changes
mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(1,6,7,8,9,10)]]='Immune_Signaling'
mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(2,3,4,5)]]='ECM'
#mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(10)]]='Neut'
mc$Desc_Group=factor(mc$Desc_Group)

for (g in levels(mc$Desc_Group)){
   temp=mc%>% filter(mc$Desc_Group %in% g)
   tiff(paste("..//MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/enr_",g,'.tiff',sep=""), units="in", width=11, height=6, res=1600)
   
p=ggplot(temp, aes(x=variable,y= -log10(value))) + geom_bar(position = 'dodge',stat='identity')+theme_light()+
   theme(axis.text.x = element_text(
      color='black',size=10, angle=90,vjust = .5, hjust=1))+facet_wrap(~Description)+theme( strip.text.x = element_text(
         size = 12, color='black'       )     )+xlab('Dataset') +ylab('-Log10(P value)')
print(p)
dev.off()
}
rm(enr_ind_meta)

#for down--- ---
require(plyr);require(dplyr)
#from indiv output : 
head(allenr_down_Ttsmall)
enr_ind_meta= merge(allenr_down_Ttsmall,menrDown[,c('Description','p.adjust')],all.y=T,all.x=T,by='Description')
colnames(enr_ind_meta)[ncol(enr_ind_meta)]='Meta-x8'


topDs=list()
k=3
dput(colnames(enr_ind_meta))
## !!remove GSE9452 bc has 0 enriched
for (k in c(2,4:10)){
   x= enr_ind_meta[,c(1,k)]
   y= x %>% arrange(x[,2]) 
   topDs=c(topDs,y$Description[1:3])
}

c=enr_ind_meta[enr_ind_meta$Description %in% topDs,]
#shorten some descriptions so fit in plot
dput(c$Description)
#                                   check for changes 
# shorten names for plotting
c$Description=c("Biological oxidations",        #1
                "Citric acid cycle (TCA cycle)",#2
                "Fatty acid metabolism",        #3
                "Glucuronidation",              #4
                "Phase I - Functionalization of compounds",   #5
                "Phase II - Conjugation of compounds", #6
                "Respiratory electron transport", #7
                "Respiratory electron transport, ATP ...", #8
                "Response to metal ions",  #9
                "SLC-mediated transmembrane transport", #10
                "The citric acid (TCA) cycle and respiratory ...") #11

rm(enr_ind_meta)
require(reshape2)
mc=melt(c)

#group pathways for plotting
mc$Description=factor(mc$Description)
levels(mc$Description)
mc$Desc_Group= 'Other' #todo: check if still same #
dput(levels(mc$Description))
#                                   check for changes 
mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(1,2,3,7,8,11)]]='Metabolism'
mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(4,5,6,9,10)]]='Transport and modification'


#mc$Desc_Group[mc$Description %in% levels(mc$Description)[c(10)]]='Neut'
mc$Desc_Group=factor(mc$Desc_Group)


for (g in levels(mc$Desc_Group)){
   temp=mc%>% filter(mc$Desc_Group %in% g)
   tiff(paste("..//MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/enr_",g,'.tiff',sep=""), units="in", width=11, height=6, res=1600)
   
   p=ggplot(temp, aes(x=variable,y= -log10(value))) + geom_bar(position = 'dodge',stat='identity')+theme_light()+
      theme(axis.text.x = element_text(
         color='black',size=10, angle=90,vjust = .5, hjust=1))+facet_wrap(~Description)+theme( strip.text.x = element_text(
            size = 12, color='black'       )     ) +xlab('Dataset')+ylab('-Log10(P value)') 
   print(p)
   dev.off()
}

#get genes in meta enriched to give examples of genes enriched in list ( in the meta deg not the indiv)
for(i in 1:20){
   print(menrUp$Description[i])
   idsplit=unlist(strsplit(menrUp$geneID[i],'/'))
   out=bitr(idsplit,fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb="org.Hs.eg.db")
   print(out)
}

for(i in 1:11){
   print(menrDown$Description[i])
   idsplit=unlist(strsplit(menrDown$geneID[i],'/'))
   out=bitr(idsplit,fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb="org.Hs.eg.db")
   print(out)
}
# 
# PCA of input datasets ####
syms=as.character(unique(mrinput$SYMBOL))

#use just fc
gsett_lfc=list()
for (k in 1:length(gsett)){
   temp=gsett[[k]]
   t2= temp[ match(syms , temp$SYMBOL),c("SYMBOL",'logFC')]
   gsett_lfc[[k]]=t2[!is.na(t2$SYMBOL),]
   colnames(gsett_lfc[[k]])[2]=paste('logFC',k,sep="")
}

merged_lfc <- Reduce(function(x,y) merge(x,y,by='SYMBOL',all.y=F), gsett_lfc)
colnames(merged_lfc)[2:ncol(merged_lfc)]=names(gse)
#add meta 
mrtemp=mr[mr$SYMBOL %in% syms,c('randomSummary','SYMBOL')]
merged_lfc_meta= merge(merged_lfc,mrtemp ,all.y = F, by='SYMBOL')
colnames(merged_lfc_meta)[ncol(merged_lfc_meta)]='Meta-x8'

data=data.frame(t(merged_lfc_meta[,2:ncol(merged_lfc_meta)])) # pick dataset
colnames(data)=as.character(merged_lfc_meta$SYMBOL)
#PCA
#convert NA to 0
data[is.na(data)] <- 0
#summary(data)
# metadata

pca_data=prcomp(data,scale. = F)
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2],
                         sample = rownames(data)
                         )#!!!!! adjust var names

head(df_pca_data)
require(pcaExplorer)
tiff("../MetaProjects/MetaDEG_v01/Output/Manuscript_doc/Figs-highqual/pca.tiff", units="in", width=7, height=7, res=1600)

ggplot(df_pca_data,aes(PC1,PC2))+geom_point(size=3)+#!!!!! adjust var names
   labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+ggrepel::geom_label_repel(aes(label = sample))+theme_bw()
dev.off()

# Correlation by selected genes (genes defined above (OSM, TREM etc)) ####
# # #rev, by gene
# 
cormr_gen_sp=data.frame(cor(data,method = 'spearman'))
# # #genes from above
c2=cormr_gen_sp[rownames(cormr_gen_sp)%in% genes,colnames(cormr_gen_sp) %in% genes]
c2=(sapply(c2,as.numeric))
rownames(c2)=colnames(c2)

require(gplots)
heatmap.2(c2,margins=c(12,9),col=colorRampPalette(c("blue", "white", "red"))(1024),scale='none', notecol="black", cellnote = round(c2,2),trace='none')

save.image('../MetaProjects/MetaDEG_v01/Output/mremx8_aug20_2020v1.Rdata')
