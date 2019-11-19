output_directory     <- "/project/histone_marks_miRNA/histone_marks/data_analysis/feature_correlation"
#input_directory = '/project/histone_marks_miRNA/histone_marks/data/HELA/'
input_directory = '/Users/Lisa/Dropbox/uni/Masterarbeit/data/HELA/'

library(RColorBrewer)
colors <- brewer.pal(9, "Set1")


#mapping function
get_ID <- function(info,indexOfID,sep){
  info = unlist(strsplit(info,";"))[indexOfID]
  info = unlist(strsplit(info,sep))[2]
  return(info)
}

#load HM data

#HM_dir = paste(input_directory,'histone_marks/pre_miRNA_500bp_HM_coverage',sep='')
HM_dir = paste(input_directory,'histone_marks/promoter_all_exp_filtered_merged_dis10_1000bp_HM_coverage',sep='')

HM_subdirs = dir(path=HM_dir,full.names=T)

table_HM = read.table(file=list.files(path = HM_subdirs[1],pattern="*.bed",full.names=T),header=F,sep="\t")
ID = as.character(table_HM$V9)

data_HM = data.frame(name=unlist(lapply(ID,get_ID,indexOfID=1,sep=':')),start=table_HM$V4)#promoter
#data_HM = data.frame(ID=unlist(lapply(ID,get_ID,indexOfID=1,sep='=')),start=table_HM$V4)#miRNA

for(i in 1:length(HM_subdirs)){
  table_HM = read.table(file=list.files(path = HM_subdirs[i],pattern="*.bed",full.names=T),header=F,sep="\t")
  HM = unlist(strsplit(HM_subdirs[i],"/"))
  HM = HM[length(HM)]
  data_HM = cbind(data_HM,log(table_HM$V10+1))
  colnames(data_HM)[i+2] = HM
}

#load DNA methylation data

#file = paste(input_directory,'DNA_methylation/pre_miRNA_500bp_DNAmeth_coverage/HAIB_HeLa-S3_fixed_Rep_1-hsa_pre_miRNA_500bp.gff',sep='')
file = paste(input_directory,'DNA_methylation/promoter_all_exp_filtered_merged_dis10_consensus_DNAmeth_coverage/HAIB_HeLa-S3_fixed_Rep_1-promoters_1000bp_miRNA_all_expressed_filtered_merged_dis10.gff',sep='')

table_DNA_meth = read.table(file=file,header=F,sep="\t")
ID = as.character(table_DNA_meth$V9)

data_DNA_meth = data.frame(name=unlist(lapply(ID,get_ID,indexOfID=1,sep=':')),start=table_DNA_meth$V4)#promoter
#data_DNA_meth = data.frame(ID=unlist(lapply(ID,get_ID,indexOfID=1,sep='=')),start=table_DNA_meth$V4)#miRNA

data_DNA_meth = cbind(data_DNA_meth, table_DNA_meth$V11)
data_DNA_meth = cbind(data_DNA_meth, table_DNA_meth$V10)
#colnames(data_DNA_meth) = c('ID','start','CpG','beta_value')#miRNA
colnames(data_DNA_meth) = c('name','start','CpG','beta_value')#promoter



#merge data frames
#data_merged = merge(data_DNA_meth,data_HM,by=c("start",'ID'))#miRNA

data_merged = merge(data_DNA_meth,data_HM,by=c("start",'name'))#promoter

data_filter = data_merged[data_merged$CpG != 'unknown',]
print(paste('data set with miRNAs having ß-values:', as.character(nrow(data_filter))))

data_set = data_filter[,4:ncol(data_filter)]
data_matrix = as.matrix(data_set)


###################Correlation########################

n <- ncol(data_matrix)
corResults <- matrix(ncol=n,nrow=n)
rownames(corResults) <- colnames(data_matrix)
colnames(corResults) <- colnames(data_matrix)

for(col in 1:ncol(corResults)){
  for(row in 1:nrow(corResults)){
    hm1 <- data_matrix[,col]
    hm2 <- data_matrix[,row]
    corResults[row,col] <- round(cor(hm2,hm1),1)
  }
}




z <- seq(0,1,length=n)
x <- c()
for(i in 1:n){x <- c(x,rep(z[i],n))}
y <- rep(z,n)

file_name = "correlation.pdf"
pdf(paste(output_directory,file_name,sep="/"))

par(oma=c(2,2,0,0),mar=c(5,5,3,1))
image(corResults,xaxt="n",yaxt="n",col=rev(terrain.colors(20)),main="correlation between features")
axis(1,at=z,labels=colnames(corResults),las=2)
axis(2,at=z,labels=colnames(corResults),las=2)
text(y,x,round(unlist(as.list(corResults)),1),cex=0.5)

dev.off()
  

  
  
  
