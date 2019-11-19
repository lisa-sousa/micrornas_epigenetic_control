###############################
#Directories                  #
###############################

input_directory = '/project/histone_marks_miRNA/data/HELA/'
#input_directory = '/project/histone_marks_miRNA/data/IMR90/'

###############################
#Functions                    #
###############################

#mapping function
get_ID <- function(info,indexOfID,sep){
  info = unlist(strsplit(info,";"))[indexOfID]
  info = unlist(strsplit(info,sep))[2]
  return(info)
}

################################
#Load Data                     #
################################

#load pre-miRNA HM data
HM_dir = paste(input_directory,'histone_marks/pre_miRNA_500bp_HM_coverage',sep='')
HM_subdirs = dir(path=HM_dir,full.names=T)

table_HM = read.table(file=list.files(path = HM_subdirs[1],pattern="*.bed",full.names=T),header=F,sep="\t")
ID = unlist(lapply(as.character(table_HM$V9),get_ID,indexOfID=1,sep='='))
name = unlist(lapply(as.character(table_HM$V9),get_ID,indexOfID=3,sep='='))

data_HM = data.frame(miRNA_name = name, miRNA_ID = ID)

for(i in 1:length(HM_subdirs)){
  table_HM = read.table(file=list.files(path = HM_subdirs[i],pattern="*.bed",full.names=T),header=F,sep="\t")
  HM = unlist(strsplit(HM_subdirs[i],"/"))
  HM = HM[length(HM)]
  data_HM = cbind.data.frame(data_HM,log(table_HM$V10+1),stringsAsFactors = FALSE)
  colnames(data_HM)[i+2] = paste('m',HM,sep='_')
}
data_HM_miRNA = data_HM
head(data_HM_miRNA)


#load promoter HM data                                           
HM_dir = paste(input_directory,'histone_marks/promoter_encode_all_exp_filtered_merged_dis10_1000bp_HM_coverage',sep='')
HM_subdirs = dir(path=HM_dir,full.names=T)

table_HM = read.table(file=list.files(path = HM_subdirs[1],pattern="*.bed",full.names=T),header=F,sep="\t")
name = unlist(lapply(as.character(table_HM$V9),get_ID,indexOfID=1,sep=':'))

data_HM = data.frame(miRNA_name = name,start=table_HM$V4)

for(i in 1:length(HM_subdirs)){
  table_HM = read.table(file=list.files(path = HM_subdirs[i],pattern="*.bed",full.names=T),header=F,sep="\t")
  HM = unlist(strsplit(HM_subdirs[i],"/"))
  HM = HM[length(HM)]
  data_HM = cbind.data.frame(data_HM,log(table_HM$V10+1),stringsAsFactors = FALSE)
  colnames(data_HM)[i+2] = paste('p',HM,sep='_')
}

data_HM_promoter = data_HM
head(data_HM_promoter)

#load pre-miRNA DNA methylation data
file = paste(input_directory,'DNA_methylation/pre_miRNA_500bp_DNAmeth_coverage/HAIB_HeLa-S3_fixed_Rep_1-hsa_pre_miRNA_500bp.gff',sep='')
table_DNA_meth = read.table(file=file,header=F,sep="\t")
ID = unlist(lapply(as.character(table_DNA_meth$V9),get_ID,indexOfID=1,sep='='))
name = unlist(lapply(as.character(table_DNA_meth$V9),get_ID,indexOfID=3,sep='='))

data_DNA_meth = data.frame(miRNA_name = name, miRNA_ID = ID)

data_DNA_meth = cbind.data.frame(data_DNA_meth, as.character(table_DNA_meth$V11),stringsAsFactors = FALSE)
data_DNA_meth = cbind.data.frame(data_DNA_meth, table_DNA_meth$V10,stringsAsFactors = FALSE)
colnames(data_DNA_meth) = c('miRNA_name','miRNA_ID','m_CpG','m_beta_value')

data_DNA_meth$m_CpG[data_DNA_meth$m_CpG == 'CpG'] = 1
data_DNA_meth$m_CpG[data_DNA_meth$m_CpG == 'unknown'] = 0
data_DNA_meth$m_CpG = as.factor(data_DNA_meth$m_CpG)
x = sum(data_DNA_meth$m_CpG == 0)
sampled_beta_values = runif(x,min=0.2,max=0.6)
data_DNA_meth$m_beta_value[data_DNA_meth$m_CpG == 0] = sampled_beta_values

data_DNA_meth_miRNA = data_DNA_meth
head(data_DNA_meth_miRNA)

#load promoter DNA methylation data                              
file = paste(input_directory,'DNA_methylation/promoter_encode_all_exp_filtered_merged_dis10_consensus_DNAmeth_coverage/HAIB_HeLa-S3_fixed_Rep_1-promoters_1000bp_miRNA_all_expressed_filtered_merged_dis10.gff',sep='')
table_DNA_meth = read.table(file=file,header=F,sep="\t")
name = unlist(lapply(as.character(table_DNA_meth$V9),get_ID,indexOfID=1,sep=':'))

data_DNA_meth = data.frame(miRNA_name = name,start=table_DNA_meth$V4)

data_DNA_meth = cbind.data.frame(data_DNA_meth, as.character(table_DNA_meth$V11),stringsAsFactors = FALSE)
data_DNA_meth = cbind.data.frame(data_DNA_meth, table_DNA_meth$V10,stringsAsFactors = FALSE)
colnames(data_DNA_meth) = c('miRNA_name','start','p_CpG','p_beta_value')

data_DNA_meth$p_CpG[data_DNA_meth$p_CpG == 'CpG'] = 1
data_DNA_meth$p_CpG[data_DNA_meth$p_CpG == 'unknown'] = 0
data_DNA_meth$p_CpG = as.factor(data_DNA_meth$p_CpG)
x = sum(data_DNA_meth$p_CpG == 0)
sampled_beta_values = runif(x,min=0.2,max=0.6)
data_DNA_meth$p_beta_value[data_DNA_meth$p_CpG == 0] = sampled_beta_values

data_DNA_meth_promoter = data_DNA_meth
head(data_DNA_meth_promoter)

#load mature-miRNA expression data
file = paste(input_directory,'encode_smallRNA_seq/miRNA_expression_RPM/expression_values.txt',sep='')   

table_expression = read.table(file=file,header=T,sep="\t")
ID = unlist(lapply(as.character(table_expression$group),get_ID,indexOfID=2,sep='='))                             
name = unlist(lapply(as.character(table_expression$group),get_ID,indexOfID=1,sep='='))                           

data_expression = data.frame(miRNA_name = name, miRNA_ID = ID)
data_expression = cbind.data.frame(data_expression, table_expression[,10:11],stringsAsFactors = FALSE)

data_expression = data_expression[data_expression$expression != 1, ]
print(paste('data set without low expressed miRNAs:', nrow(data_expression)))

data_expression$expression[data_expression$expression == 2] = 1
data_expression$expression[data_expression$expression == 0] = 0
data_expression$expression = as.factor(data_expression$expression)

head(data_expression)

#load mature-miRNA genomic classification
file_intergenic = '/project/histone_marks_miRNA/data/raw_data/miRNA/miRNA_classification/miriad_human_intergenic_v.2014.txt'
table_intergenic = read.table(file_intergenic,sep='\t',header=T)
table_intergenic = data.frame(miRNA_name = table_intergenic$precursor, miRNA_ID = table_intergenic$mirbase_id, class = rep('intergenic',nrow(table_intergenic)))
table_intergenic = unique(table_intergenic)

file_intragenic = '/project/histone_marks_miRNA/data/raw_data/miRNA/miRNA_classification/miriad_human_intragenic_v.2014.txt'
table_intragenic = read.table(file_intragenic,sep='\t',header=T)
table_intragenic = data.frame(miRNA_name = table_intragenic$precursor, miRNA_ID = table_intragenic$mirbase_id, class = rep('intragenic',nrow(table_intragenic)))
table_intragenic = unique(table_intragenic)

data_class = rbind.data.frame(table_intergenic,table_intragenic)
data_class$class = as.character(data_class$class)
data_class$class[data_class$class == 'intergenic'] = 1
data_class$class[data_class$class == 'intragenic'] = 0
data_class$class = as.factor(data_class$class)

head(data_class)

################################
#Merge Data                    #
################################

#create data_set and target
get_data_set <- function(data_merged,i){
  
  data_set = data_merged[,c(i:ncol(data_merged))]
  rownames(data_set) = make.names(data_merged$miRNA_name,unique=T)
  
  #normalization -> no mormalization for dummy variables (factors)
  for(j in 1:(ncol(data_set))){ 
    if(!is.factor(data_set[,j])){
      data_set[,j] = (data_set[,j] - mean(data_set[,j])) / sd(data_set[,j])
    }
  } 
  return(data_set)
}


####Full Data Set
#promoter (HM + DNAmeth) + pre-miRNA (HM + DNAmeth) + class
data_merged_miRNA = merge(data_expression,data_DNA_meth_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged_miRNA = merge(data_merged_miRNA,data_HM_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged_miRNA = unique(data_merged_miRNA)

data_merged_promoter = merge(data_expression,data_DNA_meth_promoter,by = c('miRNA_name'))
data_merged_promoter = merge(data_merged_promoter,data_HM_promoter,by = c('miRNA_name','start'))
data_merged_promoter = unique(data_merged_promoter)

data_merged = merge(data_merged_promoter,data_merged_miRNA,by=c('miRNA_ID','miRNA_name','log_RPM','expression'))            
data_merged = merge(data_merged,data_class,by=c('miRNA_name','miRNA_ID'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,6)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/pre_miRNA_promoter_class.RData')



####Partial Data Sets
#only pre_miRNA DNA meth
data_merged = merge(data_expression,data_DNA_meth_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,5)
save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/pre_miRNA_DNAmeth.RData')



#only promoter DNA meth
data_merged = merge(data_expression,data_DNA_meth_promoter,by = c('miRNA_name'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,6)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/promoter_DNAmeth.RData')



#only pre-miRNA HM
data_merged = merge(data_expression,data_HM_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,5)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/pre_miRNA_HM.RData')



#only promoter HM
data_merged = merge(data_expression,data_HM_promoter,by = c('miRNA_name'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,6)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/promoter_HM.RData')



#only pre-miRNA DNAmeth + HM
data_merged = merge(data_expression,data_DNA_meth_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged = merge(data_merged,data_HM_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,5)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/pre_miRNA_DNAmeth_HM.RData')



#only promoter DNAmeth + HM
data_merged = merge(data_expression,data_DNA_meth_promoter,by = c('miRNA_name'))
data_merged = merge(data_merged,data_HM_promoter,by = c('miRNA_name','start'))
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,6)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/promoter_DNAmeth_HM.RData')



#promoter (HM + DNAmeth) + pre-miRNA (HM + DNAmeth)
data_merged_miRNA = merge(data_expression,data_DNA_meth_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged_miRNA = merge(data_merged_miRNA,data_HM_miRNA,by = c('miRNA_ID','miRNA_name'))
data_merged_miRNA = unique(data_merged_miRNA)

data_merged_promoter = merge(data_expression,data_DNA_meth_promoter,by = c('miRNA_name'))
data_merged_promoter = merge(data_merged_promoter,data_HM_promoter,by = c('miRNA_name','start'))
data_merged_promoter = unique(data_merged_promoter)

data_merged = merge(data_merged_promoter,data_merged_miRNA,by=c('miRNA_ID','miRNA_name','log_RPM','expression'))           
data_merged = unique(data_merged)

target = data_merged$expression
data_set = get_data_set(data_merged,6)

save(data_set,target,file = '/project/histone_marks_miRNA/data/RData_classification/HELA/pre_miRNA_promoter.RData')

