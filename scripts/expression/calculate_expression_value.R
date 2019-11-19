##################################
#input files                     #
##################################

library(Cairo)
library(RColorBrewer)
library(mixtools)

#set PATH environment variable for bedTools: 
#PATH.save = Sys.getenv('PATH')
#Sys.setenv(PATH=paste('/home/lisasous/programs_libs/programs/bedtools2/bin',PATH.save,sep=':'))

output_dir = '/project/histone_marks_miRNA/pipeline/calculate_expression/plots/'
cell_line = 'hela'

bam_files = c('/project/histone_marks_miRNA/data/IMR90/encode_smallRNA_seq/bam_files/shortRnaSeq_Imr90_Rep1.bam',
              '/project/histone_marks_miRNA/data/IMR90/encode_smallRNA_seq/bam_files/shortRnaSeq_Imr90_Rep2.bam')

bam_files = c('/project/histone_marks_miRNA/data/HELA/encode_smallRNA_seq/bam_files/shortRnaSeq_Hela_Tap_Rep1.bam',
              '/project/histone_marks_miRNA/data/HELA/encode_smallRNA_seq/bam_files/shortRnaSeq_Hela_Tap_Rep2.bam')


totalNofReads <- c(162782772,193436115)#library sizes untreated IMR90
totalNofReads <- c(41301918,40798294)#library sizes TAP HELA

gff_files = c('/project/histone_marks_miRNA/data/raw_data/miRNA/miRNA_region/hsa_mature_miRNA.gff')

pre_miRNA_gff_files = c('/project/histone_marks_miRNA/data/raw_data/miRNA/miRNA_region/hsa_pre_miRNA.gff')

#set colors
if(length(gff_files) > 1){
  numberOfColors = length(gff_files)
}else{
  numberOfColors = length(bam_files)
}
colors = brewer.pal(max(numberOfColors,5), "Set1")


##################################
#bed file analysis               #
##################################
bed_output = data.frame(dir=character(), feature=character(), file=character())

for(gff_file in gff_files){
  gff = head(unlist(strsplit( tail(unlist(strsplit(gff_file,'/')),1), '\\.')),1)
  print(paste('Directory created:',gff))
  system(paste('mkdir',gff))
  for(bam_file in bam_files){
    bam = head(unlist(strsplit( tail(unlist(strsplit(bam_file,'/')),1), '\\.')),1)
    bed_analysis_file = paste('.', gff, paste( paste(bam ,gff , sep='-'), 'bed', sep='.'), sep='/')
    #uses bedTools v 2.25.0 -> -a and -b are swapped!!
    bedCmd = paste("coverageBed -split -s -a" ,gff_files , "-b", bam_file, ">", bed_analysis_file)
    bedCmd
    system(bedCmd)
    bed_output = rbind(bed_output,t(c(gff,bam,bed_analysis_file)))
  }
}
colnames(bed_output) = c('dir','feature','file')
bed_output <- data.frame(lapply(bed_output, as.character), stringsAsFactors=FALSE)


##################################
#functions                       #
##################################

read_raw_counts <- function(bed_output_sub){
  
  raw_counts = read.table(bed_output_sub[1,3],header=FALSE)[,c(1:9,12)]
  colnames(raw_counts) = c('chr','source','feature','start','end','score','strand','frame','group','readlength')
  
  for(i in 1:nrow(bed_output_sub)){
    data = read.table(bed_output_sub[i,3],header=FALSE)
    raw_counts = cbind(raw_counts,data$V10)
    colnames(raw_counts)[10+i] = paste('rep',i,sep='')
  }
  
  return(raw_counts)
}

extract_info <- function(group,indexOfID){
  group = unlist(strsplit(as.character(group),';'))[indexOfID]
  group = unlist(strsplit(group,'='))[2]
  return(group)
}

get_mapping_name_id <- function(pre_miRNA_gff_file){
  table = read.table(pre_miRNA_gff_file,sep='\t',header=F,comment.char='#')
  ID = unlist(lapply(table$V9,extract_info,indexOfID=1))
  name = unlist(lapply(table$V9,extract_info,indexOfID=3))
  mapping = data.frame(start = table$V4,end = table$V5,ID = ID, name = name)
  return(mapping)
}

remove_arm_description <- function(group,arm){
  if(arm != ''){
    group = unlist(strsplit(as.character(group),arm))
    group = paste(group[1],group[2], sep="")
  }
  group = unlist(strsplit(as.character(group),';'))
  group = paste(group[3],group[4],sep=';')
  return(group)
}

sum_counts <- function(readcount.x,readcount.y){
  readcount.x[is.na(readcount.x)] = 0
  readcount.y[is.na(readcount.y)] = 0
  readcount.all = readcount.x + readcount.y
  return(readcount.all)
}

get_readlength <- function(readlength.x,readlength.y){
  readlength.all = c()
  readlength.merged = cbind(readlength.x,readlength.y)
  for(i in c(1:nrow(readlength.merged))){
    if(is.na(readlength.merged[i,1])){
      readlength.all = append(readlength.all,readlength.merged[i,2])
    }
    else if(is.na(readlength.merged[i,2])){
      readlength.all = append(readlength.all,readlength.merged[i,1])
    }else{
      readlength.all = append(readlength.all,sum(readlength.merged[i,])/2)
    }
  }
  return(readlength.all)
}

merge_arm_counts <- function(raw_counts,mapping){
  
  five_arm = raw_counts[grep('-5p',raw_counts$group),]
  five_arm$group = unlist(lapply(five_arm$group,remove_arm_description,arm='-5p'))
  
  ID = unlist(lapply(five_arm$group,extract_info,indexOfID=2))
  five_arm$ID = ID
  five_arm_merged = merge(five_arm,mapping,by='ID')
  five_arm_merged = five_arm_merged[c(2:ncol(five_arm),1,(ncol(five_arm)+1):ncol(five_arm_merged))]
  colnames(five_arm_merged)[c(4,5,(ncol(five_arm)+1),(ncol(five_arm)+2))] = c('mm_start','mm_end','pm_start','pm_end')
  
  
  three_arm = raw_counts[grep('-3p',raw_counts$group),]
  three_arm$group = unlist(lapply(three_arm$group,remove_arm_description,arm='-3p'))
  
  ID = unlist(lapply(three_arm$group,extract_info,indexOfID=2))
  three_arm$ID = ID
  three_arm_merged = merge(three_arm,mapping,by='ID')
  three_arm_merged = three_arm_merged[c(2:ncol(three_arm),1,(ncol(three_arm)+1):ncol(three_arm_merged))]
  colnames(three_arm_merged)[c(4,5,(ncol(three_arm)+1),(ncol(three_arm)+2))] = c('mm_start','mm_end','pm_start','pm_end')
  
  
  no_arm = raw_counts[intersect(grep("-5p",raw_counts$group,invert=TRUE),grep("-3p",raw_counts$group,invert=TRUE)),]
  no_arm$group = unlist(lapply(no_arm$group,remove_arm_description,arm=''))
  ID = unlist(lapply(no_arm$group,extract_info,indexOfID=2))
  no_arm$ID = ID
  no_arm_merged = merge(no_arm,mapping,by='ID')
  no_arm_merged = no_arm_merged[c(2:ncol(no_arm),1,(ncol(no_arm)+1):ncol(no_arm_merged))]
  colnames(no_arm_merged)[c(4,5,(ncol(no_arm)+1),(ncol(no_arm)+2))] = c('mm_start','mm_end','pm_start','pm_end')
  no_arm_final = data.frame(chr = no_arm_merged$chr, source = no_arm_merged$source, feature = no_arm_merged$feature, 
                            start = no_arm_merged$pm_start, end = no_arm_merged$pm_end,score = no_arm_merged$score, 
                            strand = no_arm_merged$strand, frame = no_arm_merged$frame, 
                            group = paste('Name=',no_arm_merged$name,';ID=',no_arm_merged$ID,sep=''), 
                            readlength = no_arm_merged$readlength,no_arm_merged[11:(ncol(no_arm)-1)])
  
  
  fmt = merge(five_arm_merged,three_arm_merged,by = c('chr','source','feature','pm_start','pm_end', 'score','strand','frame','name','ID'),all=T)
  read_counts_fmt = sum_counts(fmt[15:ncol(five_arm_merged)],fmt[(ncol(five_arm_merged)+5):ncol(fmt)])
  
  raw_counts_merged = data.frame(chr = fmt$chr, source = fmt$source, feature = fmt$feature, start = fmt$pm_start, 
                                 end = fmt$pm_end,score = fmt$score, strand = fmt$strand, frame = fmt$frame, 
                                 group = paste('Name=',fmt$name,';ID=',fmt$ID,sep=''), 
                                 readlength = round(get_readlength(fmt$readlength.x,fmt$readlength.y)),read_counts_fmt)
  
  
  
  for(i in 11:ncol(raw_counts_merged)){
    colnames(raw_counts_merged)[i] = paste('rep',i-10,sep='')
  }
  
  raw_counts_merged = rbind(raw_counts_merged,no_arm_final)
  
  return(raw_counts_merged)
}

plot_density_df <- function(data_frame,main,isLog){
  
  if(!isLog){
    data_frame = log(data_frame + 1)
  }
  limits = data.frame(yMin=numeric(),yMax=numeric(),xMin=numeric(),yMin=numeric())
  for(i in 1:ncol(data_frame)){
    dens = density(data_frame[,i],bw=0.4, na.rm=TRUE)
    limits = rbind(limits,t(c(min(dens$y),max(dens$y),min(dens$x),max(dens$x))))
  } 
  colnames(limits) = c('yMin','yMax','xMin','xMax')
  
  dens = density(data_frame[,1],bw=0.4, na.rm=TRUE)
  plot(dens,type="l",col=colors[1],xlim=c(min(limits$xMin),max(limits$xMax)),ylim=c(min(limits$yMin),max(limits$yMax)),ylab="density",main=main)
  
  for(i in 2:ncol(data_frame)){
    dens = density(data_frame[,i],bw=0.4, na.rm=TRUE)
    lines(dens,type="l",col=colors[i])
  } 
  legend("topright",legend=colnames(data_frame),lty=c(1,1),col=colors)
  
}

##################################
#calculate expression value      #
##################################

dirs = unique(bed_output$dir)

#plot
estWidth <- 10
estHeight <- 10
CairoPDF(file = paste(output_dir,cell_line,'_expression.pdf',sep=''), width = estWidth, height = estHeight)

for(dir in dirs){
  bed_output_sub = bed_output[bed_output$dir == dir,]
  
  #load raw counts
  raw_counts = read_raw_counts(bed_output_sub)
  
  
  #plot raw counts
  main = "Density of log(raw counts) of short RNA Seq"
  plot_density_df(raw_counts[11:ncol(raw_counts)], main,FALSE)
  
  #get mapping for pre-miRNA from ID to name
  pre_miRNA_gff_file = pre_miRNA_gff_files[which(dirs == dir)[1]]
  mapping = get_mapping_name_id(pre_miRNA_gff_file)
  
  #megre counts from 5' arm with counts from 3' arm for miRNAs
  raw_counts = merge_arm_counts(raw_counts,mapping)
  
  #plot mixture model
  mix_counts = log(raw_counts[11:ncol(raw_counts)] + 1)
  mix_mod = npEM(mix_counts, 2, blockid=rep(1,ncol(mix_counts)), verb=F)
  plot(mix_mod, breaks = 30,xlab='log(raw counts)')
  
  lower_thr = 1
  upper_thr = 3
  abline(v=lower_thr,col=colors[4])
  abline(v=upper_thr,col=colors[5])
  
  #mean over repeats
  mean_counts = data.frame(mean_counts = round(rowMeans(raw_counts[11:ncol(raw_counts)])))
  raw_counts_mean = cbind(raw_counts[1:10],mean_counts)
  
  #assign expression levels -> 0 - non expressed, 1 - low expressed, 2 - fully expressed
  raw_counts_mean$expression = rep(1,nrow(raw_counts_mean))
  raw_counts_mean[log(raw_counts_mean$mean_counts+1) < lower_thr,ncol(raw_counts_mean)] = 0
  raw_counts_mean[log(raw_counts_mean$mean_counts+1) > upper_thr,ncol(raw_counts_mean)] = 2
  
  #plot raw counts grouped by expression level
  max_len = max(as.data.frame(table(raw_counts_mean$expression))$Freq)
  non_exp = raw_counts_mean$mean_counts[raw_counts_mean$expression == 0]
  low_exp = raw_counts_mean$mean_counts[raw_counts_mean$expression == 1]
  ful_exp = raw_counts_mean$mean_counts[raw_counts_mean$expression == 2]
  df = data.frame(non_exp = append(non_exp,rep(NA,max_len - length(non_exp))), low_exp = append(low_exp,rep(NA,max_len - length(low_exp))), 
                  ful_exp = append(ful_exp,rep(NA,max_len - length(ful_exp))))
  plot_density_df(df,main='Density of log(raw counts) of short RNA Seq',FALSE)
  
  #calculate RPM and take the mean
  raw_counts_RPM <- raw_counts[,11:ncol(raw_counts)]
  for(i in 1:ncol(raw_counts_RPM)){
    for(j in 1:nrow(raw_counts_RPM)){
      raw_counts_RPM[j,i] = (raw_counts_RPM[j,i] + 1) / totalNofReads[i] * 1000000
    }
  }
  raw_counts_mean$log_RPM = log(rowMeans(raw_counts_RPM[,1:2]))
  
  #plot RPM values by expression level
  non_exp = raw_counts_mean$log_RPM[raw_counts_mean$expression == 0]
  low_exp = raw_counts_mean$log_RPM[raw_counts_mean$expression == 1]
  ful_exp = raw_counts_mean$log_RPM[raw_counts_mean$expression == 2]
  df = data.frame(non_exp = append(non_exp,rep(NA,max_len - length(non_exp))), low_exp = append(low_exp,rep(NA,max_len - length(low_exp))), 
                  ful_exp = append(ful_exp,rep(NA,max_len - length(ful_exp))))
  plot_density_df(df,main='Density of log(RPM) of short RNA Seq', TRUE)
  
  #write data to file
  write.table(raw_counts_mean[,c(1:9,13,12)],file="./expression_values.txt",row.names=F, col.names=T, sep= "\t", quote=F)
  
  non_exp = raw_counts_mean[raw_counts_mean$expression == 0,]
  ful_exp = raw_counts_mean[raw_counts_mean$expression == 2,]
  write.table(non_exp[,c(1:9,13,12)],file="./expression_values_non_expressed_miRNAs.txt",row.names=F, col.names=F, sep= "\t", quote=F)
  write.table(ful_exp[,c(1:9,13,12)],file="./expression_values_expressed_miRNAs.txt",row.names=F, col.names=F, sep= "\t", quote=F)
}

dev.off()

