plot_density <- function(feature,listOfFiles,colors){
  
  plotFile <- paste(feature,"_density.pdf",sep="")
  pdf(plotFile)
  
  data = read.table(listOfFiles[1],header=FALSE)
  data$V10[data$V10 == 0 & data$V11 == 'unknown'] = NA #mark missing values
  DNA_meth_data <- data$V10
  dens <- density(DNA_meth_data,na.rm=TRUE)
  
  plot(dens,type="l",col=colors[1],ylab="density",main=feature)
  abline(v=0.2,col=colors[2])
  abline(v=0.6,col=colors[2])
  
  legend("topright",legend="beta-Values",lty=c(1,1),col=colors)
  
  dev.off()
}


args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] #inputfile directory
feature <- args[2] #considered histonemodification for plot

setwd(dir)
listOfFiles <- list.files(path = ".", pattern = "*.gff", full.names=TRUE)
listOfFiles <- listOfFiles[grep(feature,listOfFiles)]

print(paste("Input directory: ",dir))
print(paste("feature: ",feature))
print("list of files in input directory: ")
print(listOfFiles)


library(RColorBrewer)
numberOfColors = length(listOfFiles)
colors <- brewer.pal(max(numberOfColors,3), "Set1")

plot_density(feature,listOfFiles,colors)






