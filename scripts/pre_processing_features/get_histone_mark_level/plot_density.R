cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

readData <- function(dir){

	dirs <- dir(path=dir,full.names=T)

	listOfFiles <- list.files(path = dirs[1], pattern = "*.bed", full.names=TRUE)
	print(listOfFiles)	
	coverage = read.table(listOfFiles[1],header=FALSE,sep="\t")
	
	HM = unlist(strsplit(dirs[1],"/"))
	HM = HM[length(HM)]
	print(HM)

	coverage_dataFrame <- data.frame(HM = log(coverage$V10 + 1))
	colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM

	feature = unlist(strsplit(listOfFiles[1],"/"))
	feature = unlist(strsplit(feature[length(feature)],"-"))
	feature = unlist(strsplit(feature[length(feature)],".bed"))
	feature = feature[1]
	legend = c(feature)

	if(length(listOfFiles) > 1){
		for(i in c(2:length(listOfFiles))){
			coverage <- read.table(listOfFiles[i],header=FALSE,sep="\t")
			coverage_dataFrame <- cbind.fill(coverage_dataFrame,data.frame(HM = log(coverage$V10 + 1)))
		  colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM

			feature = unlist(strsplit(listOfFiles[i],"/"))
			feature = unlist(strsplit(feature[length(feature)],"-"))
			feature = unlist(strsplit(feature[length(feature)],".bed"))
			feature = feature[1]
			legend = append(legend,feature)
		} 
	}
	


	for(j in c(2:length(dirs))){
		listOfFiles <- list.files(path = dirs[j], pattern = "*.bed", full.names=TRUE)
		print(listOfFiles)
		coverage = read.table(listOfFiles[1],header=FALSE,sep="\t")
		coverage_dataFrame <- cbind(coverage_dataFrame,log(coverage$V10 + 1))
		

		HM = unlist(strsplit(dirs[j],"/"))
		HM = HM[length(HM)]
		print(HM)
		colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM
		
		if(length(listOfFiles) > 1){
			for(k in c(2:length(listOfFiles))){
				coverage <- read.table(listOfFiles[k],header=FALSE,sep="\t")
				coverage_dataFrame <- cbind.fill(coverage_dataFrame,data.frame(HM = log(coverage$V10 + 1)))
				colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM
			} 
		}
	}
	return(list(coverage_dataFrame,legend))
}






plotCoverage <- function(dir,coverage_dataFrame,legendValues){
	library(RColorBrewer)
	setwd(dir)

	yMin <- c()
	yMax <- c()
	xMin <- c()
	xMax <- c()
	for(i in c(1:ncol(coverage_dataFrame))){
		dens <- density(coverage_dataFrame[,i],na.rm=T)
		yMin[i] <- min(dens$y)
		yMax[i] <- max(dens$y)
		xMin[i] <- min(dens$x)
		xMax[i] <- max(dens$x)
	
	} 
	yMin <- min(yMin)
	yMax <- max(yMax)
	xMin <- min(xMin)
	xMax <- max(xMax)
		
		
	pdf("allDensities_plot.pdf")
	listOfHM <- unique(colnames(coverage_dataFrame))
	numbOfCol <- length(legendValues)
	colors <- brewer.pal(numbOfCol, "Set1")
		
	for(i in c(1:length(listOfHM))){
		subDataFrame <- as.data.frame(coverage_dataFrame[,colnames(coverage_dataFrame)==listOfHM[i]])
		plot(density(subDataFrame[,1],na.rm=T),type="l",col=colors[1],xlim=c(xMin,xMax),ylim=c(yMin,yMax),ylab="density of log(data)",main=listOfHM[i])
		if(ncol(subDataFrame) > 1){
			for(i in c(2:ncol(subDataFrame))){
				lines(density(subDataFrame[,i],na.rm=T),type="l",col=colors[i])
			} 
		}
		legend("topleft",legend=legendValues,lty=c(1,1),col=colors)
	}

	dev.off()

}






args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] #inputfile directory

list <- readData(dir)

coverage_dataFrame = list[[1]]
legendValues = list[[2]]
print(head(coverage_dataFrame))

plotCoverage(dir,coverage_dataFrame,legendValues)


