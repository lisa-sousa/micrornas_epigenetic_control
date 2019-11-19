readData <- function(dir){

	dirs <- dir(path=dir,full.names=T)
	print("DIR")
	print(dir)
	print("SUBDIRS")
	print(dirs)

	listOfFiles <- list.files(path = dirs[1], pattern = "*.coverageOverRegion", full.names=TRUE)
	print("FILES")
	print(listOfFiles)	
	coverage = unlist(read.delim(listOfFiles[1],header=FALSE))
	
	HM = unlist(strsplit(dirs[1],"/"))
	HM = HM[length(HM)]

	coverage_dataFrame <- data.frame(HM = coverage)
	colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM

	feature = unlist(strsplit(listOfFiles[1],"/"))
	feature = unlist(strsplit(feature[length(feature)],"-"))
	feature = unlist(strsplit(feature[length(feature)],".coverage"))
	feature = feature[1]
	legend = c(feature)

	if(length(listOfFiles) > 1){
		for(i in c(2:length(listOfFiles))){
			coverage <- unlist(read.delim(listOfFiles[i],header=FALSE))
			coverage_dataFrame <- cbind(coverage_dataFrame,coverage)
			colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM

			feature = unlist(strsplit(listOfFiles[i],"/"))
			feature = unlist(strsplit(feature[length(feature)],"-"))
			feature = unlist(strsplit(feature[length(feature)],".coverage"))
			feature = feature[1]
			legend = append(legend,feature)
		}
	} 

	if(length(dirs) > 1){
		for(j in c(2:length(dirs))){
			listOfFiles <- list.files(path = dirs[j], pattern = "*.coverageOverRegion", full.names=TRUE)
			print("LIST IN LOOP")
			print(listOfFiles)
			coverage = unlist(read.delim(listOfFiles[1],header=FALSE))
			coverage_dataFrame <- cbind(coverage_dataFrame,coverage)
		

			HM = unlist(strsplit(dirs[j],"/"))
			HM = HM[length(HM)]
			colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM
		
			if(length(listOfFiles) > 1){
				for(k in c(2:length(listOfFiles))){
					coverage <- unlist(read.delim(listOfFiles[k],header=FALSE))
					coverage_dataFrame <- cbind(coverage_dataFrame,coverage)
					colnames(coverage_dataFrame)[ncol(coverage_dataFrame)] <- HM
				} 
			}
		}
	}
	return(list(coverage_dataFrame,legend))
}






plotCoverage <- function(dir,coverage_dataFrame,legendValues){
	library(RColorBrewer)
	setwd(dir)


	xMax <- round(nrow(coverage_dataFrame)/2)
	xMin <- (-1) * xMax
	if(nrow(coverage_dataFrame) %% 2 == 0){
		xMin <- xMin + 1
	}
	basePair <- seq(xMin,xMax)
	yMin = min(coverage_dataFrame)
	yMax = max(coverage_dataFrame)
		
		
	pdf("allHistoneModification_plot.pdf")
	listOfHM <- unique(colnames(coverage_dataFrame))
	numbOfCol <- length(colnames(coverage_dataFrame))/length(listOfHM)
	colors <- brewer.pal(numbOfCol, "Set1")
	colMarkerLine = brewer.pal(2, "Set2")

	for(i in c(1:length(listOfHM))){
		subDataFrame <- as.data.frame(coverage_dataFrame[,colnames(coverage_dataFrame)==listOfHM[i]])
		plot(basePair,subDataFrame[,1],type="l",col=colors[1],xlim=c(xMin,xMax),ylim=c(yMin,yMax),ylab="coverage",main=listOfHM[i])
		if(ncol(subDataFrame) > 1){
			for(i in c(2:ncol(subDataFrame))){
				lines(basePair,subDataFrame[,i],type="l",col=colors[i])
			} 
		}
		#abline(v = -300,col=colMarkerLine[1])
		#abline(v = 100,col=colMarkerLine[1])
		#abline(v = -35,col=colMarkerLine[2])
		#abline(v = 35,col=colMarkerLine[2])
		legend("topleft",legend=legendValues,lty=c(1,1),col=colors)
	}

	dev.off()

}






args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] #inputfile directory

list <- readData(dir)
coverage_dataFrame = list[[1]]
legendValues = list[[2]]

plotCoverage(dir,coverage_dataFrame,legendValues)


