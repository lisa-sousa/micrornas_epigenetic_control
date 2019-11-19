library(RColorBrewer)
color = brewer.pal(9,'Set1')

file = '/project/histone_marks_miRNA/histone_marks/data_analysis/DNA_methylation/HELA_pre_miRNA_1000bp_DNAmethylation/HAIB_HeLa-S3_fixed_Rep_1-hsa_pre_miRNA_1000bp.coverageOverRegion'
DNAcoverage = read.table(file,header=F)
x = 1:1000
lo = loess(as.integer(DNAcoverage$V1)~x,span=0.1)

pdf('/project/histone_marks_miRNA/histone_marks/data_analysis/DNA_methylation_profile/DNA_methylation_profile.pdf')
par(mfrow=c(1,1),mar=c(7,5,3,1),oma=c(1,1,1,1))
plot(predict(lo),type='l',xaxt = 'n',xlab='base pair',ylab='# of measured CpG sites',main='Distribution of measured CpG site across pre-miRNA region',col=color[2],lwd=2,ylim=c(0,10))
axis(1, at=seq(100,1000,200),labels=c(-400,-200,0,200,400), las=1)
dev.off()
