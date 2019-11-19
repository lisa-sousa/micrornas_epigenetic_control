import os.path
import argparse

#===============================================================
#read terminal parameters
#===============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-i","--biomartFile",help="biomart transcript file",type=str,required=True)
  args_parser.add_argument("-o","--outputFile",help="modified output file",type=str,required=True)
  return args_parser.parse_args()
 
#===============================================================

def proofFile(filename):
  if not (os.path.exists(filename)):
    exit("{0} does not exist!".format(filename))
    return False
  if not (os.path.isfile(filename)):
    exit("{0} is not a file!".format(filename))
    return False
  return True

#===============================================================
#Methods
#===============================================================

def processBiomartFile(biomartFile,outputFile,filterList):
	readhandle = open(biomartFile)
	fieldsDict = {}
	readhandle.readline()
	with open(outputFile,"w") as writehandle:
		for line in readhandle:
			line = line.rstrip().split(",")
			chrom = line[1]
			strand = line[2]
			start = int(line[3]) - 100
			end = int(line[4]) + 100
			biotype = line[6]
			if biotype in filterList:
				if len(chrom) <= 2:
					chrom = 'chr' + str(chrom)
					if int(strand) == 1:
						strand = '+'
					else:
						strand = '-'
				
					writehandle.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(chrom,start,end,strand,biotype))
			
#===============================================================    
#File Parameters   
#===============================================================     
      
parameters = args()

biomartFile = parameters.biomartFile
proofFile(biomartFile)

outputFile = parameters.outputFile

#===============================================================
#executive part
#===============================================================

#ensemble biotypes
long_noncoding = ['3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense', 'antisense_RNA', 'lincRNA', 'ncrna_host', 
					'non_coding', 'processed_transcript', 'retained_intron', 'sense_intronic', 'sense_overlapping']
short_noncoding = ['miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA', 'Mt_tRNA_pseudogene', 
					'ncRNA', 'ncRNA_pseudogene', 'rRNA', 'rRNA_pseudogene', 'scRNA', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 
					'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'tRNA', 'tRNA_pseudogene']
filterList = long_noncoding + short_noncoding

processBiomartFile(biomartFile,outputFile,filterList)