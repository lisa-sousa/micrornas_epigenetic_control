import os.path
import argparse
import numpy

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-p","--promoterFile",help="Text file with predicted promoters",type=str,required=True)
  args_parser.add_argument("-csp","--cellLineSpecificPromoterFile",help="Text file with cell line specific promoter predictions",type=str,required=True)
  args_parser.add_argument("-t","--tagCount",help="tag count threshold for cell line specific promoters that should not be used for filtering; default: 4",type=int, default=4)
  args_parser.add_argument("-o","--outputFile",help="outputfile",type=str,required=True)
  return args_parser.parse_args()
 
#============================================================

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

def filterPredictedPromoters(cellLineSpecificPromoterFile,tagCount):
	filteredPromotersFile = './filteredCellLineSpecificPromoters.gff'
	readhandle = open(cellLineSpecificPromoterFile)
	with open(filteredPromotersFile,'w') as writehandle:
		for line in readhandle:
		 	line = line.rstrip().split('\t')
		 	tag_count = float(line[8].split(';')[1].split(':')[1])
		 	if tag_count > tagCount:
			 	writehandle.write('\t'.join(line) + '\n')
	return filteredPromotersFile

#============================================================

def intersectDifferentPromFiles(promoterFile,filteredPromotersFile,outputFile):
	bedCmd =  'intersectBed -s -v -a ' + promoterFile + ' -b ' + filteredPromotersFile + ' > ' + outputFile
	os.system(bedCmd)

#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterFile = parameters.promoterFile
proofFile(promoterFile)

cellLineSpecificPromoterFile = parameters.cellLineSpecificPromoterFile
proofFile(cellLineSpecificPromoterFile)

tagCount = parameters.tagCount

outputFile = parameters.outputFile

#=============================================================
#executive part
#=============================================================

filteredPromotersFile = filterPredictedPromoters(cellLineSpecificPromoterFile,tagCount)
promoterOverlapFile = intersectDifferentPromFiles(promoterFile,filteredPromotersFile,outputFile)

