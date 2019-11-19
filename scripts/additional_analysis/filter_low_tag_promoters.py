import os.path
import argparse

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-p","--promoterFile",help="Text file with predicted promoters",type=str,required=True)
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

def filterPredictedPromoters(promoterFile,outputFile):
  readhandle = open(promoterFile)
  with open(outputFile,'w') as writehandle:
	  for line in readhandle:
		  line = line.rstrip().split('\t')
		  tag_count = float(line[8].split(';')[1].split(':')[1])
		  if tag_count > 4:
			  writehandle.write('\t'.join(line) + '\n')
		
#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterFile = parameters.promoterFile
proofFile(promoterFile)

outputFile = parameters.outputFile

#=============================================================
#executive part
#=============================================================

filterPredictedPromoters(promoterFile,outputFile)


