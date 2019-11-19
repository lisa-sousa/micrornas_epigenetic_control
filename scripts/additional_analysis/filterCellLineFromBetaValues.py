import os.path
import argparse

#===============================================================
#read terminal parameters
#===============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-i","--betaValuesFile",help="file with beta values from DNA methylation",type=str,required=True)
  args_parser.add_argument("-c","--cellLine",help="cell line",type=str,required=True)
  args_parser.add_argument("-o","--outputFile",help="file with DNA methylations filtered by cell line",type=str,required=True)
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

def processBeteValuesFile(betaValuesFile,outputFile,cellLine):
	readhandle = open(betaValuesFile)
	header = readhandle.readline()
	#print header
	header = header.rstrip().split("\t")
	index = [i for i, elem in enumerate(header) if cellLine.lower() in elem.lower()][0]
	#print index
	with open(outputFile,"w") as writehandle:
		for line in readhandle:
			line = line.rstrip().split("\t")
			if len(line) > index:
				writehandle.write(line[0] + '\t' + line[index] + '\n')
			
#===============================================================    
#File Parameters   
#===============================================================     
      
parameters = args()

betaValuesFile = parameters.betaValuesFile
proofFile(betaValuesFile)

cellLine = parameters.cellLine


outputFile = parameters.outputFile

#===============================================================
#executive part
#===============================================================

processBeteValuesFile(betaValuesFile,outputFile,cellLine)
