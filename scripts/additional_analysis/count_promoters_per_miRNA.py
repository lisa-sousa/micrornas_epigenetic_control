import os.path
import argparse

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-p","--promoterList",help="Text file with predicted promoters",type=str,required=True)
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

def getPredictedPromoters(promoterFile):
  readhandle = open(promoterFile)
  miRNA_promoter_mapping = {}
  for line in readhandle:
    line = line.rstrip().split('\t')
    
    miRNA = str(line[8].split(';')[0].split(':')[1])
    
    if miRNA_promoter_mapping.has_key(miRNA):
      miRNA_promoter_mapping[miRNA].append(line[3])
    else:
      miRNA_promoter_mapping[miRNA] = [line[3]]
  return miRNA_promoter_mapping
		
#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterList = parameters.promoterList
proofFile(promoterList)

outputFile = parameters.outputFile

#=============================================================
#executive part
#=============================================================

miRNA_promoter_mapping = getPredictedPromoters(promoterList)

for miRNA, promoterList in miRNA_promoter_mapping.iteritems():
	miRNA_promoter_mapping[miRNA] = len(promoterList)

sorting = sorted(miRNA_promoter_mapping.items(), key=lambda x:x[1])

with open(outputFile,'w') as writehandle:
	for tuple in sorting:
		writehandle.write(tuple[0] + '\t' + str(tuple[1]) + '\n')


