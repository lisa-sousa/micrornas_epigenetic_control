import os.path
import argparse
import numpy

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-i","--promoterFile",help="Text file with predicted promoters that should be merged if overlap",type=str,required=True)
  args_parser.add_argument("-d","--distance",help="Maximal distance between promoter that should be merged. Default: 0",type=int, default=0)
  args_parser.add_argument("-o","--outputFile",help="outputfile",type=str,required=True)
  return args_parser.parse_args()
 
#===============================================================
#Methods
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

def get_predicted_promoters(promoterFile):
  readhandle = open(promoterFile)
  miRNA_promoter_mapping = {}
  promoter_mapping = {}
  for line in readhandle:
    line = line.rstrip().split('\t')
    miRNA = str(line[8].split(';')[0].split(':')[1])
    prom_start = line[3]

    promoter_mapping[prom_start] = line[0:9]
    if miRNA in miRNA_promoter_mapping:
      miRNA_promoter_mapping[miRNA].append(prom_start)
    else:
      miRNA_promoter_mapping[miRNA] = [prom_start]
  return miRNA_promoter_mapping,promoter_mapping

#===============================================================

def calculateXbpRegion(start,end,x):
  middle = int(round((end-start)/2))
  new_start = start - (x - middle)
  new_end = end + (x - middle)

  return str(max(1,new_start)),str(new_end)

#============================================================

def merge_promoters(miRNA_promoter_mapping,promoter_mapping,distance,outputFile):
  bed_file = './merge_prom.gff'
  bed_results = './merge_prom.bed'

  with open(outputFile,'w') as writehandle_res:
    for miRNA, promoters_start in miRNA_promoter_mapping.items():

      with open(bed_file,'w') as writehandle_tmp:
        for prom_start in promoters_start:
          promoter = promoter_mapping[prom_start]
          prom_in_start = int(promoter[3])
          prom_in_end = int(promoter[4])
          promoter[3],promoter[4] = calculateXbpRegion(prom_in_start,prom_in_end,250)
          writehandle_tmp.write('\t'.join(promoter) + '\n')
          source = promoter[1]
          feature = promoter[2]

      bedCmd = "mergeBed -s -d " + str(distance) + " -i " + bed_file + " > " + bed_results
      os.system(bedCmd)


      miRNA_promoter_mapping[miRNA] = []
      readhandle = open(bed_results)
      for line in readhandle:
      	merged_prom = []      
      	line = line.rstrip().split('\t')
      	prom_start = int(line[1])
      	prom_end = int(line[2])
      	new_start,new_end = calculateXbpRegion(prom_start,prom_end,500)

      	new_prom = [line[0],source,feature,new_start,new_end,'.',line[3],'.','mirna:'+miRNA]
      	writehandle_res.write('\t'.join(new_prom) + '\n')

  rmCmd = 'rm -f ' + bed_file
  os.system(rmCmd)
  rmCmd = 'rm -f ' + bed_results
  os.system(rmCmd)

#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterFile = parameters.promoterFile
proofFile(promoterFile)

distance = parameters.distance

outputFile = parameters.outputFile

#=============================================================
#executive part
#=============================================================

miRNA_promoter_mapping,promoter_mapping = get_predicted_promoters(promoterFile)
merge_promoters(miRNA_promoter_mapping,promoter_mapping,distance,outputFile)
