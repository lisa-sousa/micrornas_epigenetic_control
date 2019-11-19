import os.path
import argparse

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-p","--promoterList",help="Text file with predicted promoters",type=str,required=True)
  args_parser.add_argument("-m","--miRNAlist",help="Text file with non-expressed or expressed miRNAs",type=str,required=True)
  args_parser.add_argument("-f","--format",help="0: in gff format; 1: in PROmiRNA format",type=int,required=True)
  args_parser.add_argument("-i","--ID",help="position of miRNA name in group field (start counting from 0)",type=int,required=True)
  args_parser.add_argument("-o","--outputfileGFF",help="Outputfile with promoters from non-expressed/expressed miRNAs",type=str,required=True)
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

def getMiRNAlist(miRNAlist,ID):
	readhandle = open(miRNAlist)
	miRNAs = []
	for line in readhandle:
		if not line.startswith("#"):
			miRNA = line.rstrip().split(";")[ID].split("=")[1]
			miRNAs.append(miRNA)
	return(miRNAs)
			
#============================================================

def filterPromoters(miRNAs,promoterList,outputfileGFF):
	readhandle = open(promoterList)
	with open(outputfileGFF,"w") as writehandle:	
		for line in readhandle:		
			line = line.rstrip().split("\t")
			miRNA = line[0]
			if miRNA in miRNAs:
				output = "chr" + line[1] + "\t" + "PROmiRNA\tpromoter" + "\t" + line[2] + "\t" + line[3] + "\t.\t" + line[6] + "\t.\t" + "mirna:" + line[0] + ";" + line[4] + ";distance:" + line[5] + "\n"
				writehandle.write(output)

#============================================================

def filterPromotersinGFFformat(miRNAs,promoterList,outputfileGFF):
	readhandle = open(promoterList)
	with open(outputfileGFF,"w") as writehandle:	
		for line in readhandle:		
			line = line.rstrip().split("\t")
			miRNA = line[8].split(";")[0].split(":")[1]
			if miRNA in miRNAs:
				writehandle.write("\t".join(line) + "\n")
		
#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterList = parameters.promoterList
proofFile(promoterList)

miRNAlist = parameters.miRNAlist
proofFile(miRNAlist)

format = parameters.format
ID = parameters.ID

outputfileGFF = parameters.outputfileGFF

#=============================================================
#executive part
#=============================================================

miRNAs = getMiRNAlist(miRNAlist,ID)

if format == 0:
	filterPromotersinGFFformat(miRNAs,promoterList,outputfileGFF)
else:
	filterPromoters(miRNAs,promoterList,outputfileGFF)



