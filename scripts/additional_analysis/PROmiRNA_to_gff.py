import os.path
import argparse

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-p","--promoterList",help="Text file with predicted promoters",type=str,required=True)
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

def formatPromoters(promoterList,outputfileGFF):
        readhandle = open(promoterList)
        with open(outputfileGFF,"w") as writehandle:    
                for line in readhandle:         
                        line = line.rstrip().split("\t")
                        output = "chr" + line[1] + "\t" + "PROmiRNA\tpromoter" + "\t" + line[2] + "\t" + line[3] + "\t.\t" + line[6] + "\t.\t" + "mirna:" + line[0] + ";" + line[4] + ";distance:" + line[5] + "\n"
                        writehandle.write(output)


#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

promoterList = parameters.promoterList
proofFile(promoterList)

outputfileGFF = parameters.outputfileGFF

#=============================================================
#executive part
#=============================================================


formatPromoters(promoterList,outputfileGFF)
