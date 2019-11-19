import os.path
import argparse

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-i","--inputAnnotationMIRNA",help="Text file with miRNA annotation",type=str,required=True)
  args_parser.add_argument("-x","--x",help="number of bp upstream of feature",type=int,required=True)
  args_parser.add_argument("-y","--y",help="number of bp downstream of feature",type=int,required=True)
  args_parser.add_argument("-o","--outputfileGFF",help="Outputfile with X bp region around miRNA",type=str,required=True)
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

def calculateXbpRegion(inputAnnotationMIRNA,x,y,outputfileGFF):
	readhandle = open(inputAnnotationMIRNA)
	with open(outputfileGFF,"w") as writehandle:
		for line in readhandle:
			if line.startswith("#"):
				writehandle.write(line)
			else:
				line = line.split()
				start = int(line[3])
				end = int(line[4])
				strand = line[6]
				middle = int(round((end-start)/2))
				if strand == "+":
					newStart = start - (x - middle)
					newEnd = end + (y - middle)
					#newStart = start - x
					#newEnd = start
				else:
					newStart = start - (y - middle)
					newEnd = end + (x - middle)
					#newStart = end
					#newEnd = end + x
				line[3] = str(max(1,newStart))
				line[4] = str(newEnd)
				writehandle.write("\t".join(line) + "\n")

#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

inputAnnotationMIRNA = parameters.inputAnnotationMIRNA
proofFile(inputAnnotationMIRNA)

x = parameters.x
y = parameters.y

outputfileGFF = parameters.outputfileGFF

#=============================================================
#executive part
#=============================================================

calculateXbpRegion(inputAnnotationMIRNA,x,y,outputfileGFF)
