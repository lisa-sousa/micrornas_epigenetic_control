import os.path
import argparse
import math
import numpy as np


#===============================================================
#read terminal parameters
#===============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-b","--bamFile",help="File with coverageBed data for each bp",type=str,required=True)
  args_parser.add_argument("-g","--gffFile",help="Outputfile as matrix bp*region",type=str,required=True)
  args_parser.add_argument("-l","--lengthOfFeature",help="length of the features in GFF file, representing the regions",type=int,required=True)
  args_parser.add_argument("-m","--histoneModification",help="the considered histone modification",type=str,required=True)
  args_parser.add_argument("-d","--outputDir",help="output directory for bed analysis files (not required): default: ./",type=str,default="./")
  args_parser.add_argument("-s","--swapMinusStrand",help="Swap read position on minus strand (default:false)",action='store_false')
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

def setOutputFilenames(bamFile,gffFile,outputDir):
  bam = bamFile.split("/")
  bam = bam[len(bam)-1].split(".")
  bam = "_".join(bam[0:(len(bam)-1)])

  gff = gffFile.split("/")
  gff = gff[len(gff)-1].split(".")
  gff = "_".join(gff[0:(len(gff)-1)])

  bedAnalysisFileName = bam + "-" + gff 
  global bedAnalysisFile
  bedAnalysisFile = outputDir + bedAnalysisFileName + ".bed"
  
  global coverageFile
  coverageFile = outputDir + bedAnalysisFileName + ".coverageOverRegion"

#===============================================================

def runBedAnalysis(bamFile,gffFile,bedAnalysisFile):
  bedCmd = "coverageBed -d -abam " + bamFile + " -b " + gffFile + " > " + bedAnalysisFile
  print "Bed Tool command: ",bedCmd
  os.system(bedCmd)

#===============================================================

def runBedAnalysisParallel(bamFile,gffFile,bedAnalysisFile):
  bedCmd = "coverageBed -d -abam " + bamFile + " -b " + gffFile + " > " + bedAnalysisFile
  print "Bed Tool command: ",bedCmd
  print os.popen(bedCmd).read()
 
#===============================================================

def calculateCoverageOverRegionShift(bedAnalysisFile,lengthOfFeature,coverageFile,swapMinusStrand):
  readhandle = open(bedAnalysisFile)
  countsForFeature_plus = []
  countsForFeature_minus = []
  tmp = []
  maxima_plus = []
  maxima_minus = []
  strand = True #plus strand
  for line in readhandle:
    columns = line.rstrip("\n").split("\t")
    bp = int(columns[9])
    if bp == 1:
      if len(tmp) == lengthOfFeature:
        if(strand):
          maxima_plus.append(np.argmax(tmp))
          countsForFeature_plus.append(tmp)
        else:
          maxima_minus.append(np.argmax(tmp))
          countsForFeature_minus.append(tmp)
      tmp = []
    if bp <= lengthOfFeature:
      if columns[6] == "-" and swapMinusStrand:
        strand = False
        tmp.append(int(columns[10]))
      else:
        strand = True
        tmp.append(int(columns[10]))
  with open(coverageFile,"w") as writehandle:
    d = int(abs(np.mean(maxima_plus) - np.mean(maxima_minus)))

    countsForFeature_plus = np.array(countsForFeature_plus)
    tmp_plus = countsForFeature_plus[:,range(d, lengthOfFeature)]
    nrow = countsForFeature_plus.shape[0]
    z = np.zeros(d*nrow).reshape((nrow,d))
    countsForFeature_plus = np.concatenate((z,tmp_plus),axis=1)

    countsForFeature_minus = np.array(countsForFeature_minus)
    tmp_minus = countsForFeature_minus[:,range(0,lengthOfFeature-d)]
    nrow = countsForFeature_minus.shape[0]
    z = np.zeros(d*nrow).reshape((nrow,d))
    countsForFeature_minus = np.concatenate((tmp_minus,z),axis=1)

    countsForFeature = np.concatenate((countsForFeature_minus,countsForFeature_plus),axis=0)

    coverageOverRegion = []
    for i in range(countsForFeature.shape[1]):
      coverageOverRegion.append(str(np.mean(countsForFeature[:,i])))
    writehandle.write("\n".join(list(coverageOverRegion)) + "\n")

#===============================================================

def calculateCoverageOverRegionSwapp(bedAnalysisFile,lengthOfFeature,coverageFile,swapMinusStrand):
  readhandle = open(bedAnalysisFile)
  countsForFeature = []
  tmp = []
  for line in readhandle:
    columns = line.rstrip("\n").split("\t")
    bp = int(columns[9])
    if bp == 1:
      if len(tmp) == lengthOfFeature:
        countsForFeature.append(tmp)
      tmp = []
    if bp <= lengthOfFeature:
      if columns[6] == "-" and swapMinusStrand:
        tmp.insert(0,int(columns[10]))
      else:
        tmp.append(int(columns[10]))
  with open(coverageFile,"w") as writehandle:
    countsForFeature = np.array(countsForFeature)
    coverageOverRegion = []
    for i in range(countsForFeature.shape[1]):
      coverageOverRegion.append(str(np.mean(countsForFeature[:,i])))
    writehandle.write("\n".join(list(coverageOverRegion)) + "\n")

#===============================================================

def runPlotAnalysis(outputDir):  
  outputDir = outputDir[0:-1]
  rCmd = "R CMD BATCH --slave \"--args " + outputDir + "\" allProfilePlots.R reportAllProfilePlots.txt" 
  print rCmd
  os.system(rCmd)

#===============================================================    
#Main: Parameters and executive part
#===============================================================    

if __name__ == '__main__':
  parameters = args()

  bamFile = parameters.bamFile
  gffFile = parameters.gffFile
  outputDir = parameters.outputDir
  lengthOfFeature = parameters.lengthOfFeature
  swapMinusStrand = parameters.swapMinusStrand
  histoneModification = parameters.histoneModification
  
  proofFile(bamFile)
  proofFile(gffFile)

  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"

  analysisDir = outputDir + histoneModification + "/"
  if not (os.path.isdir(analysisDir)):
    dirCmd = "mkdir " + analysisDir
    print os.popen(dirCmd)
  print analysisDir, os.path.isdir(analysisDir)

  setOutputFilenames(bamFile,gffFile,analysisDir)



  print "\nHistone modification: ", histoneModification
  print "Minus strand swapped: ", swapMinusStrand
  print "Outpu directory: ", outputDir
  print "Path for analysis files: ", analysisDir
  print "\nBed Analysis File: ",bedAnalysisFile
  print "coverageFile: ",coverageFile


  #runBedAnalysis(bamFile,gffFile,bedAnalysisFile)
  print "\nBed Analysis done!"
  
  calculateCoverageOverRegionSwapp(bedAnalysisFile,lengthOfFeature,coverageFile,swapMinusStrand)
  print "\nCoverage calculated!"

  runPlotAnalysis(outputDir)


