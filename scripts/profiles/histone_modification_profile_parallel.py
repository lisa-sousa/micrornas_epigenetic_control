import os.path
import argparse
import pp
import histone_modification_profile_pipeline as hMPP

#===============================================================
#read terminal parameters
#===============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-c","--configFile",help="Config fileincluding information about bamFile, gffFile, lengthOfFeature, shift, stepSizeForPlot and histoneModification",type=str,required=True)
  args_parser.add_argument("-d","--outputDir",help="output directory for bed analysis files (not required): default: ./",type=str, default="./")
  args_parser.add_argument("-n","--ncpus",help="number of CPUs used for parallelization",type=int, default=0)
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
  bedAnalysisFile = outputDir + bedAnalysisFileName + ".bed"
  
  coverageFile = outputDir + bedAnalysisFileName + ".coverageOverRegion"

  return bedAnalysisFile,coverageFile

#===============================================================

def runHistoneModificationProfilePipeline(line, outputDir):
  line = line.replace("\n","").split("\t")
  
  bamFile = line[0]
  gffFile = line[1]
  lengthOfFeature = int(line[2])
  histoneModification = line[3]
  swapMinusStrand = bool(int(line[4]))

  proofFile(bamFile)
  proofFile(gffFile)

  path = outputDir + histoneModification
  if not (os.path.isdir(path)):
    dirCmd = "mkdir " + path
    print os.popen(dirCmd)

  outputDir = path + "/"

  bedAnalysisFile,coverageFile = setOutputFilenames(bamFile,gffFile,outputDir)

  print "Histone modification: ", histoneModification
  print "Length of feature: ", lengthOfFeature
  print "Minus strand swapped: ", swapMinusStrand

  print "\nBed Analysis File: ",bedAnalysisFile
  print "\nCoverage File: ",coverageFile

  histone_modification_profile_pipeline.runBedAnalysisParallel(bamFile,gffFile,bedAnalysisFile)
  print "Bed Analysis done!"

  histone_modification_profile_pipeline.calculateCoverageOverRegionSwapp(bedAnalysisFile,lengthOfFeature,coverageFile,swapMinusStrand)
  print "Coverage file created!"

  delCmd = "rm -f " + bedAnalysisFile
  print os.popen(delCmd).read()
  return True

#===============================================================

def parallelization(configFile,outputDir,ncpus = 0):
  if ncpus == 0:
    job_server = pp.Server(ppservers = ())
  else:
    job_server = pp.Server(ncpus, ppservers = ())
  jobs = []
  readhandle = open(configFile)
  print "Number of CPUs: ", job_server.get_ncpus()
  for line in readhandle:
    if not line.startswith("#"):
      jobs.append(job_server.submit(runHistoneModificationProfilePipeline, (line, outputDir), (proofFile,setOutputFilenames), ("os.path","histone_modification_profile_pipeline")))
  for job in jobs:
    result = job()
 
#===============================================================    
#Main: Parameters and executive part
#===============================================================  

if __name__ == '__main__':
  print "\nWorking directory: ", os.getcwd()

  parameters = args()

  configFile = parameters.configFile
  proofFile(configFile)

  outputDir = parameters.outputDir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"

  ncpus = parameters.ncpus

  parallelization(configFile,outputDir,ncpus)

  rCmd = "R CMD BATCH --slave \"--args " + outputDir + "\"  plot_profile.R report_plot_profile.txt"
  print os.popen(rCmd).read()
