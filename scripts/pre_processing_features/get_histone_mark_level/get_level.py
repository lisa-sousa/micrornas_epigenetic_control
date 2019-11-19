import os.path
import argparse
import pp

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
      jobs.append(job_server.submit(runBedAnalysis, (outputDir, line), (proofFile,plotDistribution), ("os.path",)))
  for job in jobs:
    result = job()

#===============================================================

def runBedAnalysis(outputDir,line):

  line = line.replace("\n","").split("\t")
  bamFile = line[0]
  gffFile = line[1]
  histoneModification = line[2]
  path = outputDir + histoneModification
      
  bam = bamFile.split("/")
  bam = bam[len(bam)-1].split(".")
  bam = "_".join(bam[0:(len(bam)-1)])
  gff = gffFile.split("/")
  gff = gff[len(gff)-1].split(".")
  gff = "_".join(gff[0:(len(gff)-1)])
  bedAnalysisFile = path + "/" + bam + "-" + gff + ".bed"
      
  if not (os.path.isdir(path)):
    dirCmd = "mkdir " + path
    print os.popen(dirCmd).read()
      
  bedCmd = "coverageBed -abam " + bamFile + " -b " + gffFile + " > " + bedAnalysisFile
  print "\nBed Tool command: ",bedCmd
  print os.popen(bedCmd).read()
  
#===============================================================    
#Parameters   
#===============================================================    

parameters = args()

configFile = parameters.configFile
proofFile(configFile)

outputDir = parameters.outputDir
ncpus = parameters.ncpus

#===============================================================
#executive part
#===============================================================

print "\nWorking directory: ", os.getcwd()

if not outputDir.endswith("/"):
  outputDir = outputDir + "/"

parallelization(configFile,outputDir,ncpus)

#rCmd = "R CMD BATCH --slave \"--args " + outputDir + "\"  plot_density.R reportAllDensities.txt" 
#print os.popen(rCmd).read()
#print("Combined Plot done!")
