import os.path
import argparse
import numpy

#===============================================================
#read terminal parameters
#===============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-s","--fasta_file",help="file with promoter sequence",type=str,required=True)
  args_parser.add_argument("-m","--meth_file",help="file with methylation coverage",type=str,required=True)
  args_parser.add_argument("-b","--bed_file",help="mapping file for varying start sites",type=str,default='.')
  args_parser.add_argument("-o","--output_dir",help="output directory for matrix files",type=str,default="./")
  return args_parser.parse_args()

#===============================================================
#Methods
#===============================================================

#proof if file exists
def proof_file(filename):
  if not (os.path.exists(filename)):
    exit("{0} does not exist!".format(filename))
    return False
  if not (os.path.isfile(filename)):
    exit("{0} is not a file!".format(filename))
    return False
  return True

#===============================================================


def get_mapping_start_sites(bed_file):

  readhandle = open(bed_file)
  mapping_start_sites = {}
  for line in readhandle:
    line = line.rstrip().split('\t')
    original_start = int(line[5])
    new_start = int(line[3])
    mapping_start_sites[new_start] = original_start
  return mapping_start_sites

#===============================================================

def get_number_of_CpG_sites(fasta_file,mapping,mapping_start_sites):
  
  readhandle = open(fasta_file)
  CpG_counts = {}
  for line in readhandle:
    if line.startswith('>'):
      header = line.rstrip().split(':')
      chrom = header[0].replace('>','')
      region = header[1].split('-')
      start = int(region[0]) + 1
      end = int(region[1].split('(')[0])
      strand = region[1].split('(')[1].replace(')','')
    else:
      count1 = line.count('CG')
      count2 = line.count('cg')
      if mapping:
        original_start = mapping_start_sites[start]
        CpG_counts[original_start] = count1 + count2
      else:
        CpG_counts[start] = count1 + count2
  return CpG_counts

#===============================================================

def get_450K_coverage(meth_file):
  readhandle = open(meth_file)
  meth_counts = {}
  for line in readhandle:
    line = line.rstrip().split('\t')
    start = int(line[3])
    end = int(line[4])
    count = line[13]
    meth_counts[start] = [int(count),line[0:9]]
  return meth_counts

#===============================================================

def get_coverage(CpG_counts,meth_counts):

  coverage_dict = {}
  not_covered = 0
  no_CpG_site = 0
  with open('coverage.txt','w') as writehandle:
    for prom_start, counts_CpG in CpG_counts.iteritems():
      if meth_counts.has_key(prom_start):
        counts_meth = meth_counts[prom_start][0]
        if counts_CpG == 0:
          coverage = 100
          no_CpG_site = no_CpG_site + 1
          #coverage_dict[prom_start] = coverage
        elif counts_meth == 0:
          not_covered = not_covered + 1
          coverage = 0
        else:
          coverage = float(counts_meth)/counts_CpG * 100
          coverage_dict[prom_start] = coverage
        writehandle.write('\t'.join(meth_counts[prom_start][1]) + '\t' + str(counts_meth) + '\t' + str(counts_CpG) + '\t' + str(coverage) + '\n')
      else:
        print prom_start
  
  coverage_list = coverage_dict.values()
  mean_cov = numpy.mean(coverage_list)
  print mean_cov, len(coverage_list)
  print not_covered, no_CpG_site, not_covered+no_CpG_site

#===============================================================    
#Parameters   
#===============================================================    

parameters = args()

fasta_file = parameters.fasta_file
proof_file(fasta_file)

meth_file = parameters.meth_file
proof_file(meth_file)

bed_file = parameters.bed_file
if bed_file == '.':
  mapping = False
else:
  mapping = True
  proof_file(bed_file)

output_dir = parameters.output_dir


#===============================================================
#executive part
#===============================================================

if not output_dir.endswith("/"):
  output_dir = output_dir + "/"

if mapping:
  mapping_start_sites = get_mapping_start_sites(bed_file)
else:
  mapping_start_sites = {}
print mapping_start_sites

CpG_counts = get_number_of_CpG_sites(fasta_file,mapping,mapping_start_sites)
#print CpG_counts

meth_counts = get_450K_coverage(meth_file)
#print meth_counts

get_coverage(CpG_counts,meth_counts)
