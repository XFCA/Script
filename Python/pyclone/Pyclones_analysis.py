#coding=utf-8
import os,commands,sys,argparse,HTSeq 
from argparse import RawTextHelpFormatter
import threading
import Queue
import time
import math
parser = argparse.ArgumentParser(description="Pyclone modules.",formatter_class=RawTextHelpFormatter)
parser.add_argument('-i','--input',help="A list file contains maf and cnv files, one sample per line, seperated by tab.\n   #sampleID\tmaf\tcnv\n",required=True)
parser.add_argument('-d','--depth', help="Min depth for varitions", default=20, type=int)
parser.add_argument('-o','--outdir',help="Output directory",required=True)
parser.add_argument('-c','--chr',type=int,default=0,help="For CNV, chr index number, default: 0")
parser.add_argument('-s','--start',type=int,default=1,help="For CNV, start index number, default: 1")
parser.add_argument('-e','--end',type=int,default=2,help="For CNV, end index number, default: 2")
parser.add_argument('-n','--cn',type=int,default=4,help="For CNV, copy number index number, default: 4")
parser.add_argument('-g','--genotype',type=int,default=7,help="For CNV, copy number\' genotype(AAB, AAA ...), default: 7")
parser.add_argument('-v','--cnvtype',help="CNV column type, default: cnv", choices=['cnv','ratio','log2ratio'], default='cnv')
parser.add_argument('-p','--purity',type=float,help="the cell purity")
parser.add_argument('-r','--region',help="Region for filter mutations, default: gene", choices=['all','gene','exon','coding','non-synonymous'],default='gene')
parser.add_argument('-t','--thread',type=int,help="Thread number",default=6)
parser.add_argument('-m','--testmode',help="TEST Model",action='store_true',default=False)
parser.add_argument('-l','--indel',help="Analysis INDEL variants or not",action='store_true',default=False)
parser.add_argument('-f',"--filter",help="filter mutation left number",default="400")
argv=parser.parse_args()

exitFlag = 0
class myThread ( threading.Thread ):
	def __init__ (self, threadID, name, q):
		threading.Thread.__init__ (self)
		self.threadID = threadID
		self.name = name
		self.q = q
	def run (self):
		print 'Starting ' + self.name
		doJob(self.name,self.q)
		print 'Exiting: ' + self.name

def doJob(threadName, q):
	while not exitFlag:
		queueLock.acquire()
		if not workQueue.empty():
			data = q.get()
			queueLock.release()
			os.system(data['script'])
			print '%s process %s' % (threadName, data['name'])
		else:
			queueLock.release()
		time.sleep(1)
#<bo>
'''
def filter_region(func, type='coding'):
	coding_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Silent','Splice_Site']
	nonsyn_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation','Splice_Site']
	exon_funcs    = ['3\'UTR','5\'UTR','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site']
	notgene_funcs = ['3\'Flank','5\'Flank','IGR']
	if type == 'coding':
		return func in coding_funcs
	if type == 'non-synonymous':
		return func in nonsyn_funcs
	if type == 'exon':
		return func in exon_funcs
	if type == 'gene':
		return func not in notgene_funcs
	return True
'''
#</bo>
def cal_cnv(cn,cnvtype):
	if cnvtype == 'cnv':
		return cn
	elif cnvtype == 'ratio':
		return cn*2
	elif cnvtype == 'log2ratio':
		return 2**(cn+1)

def get_purity(file,name) :
	purity_dict={}
	for each in open(file,'r') :
		array=each.strip().split('\t')
		purity_dict[array[0]]=array[2]
	if  purity_dict.get(name,'1')!='NA' : 
		print 'the purity of %s is %s '%(name,purity_dict.get(name,'1'))
#	print name,purity_dict.get(name,'1')
		return purity_dict.get(name,'1')
	else :
		print "the Absolute  result of %s is 'NA',defult purity=1"%(name)
		return '1'

def constructGA(cnvfile,chr,start,end,cn,cnvtype,genotype=None):#chr 0 start 1 end 2 cn 4 cnvtype cnv cnvfile *.somatic.CNV.txt
	cnvGA = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	i = 1
	for line in open(cnvfile):
		if i == 1:
			i += 1
			continue
		array = line.strip().split('\t')#行列表
		cn_tmp = cal_cnv(float(array[cn]),cnvtype=cnvtype)#cn_tmp拷贝数
		iv = HTSeq.GenomicInterval(array[chr].replace('chr',''),int(array[start]),int(array[end]))
		if genotype:
			cns = [array[genotype].count('A'),array[genotype].count('B')]
			cnvGA[iv] = [str(int(cn_tmp)),min(cns),max(cns)]
		else:
			cnvGA[iv] = [str(int(cn_tmp)),0,int(cn_tmp)]
	return cnvGA
		
#<bo/>def formatTSV(maf_file,tsv_file,cnvGA,sample,region='coding',depth=60):
def formatTSV(maf_file,tsv_file,cnvGA,sample,depth=60):#maf_file *.muTect.somatic.snv.maf,tsv_file ,cnvGA,sample,depth60
	tsv = open(tsv_file,'w')
	tsv.write('mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\tvariant_case\tvariant_freq\tgenotype\tRef\tAlt\tfuncs\n')
	i = 1
	for line in open(maf_file):
		if i == 1:
			i += 1
			continue
		array=line.strip().split('\t')
		## Only SNV 
		if not array[9] == 'SNP' and not argv.indel:
			continue
		## defined sample
		if sample and array[15] != sample:
			continue
		## region filter
#<bo/>		if not filter_region(array[8], type=region):
#<bo/>			continue
		## filter variations on chrY or other
		if array[4] in ['X','Y','chrX','chrY','MT','chrM']:
			continue
		## depth filter
		if int(array[32])+int(array[33]) < depth:
			continue
		id = ':'.join([array[0],array[4],array[5]])#id
		if array[0] == ".":
			id = ':'.join([array[4],array[5]])
		vaf = '%0.4f' % (float(array[33])/(int(array[32])+int(array[33])))#vaf
		normal_cn = '2'#normal_cn
		gp = HTSeq.GenomicPosition(array[4].replace('chr',''),int(array[5]))
		gt = 'A'*cnvGA[gp][1]+'B'*cnvGA[gp][2]#最小拷贝，最大拷贝
		if gt == "":
			continue
		tsv.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
			(id, array[32], array[33], normal_cn, cnvGA[gp][1], cnvGA[gp][2], array[15], vaf, gt, array[10], array[12], array[8]))
	tsv.close()


config ='''# Specifies working directory for analysis. All paths in the rest of the file are relative to this.
working_dir: %s

# Where the trace (output) from the PyClone MCMC analysis will be written.
trace_dir: %s

# Specifies which density will be used to model read counts. Most people will want pyclone_beta_binomial or pyclone_binomial
density: pyclone_binomial

# Number of iterations of the MCMC chain.#<bo 10000->50000/>
num_iters: 30000 

# Specifies parameters in Beta base measure for DP. Most people will want the values below.
base_measure_params:
  alpha: 1
  beta: 1

# Specifies initial values and prior parameters for the prior on the concentration (alpha) parameter in the DP. If the prior node is not set the concentration will not be estimated and the specified value will be used.
concentration:
  # Initial value if prior is set, or fixed value otherwise for concentration parameter.
  value: 1.0

  # Specifies the parameters in the Gamma prior over the concentration parameter.
  prior:
    shape: 1.0
    rate: 0.001

samples:
  # Unique sample ID
  %s:
    # Path where YAML formatted mutations file for the sample is placed.
    mutations_file: %s

    tumour_content:
      # The predicted tumour content for the sample. If you have no estimate set this to 1.0.
      value: %s

    # Expected sequencing error rate for sample
    error_rate: 0.001
'''
maf2xls_py = '''#!/usr/bin/python
import os,sys
if not len(sys.argv) == 4:
	print '\\tpython %s [cluster.txt] [tsv_file] [cluster.xls]\\n' % sys.argv[0]
	sys.exit(1)

pyclone = {}
i = 1
for line in open(sys.argv[1]):
	if i == 1:
		i += 1
		continue
	array = line.strip().split('\\t')
	pyclone[array[0]] = '\\t'.join(array[1:])

out = open(sys.argv[3],'w')
i = 1
for line in open(sys.argv[2]):
	if i == 1:
		i += 1
		out.write('Gene\\tChr\\tPos\\tRef\\tAlt\\tVariant_freq\\tCopy_number\\tCn_genotype\\tSample\\tFunction\\tCellular_frequencies\\tCf_std\\tCluster_id\\n')
		continue
	array = line.strip().split('\\t')
	if not array[0] in pyclone:
		continue
	out.write('%s\\t%s\\t%s\\t%s\\t%d\\t%s\\t%s\\t%s\\t%s\\n' % \\
		(array[0].replace(":","\\t"),array[9],array[10],array[7],len(array[8]),array[8],array[6],array[11],pyclone[array[0]]))
out.close()
'''
#<bo>
custom_mution_py='''#!/PROJ/HUMAN/share/Cancer/python/Python-2.7.11/bin/python
#coding=utf-8
###########################################
#                                         #
#       date: 2017/02/21                  #
#       author: Li Shaobo                 #
#                                         #
###########################################
import os,sys,random,argparse,commands
from collections import defaultdict
orgenefile="/PUBLIC/software/CANCER/Module/CancerGenome/Clone/or_gene"
parser = argparse.ArgumentParser(description="Pyclone custom_mution")
parser.add_argument('-c','--cluster',help="*cluster.txt,header=T",required=True)
parser.add_argument('-t','--tsv',help="*.tsv,header=T",required=True)
parser.add_argument('-o','--outfile',help="output_file name default:custom_mution_cluster.txt",default="custom_mution_cluster.txt")
parser.add_argument('-r','--region',help="Region for filter mutations, default: gene", choices=['all','gene','exon','coding','non-synonymous','random'],default='coding')
parser.add_argument('-n','--number',help="random numder", default="50")
argv=parser.parse_args()
ids=[]
orgene=[]

with open(argv.tsv,"r") as tsv1:
        tsv1.readline()
        for e in tsv1:
                id=e.strip().split('\\t')[0]
                if id not in ids:
                        ids.append(id)
        pass
pass

random_ids=random.sample(ids,min(int(argv.number),len(ids)))

try:
	with open(orgenefile,"r") as org :
			for e in org:
					if e.strip() not in orgene:
							orgene.append(e.strip())
			pass
except:
	print 'orgenefile not exeit.'

def filter_region(func,id, type='coding'):
        coding_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Silent','Splice_Site']
        nonsyn_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation','Splice_Site']
        exon_funcs    = ['3\\'UTR','5\\'UTR','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site']
        notgene_funcs = ['3\\'Flank','5\\'Flank','IGR']
        if type == 'coding':
                return func in coding_funcs
        if type == 'non-synonymous':
                return func in nonsyn_funcs
        if type == 'exon':
                return func in exon_funcs
        if type == 'gene':
                return func not in notgene_funcs
        if type == 'random' :
                return id in random_ids
        return True
pass
def filter_orgene(id):
        list=id.split(":")
        if len(list)==3:
                return (list[0] not in orgene)
        return True
pass

outf=open(argv.outfile,"w")
getids=[]
with open(argv.tsv,"r") as tsv:
        outf.write(tsv.readline().strip()+"\\n")
        for e in tsv:
                array=e.strip().split('\\t')
                if (filter_region(array[-1],array[0],argv.region) and filter_orgene(array[0])) :
                        #cmd="grep \\"^"+array[0]+"\\" "+argv.cluster
                        if array[0] not in getids :
                                getids.append(array[0])
                        #return_code,cluster=commands.getstatusoutput(cmd)
                pass
        pass
pass
finalids=random.sample(getids,min(int(argv.number),len(getids)))
for e in finalids :
        cmd="grep \\"^"+e+"\\" "+argv.cluster
        return_code,cluster=commands.getstatusoutput(cmd)
        outf.write(cluster.strip()+"\\n")
outf.close()
'''

pyclone_plot_py='''#!/PUBLIC/software/public/Python-2.7.6/bin/python
#coding=utf-8
###########################################
#                                         #
#       date: 2017/02/21                  #
#       author: Li Shaobo                 #
#                                         #
###########################################
import os,sys,random,argparse
import numpy as np
from collections import defaultdict
from pyclone.post_process.utils import load_cellular_frequencies_trace, load_cluster_labels_trace
from pyclone.post_process.plot.cellular_frequencies import CellularFrequencyPlot
from pyclone.post_process.plot.similarity_matrix import SimilarityMatrixPlot
from pyclone.post_process.plot.densities import PosteriorDensity
import brewer2mpl

bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
colors = bmap.mpl_colormap
parser = argparse.ArgumentParser(description="Pyclone plot")
parser.add_argument('--trace_file',help="trace_file<.bz2>",required=True)
parser.add_argument('--labels_file',help="labels_file<.bz2>",required=True)
parser.add_argument('--filter_id_file',help="filter ids of mutions,header=T",required=True)
parser.add_argument('--sample_name',help="sample_name",required=True)
parser.add_argument('--burnin',help="burnin",default='10000')
parser.add_argument('--thin',help="thin",default='1')
parser.add_argument('--outdir',help="outdir",required=True)
argv=parser.parse_args()
#out_cf_plot=argv.outdir+"/"+argv.sample_name+".cellular_frequencies"
#out_sm_plot=argv.outdir+"/"+argv.sample_name+".similarity_matrix"
out_cf_plot=os.path.join(argv.outdir,argv.sample_name+".cellular_frequencies.random50")
out_sm_plot=os.path.join(argv.outdir,argv.sample_name+".similarity_matrix.random50")

ids=[]
# get custom mution ids 
with open(argv.filter_id_file,"r") as filter_file:
        filter_file.readline()
        for e in filter_file:
                id=e.strip().split("\\t")[0]
                if id not in ids :
                        ids.append(id)
        pass
pass
#plot cellular_frequencies 
trace_cf_temp = load_cellular_frequencies_trace(argv.trace_file, int(argv.burnin), int(argv.thin))
trace_cf = defaultdict(list)
for id in trace_cf_temp:
        if id in ids:
            trace_cf[id] = trace_cf_temp[id]
pass
plotter = CellularFrequencyPlot(trace_cf, cmap=colors)
plotter.plot()
plotter.save(out_cf_plot+'.png')
plotter.save(out_cf_plot+'.pdf')
#plot similarity_matrix
trace_sm_temp = load_cluster_labels_trace(argv.labels_file, int(argv.burnin), int(argv.thin))
trace_sm = defaultdict(list)
for id in trace_sm_temp:
        if id in ids:
            trace_sm[id] = trace_sm_temp[id]
pass
plotter = SimilarityMatrixPlot(trace_sm)
plotter.plot()
plotter.save(out_sm_plot+'.png')
plotter.save(out_sm_plot+'.pdf')
'''
filter_tsv_py='''#!/PROJ/HUMAN/share/Cancer/python/Python-2.7.11/bin/python
#coding=utf-8
###########################################
#                                         #
#       date: 2017/02/21                  #
#       author: Li Shaobo                 #
#                                         #
###########################################
import os,sys,random,argparse,commands
parser = argparse.ArgumentParser(description="Pyclone custom_mution")
parser.add_argument('-t','--tsv',help="*.tsv,header=T",required=True)
parser.add_argument('-o','--outfile',help="output filename",required=True)
parser.add_argument('-n','--number',help="left numder", default="400")
argv=parser.parse_args()

def filter_region(func, type='coding'):
        coding_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Silent','Splice_Site']
        nonsyn_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation','Splice_Site']
        exon_funcs    = ['3\\'UTR','5\\'UTR','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site']
        notgene_funcs = ['3\\'Flank','5\\'Flank','IGR']
        if type == 'coding':
                return func in coding_funcs
        if type == 'non-synonymous':
                return func in nonsyn_funcs
        if type == 'exon':
                return func in exon_funcs
        if type == 'gene':
                return func not in notgene_funcs
        return True
pass
num={}.fromkeys(('all','coding', 'non-synonymous','exon','gene'), 0)
con={}.fromkeys(('all','coding', 'non-synonymous','exon','gene','random_non-synonymous'), "")
title=""
randomids=[]
with open(argv.tsv,"r") as tsv:
        title=tsv.readline().strip()
        for e in tsv :
                array=e.strip().split("\\t")
                if filter_region(array[-1],"coding"):
                        num['coding'] +=1
                        con['coding'] +=(e.strip()+"\\n")
                if filter_region(array[-1],"non-synonymous"):
                        num['non-synonymous'] +=1
                        con['non-synonymous'] +=(e.strip()+"\\n")
                        if array[0] not in randomids :
                                randomids.append(array[0])
                if filter_region(array[-1],"exon"):
                        num['exon'] +=1
                        con['exon'] +=(e.strip()+"\\n")
                if filter_region(array[-1],"gene"):
                        num['gene'] +=1
                        con['gene'] +=(e.strip()+"\\n")
                num['all'] +=1
                con['all'] +=(e.strip()+"\\n")
        pass
pass
stat="coding: %d,synonymous: %d,exon: %d,gene: %d,all: %d" %(num['coding'],num['non-synonymous'], num['exon'],num['gene'],num['all'])
noncon=con['non-synonymous'].split("\\n")
randomnon=random.sample(randomids,min(int(argv.number),len(randomids)))
for e in noncon :
        id=e.strip().split("\\t")[0]
        if id in randomnon:
                con['random_non-synonymous'] +=(e.strip()+"\\n")
with open(argv.outfile,"w") as outf:
        if (num['all'] <= int(argv.number)):
                outf.write(title+"\\t#filter by all"+stat+"\\n")
                outf.write(con['all'])
        elif (num['gene'] <= int(argv.number)):
                outf.write(title+"\\t#filter by gene"+stat+"\\n")
                outf.write(con['gene'])
        elif (num['exon'] <= int(argv.number)):
                outf.write(title+"\\t#filter by exon"+stat+"\\n")
                outf.write(con['exon'])
        elif (num['coding'] <=int(argv.number)) :
                outf.write(title+"\\t#filter by coding"+stat+"\\n")
                outf.write(con['coding'])
        elif (num['non-synonymous'] <= int(argv.number)) :
                outf.write(title+"\\t#filter by non-synonymous"+stat+"\\n")
                outf.write(con['non-synonymous'])
        else :
                outf.write(title+"\\t#filter by random non-synonymous by "+argv.number+" "+stat+"\\n")
                outf.write(con['random_non-synonymous'])
        pass
pass
'''
#</bo>
def pyclone_cmd(tsv_file,sample,resultdir):
	cfg = os.path.join(resultdir,sample+'.pyclone.config')
	yaml = os.path.join(resultdir,sample+'.yaml')
	sim = os.path.join(resultdir,sample+'.similarity_matrix')
	cluster = os.path.join(resultdir,sample+'.pyclone_cluster.txt')
	newcluster=os.path.join(resultdir,sample+'.pyclone_cluster_format.txt')
	out_file = os.path.join(resultdir,sample+'.pyclone.xls')
	open(os.path.join(resultdir,'table_format.py'),'w').write(maf2xls_py)
#<bo PyClone 0.13.0>
#       output custom_mution.py and pyclone_plot.py
	filter_tsv_py_file=os.path.join(resultdir,'filter_tsv.py')
	open(filter_tsv_py_file,"w").write(filter_tsv_py) #make script: filter_tsv.py
        custom_mution_file=os.path.join(resultdir,'custom_mution.py')
        open(custom_mution_file,"w").write(custom_mution_py) #make script: custom_mution.py
        pyclone_plot_file=os.path.join(resultdir,'pyclone_plot.py')
        open(pyclone_plot_file,"w").write(pyclone_plot_py) #make script: pyclone_plot.py
	trace_file=resultdir+"/"+sample+".cellular_frequencies.tsv.bz2"
	filter_tsv_file=tsv_file.strip(".tsv")+".filter.tsv"
	script = 'echo Pyclone analysis of %s\n' % (sample)
	script += 'cd %s \n' % (resultdir)
	script += 'mkdir -p %s/plot \n\n' % (resultdir) 
	script += 'python filter_tsv.py \\\n'
	script += '\t -t %s \\\n' %(tsv_file)
	script += '\t -o %s \\\n' %(filter_tsv_file)
	script += '\t -n %s \n\n' %(argv.filter)
	script += 'export PYTHONPATH=/PROJ/HUMAN/share/Cancer/Python-2.7.12/lib/python2.7:/PROJ/HUMAN/share/Cancer/Python-2.7.12/lib/python2.7/site-packages \n' 
	script += 'export PATH=/PROJ/HUMAN/share/Cancer/Python-2.7.12/bin:$PATH \n'
	script += 'PyClone build_mutations_file --prior total_copy_number --in_file %s --out_file %s && \\\n' % (filter_tsv_file,yaml)
	script += 'time PyClone run_analysis --seed 5 --config_file %s && \\\n' % (cfg)
	script += 'echo PyClone run_analysis %s done && \\\n' %(sample)
	script += 'PyClone build_table --config_file %s --out_file %s --table_type loci && \\\n' % (cfg,cluster)
	script += 'cut -f 1,3,4,5 %s > %s \n' %(cluster,newcluster)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.density.png --plot_type density && \\\n' % (cfg,resultdir,sample)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.parallel_coordinates.png --plot_type parallel_coordinates  && \\\n' % (cfg,resultdir,sample)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.scatter.png --plot_type scatter && \\\n' % (cfg,resultdir,sample)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.similarity_matrix.png --plot_type similarity_matrix  && \\\n' % (cfg,resultdir,sample)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.vaf_parallel_coordinates.png --plot_type vaf_parallel_coordinates  && \\\n' % (cfg,resultdir,sample)
	script += 'PyClone plot_loci --config_file %s --plot_file %s/plot/%s.vaf_scatter.png --plot_type vaf_scatter  \n' % (cfg,resultdir,sample)
	#use PyClone 0.12.7 to plot
	script += 'python %s/custom_mution.py \\\n' %(resultdir)
	script += '\t -c %s \\\n' %(cluster)
	script += '\t -t %s \\\n' %(filter_tsv_file)
	script += '\t -o %s/%s.filter_pyclone_cluster.txt \\\n' %(resultdir,sample)
	script += '\t -r %s \n\n' %(argv.region)
	script += 'export PATH=/PUBLIC/software/public/Python-2.7.6/bin:$PATH \n'
	script += 'export PYTHONPATH=/PUBLIC/software/public/Python-2.7.6/lib/python2.7/site-packages:/PUBLIC/software/public/Python-2.7.6/lib \n'
	script += 'python %s/pyclone_plot.py \\\n' %(resultdir)
	script += '\t --trace_file %s/trace/%s.cellular_prevalence.tsv.bz2 \\\n' % (resultdir,sample)
	script += '\t --labels_file %s/trace/labels.tsv.bz2 \\\n' % (resultdir)
	script += '\t --filter_id_file %s/%s.filter_pyclone_cluster.txt \\\n' % (resultdir,sample)
	script += '\t --sample_name %s \\\n' % (sample)
	script += '\t --outdir %s \n\n' %(resultdir) 
	script += 'python %s %s %s %s\n' % (os.path.join(resultdir,'table_format.py'), cluster, tsv_file, out_file)
	script += 'source ~/.bash_profile \n'
#</bo>
	return script

def create_dir(dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir '+dir)

outdir=argv.outdir
create_dir(outdir)

scripts = {}
samples = []

for line in open(argv.input):

	#创建样本目录和个文件名
	if line.startswith('#'):
		continue
	array = line.strip().split()#[sample,*.muTect.somatic.snv.maf,*.somatic.CNV.txt ]
	## tmpdir, tracedir, cfg_file, tsv_file
	tmpdir = os.path.join(outdir,array[0])
	create_dir(tmpdir)
	cfg_file = os.path.join(tmpdir,array[0]+'.pyclone.config')
	yaml = os.path.join(tmpdir,array[0]+'.yaml')
	tracedir = os.path.join(tmpdir,'trace')
	
	#获取纯度purity
	if argv.purity :
		purity=argv.purity
		print 'You force set %s purity = %s'%(array[0],purity)
	else :
		purity_file= '/'.join(argv.outdir.split('/')[0:-2]+['Absolute','Absolute.purity_plody.summary.xls'])
		if not os.path.exists(purity_file) :
			print 'I cannot find %s,%s defult purity=1'%(purity_file,array[0])
			purity=1
		else :	
			purity=get_purity(purity_file,array[0])
	
#	output custom_mution.py and pyclone_plot.py 创建py文件
	custom_mution_file=os.path.join(tmpdir,'custom_mution.py')
	open(custom_mution_file,"w").write(custom_mution_py)
	pyclone_plot_file=os.path.join(tmpdir,'pyclone_plot.py')
        open(pyclone_plot_file,"w").write(pyclone_plot_py)
	open(cfg_file,'w').write(config % (tmpdir, tracedir,array[0], yaml, purity))
	tsv_file = os.path.join(tmpdir,array[0]+'.tsv')
	
	cnvGA = constructGA(array[2],chr=argv.chr,start=argv.start,end=argv.end,cn=argv.cn,cnvtype=argv.cnvtype,genotype=argv.genotype)
	#得到一个HTSeq.GenomicArrayOfSets对象
#<bo>	formatTSV(array[1], tsv_file, cnvGA, sample=array[0], region=argv.region, depth=argv.depth)
	formatTSV(array[1], tsv_file, cnvGA, sample=array[0], depth=argv.depth)
#</bo>
	## script创建.sh
	script = pyclone_cmd(tsv_file, array[0], tmpdir)
	script_file = os.path.join(tmpdir,'pyclone_'+array[0]+'.sh')
	open(script_file,'w').write(script)
	scripts[array[0]] = script_file
	if array[0] not in samples:
		samples.append(array[0])

#####################################
## RUN Pyclone with threading mode ##
#####################################
if not argv.testmode:
	threadList = ['Thread-%d' % each for each in range(0,argv.thread)]
	myjobs = []
	for each in samples:
		myjobs.append({'name':each,'script':'sh %s'%scripts[each]})

	queueLock = threading.Lock()
	workQueue = Queue.Queue(0)
	threads = []
	threadID = 1
	## start threads
	for tName in threadList:
		thread = myThread(threadID, tName, workQueue)
		thread.start()
		threads.append(thread)
		threadID += 1
	## Queue and job
	queueLock.acquire()
	for eachjob in myjobs:
		workQueue.put(eachjob)
	queueLock.release()

	while not workQueue.empty():
		pass
	exitFlag = 1

	for t in threads:
		t.join()
