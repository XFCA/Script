#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import time
import argparse
import getpass
from argparse import RawTextHelpFormatter
from MaPping import Mapping
from Alignstat import AlnStat
from SmallMutationCalling import MutationCalling
from StructureVariationCalling import SV
from SomaticVariation import Somatic
from AdvanceCancer_v2 import AdvAnalysis #for cancer
from AdvanceDisease import Advance  #for disease; by sun
from QualityControl import QC #20180112


advance_cancer = {'1'   :  ' predispose_filter',
		'2'   :  ' signature_spectrum',
		'3'   :  ' drivergenes_filter',
		'4'   :  ' smg',
		'5'   :  ' gistic',
		'6'   :  ' fusion_gene',
		'7'   :  ' absolutes',
		'8.1' :  'pyclones',
		'8.2' :  'sciclones',
		'9'   :  ' phylip_evolution',
		'11'  :  'mrt_music',
		'12'  :  'oncodriveclust',
		'15'  :  'noncoding_filter',
                '17'  :  'mutationshow'}

advance_disease = {'6'  :  '   Filter',
                '6.1' :  ' filterDB',
                '7'   :  '   ModelFilter',
                '7.1' :  ' ModelD',
                '7.2' :  ' ModelR',
                '8'   :  '   Denovo',
                '8.1' :  ' Denovos_samtools',
                '8.2' :  ' Denovos_denovogear',
                '8.3' :  ' DenovoF',
                '8.4' :  ' DenovoR',
                '8.5' :  ' DenovoSV',
                '8.6' :  ' DenovoCNV',
                '9'   :  '   Linkage',
                '9.1' :  ' Merlinkage',
                '10'  :  '   Other',
                '10.1':  ' ROH'}  #by sun

parser = argparse.ArgumentParser(description="test mapping module",formatter_class=RawTextHelpFormatter)
parser.add_argument('--infile',help="a file include pacient sample ad library information. LaneID PatientID SampleID LibID NovoID Index Path Type Sex",required=True)
parser.add_argument('--seqstrag',help="WGS, WES_illu, WES_agi or TS",default='WGS',choices=['WGS','WES_agi','WES_illu','TS'])
parser.add_argument('--bwath',help="the CPU number bwa-mem will use",default='8')
parser.add_argument('--analydir',help="a path to store analysis results(defalut: ./)",default=os.path.abspath('./'))
parser.add_argument('--TR',help="a bed file for sequence region",default='/PUBLIC/source/HW/CANCER/HW_v1.0/Data/TR_bed/b37.chr25Region.bed')
parser.add_argument('--analy_array',help="which steps of ananlysis will do. Default (1)\n	"+"\n	".join(["1    quality_control",
			   "2    mapping"," 2.1 mapping_noGATK"," 2.2 mapping_withGATK_bestPrac (BQSR)"," 2.3 mapping_withGATK (BQSR+realignment)","3    snpindl_call"," 3.1 snpindl_call_samtools"," 3.2 snpindl_call_GATK",
			   "4    sv_call"," 4.1 sv_call_breakdancer"," 4.2 sv_call_crest"," 4.3 sv_call_delly"," 4.4 sv_call_lumpy","5    cnv_call"," 5.1 cnv_call_freec"," 5.2 cnv_call_cnvnator (only WGS)","5.3 cnv_call_freec_loh","6    somatic"," 6.1 somatic_snvindl_samtools"," 6.2 somatic_snvindl_varscan"," 6.3 somatic_snvindl_muTect&Strelka","7    somatic_sv"," 7.1 somatic_sv_breakdancer"," 7.2 somatic_sv_crest"," 7.3 somatic_sv_delly", " 7.4 somatic_sv_lumpy","8    somatic_cnv"," 8.1 somatic_cnv_freec"," 8.2 somatic_cnv_exomcnv (only WES)","8.3 somatic_cnv_varscan", "8.4 somatic_cnv_freec_loh"]),default='1')
parser.add_argument('--adv_array_cancer',help="which steps of cancer advance ananlysis will do. Only 8.2 is available for 8. Default: None\n    "+"\n   ".join(['%s   %s'%(e,advance_cancer[e]) for e in ['1','2','3','4','5','6','7','8.1','8.2','9','11','12','15','17']]),default=None)
parser.add_argument('--adv_array_disease',help="which steps of disease advance ananlysis will do. Default: None\n    "+"\n   ".join(['%s   %s'%(e,advance_disease[e]) for e in ['6','6.1','7','7.1','7.2','8','8.1','8.2','8.3','8.4','8.5','8.6','9','9.1','10','10.1']]),default=None)
parser.add_argument('--jobstat',help="the job.statues file")
parser.add_argument('--startpoint',help="give a start analysis point,if need!",default=None)
parser.add_argument('--newjob',help="the new job file",default=time.strftime('%m.%d',time.localtime(time.time()))+'.job')
parser.add_argument('--singlecell',help="single cell sequence data or not",default=None)
parser.add_argument('--opts',help="Options for raw2clean_QC(except -i,-a,-o), quoted. \nAvailable options: -N, -q, -L, -p, -r, -1, -2, -n, -d",default=None)
parser.add_argument('--rawopts',help="Options for raw2clean_QC to handle with rawdata(except -i,-a,-o), quoted. \nAvailable options: -N, -q, -L, -p, -r, -1, -2, -n, -d",default=None)
#parser.add_argument('--gz',help="If gzip clean data or not, apply it for only QC projects",action='store_true')
parser.add_argument('--fqcheck',help="If fqcheck for rawdata.",action='store_true')
parser.add_argument('--rmadapter',help="Wether remove adapter from raw fastqs.",action='store_true')
parser.add_argument('--cutadapter',help="Whether cut adapter from raw fastqs.",action='store_true')
parser.add_argument('--rmdup',help="Wether remove dup from raw fastqs.",action='store_true')
parser.add_argument('--duplevel',help="the duplevel you set ",default='15',type=int)
parser.add_argument('--pollution',help="If do pollution analysis or not.",action='store_true')
parser.add_argument('--TumorGermline',help="If do germline variation analysis for paired tumor samples,default yes.",action='store_true')
parser.add_argument('--nonflanking',help="Detect SNP/INDEL not on flanking region.",action='store_true')
parser.add_argument('--order',help="Samples order, seperated by comma(,)",default=None)
#parser.add_argument('--english',help="English Report",action='store_true')
parser.add_argument('--genome',help="Genome for analysis, default human_B37",default='human_B37',choices=['human_B37','human_hg19','human_hg38','human_B38','mm10','mm9'])
parser.add_argument('--germcallbychrom',help="gatk calling by chromosome. used in WGS projects to accelerate the germline SNPIndel calling step",action='store_true')
parser.add_argument('--socallbychrom',help="gatk calling by chromosome. used in WGS projects to accelerate the somatic SNP calling step",action='store_true')
parser.add_argument('--freec-step',help="Set the window and step size for freec, window = step * 2, default: 250 for WES/TS; 1000 for WGS", type=int, default=None)
parser.add_argument('--ffpe',help="If the sample is FFPE,Filter Cover when using Crest to detect SVs",action='store_true')#20180112
parser.add_argument('--moduledir',help="The directory of the module stored\n"
					"(default=/PUBLIC/source/HW/Disease/moduledir/)\n",default='/PUBLIC/source/HW/Disease/moduledir/')
parser.add_argument('--splitfq',help="For WGS sample sequenced on one single lane, split the clean.fq.gz file into several parts in order to speed up the mapping jobs. Note, all lanes listed in the sample_list file will be splitfq. Defaule=N", choices=['N', 'Y'], default='N')
parser.add_argument('--splitfq_num',help="the number of fq files to be splited. Use together with --splitfq. The value could be any positive int number. Default=8", default='8')
parser.add_argument('--PCRFree',help="Libraries were generated using PCR-Free", choices=['N', 'Y'], default='N')
parser.add_argument('--DNM',help="The file do denovo mutation need, 5 tab delimited fields: "
                                "#familyID   childID fatherID    motherID    ChildSEX(F or M)",default='')
parser.add_argument('--mlin',help="The config file do linkage analysis(merlin) need.(default=None)\n"  #by sun
                                "Format(tab separated): First line must startswith \"#FamilyID\" and include \"samples\"(\",\" separated),\"ped\"(file) and \"ws\"(file) at least,\"pop\",\"tbl\",\"lod\",\"prog\" can ber selectly included.",default='')  #by sun
parser.add_argument('--samp_info',help="The file of sample's information. 5 tab delimited fields:\n"
                                        "#FamilyID   SampleID    SEX Normal/Patient ProjectNo",default=None) #by sun
parser.add_argument('--ROH_soft',help="the software to detect ROH in disease advance analysis",choices=["PLINK","H3M2"],default="PLINK")
argv = vars(parser.parse_args())

# int memory
memory = {
  'ln' : '1M',
  'qc' : '5G',
  'md5' : '1G',
  'split_fq_read1' : '1G',
  'split_fq_read2' : '1G',
  'bamMD5':'1G',
  'bwa_mem' : '7G',
  'gzip_fq' : '1G',
  'pollution' : '1G',
  'samtools_sort' : '7G',
  'picard_mergebam' : '8G',
  'picard_rmdupBam' : '15G',
  'gatk_realign' : '12G',
  'gatk_bqsr' : '12G',
  'finalbam' : '1M',
  'mpileup' : '4G',
  'extractSoftClip_sub' : '1G',
  'catSoftClip' : '100M',
  'remove_rmdup' : '1M',
  'remove_sortbam' : '1M',
  'depthStat' : '1G',
  'covXdepth' : '2G',
  'flagstat' : '1G',
  'combine' : '1M',
  'samtoolsMpileup' : '2G',
  'samtoolsMpileup_sub' : '1G',
  'samtoolsMpileup_chr' : '300M',
  'cat_sub_all' : '1G',
  'filtersamtoolsCalling' : '2G',
  'seperateSnpIndel' : '1G',
  'gatk_calling' : '13G',
  'gatk_variantion_filter' : '8G',
  'annotatVcf' : '4G',
  'crestSV' : '8G',
  'crestSVann' : '1G',
  'lumpySV' : '20G',
  'lumpySVann' : '1G',
  'bam2cfg' : '1G',
  'breakdancerSV' : '1G',
  'breakdancerSvAnno' : '1G',
  'cnvnatorCNV' : '12G',
  'cnvnatorCNVann' : '1G',
  'freec_cnv' : '4G',
  'freec_cnv_loh' : '9G',
  'info_snp_circos': '1G',
  'info_indel_circos': '1G',
  'info_depth_circos': '5G',
  'info_crest_circos': '10M',
  'info_delly_circos': '10M',
  'info_breakdancer_circos': '10M',
  'info_lumpy_circos': '10M',
  'info_freec_circos': '10M',
  'conf_circos': '100M',
  'Circos': '100M',
  'varscan' : '5G',
  'somaticvarscanannovar' : '1G',
  'somaticsamtoolsMpileup' : '5G',
  'somaticsamtoolscalling' : '5G',
  'somaticsamtoolsannovar' : '1G',
  'somatic_muTect' : '12G',
  'somatic_Strelka' : '5G',
  'diff_sclip' : '1G',
  'crest_somaticSV' : '10G',
  'crest_somaticSVann' : '1G',
  'lumpy_somaticSV' : '5G',
  'lumpy_somaticSVann' : '1G',
  'sobam2cfg' : '1G',
  'sobreakdancerSV' : '1G',
  'sobreakdancerSvAnno' : '1G',
  'freec_somaticcnv' : '4G',
  'freec_somaticcnv_loh' : '9G',
  'freec_somaticCNV' : '1G',
  'freec_somaticCNVann' : '1G',
  'varscan_somaticCNV' : '5G',
  'varscan_somaticCNVann' : '1G',
  'exomeCNV_Ndepth' : '5G',
  'exomeCNV_Tdepth' : '5G',
  'exomeCNV_CNV' : '5G',
  'spectrum' : '500M',
  'merge_vcf' : '100M',
  'qc_report' : '100M',
  'mapping_report' : '100M',
  'result' : '500M',
  'primary_report' : '100M',
  'combine_maf' : '100M',
  'drivergenes_filter' : '100M',
  'predispose_filter' : '100M',
  'signature_spectrum' : '100M',
  'driverGene' : '2G',
  'smg' : '1G',
  'gistic' : '2G',
  'fusion_gene' : '100M',
  'absolutes' : '1G',
  'pyclones' : '1G',
  'pyclone_connect' : '2G',
  'sciclones' : '1G',
  'expands' : '2G',
  'phylip_evolution' : '1G',
  'circos' : '2G',
  'mrt_music':'500M',
  'oncodriveclust':'1G',
  'drug_target':'500M',
  'drug_resistance':'500M',
  'noncoding_filter':'500M',
  'bubbletree':'1G',
  'mutationshow' : '500M',
  'merged_vcf' : '10G',  #by sun <
  'samtools_DNM' : '5G',
  'denovogear_DNM':'5G',
  'catchr' : '500M',
  'samtools_DNM_snpanno' : '1G',
  'samtools_DNM_indelanno' : '1G',
  'denovogear_DNM_snpanno' : '1G',
  'denovogear_DNM_indelanno' : '1G',
  'mkdenovo': '200M',
  'Denovo':'1G',
  'ModelF': '200M',
  'ROH': '2G',
  'linkdatagen':'1G',
  'merlin':'1G',
  'merlin2R':'9G',
  'merlin2Excel':'200M',
  'DenovoSV':'2G',
  'DenovoCNV':'2G',  #by sun >
  'adv_result' : '500M',
  'adv_report' : '100M',
  'dellySV_del' : '1G',
  'dellySV_dup' : '1G',
  'dellySV_ins' : '1G',
  'dellySV_inv' : '1G',
  'dellySV_tra' : '1G',
  'dellySVanno' : '3G',
  'delly_somatic_sv_annot' : '3G',
  'delly_somatic_sv_del' : '2G',
  'delly_somatic_sv_dup' : '2G',
  'delly_somatic_sv_ins' : '2G',
  'delly_somatic_sv_inv' : '2G',
  'delly_somatic_sv_tra' : '2G',
  'gatk_calling_chr1to5' : '15G',
  'gatk_calling_chr6to12' : '15G',
  'gatk_calling_chr13toMT' : '15G',
  'vcf_indexing' : '300M',
  'somatic_muTect_chr1to5': '12G',
  'somatic_muTect_chr6to12': '12G',
  'somatic_muTect_chr13toMT': '12G',
  'somatic_muTect_annot':'3G'}

if argv['rawopts']:
	memory['ln'] = '500M'

threads = {
  'combine_maf' : '1',
  'drivergenes_filter' : '1',
  'predispose_filter' : '1',
  'signature_spectrum' : '1',
  'driverGene' : '1',
  'smg' : '6',
  'gistic' : '1',
  'fusion_gene' : '1',
  'absolutes' : '6',
  'pyclones' : '6',
  'sciclones' : '6',
  'expands' : '6',
  'phylip_evolution' : '1',
  'circos' : '6',
  'mrt_music':'1',
  'oncodriveclust':'1',
  'drug_target':'1',
  'drug_resistance':'1',
  'noncoding_filter':'1',
  'mutationshow' : '1',
  'bubbletree':'1'}
	

if os.path.isfile(argv['infile']):
	infile = argv['infile']
else:
	raise IOError('%s is not exsit' % argv['infile'])
for row in open(infile):
	array=row.strip().split("\t")
	if "-" in array[2]:
		raise IOError("invalid sample name: "+array[2])

seqstrag = argv['seqstrag']
th = argv['bwath']

if os.path.exists(argv['analydir']):
	analydir = argv['analydir']
else:
	raise IOError('%s is not exsit' % argv['analydir'])
tFreq = 0.1
nFreq = 0.05
jobstat = argv['jobstat']

if argv['ROH_soft']:   #by sun
        ROH_soft = argv['ROH_soft']
        

## freec step and window
freec_step = 250
if 'WGS' == seqstrag:
	freec_step = 1000
if argv['freec_step']:
	freec_step = argv['freec_step']
freec_window = freec_step * 2


if jobstat:
	if not jobstat.endswith('status'):
		print '\033[1;32;40m\nJob status file ends with \".status\", not the \".status.bak\" one.\nIf you DON\'T understand, please consult others.\n\033[0m\n'
		sys.exit(1)

startpoint = argv['startpoint']
startpoints = []
if startpoint:
	startpoints = [str(x) for x in argv['startpoint'].strip().split(',')]
singlecell = argv['singlecell']
#crest = argv['crest']
## qc
qc_opts=''
qc_opts2=""
opts=['-N','-q','-L','-p','-r','-1','-2','-n','-d']
clean = True
if argv['opts']:
	flag=False
	myopts = argv['opts'].strip().split()
	flags = [x for x in myopts if x.startswith('-')]
	myopts = myopts + ['']
	tmp = [(each,'') for each in flags if myopts[myopts.index(each)+1].startswith('-')] + [(each,myopts[myopts.index(each)+1]) for each in flags if not myopts[myopts.index(each)+1].startswith('-')]
	opts_dict = dict(tmp)
	for each in opts_dict:
		if each not in opts:
			flag = True
#		if each == '-n' and argv['gz']:
#			print '--gzip conflicts with -n options in raw2clean_QC\n'
#			sys.exit(1)
#		if each == '-n':
#			clean = False
	if flag:
		print '--opts is wrong!\nAvailable options for raw2clean_QC are:\n  ' + \
		', '.join(opts)
	qc_opts = '-option "%s"' % ' '.join(' '.join(each) for each in tmp)
	qc_opts2=' '.join(' '.join(each) for each in tmp)
	
#reAnaly = argv['reAnaly']
newjob = argv['newjob']

if 'WGS' not in seqstrag :
	memory['qc'] = '3G'
	memory['gatk_realign'] = '5G'
	memory['gatk_calling'] = '4G'
	memory['annotatVcf'] = '1G'
        memory['dnm'] = '500M'
        memory['merged_vcf'] = '2G' 

def safe_open(file_name,mode='r'):
	try:
		if not file_name.endswith('.gz'):
			return open(file_name,mode)
		else:
			import gzip
			return gzip.open(file_name,mode)
	except IOError:
		print file_name + ' do not exist!'

def parse_jobstatus (jobstatus):
	logdir = ''
	orders = []
	jobs = {}
	status_h = safe_open(jobstatus,'r')
	jobtxt=''
	jobtmp = {}
	for line in status_h:
		if line.startswith('log_dir'):
			logdir = line.strip().split()[1]
			continue
		elif line.startswith('order'):
			orders.append(line.strip())
			continue
		array = line.strip().split()
		jobtxt += line
		if len(array) == 2:
			jobtmp[array[0]] = array[1]
		if line.startswith('job_end'):
			if jobtmp['status'] == 'running':
				jobtmp['status'] = 'waiting'
				print 'JOB \"'+jobtmp['name']+'\" status: running, changed to: waiting.'
				jobtxt = jobtxt.replace('status running','status waiting')
			jobs[jobtmp['name']] = {'status':jobtmp['status'],'txt':jobtxt.strip()}
			jobtxt = ''
			jobtmp = {}
	status_h.close()
	return logdir,jobs,orders


def add_items(a,b):
	if type(b) == str:
		a.append(b)
	else:
		for mm in b:
			a.append(mm)

def create_dir (dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir %s' % dir)

def mv_0(x):
	if int(x) == x:
		return int(x)
	else:
		return x

## genome related
database_dir = '/PUBLIC/source/HW/Disease/database'
moduledir = argv['moduledir'].strip()
genome_version = argv['genome'].strip()
print 'The version of the reference genome is %s'% genome_version

ref=""
TR=""
if genome_version == 'human_B37':
    ref = "b37"
    if seqstrag=="WGS":
        TR="/PUBLIC/source/HW/CANCER/HW_v1.0/Data/TR_bed/b37.chr25Region.bed"
    elif seqstrag=="WES_agi":
        TR="/PUBLIC/source/HW/Disease/database/Cap/Exome_bed/Agilent/SureSelect.Human.All.Exon.V6.r2/agilent_region.bed"
elif genome_version == 'human_hg19':
    ref = "hg19"
    if seqstrag=="WGS":
        TR="/PUBLIC/source/HW/CANCER/HW_v1.0/Data/TR_bed/b37.chr25Region.bed"
    elif seqstrag=="WES_agi":
        TR="/PUBLIC/source/HW/Disease/database/Cap/Exome_bed/Agilent/SureSelect.Human.All.Exon.V6.r2/agilent_region.hg19.bed"
elif genome_version=="human_hg38":
    ref = "hg38"
    if seqstrag=="WGS":
        TR="%s/genome/human/GRCh38/GRCh38.chr25Region.bed" % database_dir
    elif seqstrag=="WES_agi":
        TR="/PUBLIC/source/HW/Disease/database/Cap/Exome_bed/Agilent/SureSelect.Human.All.Exon.V6.r2/agilent_region.liftover.hg38.bed"
elif genome_version=="human_B38":
    ref = "hg38"
    if seqstrag=="WGS":
        TR="%s/genome/human/GRCh38/GRCh38.chr25Region.bed" % database_dir
elif genome_version == 'mm9':
    ref = "mm9"
    database_dir='/PUBLIC/source/HW/CANCER/HW_v2.2/mm9'
    if seqstrag=="WGS":
        TR="%s/GRCm37.chrom_bed" % database_dir
    elif seqstrag=="WES_agi":
        TR="/PUBLIC/source/HW/Disease/database/genome/mouse/mm9/agilent_mouse_region_v1.bed"
elif genome_version == 'mm10':
    ref = "mm10"
    database_dir='/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/Mouse/mm10'
    if seqstrag=="WGS":
        TR="%s/GRCm38.chrom_bed" % database_dir
    elif seqstrag=="WES_agi":
        TR="/PUBLIC/source/HW/Disease/database/genome/mouse/mm10/agilent_region_mm10.bed"
else:
    raise IOError("wrong genome")
if argv["TR"].strip()!="":
    TR=argv["TR"].strip()
    if not os.path.isfile(TR):
        raise IOError('%s does not exsit' % argv['TR'])

genome_files = {
	'human_B37':{
#		'fasta':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/human_g1k_v37_decoy.fasta',
		'fasta':'%s/genome/human/b37_gatk/human_g1k_v37_decoy.fasta' % database_dir,
#		'fa2bit':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/human_g1k_v37_decoy.2bit', ## crest  -bit
		'fa2bit':'%s/genome/human/b37_gatk/human_g1k_v37_decoy.2bit' % database_dir,
#		'bychr':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/byChr', ## freec config
		'bychr':'%s/genome/human/b37_gatk/byChr' % database_dir,
		'annovardb':'/PUBLIC/database/HW/CANCER/HW_v2.0/Module/Annotation/humandb2', ## annovar database
#		'chrbed':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/b37.chr25Region.bed',  ## crest -regionList
		'chrbed':'%s/genome/human/human_b37/b37.chr25Region.bed' % database_dir,
#		'chrlen':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/freec/chr.24.length',  ## freec config
		'chrlen':'%s/genome/human/human_b37/freec/chr.24.length' % database_dir,
#		'Nblock':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/freec/all.NBlock.larger1000bp.bed', ## cnvnator
		'Nblock':'%s/genome/human/human_b37/freec/all.NBlock.larger1000bp.bed' % database_dir,
#		'indel1000':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf', ## bqsr; realn
		'indel1000':'%s/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf' % database_dir,
#		'dbsnp':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/dbSNP/138/VCF/00-All.vcf', ## bqsr
		'dbsnp':'%s/HEALTH/dbSNP/138/VCF/00-All.vcf' % database_dir,
		'lumpy-exclude':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_GRCh37/ceph18.b37.lumpy.exclude.2014-01-15.bed', ## lumpy
#		'mappability':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/freec/human_g1k_v37_decoy.mappability.100bp.out.mappability', ## freec config
		'mappability':'%s/genome/human/human_b37/freec/human_g1k_v37_decoy.mappability.100bp.out.mappability' % database_dir,
#		'dbsnp.txt':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/freec/b37_snp137.SingleDiNucl.1based.txt',  ## freec config
		'dbsnp.txt':'/PUBLIC/source/HW/Disease/database/genome/human/hg19/bwa_index/freec/b37_snp142.SingleDiNucl.1based.txt',
		'dbsnp.mutect':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_GRCh37/muTect/dbsnp_132_b37.leftAligned.vcf',  ## mutect
		'cosmic.mutect':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_GRCh37/muTect/b37_cosmic_v54_120711.vcf',  ## mutect
		'armsize':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_GRCh37/Varscan/ref-arm-sizes', ## varscan
		'annovarbuild':'hg19', ## annovar database build version
		'genomesplit_dir':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/GenomeSplit'
	},
	'human_hg19':{
#		'fasta':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/human_hg19_whole/ucsc.hg19.fasta',
		'fasta':'%s/genome/human/hg19/broad_Y_PAR_masked/ucsc.hg19.fasta' % database_dir,
#		'fa2bit':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/human_hg19_whole/ucsc.hg19.fasta.dict',
		'fa2bit':'%s/genome/human/hg19/broad_Y_PAR_masked/ucsc.hg19.2bit' % database_dir,
#		'bychr':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/byChr',
		'bychr':'%s/genome/human/hg19/byChr_Y_PAR_masked' % database_dir,
		'annovardb':'/PUBLIC/database/HW/CANCER/HW_v2.0/Module/Annotation/humandb2',
#		'chrbed':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/b37.chr25Region.bed',
		'chrbed':'%s/genome/human/hg19/bwa_index/hg19.chr25Region.bed' % database_dir,
#		'chrlen':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/chr.24.length',
		'chrlen':'%s/genome/human/hg19/bwa_index/freec/hg19.len' % database_dir,
#		'Nblock':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/all.NBlock.larger1000bp.bed',
		'Nblock':'%s/genome/human/hg19/bwa_index/freec/hg19.all.NBlock.larger1000bp.bed' % database_dir,
#		'indel1000':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/human_B37_hs37d5/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf',
		'indel1000':'%s/GATK_hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf' % database_dir,
#		'dbsnp':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/dbSNP/138/VCF/00-All.vcf',
		'dbsnp':'%s/GATK_hg19/dbsnp_138.hg19.vcf' % database_dir,
		'lumpy-exclude':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/ceph18.hg19.lumpy.exclude.2014-01-15.bed',
		'mappability':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/freec/human_g1k_v37_decoy.mappability.100bp.out.mappability',
#		'dbsnp.txt':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_B37/freec/hg19_snp137.SingleDiNucl.1based.txt',
		'dbsnp.txt':'/PUBLIC/source/HW/Disease/database/genome/human/hg19/bwa_index/freec/hg19_snp142.SingleDiNucl.1based.txt.gz',
		'dbsnp.mutect':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/muTect/dbsnp_132_hg19.leftAligned.vcf',
		'cosmic.mutect':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/muTect/hg19_cosmic_v54_120711.vcf',
		'armsize':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/human_hg19/ref-arm-sizes',
		'annovarbuild':'hg19',
		'genomesplit_dir':'%s/genome/human/hg19/GenomeSplit' % database_dir
	},
	'human_hg38':{
		'fasta':'%s/genome/human/hg38/bwa_index_whole/hg38.fa' % database_dir,
		'fa2bit':'%s/genome/human/hg38/bwa_index_whole/hg38.2bit' % database_dir,
		'bychr':'%s/genome/human/hg38/byChr/' % database_dir,
		'annovardb':'/PUBLIC/source/HW/Disease/database/genome/human/hg38/ANNOVAR/humandb/' ,
		'chrbed':'%s/genome/human/hg38/bwa_index/hg38.chr25Region.bed' % database_dir,
		'chrlen':'%s/genome/human/hg38/freec/hg38.chr.length' % database_dir,
		'Nblock':'%s/genome/human/hg38/freec/hg38.NBlock.larger5000bp.bed' % database_dir,
		'indel1000':'%s/genome/human/hg38/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf' % database_dir,
		'dbsnp':'%s/genome/human/hg38/GATK/dbsnp_144.hg38.addchr.vcf' % database_dir,
		'lumpy-exclude':'',
		'mappability':'/PUBLIC/software/CANCER/Database/Genome/human/hg38/human_hg38.mappability.100bp.out.mappability',#20180112
		'dbsnp.txt':'/PUBLIC/source/HW/Disease/database/genome/human/hg38/GATK/dbsnp147_00-All.hg38.vcf',
		'dbsnp.mutect':'/PUBLIC/source/HW/Disease/database/genome/human/hg38/GATK/dbsnp147_00-All.hg38.addchr.vcf',
#		'cosmic.mutect':'/PUBLIC/software/CANCER/Database/Genome/human/hg38/b38_cosmic_v54_120711.vcf',#20180112
		'armsize':'/PUBLIC/software/CANCER/Database/Genome/human/hg38/ref-arm-sizes',#20180112
		'annovarbuild':'hg38',
		'genomesplit_dir':'%s/genome/human/hg38/GenomeSplit' % database_dir
	},
	'human_B38':{#20180112
		'fasta':'/PUBLIC/software/CANCER/Database/Genome/human/b38/human_B38.fa',
		'fa2bit':'/PUBLIC/software/CANCER/Database/Genome/human/b38/human_B38.fa.2bit',
		'bychr':'/PUBLIC/software/CANCER/Database/Genome/human/b38/bychr',
		'annovardb':'/PUBLIC/software/CANCER/Database/ANNOVAR/humandb/humandb_B38',
		'chrbed':'/PUBLIC/software/CANCER/Database/Genome/human/b38/B38.chr25Region.bed',
		'chrlen':'/PUBLIC/software/CANCER/Database/Genome/human/b38/chr.24.length',
		'Nblock':'/PUBLIC/software/CANCER/Database/Genome/human/b38/B38.Nblock.larger1000bp.bed',
		'indel1000':'/PUBLIC/software/CANCER/Database/Genome/human/b38/Mills_and_1000G_gold_standard.indels.B38.vcf',
		'dbsnp':'/PUBLIC/software/CANCER/Database/Genome/human/b38/00-All.vcf',
		'lumpy-exclude':'',
		'mappability':'/PUBLIC/software/CANCER/Database/Genome/human/b38/human_B38.mappability.100bp.out.mappability',
		'dbsnp.txt':'/PUBLIC/software/CANCER/Database/Genome/human/b38/B38_snp137.SingleDiNucl.1based.txt',
		'dbsnp.sorted.txt':'/PUBLIC/software/CANCER/Database/Genome/human/b38/dbsnp_144_B38.leftAligned.vcf',
		'dbsnp.mutect':'/PUBLIC/software/CANCER/Database/Genome/human/b38/common_all_20160527.vcf',
		'cosmic.mutect':'/PUBLIC/software/CANCER/Database/Genome/human/b38/b38_cosmic_new.vcf',
		'armsize':'/PUBLIC/software/CANCER/Database/Genome/human/b38/ref-arm-sizes',
		'annovarbuild':'hg38'
	},
	## no indel1000, dbsnp, cosmic.mutect
	'mm9':{#20180112
		'fasta':os.path.join(database_dir,'mm9.fa'),
		'fa2bit':os.path.join(database_dir,'mm9.fa.2bit'),
		'bychr':os.path.join(database_dir,'byChr'),
		'annovardb':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/Mouse/annovar',
		'chrbed':os.path.join(database_dir,'mm9.chrom_bed'),
		'chrlen':os.path.join(database_dir,'mm9.chrom_length'),
		'Nblock':os.path.join(database_dir,'freec/mm9.Nblock.larger1000bp.bed'),
		'indel1000':'',
		'dbsnp':'',
		'lumpy-exclude':'',
		'mappability':os.path.join(database_dir,'freec/out50m2_mm9.gem'),
		'dbsnp.txt':os.path.join(database_dir,'freec/mm9_snp128'),
		'dbsnp.mutect':os.path.join(database_dir,'snp128_sorted.vcf'),
		'armsize':'',
		'annovarbuild':'mm9',
		'genomesplit_dir':os.path.join(database_dir,'GenomeSplit')
	},
	'mm10':{
		'fasta':os.path.join(database_dir,'mm10.fa'),
		'fa2bit':os.path.join(database_dir,'mm10.fa.2bit'),
		'bychr':os.path.join(database_dir,'byChr'),
		'annovardb':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Genome/Mouse/annovar',
		'chrbed':os.path.join(database_dir,'mm10.chrom_bed'),
		'chrlen':os.path.join(database_dir,'mm10.chrom_length'),
		'Nblock':os.path.join(database_dir,'freec/mm10.Nblock.larger1000bp.bed'),
		'indel1000':'',
		'dbsnp':'',
		'lumpy-exclude':'',
		'mappability':os.path.join(database_dir,'freec/GEM_mapp_GRCm38_68_mm10.gem'),
		'dbsnp.txt':os.path.join(database_dir,'freec/mm10_snp138.SingleDiNucl.1based.txt'),
		'dbsnp.mutect':os.path.join(database_dir,'mm10_snp138.vcf'),
		'armsize':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10/varscan/ref-arm-sizes',
		'annovarbuild':'mm10',
		'genomesplit_dir':os.path.join(database_dir,'GenomeSplit')
	}
}
genome_info = genome_files[genome_version]

softwares = {
	'java6':'/PUBLIC/software/public/System/jre1.6.0_33/bin/java',
	'java7':'/PUBLIC/software/public/System/jdk1.7.0_55/bin/java',
        'java8':'/PUBLIC/software/public/System/jre1.8.0_25/bin/java',
	'r3.2.1':'/PUBLIC/software/public/System/R-3.2.1/bin/Rscript',
	'r3.0.3':'/PUBLIC/software/public/R/v3.0.3/bin/Rscript',
	'gatk':'/PUBLIC/source/HW/Disease/software/GATK3.8/GenomeAnalysisTK-3.8-0-ge9d806836',
	'crest':'/PUBLIC/software/HW/CANCER/HW_v2.0/CREST',
	'breakdancer':'/PUBLIC/software/HW/CANCER/HW_v2.0/breakdancer-max1.4.4',
	'annovar':'/PUBLIC/database/HW/CANCER/HW_v2.0/Module/Annotation/Var_annotation_cancer_ANNOVAR2015Dec14_HW_v2.2.sh',
	'picard':'/PUBLIC/software/public/VarCall/picard/current',
	'bin':'/PUBLIC/software/HUMAN/bin',
	'exomecnv':'/PUBLIC/software/HW/CANCER/HW_v2.0/ExomeCNV',
	'mutect':'/PUBLIC/software/HW/CANCER/HW_v2.0/muTect_v1.1.4',
	'strelka':'/PUBLIC/software/HW/CANCER/HW_v2.0/strelka_v1.0.13',
	'freec':'/PUBLIC/software/HW/CANCER/HW_v2.0/control_freec_v7.0',
	'cnvnator':'/PUBLIC/software/HW/CANCER/HW_v2.0/CNVnator_v0.3',
	'lumpy':'/PUBLIC/software/HW/CANCER/HW_v2.0/speedseq_v0.0.3b/bin',
	'advbin':'/PUBLIC/source/HW/CANCER/HW_v2.2/Advance',
	'varscan':'/PUBLIC/software/HW/CANCER/HW_v2.0/varScan',
	'somatic_varScan':'/PUBLIC/database/HW/CANCER/HW_v2.0/Module/varScan/somatic_varScan.sh',
	'speedseq':'/PUBLIC/software/HW/CANCER/HW_v2.0/speedseq_v0.0.3b',
	'samtools.0.1.18':'/PUBLIC/software/HW/CANCER/HW_v2.0/samtools_v0.1.18',
	'samtools_v1.0':'/PUBLIC/source/HW/Disease/bin/',
	'sambamba_v0.5.9':'/PUBLIC/software/HW/CANCER/HW_v2.0/sambamba',
	'bcftools':'/PUBLIC/source/HW/Disease/bin/',
	'bwa':'/PUBLIC/source/HW/Disease/bin/',
	'igvtools':'/PUBLIC/software/HW/CANCER/src/igvtools/igvtools_2.3.67/igvtools',
	'DELLY':'/PUBLIC/software/HW/CANCER/HW_v2.0/DELLY_v0.7.3',
	'pigz':'/PUBLIC/software/HW/CANCER/HW_v2.0/pigz',
	'FFPE':'/PUBLIC/software/CANCER/Module/CancerGenome/FFPE',#20180112
	'QC_report':'/PUBLIC/source/HW/CANCER/HW_v2.2/QC_report',
	'Script':'/PUBLIC/source/HW/CANCER/HW_v2.2/Script'
}
#if ref == "hg38":
#	softwares['annovar'] = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/Annotation/Var_annotation_cancer_ANNOVAR2015Dec14_HW_hg38.sh'
## Other script or path used
generate_qc = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/QC_report/Human_reseq_qc.v1.1.pl'
pollution_script = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/QC_report/pollution/generate_contamination_QC_sh.v2.1.pl'
report_script ='/HWPROJ2/HW/wangjinhao/software/HW_v2.3/primary/cancer_reporter_v2.2.py'  #20190715
#report_script ='/PUBLIC/source/HW/CANCER/HW_v2.2/primary/cancer_reporter_v2.2.py'  #20180226
release_script ='/HWPROJ2/HW/wangjinhao/software/HW_v2.3/data_release_v2.2.1.py'  #20190715
#release_script ='/PUBLIC/source/HW/CANCER/HW_v2.2/data_release_v2.2.py'
disease_adv_release_script ='/PUBLIC/source/HW/Disease/software/Data_Release_Advance/disease_adv_release.py'
READMEdir = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/README/README_HW_v2.2'

advreport_script = '/HWPROJ2/HW/wangjinhao/software/HW_v2.3/Report/adv_cancer_report/cancer_reporter_adv_HW_v2.2.1.py' #by sun
#advreport_script = '/PUBLIC/source/HW/CANCER/HW_v2.2/Report/adv_cancer_report/cancer_reporter_adv_HW_v2.2.py' #by sun
advreport_disease = '/PUBLIC/source/HW/Disease/software/Disease_Report_Advance/disease_adv_report.py'
project_cost_script = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/project_cost_stat/project_cost_stat.py'
SV_pipeline_dir = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/SV_pipeline'
backup_script = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/project_info/Extract_samp_infov8.py'
backup_dir = {
	'WES':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Cancer/WES', 
	'WGS':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Cancer/WGS', 
	'NonSomatic':'/PUBLIC/database/HW/CANCER/HW_v2.0/Database/Cancer/NonSomatic'
}
data_relase_READMEFILE = '/PUBLIC/database/HW/CANCER/HW_v2.0/Module/Release/HW_v2.2/README.txt'  

times = time.strftime('%m.%d',time.localtime(time.time()))

mut_method = 'samtools'
sv_method = 'breakdancer'
cnv_method = 'freec'
sosnp_method = 'muTect'
soindel_method = 'Strelka'
sosv_method = 'breakdancer'
socnv_method = 'freec'


#give a value to analysis steps
symbol_code =  {1:'quality_control',
		2:'mapping',
		2.1:'mapping_noGATK',
		2.2:'mapping_withGATK_bestPrac',
                2.3:'mapping_withGATK',
		3:'snpindl_call',
		3.1:'snpindl_call_samtools',
		3.2:'snpindl_call_GATK',
		4:'sv_call',
		4.1:'sv_call_breakdancer',
		4.2:'sv_call_crest',
		4.3:'sv_call_delly',
		4.4:'sv_call_lumpy',
		5:'cnv_call',
		5.1:'cnv_call_freec',
		5.2:'cnv_call_cnvnator',
		5.3:'cnv_call_freec_loh',
		6:'somatic_snvindl',
		6.1:'somatic_snvindl_samtools',
		6.2:'somatic_snvindl_varscan',
		6.3:'somatic_snvindl_muTect',
		7:'somatic_sv',
		7.1:'somatic_sv_breakdancer',
		7.2:'somatic_sv_crest',
		7.3:'somatic_sv_delly',
		7.4:'somatic_sv_lumpy',	
		8:'somatic_cnv',
		8.1:'somatic_cnv_freec',
		8.2:'somatic_cnv_exomcnv',
		8.3:'somatic_cnv_varscan',
		8.4:'somatic_cnv_freec_loh'}

user = getpass.getuser()
#default_queues = set(['cancer.q','cancer1.q','disease.q','all.q'])
#user_queues = set([each.split('@')[0] for each in os.popen('qselect -U %s'%user).read().strip().split('\n')])
#useful_queues = user_queues & default_queues
user_queues= ','.join([ x.strip("\n") for x in os.popen("qselect -U $(whoami) | grep -E 'hw|hweu|hwas|hwus|hwsg|novo|all|joyce' | awk -F \"@\" '{print $1}' | sort -u").readlines() ])
queue_list='-q %s' % user_queues

#queue_list = ' '.join(['-q %s'%q for q in useful_queues])
#queue_list = '-q novo.q -q hw1.q -q hw2.q -q all.q'
#for i in [1, 2, 2.1, 2.2, 3, 3.1, 3.2, 4, 4.1, 4.2, 5, 5.1, 5.2, 6, 6.1, 6.2, 6.3, 7, 7.1, 7.2, 8, 8.1, 8.2]:
#	vars()[symbol_code[i]] = False
analy_array = [float(x) for x in argv['analy_array'].strip().split(',')]
analy_array = map(mv_0,analy_array)
analy_symbol = [symbol_code[each] for each in analy_array]
analysis = dict([(int(each),each) for each in analy_array])
includes = set(analysis.keys())
include_symbol = [symbol_code[each] for each in includes]
## for samtools call SNP/INDEL and crest softclip
chr_subs = [['1','2','3'],['4','5','6'],['7','8','9'],['10','11','12'],['13','14','15'],['16','17','18'],['19','20','21','22'],['X','Y','MT']]
chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
if ref == "hg19" or ref=="hg38":
	chr_subs = [['chr1','chr2','chr3'],['chr4','chr5','chr6'],['chr7','chr8','chr9'],['chr10','chr11','chr12'],['chr13','chr14','chr15'],['chr16','chr17','chr18'],['chr19','chr20','chr21','chr22'],['chrX','chrY','chrM']]
	chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
elif ref=="mm9" or ref=="mm10":
	chr_subs = [['chr1','chr2','chr3'],['chr4','chr5','chr6'],['chr7','chr8','chr9'],['chr10','chr11','chr12'],['chr13','chr14','chr15'],['chr16','chr17','chr18'],['chr19','chrX','chrY','chrM']]
	chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM']

# for disease_adv 
if genome_version == 'human_hg19' or genome_version == 'human_B38' or genome_version == 'human_hg38':
        chrlist = ['chr'+str(i) for i in  range(1,23)]+['chrX','chrY']
elif genome_version=="mm9" or genome_version=="mm10":
        chrlist = ['chr'+str(i) for i in  range(1,20)]+['chrX','chrY']
else:
        chrlist = [str(i) for i in range(1,23)]+['X','Y']

germcallbychrom = argv['germcallbychrom']
socallbychrom = argv['socallbychrom']
flank = not argv['nonflanking']
splitfq = argv['splitfq'].strip()
splitfq_num =argv['splitfq_num'].strip()
pcrfree=argv['PCRFree'].strip()
dnmf = argv['DNM'].strip()  #by sun
mlin = argv['mlin'].strip()  #by sun

###add DNM
#DNM = False
if dnmf:
    dnmf = os.path.abspath(dnmf)
    dnminfo = {}
    #dnmID = open(dnmf,'r')
    for line in safe_open(dnmf,'r'):
        if line.startswith('#'):continue
        if len(line.split('\t')) == 5:
        #DNM = True
            if not dnminfo.has_key(line.split('\t')[0].strip()):
                dnminfo[line.split('\t')[0].strip()] = {}

            dnminfo[line.split('\t')[0].strip()][line.split('\t')[1].strip()] = []
            #dnminfo[line.split('\t')[0].strip()]=[line.split('\t')[1].strip(),line.split('\t')[2].strip(),line.split('\t')[3].strip(),line.split('\t')[4].strip()]
            dnminfo[line.split('\t')[0].strip()][line.split('\t')[1].strip()]=[line.split('\t')[2].strip(),line.split('\t')[3].strip(),line.split('\t')[4].strip()]
        else:
            exit('Error:\tThere must be 5 ID for family,child ,father , mother and child_sex in every family.')



for i in analy_array:
	if i not in symbol_code:
		print str(i) + ' is not in analysis'
		sys.exit(1)
#	vars()[symbol_code[i]] = True
#	vars()[symbol_code[int(i)]] = True
somatic = False
if 'somatic_snvindl' in include_symbol or 'somatic_sv' in include_symbol or 'somatic_cnv' in include_symbol:
	somatic = True

sclip = False
#if seqstrag == 'WGS':
#	sclip = True
if 'sv_call_crest' in analy_symbol or 'somatic_sv_crest' in analy_symbol:
	sclip = True

mpileup = False
raw_mpileup = mpileup
# if 'somatic_cnv_varscan' in analy_symbol:
	# mpileup = True
if 'cnv_call_freec' in analy_symbol or 'somatic_cnv_freec' in analy_symbol or 'somatic_cnv_varscan' in analy_symbol or argv['ffpe']:
	mpileup = True
	
if ('TS' in seqstrag or 'WES' in seqstrag) and 'cnv_call' in include_symbol:
	print 'Germline CNV analysis is not applicable for WES or target sequencing!\n'
	sys.exit(1)
if 'WGS' in seqstrag and 'somatic_cnv' in include_symbol and 'somatic_cnv_exomcnv' in analy_symbol:
	print 'ExomeCNV only for WES, not WGS. Please choose other tools for Somatic CNV analysis.\n'
	sys.exit(1)

#give a relationship of analysis steps
relationship = {1:[],
		2:[1],
		2.1:[1],
		2.2:[1],
                2.3:[1],
		3:[1,2],
		3.1:[1,2],
		3.2:[1,2],
		4:[1,2],
		4.1:[1,2],
		4.2:[1,2],
		4.3:[1,2],
		4.4:[1,2],		
		5:[1,2],
		5.1:[1,2],
		5.2:[1,2],
		5.3:[1,2],
		6:[1,2],
		6.1:[1,2],
		6.2:[1,2],
		6.3:[1,2],
		7:[1,2],
		7.1:[1,2],
		7.2:[1,2],
		7.3:[1,2],
		7.4:[1,2],
		8:[1,2],
		8.1:[1,2],
		8.2:[1,2],
		8.3:[1,2],
		8.4:[1,2]}
## n opt in raw2clean_QC and 2
analys = [int(x) for x in analy_array]
if not clean and 2 in analys:
	print '-n option for raw2clean_QC is conflict with --analy_array.'
	sys.exit(1)

will_quit = False
for i in analy_array:
	s = set(relationship[i])
	if not s.issubset(set(analys)):
		print 'you need do %s before %f' %(', '.join(relationship[i]), i)
		will_quit = True
	else:pass
if will_quit:
	sys.exit(0)
if argv['TumorGermline'] :
	TumorGermline=True
else:
	TumorGermline=False
#rawdata and QC
print "QC ...\n"
qc_script = ''
if argv['fqcheck']:
	qc_opts += ' -fqcheck '
if singlecell:
	if argv['rawopts']:
		qc_script = '/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s -single %s -rawopts %s 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts, argv['rawopts'])
		assert not os.system('/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s -single %s -rawopts "%s" 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts, argv['rawopts']))
	else:
		qc_script = '/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s -single %s 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts)
		assert not os.system('/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s -single %s 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts))
else:
	if argv['rawopts']:
		qc_script = '/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s %s -rawopts "%s" 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts,argv['rawopts'])
		assert not os.system('/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s %s -rawopts "%s" 2>/dev/null ' % (generate_qc, analydir, infile, times, qc_opts, argv['rawopts']))
	else:
		qc_script = '/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s %s 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts)
		assert not os.system('/PUBLIC/software/public/System/Perl-5.18.2/bin/perl %s -pwd %s -list %s -date %s %s 2>/dev/null' % (generate_qc, analydir, infile, times, qc_opts))

open('%s/QC/generate_QCscript.sh' % analydir,'w').write(qc_script+'\n')

## parse qc_list and sample config
qclist = safe_open(os.path.join(analydir,'qc_list'),'r')
samplelist = safe_open(infile,'r')

list_in_sample = {}
list_in_qc = {}
sam2pid = {}
patientInfo = {}
sam2sex = {}
for line in qclist:
	if line.startswith('#'):
		continue
	array = line.strip().split('\t')
	if not array[2] in list_in_qc:
		list_in_qc[array[2]] = {}
	tmpid = array[-1].strip('"')
	list_in_qc[array[2]][tmpid] = '_'.join([array[4],array[0]])  ## pid,sampleID,libID,LaneID
	if array[2] not in sam2pid:
		sam2pid[array[2]] = []
	if array[1] not in sam2pid[array[2]]:
		sam2pid[array[2]].append(array[1])

	array[7] = array[7].strip()
	if array[1] not in patientInfo:
		patientInfo[array[1]] = {'N':'','T':[]}
	if array[7] == 'N':
		patientInfo[array[1]]['N'] = array[2]
	elif array[7] == 'T':
		if array[2] not in patientInfo[array[1]]['T']:
			patientInfo[array[1]]['T'].append(array[2])
	sam2sex[array[2]] = array[8]
Nsamples = []
for line in samplelist:
	if line.startswith('#') or line.strip() == '':
		continue
	array = line.strip().split('\t')
	tmpid = '_'.join([array[2],array[4],array[6],array[0]])
	sampleID = array[2]
	libid = list_in_qc[sampleID][tmpid]
	assert re.search(u'(\d+)',array[0])
	laneid = re.search(u'(\d+)',array[0]).group(1)
	tmp_fq1 = os.path.join(array[6],array[3],'%s_L%s_1.fq.gz'%(array[3],laneid))
	tmp_fq2 = os.path.join(array[6],array[3],'%s_L%s_2.fq.gz'%(array[3],laneid))
	if ',' in array[6] and array[6].endswith('.gz'):
		tmp_fq1,tmp_fq2 = array[6].split(',')
	if sampleID not in list_in_sample:
		list_in_sample[sampleID] = {}
	list_in_sample[sampleID][libid] = [tmp_fq1,tmp_fq2]
	if array[7] == 'N' and array[2] not in Nsamples:
		Nsamples.append(array[2])
	sam2sex[array[2]] = array[8]

### sampleorder
sample_order = list_in_qc.keys()
if argv['order']:
	if ',' in argv['order']:
		sample_order = argv['order'].strip(',').split(',')
	else:
		sample_order = [each.strip() for each in open(argv['order'])]

pairinfo = {}
for eachP in patientInfo:
	for each in patientInfo[eachP]['T']:
		pairinfo[each] = patientInfo[eachP]['N']

###samID for disease_adv #by sun
if argv['samp_info']:
    samp_info = argv['samp_info'].strip()
    gender = {}
    FamLi = []

    if os.path.isfile(samp_info):
        samp_info = os.path.abspath(samp_info)
        for i in open(samp_info):
            if i.startswith('#'):continue
            l = i.strip().split('\t')
            gender[l[1]] = l[2]
            FamLi.append(l[0])
        FamLi = list(set(FamLi))
    else:
        samp_info = 'Null'

    patient = {}
    QC_LIST = safe_open(os.path.join(analydir,'qc_list'),'r')
    samID=[]
    for line in QC_LIST:
        if line.startswith('#'):
            continue
        l = line.strip().split()
        if l[2] not in samID:
            samID.append(l[2])
        if samp_info == 'Null':
            gender[l[2]] = '/'
        #sample=lib_flowcell_lane
        if l[1] not in patient:
            patient[l[1]] = {}
            patient[l[1]][l[2]] = [l[4]+'_'+l[0]]  #sample=lib_flowcell_lane
        elif l[2] not in patient[l[1]]:
            patient[l[1]][l[2]] = [l[4]+'_'+l[0]]
        elif l[4]+'_'+l[0] not in patient[l[1]][l[2]]:
            patient[l[1]][l[2]].append(l[4]+'_'+l[0])
        else:pass
    QC_LIST.close()
    SamID=samID


### jobstatus file
logdir = ''
lastorders = []
lastjobs = {}
if jobstat:
	logdir,lastjobs,lastorders = parse_jobstatus(jobstat)
logdir = os.path.join(analydir,'log')
create_dir(logdir)
## dir
rawdir = os.path.join(analydir,'RawData')
qcdir = os.path.join(analydir,'QC')
mapdir = os.path.join(analydir,'Mapping')
alndir = os.path.join(analydir,'Alnstat')
mutdir = os.path.join(analydir,'Mutation')
svdir = os.path.join(analydir,'SV')
cnvdir = os.path.join(analydir,'SV')
somaticdir = os.path.join(analydir,'Somatic')
spectrumdir = os.path.join(analydir,'So_mut_charac')
reportdir = os.path.join(analydir,'Report')
resultdir = os.path.join(analydir,'Result')
adv_resultdir = os.path.join(analydir,'Result')
#advdir_cancer = os.path.join(analydir,'AdvanceAnalysis')
advdir_cancer = os.path.join(analydir,'AdvanceCancer') #by sun
if argv['adv_array_cancer']:
	os.system("mkdir -p "+advdir_cancer)
advdir_disease = os.path.join(analydir,'AdvanceDisease') #by sun
if argv['adv_array_disease']:
	os.system("mkdir -p "+advdir_disease)
adv_reportdir = os.path.join(reportdir,'advance')


new_jobs = {}
orders = []

add_jobs = {}  ## for job status
job_points = {} ## for startpoints

if set([1]).issubset(includes):
	create_dir(rawdir)
	create_dir(qcdir)
	job_points['ln'] = []
	job_points['qc'] = []
	job_points['md5'] = []
#	if argv['gz'] and clean:
	if clean:
		job_points['gzip_fq'] = []
	if argv['pollution']:
		job_points['pollution'] = []
	for eachsample in list_in_sample:
		aQC = QC(eachsample, rawdir, qcdir, list_in_sample[eachsample],softwares)
		create_dir(os.path.join(rawdir,eachsample))
		create_dir(os.path.join(qcdir,eachsample))
		for each in list_in_sample[eachsample]:
			## ln
			ln_cmd = aQC.ln(each,rawopts=argv['rawopts'],fqcheck=argv['fqcheck'],cutadapter=argv['cutadapter'],rmadapter=argv['rmadapter'],rmdup=argv['rmdup'],duplevel=argv['duplevel'])
			new_jobs['_'.join(['ln',eachsample,each])] = \
				{'name' : '_'.join(['ln',eachsample,each]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['ln'],
				'cmd' : 'sh %s' % os.path.join(rawdir,eachsample+'/'+'_'.join(['ln',eachsample,each])+'.sh')}
			safe_open(os.path.join(rawdir,eachsample+'/'+'_'.join(['ln',eachsample,each])+'.sh'),'w').write(ln_cmd)
			job_points['ln'].append('_'.join(['ln',eachsample,each]))
			## QC
			qc_cmd,order1 = aQC.qc(each,opts=qc_opts2,fqcheck=argv['fqcheck'],rmadapter=argv['rmadapter'],singlecell=argv['singlecell'],cutadapter=argv['cutadapter'])
			new_jobs['_'.join(['qc',eachsample,each])] = \
				{'name' : '_'.join(['qc',eachsample,each]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['qc'],
				'cmd' : 'sh %s' % os.path.join(qcdir,eachsample+'/'+'_'.join(['qc',eachsample,each])+'.sh')}
			#order1 = 'order %s before %s' %('_'.join(['ln',eachsample,each]), '_'.join(['qc',eachsample,each]))
			safe_open(os.path.join(qcdir,eachsample+'/'+'_'.join(['qc',eachsample,each])+'.sh'),'w').write(qc_cmd)
			#order1 = 'order %s after %s' %('_'.join(['qc',eachsample,each]),'_'.join(['ln',eachsample,each]))
			add_items(orders,order1)
			job_points['qc'].append('_'.join(['qc',eachsample,each]))
			## md5sum raw
			md5_cmd,order1 = aQC.md5(each)
			new_jobs['_'.join(['md5',eachsample,each])] = \
				{'name' : '_'.join(['md5',eachsample,each]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['md5'],
				'cmd' : 'sh %s' % os.path.join(rawdir,eachsample+'/'+'_'.join(['md5',eachsample,each])+'.sh')}
			#order1 = 'order %s before %s' %('_'.join(['qc',eachsample,each]), '_'.join(['md5',eachsample,each]))
			#order1 = 'order %s after %s' %('_'.join(['md5',eachsample,each]),'_'.join(['ln',eachsample,each]))
			add_items(orders,order1)
			safe_open(os.path.join(rawdir,eachsample+'/'+'_'.join(['md5',eachsample,each])+'.sh'),'w').write(md5_cmd)
			job_points['md5'].append('_'.join(['md5',eachsample,each]))
			## gzip clean
#			if argv['gz'] and clean:
			if clean:
				gz_cmd,order1 = aQC.gzip_fq(each)
				new_jobs['_'.join(['gzip_fq',eachsample,each])] = \
					{'name' : '_'.join(['gzip_fq',eachsample,each]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['gzip_fq'],
					'cmd' : 'sh %s' % os.path.join(qcdir,eachsample+'/'+'_'.join(['gzip_fq',eachsample,each])+'.sh')}
				#order1 = 'order %s after %s' %('_'.join(['gzip_fq',eachsample,each]),'_'.join(['qc',eachsample,each]))
				add_items(orders,order1)
				safe_open(os.path.join(qcdir,eachsample+'/'+'_'.join(['gzip_fq',eachsample,each])+'.sh'),'w').write(gz_cmd)
				job_points['gzip_fq'].append('_'.join(['gzip_fq',eachsample,each]))
		## pollution
		if argv['pollution']:
			pol_cmd,order1 = aQC.pollution(each)
			new_jobs['_'.join(['pollution',eachsample])] = \
				{'name' : '_'.join(['pollution',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['pollution'],
				'cmd' : 'sh %s' % os.path.join(qcdir,eachsample+'/'+'_'.join(['pollution',eachsample])+'.sh')}
			#order1 = 'order %s after %s' % ('_'.join(['pollution',eachsample]),'_'.join(['ln',eachsample,each]))
			safe_open(os.path.join(qcdir,eachsample+'/'+'_'.join(['pollution',eachsample])+'.sh'),'w').write(pol_cmd)
			add_items(orders,order1)
			## scripts
			# step0 = 'cd %s\nif [ ! -f "%s.pollution.xls" ];then' % (os.path.join(qcdir,eachsample),eachsample)
			# step1 = '\\\n\t'.join(['/PUBLIC/software/public/System/Perl-5.18.2/bin/perl '+pollution_script,
				# '%s/%s_%s_1.fq.gz' %(os.path.join(rawdir,eachsample),eachsample,each),
				# '500 1e-1 random.record',
				# '>%s/pollution.sh' % os.path.join(qcdir,eachsample)])
			# step2 = 'sh %s/pollution.sh' % os.path.join(qcdir,eachsample)
			# step3 = 'mv summarized.parsed.extracted.%s_%s_1.fasta.xls %s.pollution.xls\nfi\n' % (eachsample,each,eachsample)
			# safe_open(os.path.join(qcdir,eachsample+'/'+'_'.join(['pollution',eachsample])+'.sh'),'w').write(step0+'\n'+' && \\\n\t'.join([step1,step2,step3]))
			safe_open(os.path.join(qcdir,eachsample+'/'+'_'.join(['pollution',eachsample])+'.sh'),'w').write(pol_cmd)
			job_points['pollution'].append('_'.join(['pollution',eachsample]))
		
	## jobstatus
	i =1
	if lastjobs:
		for eachsample in list_in_qc:
			for each in list_in_qc[eachsample]:
				lnjob = '_'.join(['ln',eachsample,list_in_qc[eachsample][each]])
				qcjob = '_'.join(['qc',eachsample,list_in_qc[eachsample][each]])
				md5job = '_'.join(['md5',eachsample,list_in_qc[eachsample][each]])
				
				if lnjob not in new_jobs and lnjob in lastjobs:
					add_jobs[lnjob] = lastjobs[lnjob]['txt']
				if qcjob not in new_jobs and qcjob in lastjobs:
					add_jobs[qcjob] = lastjobs[qcjob]['txt']
				if md5job not in new_jobs and md5job in lastjobs:
					add_jobs[md5job] = lastjobs[md5job]['txt']
#				if argv['gz'] and clean:
				if clean:
					gzipjob = '_'.join(['gzip_fq',eachsample,each])
					if gzipjob not in new_jobs and gzipjob in lastjobs:
						add_jobs[gzipjob] = lastjobs[gzipjob]['txt']
			if argv['pollution']:
				pollutejob = '_'.join(['pollution',eachsample])
## wether run qc or not 
qc_flag = True
if startpoint and ('qc' not in startpoints or 'ln' not in startpoints):
	qc_flag = False

if set([1,2]).issubset(includes):
	print "Mapping ...\n"
	create_dir(mapdir)
	create_dir(alndir)

	job_points['split_fq_read1'] = []
	job_points['split_fq_read2'] = []
	job_points['bwa_mem'] = []
#	if not argv['gz']:
#		job_points['gzip_fq'] = []
	job_points['samtools_sort'] = []
	job_points['picard_mergebam'] = []
	job_points['picard_rmdupBam'] = []
	if analysis[2] == 2.2:
		job_points['gatk_bqsr'] = []
        if analysis[2] == 2.3:
                job_points['gatk_realign'] = []
                job_points['gatk_bqsr'] = []
	job_points['finalbam'] = []
	job_points['bamMD5'] = []
	job_points['depthStat'] = []
	job_points['flagstat'] = []
	job_points['combine'] = []
	job_points['remove_rmdup'] = []
	job_points['remove_sortbam'] = []
	if sclip:
		job_points['extractSoftClip_sub'] = []
		job_points['catSoftClip'] = []
	if mpileup:
		job_points['mpileup'] = []

	for eachsample in list_in_sample:
		mymapdir = os.path.join(mapdir,eachsample)
		create_dir(mymapdir)
		merge_bams = []
		## mapping
		#__init__(self,pacientID,sampleID,qcDir,alignDir,th,lib,genome_info,softwares)
		#print list_in_sample[eachsample]
		aMapping = Mapping(sam2pid[eachsample],eachsample,qcdir,mymapdir,th,list_in_sample[eachsample],genome_files[argv['genome']],softwares)
		for eachlib in aMapping.lib:
			if splitfq == 'Y':
				splitfq1_cmd,order1,splitfq2_cmd,order2 = aMapping.split_fq(eachlib,splitfq_num)
				add_items(orders,order1)
				add_items(orders,order2)
				new_jobs['_'.join(['split_fq_read1',eachsample,eachlib])] = \
					{'name' : '_'.join(['split_fq_read1',eachsample,eachlib]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,th),
					'memory' : memory['split_fq_read1'],
					'cmd' : 'sh '+os.path.join(qcdir+'/'+eachsample,'_'.join(['split_fq_read1',eachsample,eachlib])+'.sh')}
				new_jobs['_'.join(['split_fq_read2',eachsample,eachlib])] = \
					{'name' : '_'.join(['split_fq_read2',eachsample,eachlib]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,th),
					'memory' : memory['split_fq_read2'],
					'cmd' : 'sh '+os.path.join(qcdir+'/'+eachsample,'_'.join(['split_fq_read2',eachsample,eachlib])+'.sh')}
				safe_open(os.path.join(qcdir+'/'+eachsample,'_'.join(['split_fq_read1',eachsample,eachlib])+'.sh'),'w').write(splitfq1_cmd)
				safe_open(os.path.join(qcdir+'/'+eachsample,'_'.join(['split_fq_read2',eachsample,eachlib])+'.sh'),'w').write(splitfq2_cmd)
				job_points['split_fq_read1'].append('_'.join(['split_fq_read1',eachsample,eachlib]))
				job_points['split_fq_read2'].append('_'.join(['split_fq_read2',eachsample,eachlib]))
				for fen in range(1,int(splitfq_num)+1):
					n1 = eachlib+'-'+str(fen)
					merge_bams.append(n1)
					## bwa mem
					mem_cmd,order1 = aMapping.bwa_mem(n1,splitfq,qc=qc_flag,FFPE=argv['ffpe'])#20180116
					add_items(orders,order1)
					new_jobs['_'.join(['bwa_mem',eachsample,n1])] = \
						{'name' : '_'.join(['bwa_mem',eachsample,n1]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s -l p=%s' % (queue_list,th),
						'memory' : memory['bwa_mem'],
						'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['bwa_mem',eachsample,n1])+'.sh')}
					safe_open(os.path.join(mymapdir,'_'.join(['bwa_mem',eachsample,n1])+'.sh'),'w').write(mem_cmd)
					job_points['bwa_mem'].append('_'.join(['bwa_mem',eachsample,n1]))
					## sort bam
					sort_cmd,order1 = aMapping.samtools_sort(n1)
					add_items(orders,order1)
					new_jobs['_'.join(['samtools_sort',eachsample,n1])] = \
						{'name' : '_'.join(['samtools_sort',eachsample,n1]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s -l p=8' % queue_list,
						'memory' : memory['samtools_sort'],
						'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['samtools_sort',eachsample,n1])+'.sh')}
					safe_open(os.path.join(mymapdir,'_'.join(['samtools_sort',eachsample,n1])+'.sh'),'w').write(sort_cmd)
					job_points['samtools_sort'].append('_'.join(['samtools_sort',eachsample,n1]))
			else:
				merge_bams.append(eachlib)
				## bwa mem
				mem_cmd,order1 = aMapping.bwa_mem(eachlib,splitfq,qc=qc_flag,FFPE=argv['ffpe'])#20180116
				add_items(orders,order1)
				new_jobs['_'.join(['bwa_mem',eachsample,eachlib])] = \
					{'name' : '_'.join(['bwa_mem',eachsample,eachlib]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,th),
					'memory' : memory['bwa_mem'],
					'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['bwa_mem',eachsample,eachlib])+'.sh')}
				safe_open(os.path.join(mymapdir,'_'.join(['bwa_mem',eachsample,eachlib])+'.sh'),'w').write(mem_cmd)
				job_points['bwa_mem'].append('_'.join(['bwa_mem',eachsample,eachlib]))
				## gzip clean fq
#				if not argv['gz']:
#					gzip_cmd,order1 = aMapping.gzip_fq(eachlib)
#					add_items(orders,order1)
#					new_jobs['_'.join(['gzip_fq',eachsample,eachlib])] = \
#						{'name' : '_'.join(['gzip_fq',eachsample,eachlib]),
#						'status' : 'waiting',
#						'sched' : '-V -cwd %s' % queue_list,
#						'memory' : memory['gzip_fq'],
#						'cmd' : 'sh '+os.path.join(qcdir,eachsample,'_'.join(['gzip_fq',eachsample,eachlib])+'.sh')}
#					safe_open(os.path.join(qcdir,eachsample,'_'.join(['gzip_fq',eachsample,eachlib])+'.sh'),'w').write(gzip_cmd)
#					job_points['gzip_fq'].append('_'.join(['gzip_fq',eachsample,eachlib]))
				## sort bam
				sort_cmd,order1 = aMapping.samtools_sort(eachlib)
				add_items(orders,order1)
				new_jobs['_'.join(['samtools_sort',eachsample,eachlib])] = \
					{'name' : '_'.join(['samtools_sort',eachsample,eachlib]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['samtools_sort'],
					'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['samtools_sort',eachsample,eachlib])+'.sh')}
				safe_open(os.path.join(mymapdir,'_'.join(['samtools_sort',eachsample,eachlib])+'.sh'),'w').write(sort_cmd)
				job_points['samtools_sort'].append('_'.join(['samtools_sort',eachsample,eachlib]))
		## merge_bam
		merge_job = '_'.join(['picard_mergebam',eachsample])
		if merge_job in lastjobs:
			if lastjobs[merge_job]['status'] == 'done':
				merge_bam = ''
				if os.path.exists(os.path.join(mymapdir,eachsample+'.sort.bam')):
					merge_bam = os.path.join(mymapdir,eachsample+'.sort.bam')
				elif os.path.join(mymapdir,eachsample+'.sorted.bam'):
					merge_bam = os.path.join(mymapdir,eachsample+'.sorted.bam')
				else:
					print "Check you sort.bam file: "+ eachsample
					sys.exit(1)
				merge_cmd,order1 = aMapping.picard_mergebam(merge_bams,merge_bam)
				add_items(orders,order1)
				new_jobs[merge_job] = \
					{'name' : '_'.join(['picard_mergebam',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['picard_mergebam'],
					'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh')}
				safe_open(os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh'),'w').write(merge_cmd)
			else:
				merge_cmd,order1 = aMapping.picard_mergebam(merge_bams)
				add_items(orders,order1)
				new_jobs[merge_job] = \
					{'name' : '_'.join(['picard_mergebam',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['picard_mergebam'],
					'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh')}
				safe_open(os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh'),'w').write(merge_cmd)
		else:
			merge_cmd,order1 = aMapping.picard_mergebam(merge_bams)
			add_items(orders,order1)
			new_jobs[merge_job] = \
				{'name' : '_'.join(['picard_mergebam',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['picard_mergebam'],
				'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh')}
			safe_open(os.path.join(mymapdir,'_'.join(['picard_mergebam',eachsample])+'.sh'),'w').write(merge_cmd)
		job_points['picard_mergebam'].append('_'.join(['picard_mergebam',eachsample]))
		## mark duplicates
		rmdup_cmd,order1 = aMapping.picard_rmdupBam()
		add_items(orders,order1)
		new_jobs['_'.join(['picard_rmdupBam',eachsample])] = \
			{'name' : '_'.join(['picard_rmdupBam',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['picard_rmdupBam'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['picard_rmdupBam',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['picard_rmdupBam',eachsample])+'.sh'),'w').write(rmdup_cmd)
		job_points['picard_rmdupBam'].append('_'.join(['picard_rmdupBam',eachsample]))


                ## recal
                realn_recal = False
                realn_recal_bestPrac = False
                if analysis[2] == 2.2 and 'dbsnp' in genome_files[argv['genome']] and 'indel1000' in genome_files[argv['genome']]:
                        realn_recal = True
                        realn_recal_bestPrac = True
                        recal_cmd,order1 = aMapping.gatk_bqsr_bestPrac()
                        add_items(orders,order1)
                        new_jobs['_'.join(['gatk_bqsr',eachsample])] = \
                                {'name' : '_'.join(['gatk_bqsr',eachsample]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s -l p=6' % queue_list,
                                'memory' : memory['gatk_bqsr'],
                                'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['gatk_bqsr',eachsample])+'.sh')}
                        safe_open(os.path.join(mymapdir,'_'.join(['gatk_bqsr',eachsample])+'.sh'),'w').write(recal_cmd)
                        job_points['gatk_bqsr'].append('_'.join(['gatk_bqsr',eachsample]))

		
		## realin and recal
		if analysis[2] == 2.3 and 'dbsnp' in genome_files[argv['genome']] and 'indel1000' in genome_files[argv['genome']]:
			realn_recal = True
			realn_cmd,order1 = aMapping.gatk_realign()
			add_items(orders,order1)
			new_jobs['_'.join(['gatk_realign',eachsample])] = \
				{'name' : '_'.join(['gatk_realign',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=6' % queue_list,
				'memory' : memory['gatk_realign'],
				'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['gatk_realign',eachsample])+'.sh')}
			safe_open(os.path.join(mymapdir,'_'.join(['gatk_realign',eachsample])+'.sh'),'w').write(realn_cmd)
			job_points['gatk_realign'].append('_'.join(['gatk_realign',eachsample]))

			recal_cmd,order1 = aMapping.gatk_bqsr()
			add_items(orders,order1)
			new_jobs['_'.join(['gatk_bqsr',eachsample])] = \
				{'name' : '_'.join(['gatk_bqsr',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=6' % queue_list,
				'memory' : memory['gatk_bqsr'],
				'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['gatk_bqsr',eachsample])+'.sh')}
			safe_open(os.path.join(mymapdir,'_'.join(['gatk_bqsr',eachsample])+'.sh'),'w').write(recal_cmd)
			job_points['gatk_bqsr'].append('_'.join(['gatk_bqsr',eachsample]))
		
		## bamMD5
		final_cmd,order1 = aMapping.bamMD5(realn_recal,realn_recal_bestPrac)
		add_items(orders,order1)
		new_jobs['_'.join(['bamMD5',eachsample])] = \
			{'name' : '_'.join(['bamMD5',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['bamMD5'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['bamMD5',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['bamMD5',eachsample])+'.sh'),'w').write(final_cmd)
		job_points['bamMD5'].append('_'.join(['bamMD5',eachsample]))		
		## final bam
		final_cmd,order1 = aMapping.finalbam(realn_recal,realn_recal_bestPrac)
		add_items(orders,order1)
		new_jobs['_'.join(['finalbam',eachsample])] = \
			{'name' : '_'.join(['finalbam',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['finalbam'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['finalbam',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['finalbam',eachsample])+'.sh'),'w').write(final_cmd)
		job_points['finalbam'].append('_'.join(['finalbam',eachsample]))
		## pileup
		if mpileup:
			mpileup_cmd,order1 = aMapping.mpileup()
			add_items(orders,order1)
			new_jobs['_'.join(['mpileup',eachsample])] = \
				{'name' : '_'.join(['mpileup',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['mpileup'],
				'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['mpileup',eachsample])+'.sh')}
			safe_open(os.path.join(mymapdir,'_'.join(['mpileup',eachsample])+'.sh'),'w').write(mpileup_cmd)
			job_points['mpileup'].append('_'.join(['mpileup',eachsample]))	
		## sclip
		if sclip:
			subtmp = []
			for i,sub in enumerate(chr_subs):
				i = str(i+1)
				subtmp.append(i)
				sclip_cmd,order1 = aMapping.extractSoftClip_sub(sub,i)
				add_items(orders,order1)
				new_jobs['_'.join(['extractSoftClip_sub%s' % i,eachsample])] = \
					{'name' : '_'.join(['extractSoftClip_sub%s'%i,eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['extractSoftClip_sub'],
					'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['extractSoftClip_sub%s'%i,eachsample])+'.sh')}
				safe_open(os.path.join(mymapdir,'_'.join(['extractSoftClip_sub%s'%i,eachsample])+'.sh'),'w').write(sclip_cmd)
				job_points['extractSoftClip_sub'].append('_'.join(['extractSoftClip_sub%s' % i,eachsample]))
			## sclip cat 
			if argv['ffpe']:
				wtfp = 'Y'
			else:
				wtfp = 'N'
			cat_sclip_cmd,order1 = aMapping.catSoftClip(subtmp,chrs,wtfp)
			add_items(orders,order1)
			new_jobs['_'.join(['catSoftClip',eachsample])] = \
				{'name' : '_'.join(['catSoftClip',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['catSoftClip'],
				'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['catSoftClip',eachsample])+'.sh')}
			safe_open(os.path.join(mymapdir,'_'.join(['catSoftClip',eachsample])+'.sh'),'w').write(cat_sclip_cmd)
			job_points['catSoftClip'].append('_'.join(['catSoftClip',eachsample]))
		######################################   alnstat   #####################
		myalndir = os.path.join(alndir,eachsample)
		create_dir(myalndir)
		##__init__(self, pacientID, sampleID, memDir, statDir,genome_info,softwares,seqstrag,TR='')
		aAlnStat = AlnStat(sam2pid[eachsample],eachsample,mymapdir,myalndir,genome_files[argv['genome']], softwares, seqstrag, ref, sam2sex[eachsample], TR)

		## depstat
		dep_cmd,order1 = aAlnStat.depthStat()
		add_items(orders,order1)
		new_jobs['_'.join(['depthStat',eachsample])] = \
			{'name' : '_'.join(['depthStat',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['depthStat'],
			'cmd' : 'sh '+os.path.join(myalndir,'_'.join(['depthStat',eachsample])+'.sh')}
		safe_open(os.path.join(myalndir,'_'.join(['depthStat',eachsample])+'.sh'),'w').write(dep_cmd)
		job_points['depthStat'].append('_'.join(['depthStat',eachsample]))
		## flagstat
		flagstat_cmd,order1 = aAlnStat.flagstat()
		add_items(orders,order1)
		new_jobs['_'.join(['flagstat',eachsample])] = \
			{'name' : '_'.join(['flagstat',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['flagstat'],
			'cmd' : 'sh '+os.path.join(myalndir,'_'.join(['flagstat',eachsample])+'.sh')}
		safe_open(os.path.join(myalndir,'_'.join(['flagstat',eachsample])+'.sh'),'w').write(flagstat_cmd)
		job_points['flagstat'].append('_'.join(['flagstat',eachsample]))
		## combine
		combine_cmd,order1 = aAlnStat.combine()
		add_items(orders,order1)
		new_jobs['_'.join(['combine',eachsample])] = \
			{'name' : '_'.join(['combine',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['combine'],
			'cmd' : 'sh '+os.path.join(myalndir,'_'.join(['combine',eachsample])+'.sh')}
		safe_open(os.path.join(myalndir,'_'.join(['combine',eachsample])+'.sh'),'w').write(combine_cmd)
		job_points['combine'].append('_'.join(['combine',eachsample]))
		## remove_rmdup
		removermdup_cmd,order1 = aAlnStat.remove_rmdup(realn_recal)
		add_items(orders,order1)
		new_jobs['_'.join(['remove_rmdup',eachsample])] = \
			{'name' : '_'.join(['remove_rmdup',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['remove_rmdup'],
			'cmd' : 'sh '+os.path.join(myalndir,'_'.join(['remove_rmdup',eachsample])+'.sh')}
		safe_open(os.path.join(myalndir,'_'.join(['remove_rmdup',eachsample])+'.sh'),'w').write(removermdup_cmd)
		job_points['remove_rmdup'].append('_'.join(['remove_rmdup',eachsample]))
		## remove_sort.bam
		remove_sortbam_cmd,order1 = aMapping.remove_sortbam(realn_recal)
		add_items(orders,order1)
		new_jobs['_'.join(['remove_sortbam',eachsample])] = \
			{'name' : '_'.join(['remove_sortbam',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['remove_sortbam'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['remove_sortbam',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['remove_sortbam',eachsample])+'.sh'),'w').write(remove_sortbam_cmd)
		job_points['remove_sortbam'].append('_'.join(['remove_sortbam',eachsample]))
	## jobstatus
	if lastjobs:
		for eachsample in list_in_qc:
			for each in list_in_qc[eachsample]:
				bwajob = '_'.join(['bwa_mem',eachsample,list_in_qc[eachsample][each]])
				sortjob = '_'.join(['samtools_sort',eachsample,list_in_qc[eachsample][each]])
				if bwajob not in new_jobs and bwajob in lastjobs:
					add_jobs[bwajob] = lastjobs[bwajob]['txt']
				if sortjob not in new_jobs and sortjob in lastjobs:
					add_jobs[sortjob] = lastjobs[sortjob]['txt']
#				if not argv['gz']:
#					gzipjob = '_'.join(['gzip_fq',eachsample,list_in_qc[eachsample][each]])
#					if gzipjob not in new_jobs and gzipjob in lastjobs:
#						add_jobs[gzipjob] = lastjobs[gzipjob]['txt']
			mergejob = '_'.join(['picard_mergebam',eachsample])
			if mergejob not in new_jobs and mergejob in lastjobs:
				add_jobs[mergejob] = lastjobs[mergejob]['txt']
			rmdupjob = '_'.join(['picard_rmdupBam',eachsample])
			if rmdupjob not in new_jobs and rmdupjob in lastjobs:
				add_jobs[rmdupjob] = lastjobs[rmdupjob]['txt']
			realnjob = '_'.join(['gatk_realign',eachsample])
			if realnjob not in new_jobs and realnjob in lastjobs:
				add_jobs[realnjob] = lastjobs[realnjob]['txt']
			bqsrjob = '_'.join(['gatk_bqsr',eachsample])
			if bqsrjob not in new_jobs and bqsrjob in lastjobs:
				add_jobs[bqsrjob] = lastjobs[bqsrjob]['txt']
			finaljob = '_'.join(['finalbam',eachsample])
			if finaljob not in new_jobs and finaljob in lastjobs:
				add_jobs[finaljob] = lastjobs[finaljob]['txt']
			bamMD5job = '_'.join(['bamMD5',eachsample])
			if bamMD5job not in new_jobs and bamMD5job in lastjobs:
				add_jobs[bamMD5job] = lastjobs[bamMD5job]['txt']			
			if mpileup:
				mpileupjob = '_'.join(['mpileup',eachsample])
				if mpileupjob not in new_jobs and mpileupjob in lastjobs:
					add_jobs[mpileupjob] = lastjobs[mpileupjob]['txt']
			if sclip:
				for i,sub in enumerate(chr_subs):
					i = str(i+1)
					sclipjob = '_'.join(['extractSoftClip_sub%s' % i,eachsample])
					if sclipjob not in new_jobs and sclipjob in lastjobs:
						add_jobs[sclipjob] = lastjobs[sclipjob]['txt']
					catsclipjob = '_'.join(['catSoftClip',eachsample])
					if catsclipjob not in new_jobs and catsclipjob in lastjobs:
						add_jobs[catsclipjob] = lastjobs[catsclipjob]['txt']
			## alnstat
			depjob = '_'.join(['depthStat',eachsample])
			if depjob not in new_jobs and depjob in lastjobs:
				add_jobs[depjob] = lastjobs[depjob]['txt']
			flagjob = '_'.join(['flagstat',eachsample])
			if flagjob not in new_jobs and flagjob in lastjobs:
				add_jobs[flagjob] = lastjobs[flagjob]['txt']
			rm_dupjob = '_'.join(['remove_rmdup',eachsample])
			if rm_dupjob not in new_jobs and rm_dupjob in lastjobs:
				add_jobs[rm_dupjob] = lastjobs[rm_dupjob]['txt']
			rm_sortbamjob = '_'.join(['remove_sortbam',eachsample])
			if rm_sortbamjob not in new_jobs and  rm_sortbamjob in lastjobs:
				add_jobs[rm_sortbamjob] = lastjobs[rm_sortbamjob]['txt']
			combinejob = '_'.join(['combine',eachsample])
			if combinejob not in new_jobs and combinejob in lastjobs:
				add_jobs[combinejob] = lastjobs[combinejob]['txt']

#######################################   Mutation   ##################################
if set([1,2,3]).issubset(includes):
	create_dir(mutdir)
	print "Mutation ..."
	if analysis[3] == 3.1:
                snvsoft = 'samtools'  #by sun
		print "   ... samtools\n"
		job_points['samtoolsMpileup_sub'] = []
		job_points['cat_sub_all'] = []
		job_points['filtersamtoolsCalling'] = []
	else:
                snvsoft = 'GATK'
		if not argv['germcallbychrom']:
			print "   ... GATK\n"
			job_points['gatk_calling'] = []
		if argv['germcallbychrom']:
			print "   ... GATK calling by chromosome"
			job_points['gatk_calling_chr1to5'] = []
			job_points['gatk_calling_chr6to12'] = []
			job_points['gatk_calling_chr13toMT'] = []
			job_points['vcf_indexing'] = []
		job_points['gatk_variantion_filter'] = []
	job_points['annotatVcf'] = []
	samp_in_mergevcf=[]
	for eachsample in list_in_sample:
		if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample: #for paired samples, do not call germline mutation for the tumor sample.
			continue
		samp_in_mergevcf.append(eachsample)
		mymapdir = os.path.join(mapdir,eachsample)
		myalndir = os.path.join(alndir,eachsample)
		mut_method = 'samtools'
		if analysis[3] == 3.1:
			mymutdir = os.path.join(mutdir,eachsample+'.samtools')
			create_dir(mymutdir)
			### samtools
			#__init__(self, pacientID, sampleID, alignDir, mutationDir, genome_info,TR='')
			aMutationCalling = MutationCalling(sam2pid[eachsample],eachsample,mymapdir,myalndir,mymutdir,softwares,genome_files[argv['genome']],ref,seqstrag,sam2sex[eachsample], TR,flank=flank,germcallbychrom=germcallbychrom)
			### samtools mpile
			subtmp = []
			for i,sub in enumerate(chr_subs):
				i = str(i+1)
				subtmp.append(i)
				samtoolsMpileupBysub_cmd,order1 = aMutationCalling.samtoolsMpileup_sub(sub,i)
				add_items(orders,order1)
				new_jobs['_'.join(['samtoolsMpileup_sub%s' % i,eachsample])] = \
					{'name' : '_'.join(['samtoolsMpileup_sub%s'%i,eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['samtoolsMpileup_sub'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['samtoolsMpileup_sub%s'%i,eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['samtoolsMpileup_sub%s'%i,eachsample])+'.sh'),'w').write(samtoolsMpileupBysub_cmd)
				job_points['samtoolsMpileup_sub'].append('_'.join(['samtoolsMpileup_sub%s' % i,eachsample]))
			## samtools cat 
			cat_sub_all_cmd,order1 = aMutationCalling.cat_sub_all(subtmp,chrs)
			add_items(orders,order1)
			new_jobs['_'.join(['cat_sub_all',eachsample])] = \
				{'name' : '_'.join(['cat_sub_all',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['cat_sub_all'],
				'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['cat_sub_all',eachsample])+'.sh')}
			safe_open(os.path.join(mymutdir,'_'.join(['cat_sub_all',eachsample])+'.sh'),'w').write(cat_sub_all_cmd)
			job_points['cat_sub_all'].append('_'.join(['cat_sub_all',eachsample]))
			### samtools filter
			filter_cmd,order1 = aMutationCalling.filtersamtoolsCalling()
			add_items(orders,order1)
			new_jobs['_'.join(['filtersamtoolsCalling',eachsample])] = \
				{'name' : '_'.join(['filtersamtoolsCalling',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['filtersamtoolsCalling'],
				'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['filtersamtoolsCalling',eachsample])+'.sh')}
			safe_open(os.path.join(mymutdir,'_'.join(['filtersamtoolsCalling',eachsample])+'.sh'),'w').write(filter_cmd)
			job_points['filtersamtoolsCalling'].append('_'.join(['filtersamtoolsCalling',eachsample]))
		### GATK
		else:
			mut_method = 'GATK'
			mymutdir = os.path.join(mutdir,eachsample+'.GATK')
			create_dir(mymutdir)
			#__init__(self, pacientID, sampleID, alignDir, mutationDir, softwares, genome_info,TR='')
			aMutationCalling = MutationCalling(sam2pid[eachsample],eachsample,mymapdir,myalndir,mymutdir,softwares,genome_info,ref,seqstrag,sam2sex[eachsample], TR,flank=flank,germcallbychrom=germcallbychrom)
			## gatk calling
			if not argv['germcallbychrom']:
				gatkcalling_cmd,order1 = aMutationCalling.gatk_calling()
				add_items(orders,order1)
				new_jobs['_'.join(['gatk_calling',eachsample])] = \
					{'name' : '_'.join(['gatk_calling',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['gatk_calling'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['gatk_calling',eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['gatk_calling',eachsample])+'.sh'),'w').write(gatkcalling_cmd)
				job_points['gatk_calling'].append('_'.join(['gatk_calling',eachsample]))
			if argv['germcallbychrom']:
				gatkcallingchr1to5_cmd,order1 = aMutationCalling.gatk_calling_chr1to5()
				add_items(orders,order1)
				new_jobs['_'.join(['gatk_calling_chr1to5',eachsample])] = \
					{'name' : '_'.join(['gatk_calling_chr1to5',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['gatk_calling_chr1to5'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['gatk_calling_chr1to5',eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['gatk_calling_chr1to5',eachsample])+'.sh'),'w').write(gatkcallingchr1to5_cmd)
				job_points['gatk_calling_chr1to5'].append('_'.join(['gatk_calling_chr1to5',eachsample]))
				
				gatkcallingchr6to12_cmd,order1 = aMutationCalling.gatk_calling_chr6to12()
				add_items(orders,order1)
				new_jobs['_'.join(['gatk_calling_chr6to12',eachsample])] = \
					{'name' : '_'.join(['gatk_calling_chr6to12',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['gatk_calling_chr6to12'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['gatk_calling_chr6to12',eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['gatk_calling_chr6to12',eachsample])+'.sh'),'w').write(gatkcallingchr6to12_cmd)
				job_points['gatk_calling_chr6to12'].append('_'.join(['gatk_calling_chr6to12',eachsample]))
				
				gatkcallingchr13toMT_cmd,order1 = aMutationCalling.gatk_calling_chr13toMT()
				add_items(orders,order1)
				new_jobs['_'.join(['gatk_calling_chr13toMT',eachsample])] = \
					{'name' : '_'.join(['gatk_calling_chr13toMT',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['gatk_calling_chr13toMT'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['gatk_calling_chr13toMT',eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['gatk_calling_chr13toMT',eachsample])+'.sh'),'w').write(gatkcallingchr13toMT_cmd)
				job_points['gatk_calling_chr13toMT'].append('_'.join(['gatk_calling_chr13toMT',eachsample]))
				
				vcf_indexing_cmd,order1 = aMutationCalling.vcf_indexing()
				add_items(orders,order1)
				new_jobs['_'.join(['vcf_indexing',eachsample])] = \
					{'name' : '_'.join(['vcf_indexing',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['vcf_indexing'],
					'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['vcf_indexing',eachsample])+'.sh')}
				safe_open(os.path.join(mymutdir,'_'.join(['vcf_indexing',eachsample])+'.sh'),'w').write(vcf_indexing_cmd)
				job_points['vcf_indexing'].append('_'.join(['vcf_indexing',eachsample]))
			## gatk_variantion_filter
			gatkfilter_cmd,order1 = aMutationCalling.gatk_variantion_filter()
			add_items(orders,order1)
			new_jobs['_'.join(['gatk_variantion_filter',eachsample])] = \
				{'name' : '_'.join(['gatk_variantion_filter',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['gatk_variantion_filter'],
				'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['gatk_variantion_filter',eachsample])+'.sh')}
			safe_open(os.path.join(mymutdir,'_'.join(['gatk_variantion_filter',eachsample])+'.sh'),'w').write(gatkfilter_cmd)
			job_points['gatk_variantion_filter'].append('_'.join(['gatk_variantion_filter',eachsample]))

		mymutdir = os.path.join(mutdir,eachsample+'.'+mut_method)
		## mutation annotation
		annovar_cmd,order1 = aMutationCalling.annotatVcf(mut_method)
		add_items(orders,order1)
		new_jobs['_'.join(['annotatVcf',eachsample])] = \
			{'name' : '_'.join(['annotatVcf',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['annotatVcf'],
			'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['annotatVcf',eachsample])+'.sh')}
		safe_open(os.path.join(mymutdir,'_'.join(['annotatVcf',eachsample])+'.sh'),'w').write(annovar_cmd)
		job_points['annotatVcf'].append('_'.join(['annotatVcf',eachsample]))
	
	### merge vcf and summary
	# if TumorGermline:
		# merge_cmd,order1 = aMutationCalling.merge_vcf(mut_method,samp_in_mergevcf,analydir,sample_order)
	# else:
		# merge_cmd,order1 = aMutationCalling.merge_vcf(mut_method,list_in_sample.keys(),analydir,sample_order)
	# add_items(orders,order1)
	# new_jobs['merge_vcf'] = \
		# {'name' : 'merge_vcf',
		# 'status' : 'waiting',
		# 'sched' : '-V -cwd %s' % queue_list,
		# 'memory' : memory['merge_vcf'],
		# 'cmd' : 'sh '+os.path.join(mutdir,'merge_vcf','merge_vcf.sh')}
	# safe_open(os.path.join(mutdir,'merge_vcf','merge_vcf.sh'),'w').write(merge_cmd)
	# job_points['merge_vcf'] = ['merge_vcf']
	if TumorGermline:
		merge_cmd,order1 = aMutationCalling.merge_vcf(mut_method,samp_in_mergevcf,analydir,sample_order)
	else:
		merge_cmd,order1 = aMutationCalling.merge_vcf(mut_method,list_in_sample.keys(),analydir,sample_order)
	add_items(orders,order1)
	new_jobs['merge_vcf'] = \
		{'name' : 'merge_vcf',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : memory['merge_vcf'],
		'cmd' : 'sh '+os.path.join(mutdir,'merge_vcf','merge_vcf.sh')}
	safe_open(os.path.join(mutdir,'merge_vcf','merge_vcf.sh'),'w').write(merge_cmd)
	job_points['merge_vcf'] = ['merge_vcf']
	
	# jobstatus
	if lastjobs:
		if analysis[3] == 3.1:
		## samtools
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				mpileupjob = '_'.join(['samtoolsMpileup_sub%s' % chr,eachsample])
				if mpileupjob not in new_jobs and mpileupjob in lastjobs:
					add_jobs[mpileupjob] = lastjobs[mpileupjob]['txt']
				catjob = '_'.join(['cat_sub_all',eachsample])
				if catjob not in new_jobs and catjob in lastjobs:
					add_jobs[catjob] = lastjobs[catjob]['txt']
				filterjob = '_'.join(['filtersamtoolsCalling',eachsample])
				if filterjob not in new_jobs and filterjob in lastjobs:
					add_jobs[filterjob] = lastjobs[filterjob]['txt']
		else:
		## gatk
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				gatkcalljob = '_'.join(['gatk_calling',eachsample])
				if gatkcalljob not in new_jobs and gatkcalljob in lastjobs:
					add_jobs[gatkcalljob] = lastjobs[gatkcalljob]['txt']
				filterjob = '_'.join(['gatk_variantion_filter',eachsample])
				if filterjob not in new_jobs and filterjob in lastjobs:
					add_jobs[filterjob] = lastjobs[filterjob]['txt']
		## annovar
		for eachsample in list_in_qc:
			if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
				continue
			annovarjob = '_'.join(['annotatVcf',eachsample])
			if annovarjob not in new_jobs and annovarjob in lastjobs:
				add_jobs[annovarjob] = lastjobs[annovarjob]['txt']




#########################   Germline SV
if set([1,2,4]).issubset(includes):
	print "SV ..."
	if analysis[4] == 4.1:
		print "   ... breakdancer\n"
		job_points['bam2cfg'] = []
		job_points['breakdancerSV'] = []
		job_points['breakdancerSvAnno'] = []
	elif analysis[4] == 4.2:
		print "   ... crest\n"
		job_points['crestSV'] = []
		job_points['crestSVann'] = []
	elif analysis[4] == 4.3:
		print "   ... delly\n"
		job_points['dellySV_del'] = []
		job_points['dellySV_dup'] = []
		job_points['dellySV_ins'] = []
		job_points['dellySV_inv'] = []
		job_points['dellySV_tra'] = []
		job_points['dellySVanno'] = []		
	else:
		print "   ... lumpy\n"
		job_points['lumpySV'] = []
		job_points['lumpySVann'] = []
		
	create_dir(svdir)
	sv_method = 'breakdancer'
	for eachsample in list_in_sample:
		if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
			continue
		mymapdir = os.path.join(mapdir,eachsample)
		mysvdir = os.path.join(svdir,eachsample)
		create_dir(mysvdir)
		aSV = SV(sam2pid[eachsample],eachsample,genome_files[argv['genome']],TR, analydir , mysvdir, softwares,SV_pipeline_dir,moduledir, database_dir, ref, sam2sex[eachsample])
		if analysis[4] == 4.1:
		## breakdancer
			## bam to cfg
			mybreakdancerdir = os.path.join(mysvdir,'breakdancer')
			create_dir(mybreakdancerdir)
			bam2cfg_cmd,order1 = aSV.bam2cfg()
			add_items(orders,order1)
			new_jobs['_'.join(['bam2cfg',eachsample])] = \
				{'name' : '_'.join(['bam2cfg',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['bam2cfg'],
				'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['bam2cfg',eachsample])+'.sh')}
			safe_open(os.path.join(mybreakdancerdir,'_'.join(['bam2cfg',eachsample])+'.sh'),'w').write(bam2cfg_cmd)
			job_points['bam2cfg'].append('_'.join(['bam2cfg',eachsample]))
			## breakdancer SV
			breakdancer_cmd,order1 = aSV.breakdancerSV()
			add_items(orders,order1)
			new_jobs['_'.join(['breakdancerSV',eachsample])] = \
				{'name' : '_'.join(['breakdancerSV',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['breakdancerSV'],
				'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['breakdancerSV',eachsample])+'.sh')}
			safe_open(os.path.join(mybreakdancerdir,'_'.join(['breakdancerSV',eachsample])+'.sh'),'w').write(breakdancer_cmd)
			job_points['breakdancerSV'].append('_'.join(['breakdancerSV',eachsample]))
			## breakdancer_SvAnno
			svann_cmd,order1 = aSV.breakdancerSvAnno()
			add_items(orders,order1)
			new_jobs['_'.join(['breakdancerSvAnno',eachsample])] = \
				{'name' : '_'.join(['breakdancerSvAnno',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['breakdancerSvAnno'],
				'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['breakdancerSvAnno',eachsample])+'.sh')}
			safe_open(os.path.join(mybreakdancerdir,'_'.join(['breakdancerSvAnno',eachsample])+'.sh'),'w').write(svann_cmd)
			job_points['breakdancerSvAnno'].append('_'.join(['breakdancerSvAnno',eachsample]))
		### crest
		elif analysis[4] == 4.2:
			## crest
			sv_method = 'crest'
			mycrestdir = os.path.join(mysvdir,'crest')
			create_dir(mycrestdir)
			crest_cmd,order1 = aSV.crestSV()
			add_items(orders,order1)
			new_jobs['_'.join(['crestSV',eachsample])] = \
				{'name' : '_'.join(['crestSV',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=8' % queue_list,
				'memory' : memory['crestSV'],
				'cmd' : 'sh '+os.path.join(mycrestdir,'_'.join(['crestSV',eachsample])+'.sh')}
			safe_open(os.path.join(mycrestdir,'_'.join(['crestSV',eachsample])+'.sh'),'w').write(crest_cmd)
			job_points['crestSV'].append('_'.join(['crestSV',eachsample]))
			## crest annotation
			crestann_cmd,order1 = aSV.crestSVann()
			add_items(orders,order1)
			new_jobs['_'.join(['crestSVann',eachsample])] = \
				{'name' : '_'.join(['crestSVann',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['crestSVann'],
				'cmd' : 'sh '+os.path.join(mycrestdir,'_'.join(['crestSVann',eachsample])+'.sh')}
			safe_open(os.path.join(mycrestdir,'_'.join(['crestSVann',eachsample])+'.sh'),'w').write(crestann_cmd)
			job_points['crestSVann'].append('_'.join(['crestSVann',eachsample]))
		elif analysis[4] == 4.3:
			## delly
			sv_method = 'delly'
			mydellydir = os.path.join(mysvdir,'delly')
			create_dir(mydellydir)
			# dellySV_del
			dellysvdel_cmd,order1 = aSV.dellySV_del()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySV_del',eachsample])] = \
				{'name' : '_'.join(['dellySV_del',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySV_del'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySV_del',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySV_del',eachsample])+'.sh'),'w').write(dellysvdel_cmd)
			job_points['dellySV_del'].append('_'.join(['dellySV_del',eachsample]))
			# dellySV_dup
			dellysvdup_cmd,order1 = aSV.dellySV_dup()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySV_dup',eachsample])] = \
				{'name' : '_'.join(['dellySV_dup',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySV_dup'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySV_dup',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySV_dup',eachsample])+'.sh'),'w').write(dellysvdup_cmd)
			job_points['dellySV_dup'].append('_'.join(['dellySV_dup',eachsample]))
			# dellySV_ins
			dellysvins_cmd,order1 = aSV.dellySV_ins()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySV_ins',eachsample])] = \
				{'name' : '_'.join(['dellySV_ins',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySV_ins'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySV_ins',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySV_ins',eachsample])+'.sh'),'w').write(dellysvins_cmd)
			job_points['dellySV_ins'].append('_'.join(['dellySV_ins',eachsample]))
			# dellySV_inv
			dellysvinv_cmd,order1 = aSV.dellySV_inv()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySV_inv',eachsample])] = \
				{'name' : '_'.join(['dellySV_inv',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySV_inv'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySV_inv',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySV_inv',eachsample])+'.sh'),'w').write(dellysvinv_cmd)
			job_points['dellySV_inv'].append('_'.join(['dellySV_inv',eachsample]))
			# dellySV_tra
			dellysvtra_cmd,order1 = aSV.dellySV_tra()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySV_tra',eachsample])] = \
				{'name' : '_'.join(['dellySV_tra',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySV_tra'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySV_tra',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySV_tra',eachsample])+'.sh'),'w').write(dellysvtra_cmd)
			job_points['dellySV_tra'].append('_'.join(['dellySV_tra',eachsample]))				
			## delly annotation
			dellyann_cmd,order1 = aSV.dellySVanno()
			add_items(orders,order1)
			new_jobs['_'.join(['dellySVanno',eachsample])] = \
				{'name' : '_'.join(['dellySVanno',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['dellySVanno'],
				'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['dellySVanno',eachsample])+'.sh')}
			safe_open(os.path.join(mydellydir,'_'.join(['dellySVanno',eachsample])+'.sh'),'w').write(dellyann_cmd)
			job_points['dellySVanno'].append('_'.join(['dellySVanno',eachsample]))				
			
		## lumpy
		else:
			sv_method = 'lumpy'
			mylumpydir = os.path.join(mysvdir,'lumpy')
			create_dir(mylumpydir)
			lumpy_cmd,order1 = aSV.lumpySV()
			add_items(orders,order1)
			new_jobs['_'.join(['lumpySV',eachsample])] = \
				{'name' : '_'.join(['lumpySV',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['lumpySV'],
				'cmd' : 'sh '+os.path.join(mylumpydir,'_'.join(['lumpySV',eachsample])+'.sh')}
			safe_open(os.path.join(mylumpydir,'_'.join(['lumpySV',eachsample])+'.sh'),'w').write(lumpy_cmd)
			job_points['lumpySV'].append('_'.join(['lumpySV',eachsample]))
			## lumpy annotation
			lumpyann_cmd,order1 = aSV.lumpySVann()
			add_items(orders,order1)
			new_jobs['_'.join(['lumpySVann',eachsample])] = \
				{'name' : '_'.join(['lumpySVann',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['lumpySVann'],
				'cmd' : 'sh '+os.path.join(mylumpydir,'_'.join(['lumpySVann',eachsample])+'.sh')}
			safe_open(os.path.join(mylumpydir,'_'.join(['lumpySVann',eachsample])+'.sh'),'w').write(lumpyann_cmd)

	## job status
	if lastjobs:
		if analysis[4] == 4.1:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				bam2cfgjob = '_'.join(['bam2cfg',eachsample])
				if bam2cfgjob not in new_jobs and bam2cfgjob in lastjobs:
					add_jobs[bam2cfgjob] = lastjobs[bam2cfgjob]['txt']
				breakdancerjob = '_'.join(['breakdancerSV',eachsample])
				if breakdancerjob not in new_jobs and breakdancerjob in lastjobs:
					add_jobs[breakdancerjob] = lastjobs[breakdancerjob]['txt']
				svannjob = '_'.join(['breakdancerSvAnno',eachsample])
				if svannjob not in new_jobs and svannjob in lastjobs:
					add_jobs[svannjob] = lastjobs[svannjob]['txt']
		elif analysis[4] == 4.2:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				crestjob = '_'.join(['crestSV',eachsample])
				if crestjob not in new_jobs and crestjob in lastjobs:
					add_jobs[crestjob] = lastjobs[crestjob]['txt']
				crestannjob = '_'.join(['crestSVann',eachsample])
				if crestannjob not in new_jobs and crestannjob in lastjobs:
					add_jobs[crestannjob] = lastjobs[crestannjob]['txt']
		elif analysis[4] == 4.3:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
			# dellySV_del
				dellysvdeljob = '_'.join(['dellySV_del',eachsample])
				if dellysvdeljob not in new_jobs and dellysvdeljob in lastjobs:
					add_jobs[dellysvdeljob] = lastjobs[dellysvdeljob]['txt']
			# dellySV_dup
				dellysvdupjob = '_'.join(['dellySV_dup',eachsample])
				if dellysvdupjob not in new_jobs and dellysvdupjob in lastjobs:
					add_jobs[dellysvdupjob] = lastjobs[dellysvdupjob]['txt']
			# dellySV_ins
				dellysvinsjob = '_'.join(['dellySV_ins',eachsample])
				if dellysvinsjob not in new_jobs and dellysvinsjob in lastjobs:
					add_jobs[dellysvinsjob] = lastjobs[dellysvinsjob]['txt']
			# dellySV_inv
				dellysvinvjob = '_'.join(['dellySV_inv',eachsample])
				if dellysvinvjob not in new_jobs and dellysvinvjob in lastjobs:
					add_jobs[dellysvinvjob] = lastjobs[dellysvinvjob]['txt']
			# dellySV_tra
				dellysvtrajob = '_'.join(['dellySV_tra',eachsample])
				if dellysvtrajob not in new_jobs and dellysvtrajob in lastjobs:
					add_jobs[dellysvtrajob] = lastjobs[dellysvtrajob]['txt']
			# delly annotation
				dellyannjob = '_'.join(['dellySVanno',eachsample])
				if dellyannjob not in new_jobs and dellyannjob in lastjobs:
					add_jobs[dellyannjob] = lastjobs[dellyannjob]['txt']					
		else:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				lumpyjob = '_'.join(['lumpySV',eachsample])
				if lumpyjob not in new_jobs and lumpyjob in lastjobs:
					add_jobs[lumpyjob] = lastjobs[lumpyjob]['txt']
				lumpyannjob = '_'.join(['lumpySVann',eachsample])
				if lumpyannjob not in new_jobs and lumpyannjob in lastjobs:
					add_jobs[lumpyannjob] = lastjobs[lumpyannjob]['txt']
			
#########################   Germline CNV
if set([1,2,5]).issubset(includes):
	print "CNV ..."
	create_dir(cnvdir)
	cnv_method = 'freec'
	## seqstrag
	seq = 'WES'
	if 'WGS' in seqstrag:
		seq = 'WGS'
	if analysis[5] == 5.1:
		print "   ... control-freec\n"
		job_points['freec_cnv'] = []
		#job_points['freecCNVann'] = []
		job_points['info_freec_circos'] = []
		job_points['info_snp_circos'] = []
		job_points['info_indel_circos'] = []
		job_points['info_depth_circos'] = []
		job_points['info_delly_circos'] = []
		job_points['info_crest_circos'] = []
		job_points['info_breakdancer_circos'] = []
		job_points['info_lumpy_circos'] = []
		job_points['conf_circos'] = []
		job_points['Circos'] = []
	elif analysis[5] == 5.3:
		print "   ... control-freec_loh\n"
		job_points['freec_cnv_loh'] = []
		#job_points['freecCNVann_loh'] = []
	else:
		print "   ... cnvnator\n"
		job_points['cnvnatorCNV'] = []
		job_points['cnvnatorCNVann'] = []
		
	## CNV analysis
	for eachsample in list_in_sample:
		if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
			continue
		mymapdir = os.path.join(mapdir,eachsample)
		mycnvdir = os.path.join(cnvdir,eachsample)
		create_dir(mycnvdir)
		aSV = SV(sam2pid[eachsample],eachsample,genome_files[argv['genome']],TR, analydir, mycnvdir, softwares,SV_pipeline_dir, moduledir, database_dir, ref, sam2sex[eachsample])
		if analysis[5] == 5.1:
		## freeC
			myfreecdir = os.path.join(mycnvdir,'freec')
			create_dir(myfreecdir)
			freec_cmd,order1 = aSV.freec_cnv(seq)
			add_items(orders,order1)
			new_jobs['_'.join(['freec_cnv',eachsample])] = \
				{'name' : '_'.join(['freec_cnv',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['freec_cnv'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freec_cnv',eachsample])+'.sh')}
			safe_open(os.path.join(myfreecdir,'_'.join(['freec_cnv',eachsample])+'.sh'),'w').write(freec_cmd)
			job_points['freec_cnv'].append('_'.join(['freec_cnv',eachsample]))
		## freec annotation
		#	freecann_cmd,order1 = aSV.freecCNVann()
		#	add_items(orders,order1)
		#	new_jobs['_'.join(['freecCNVann',eachsample])] = \
		#		{'name' : '_'.join(['freecCNVann',eachsample]),
		#		'status' : 'waiting',
		#		'sched' : '-V -cwd %s' % queue_list,
		#		'memory' : memory['freecCNVann'],
		#		'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freecCNVann',eachsample])+'.sh')}
		#	safe_open(os.path.join(myfreecdir,'_'.join(['freecCNVann',eachsample])+'.sh'),'w').write(freecann_cmd)
		#	job_points['freecCNVann'].append('_'.join(['freecCNVann',eachsample]))
		elif analysis[5] == 5.3:
		## freec LOH
			myfreecdir = os.path.join(mycnvdir,'freec_loh')
			create_dir(myfreecdir)
			freec_cmd,order1 = aSV.freec_cnv_loh(seq)
			add_items(orders,order1)
			new_jobs['_'.join(['freec_cnv_loh',eachsample])] = \
				{'name' : '_'.join(['freec_cnv_loh',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['freec_cnv_loh'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freec_cnv_loh',eachsample])+'.sh')}
			safe_open(os.path.join(myfreecdir,'_'.join(['freec_cnv_loh',eachsample])+'.sh'),'w').write(freec_cmd)
			job_points['freec_cnv_loh'].append('_'.join(['freec_cnv_loh',eachsample]))
		else:
		## cnvnator
			cnv_method = 'cnvnator'
			mycnvnatordir = os.path.join(mycnvdir,'cnvnator')
			create_dir(mycnvnatordir)
			cnvnator_cmd,order1 = aSV.cnvnatorCNV()
			add_items(orders,order1)
			new_jobs['_'.join(['cnvnatorCNV',eachsample])] = \
				{'name' : '_'.join(['cnvnatorCNV',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['cnvnatorCNV'],
				'cmd' : 'sh '+os.path.join(mycnvnatordir,'_'.join(['cnvnatorCNV',eachsample])+'.sh')}
			safe_open(os.path.join(mycnvnatordir,'_'.join(['cnvnatorCNV',eachsample])+'.sh'),'w').write(cnvnator_cmd)
			job_points['cnvnatorCNV'].append('_'.join(['cnvnatorCNV',eachsample]))
		# annotation
			cnvnatorann_cmd,order1 = aSV.cnvnatorCNVann()
			add_items(orders,order1)
			new_jobs['_'.join(['cnvnatorCNVann',eachsample])] = \
				{'name' : '_'.join(['cnvnatorCNVann',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['cnvnatorCNVann'],
				'cmd' : 'sh '+os.path.join(mycnvnatordir,'_'.join(['cnvnatorCNVann',eachsample])+'.sh')}
			safe_open(os.path.join(mycnvnatordir,'_'.join(['cnvnatorCNVann',eachsample])+'.sh'),'w').write(cnvnator_cmd)
			job_points['cnvnatorCNVann'].append('_'.join(['cnvnatorCNVann',eachsample]))
		#add Circos by yangjunhui 20170310
		circos = {}
		if analysis[5] == 5.1:
			cnvsoft = 'freec'
			Usefreec = os.path.join(myfreecdir,eachsample+'.final.bam_CNVs.final.%s_multianno.xls' % genome_info['annovarbuild'])
			circos['freec'] = Usefreec
			new_jobs['_'.join(['info_freec_circos',eachsample])] = \
				{'name' : '_'.join(['info_freec_circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['info_freec_circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_freec_circos',eachsample])+'.sh')}
			add_items(orders,('order info_freec_circos_%s after freec_cnv_%s' %(eachsample,eachsample)))
			add_items(orders,('order conf_circos_%s after freec_cnv_%s' %(eachsample,eachsample)))
			add_items(orders,('order Circos_%s after info_freec_circos_%s' %(eachsample,eachsample)))
			job_points['info_freec_circos'].append('_'.join(['info_freec_circos',eachsample]))
			
			#SNV/indel
			FinalBam = os.path.join(mapdir,eachsample,eachsample+'.final.bam')
			if analysis[3] == 3.1:
				snvsoft = 'samtools'
				s_mutationCalldir = os.path.join(mutdir,eachsample+'.samtools')
				UseVcfsnp = os.path.join(s_mutationCalldir,eachsample+'.samtools.snp.vcf.gz')
				UseVcfindel = os.path.join(s_mutationCalldir,eachsample+'.samtools.indel.vcf.gz')
			if analysis[3] == 3.2:
				snvsoft = 'GATK'
				s_mutationCalldir = os.path.join(mutdir,eachsample+'.GATK')
				UseVcfsnp = os.path.join(s_mutationCalldir,eachsample+'.GATK.snp.vcf.gz')
				UseVcfindel = os.path.join(s_mutationCalldir,eachsample+'.GATK.indel.vcf.gz')
			add_items(orders,('order info_snp_circos_%s after annotatVcf_%s' %(eachsample,eachsample)))
			add_items(orders,('order info_indel_circos_%s after annotatVcf_%s' %(eachsample,eachsample)))
			add_items(orders,('order info_depth_circos_%s after finalbam_%s' %(eachsample,eachsample)))
			add_items(orders,('order Circos_%s after info_snp_circos_%s' %(eachsample,eachsample)))
			add_items(orders,('order Circos_%s after info_indel_circos_%s' %(eachsample,eachsample)))
			add_items(orders,('order Circos_%s after info_depth_circos_%s' %(eachsample,eachsample)))
			circos['snp'] = UseVcfsnp
			circos['indel'] = UseVcfindel
			circos['depth'] = FinalBam
			new_jobs['_'.join(['info_snp_circos',eachsample])] = \
				{'name' : '_'.join(['info_snp_circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['info_snp_circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_snp_circos',eachsample])+'.sh')}
			new_jobs['_'.join(['info_indel_circos',eachsample])] = \
				{'name' : '_'.join(['info_indel_circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['info_indel_circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_indel_circos',eachsample])+'.sh')}
			new_jobs['_'.join(['info_depth_circos',eachsample])] = \
				{'name' : '_'.join(['info_depth_circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['info_depth_circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_depth_circos',eachsample])+'.sh')}
			job_points['info_snp_circos'].append('_'.join(['info_snp_circos',eachsample]))
			job_points['info_indel_circos'].append('_'.join(['info_indel_circos',eachsample]))
			job_points['info_depth_circos'].append('_'.join(['info_depth_circos',eachsample]))	
				
			#SV
			if analysis[4] == 4.3:
				svsoft = 'delly'
				s_dellydir = os.path.join(svdir,eachsample, 'delly')
				Usedelly = os.path.join(s_dellydir,eachsample+'.delly.sv.%s_multianno.xls' % genome_info['annovarbuild'])
				circos['delly'] = Usedelly
				new_jobs['_'.join(['info_delly_circos',eachsample])] = \
					{'name' : '_'.join(['info_delly_circos',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['info_delly_circos'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_delly_circos',eachsample])+'.sh')}
				add_items(orders,('order info_delly_circos_%s after dellySVanno_%s' %(eachsample,eachsample)))
				add_items(orders,('order Circos_%s after info_delly_circos_%s' %(eachsample,eachsample)))
				job_points['info_delly_circos'].append('_'.join(['info_delly_circos',eachsample]))
			if analysis[4] == 4.2:
				svsoft = 'crest'
				s_crestdir = os.path.join(svdir,eachsample, 'crest')
				Usecrest = os.path.join(s_crestdir,eachsample+'.crest.sv.%s_multianno.xls' % genome_info['annovarbuild'])
				circos['crest'] = Usecrest
				new_jobs['_'.join(['info_crest_circos',eachsample])] = \
					{'name' : '_'.join(['info_crest_circos',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['info_crest_circos'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_crest_circos',eachsample])+'.sh')}
				add_items(orders,('order info_crest_circos_%s after crestSVann_%s' %(eachsample,eachsample)))
				add_items(orders,('order Circos_%s after info_crest_circos_%s' %(eachsample,eachsample)))
				job_points['info_crest_circos'].append('_'.join(['info_crest_circos',eachsample]))
			if analysis[4] == 4.1:
				svsoft = 'breakdancer'
				s_breakdancerdir = os.path.join(svdir,eachsample, 'breakdancer')
				Usebreak =  os.path.join(s_breakdancerdir,eachsample+'.breakdancer.sv.%s_multianno.xls' % genome_info['annovarbuild'])
				circos['breakdancer'] = Usebreak
				new_jobs['_'.join(['info_breakdancer_circos',eachsample])] = \
					{'name' : '_'.join(['info_breakdancer_circos',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['info_breakdancer_circos'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_breakdancer_circos',eachsample])+'.sh')}
				add_items(orders,('order info_breakdancer_circos_%s after breakdancerSvAnno_%s' %(eachsample,eachsample)))
				add_items(orders,('order Circos_%s after info_breakdancer_circos_%s' %(eachsample,eachsample)))
				job_points['info_breakdancer_circos'].append('_'.join(['info_breakdancer_circos',eachsample]))
			if analysis[4] == 4.4:
				svsoft = 'lumpy'
				s_lumpydir = os.path.join(svdir,eachsample, 'lumpy')
				Uselumpy =  os.path.join(s_lumpydir,eachsample+'.lumpy.sv.%s_multianno.xls' % genome_info['annovarbuild'])
				circos['lumpy'] = Uselumpy
				new_jobs['_'.join(['info_lumpy_circos',eachsample])] = \
					{'name' : '_'.join(['info_lumpy_circos',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['info_lumpy_circos'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['info_lumpy_circos',eachsample])+'.sh')}
				add_items(orders,('order info_lumpy_circos_%s after lumpySVann_%s' %(eachsample,eachsample)))
				add_items(orders,('order Circos_%s after info_lumpy_circos_%s' %(eachsample,eachsample)))
				job_points['info_lumpy_circos'].append('_'.join(['info_lumpy_circos',eachsample]))
			CirTy = []
			CirTyv = []
			for ana_type in circos.keys():
				CirTy.append(ana_type)
				CirTyv.append(circos[ana_type])
			vtuse = ','.join(CirTy)
			inuse = ','.join(CirTyv)
			###the picture out dir will creat by circos.py
			assert not os.system('python %s/Varition/CNV/Circos_v6.py --id %s --vt %s --in %s --out %s --ref %s' % (moduledir,eachsample,vtuse,inuse,myfreecdir,ref))
			new_jobs['_'.join(['conf_circos',eachsample])] = \
				{'name' : '_'.join(['conf_circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['conf_circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['conf_circos',eachsample])+'.sh')}
			new_jobs['_'.join(['Circos',eachsample])] = \
				{'name' : '_'.join(['Circos',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : memory['Circos'],
				'cmd' : 'sh '+os.path.join(myfreecdir,'Circos','_'.join(['Circos',eachsample])+'.sh')}
			job_points['conf_circos'].append('_'.join(['conf_circos',eachsample]))
			job_points['Circos'].append('_'.join(['Circos',eachsample]))
			add_items(orders,('order Circos_%s after conf_circos_%s' %(eachsample,eachsample)))

		## job status
	if lastjobs:
		if analysis[5] == 5.1:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				freecjob = '_'.join(['freec_cnv',eachsample])
				if freecjob not in new_jobs and freecjob in lastjobs:
					add_jobs[freecjob] = lastjobs[freecjob]['txt']
				#freecannjob = '_'.join(['freecCNVann',eachsample])
				#if freecannjob not in new_jobs and freecannjob in lastjobs:
				#	add_jobs[freecannjob] = lastjobs[freecannjob]['txt']
			#Circos
				info_freecjob = '_'.join(['info_freec_circos',eachsample])
				info_snpjob = '_'.join(['info_snp_circos',eachsample])
				info_indeljob = '_'.join(['info_indel_circos',eachsample])
				info_depthjob = '_'.join(['info_depth_circos',eachsample])
				info_dellyjob = '_'.join(['info_delly_circos',eachsample])
				info_crestjob = '_'.join(['info_crest_circos',eachsample])
				info_breakdancerjob = '_'.join(['info_breakdancer_circos',eachsample])
				info_lumpyjob = '_'.join(['info_lumpy_circos',eachsample])
				confjob = '_'.join(['conf_circos',eachsample])
				Circosjob = '_'.join(['Circos',eachsample])
				joblist = [info_freecjob,info_snpjob,info_indeljob,info_depthjob,info_dellyjob,info_crestjob,info_breakdancerjob,info_lumpyjob,confjob,Circosjob]
				for jobname in joblist:
					if jobname not in new_jobs and jobname in lastjobs:
						add_jobs[jobname] = lastjobs[jobname]['txt']
		elif analysis[5] == 5.3:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				freeclohjob = '_'.join(['freec_cnv_loh',eachsample])
				if freeclohjob not in new_jobs and freeclohjob in lastjobs:
					add_jobs[freeclohjob] = lastjobs[freeclohjob]['txt']
		else:
			for eachsample in list_in_qc:
				if TumorGermline and pairinfo.has_key(eachsample) and pairinfo[eachsample] in list_in_sample:
					continue
				cnvnatorjob = '_'.join(['cnvnatorCNV',eachsample])
				if cnvnatorjob not in new_jobs and cnvnatorjob in lastjobs:
					add_jobs[cnvnatorjob] = lastjobs[cnvnatorjob]['txt']
				cnvnatorannjob = '_'.join(['cnvnatorCNVann',eachsample])
				if cnvnatorannjob not in new_jobs and cnvnatorannjob in lastjobs:
					add_jobs[cnvnatorannjob] = lastjobs[cnvnatorannjob]['txt']


samples_for_somatic = set()
for eachsample in list_in_sample:
	if eachsample in pairinfo:
		if pairinfo[eachsample] in list_in_sample:
			samples_for_somatic |= set([eachsample])

somatic_samlist = []  # somatic sample list
############################# Somatic SNP/INDEL, SV, CNV
for each in [list(samples_for_somatic)[i:min(i+4,len(samples_for_somatic))] for i in range(0,len(samples_for_somatic),4)]:
	print 'Tumor samples: %s\n'% ', '.join(each)
if somatic:

	if set([1,2,6]).issubset(includes):
		print "Somatic SNP/INDEL ..."
		if analysis[6] == 6.1:
			print "   ... samtools\n"
			job_points['somaticsamtoolsMpileup'] = []
			job_points['somaticsamtoolscalling'] = []
			job_points['somaticsamtoolsannovar'] = []
		elif analysis[6] == 6.2:
			print "   ... varscan\n"
			job_points['varscan'] = []
			job_points['somaticvarscanannovar'] = []
		else:
			if not argv['socallbychrom']:
				print "   ... muTect/Strelka\n"
				job_points['somatic_muTect'] = []
				job_points['somatic_Strelka'] = []
			if argv['socallbychrom']:
				print "   ... muTect SNP calling by chromosome"
				job_points['somatic_muTect_chr1to5'] = []
				job_points['somatic_muTect_chr6to12'] = []
				job_points['somatic_muTect_chr13toMT'] = []
				job_points['somatic_muTect_annot'] = []
				job_points['somatic_Strelka'] = []
	if set([1,2,7]).issubset(includes):
		print "Somatic SV ..."
		if analysis[7] == 7.1:
			print "   ... breakdancer\n"
			job_points['sobam2cfg'] = []
			job_points['sobreakdancerSV'] = []
			job_points['sobreakdancerSvAnno'] = []
		elif analysis[7] == 7.2:
			print "   ... crest\n"
			job_points['diff_sclip'] = []
			job_points['crest_somaticSV'] = []
			job_points['crest_somaticSVann'] = []
		elif analysis[7] == 7.3:
			print "   ... delly\n"
			job_points['delly_somatic_sv_del'] = []
			job_points['delly_somatic_sv_dup'] = []
			job_points['delly_somatic_sv_ins'] = []
			job_points['delly_somatic_sv_inv'] = []
			job_points['delly_somatic_sv_tra'] = []
			job_points['delly_somatic_sv_annot'] = []			
		elif analysis[7] == 7.4:
			print "   ... lumpy\n"
			job_points['lumpy_somaticSV'] = []
			job_points['lumpy_somaticSVann'] = []
	if set([1,2,8]).issubset(includes):
		seq = 'WES'
		if 'WGS' in seqstrag:
			seq = 'WGS'
		print "Somatic CNV ..."
		if analysis[8] == 8.1:
			print "   ... control-freec\n"
			job_points['freec_somaticcnv'] = []
		elif analysis[8] == 8.2:
			print "   ... ExomeCNV\n"
			job_points['exomeCNV_Ndepth'] = []
			job_points['exomeCNV_Tdepth'] = []
			job_points['exomeCNV_CNV'] = []
		elif analysis[8] == 8.4:
			print "   ... control-freec loh\n"
			job_points['freec_somaticcnv_loh'] = []
			#job_points['freec_somaticCNVann'] = []
		else:
			print "   ... varscan\n"
			job_points['varscan_somaticCNV'] = []
			job_points['varscan_somaticCNVann'] = []
	create_dir(somaticdir)
	sosnp_method = 'samtools'
	soindel_method = 'samtools'
	for eachsample in samples_for_somatic :
		if not pairinfo[eachsample]:
			continue
		somatic_samlist.append(eachsample)  # somatic sample list
#		Nsam = patientInfo[sam2pid[eachsample]]['N']
		Nsam = pairinfo[eachsample]
		assert not eachsample == Nsam
		mysomaticdir = os.path.join(somaticdir,eachsample)
		create_dir(mysomaticdir)
		#__init__(self, pacientID, sampleID, Nname, mapDir, genome_info, somaticDir, softwares, TR='',tFreq=0.1, nFreq=0.05)
		aSomatic = Somatic(sam2pid[eachsample], eachsample, Nsam, analydir, genome_files[argv['genome']], mysomaticdir, softwares, SV_pipeline_dir,argv['genome'],TR, seqstrag, ref, moduledir, database_dir, sam2sex[eachsample], tFreq, nFreq, flank=flank,socallbychrom=socallbychrom)

		### somatic SNP/INDEL
		if set([1,2,6]).issubset(includes):
			if analysis[6] == 6.1:
			## samtools
				mysamtoolsdir = os.path.join(mysomaticdir,'samtools')
				create_dir(mysamtoolsdir)
				##  somaticsamtoolsMpileup
				mpileup_cmd,order1 = aSomatic.somaticsamtoolsMpileup()
				add_items(orders,order1)
				new_jobs['_'.join(['somaticsamtoolsMpileup',eachsample])] = \
					{'name' : '_'.join(['somaticsamtoolsMpileup',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['somaticsamtoolsMpileup'],
					'cmd' : 'sh '+os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolsMpileup',eachsample])+'.sh')}
				safe_open(os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolsMpileup',eachsample])+'.sh'),'w').write(mpileup_cmd)
				job_points['somaticsamtoolsMpileup'].append('_'.join(['somaticsamtoolsMpileup',eachsample]))
				##  somaticsamtoolscalling
				samtools_cmd,order1 = aSomatic.somaticsamtoolscalling()
				add_items(orders,order1)
				new_jobs['_'.join(['somaticsamtoolscalling',eachsample])] = \
					{'name' : '_'.join(['somaticsamtoolscalling',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['somaticsamtoolscalling'],
					'cmd' : 'sh '+os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolscalling',eachsample])+'.sh')}
				safe_open(os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolscalling',eachsample])+'.sh'),'w').write(samtools_cmd)
				job_points['somaticsamtoolscalling'].append('_'.join(['somaticsamtoolscalling',eachsample]))
				##  somaticsamtoolsMpileup
				annovar_cmd,order1 = aSomatic.somaticsamtoolsannovar()
				add_items(orders,order1)
				new_jobs['_'.join(['somaticsamtoolsannovar',eachsample])] = \
					{'name' : '_'.join(['somaticsamtoolsannovar',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['somaticsamtoolsannovar'],
					'cmd' : 'sh '+os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolsannovar',eachsample])+'.sh')}
				safe_open(os.path.join(mysamtoolsdir,'_'.join(['somaticsamtoolsannovar',eachsample])+'.sh'),'w').write(annovar_cmd)
				job_points['somaticsamtoolsannovar'].append('_'.join(['somaticsamtoolsannovar',eachsample]))
			elif analysis[6] == 6.2:
			## varScan
				sosnp_method = 'varScan'
				soindel_method = 'varScan'
				myvarscandir = os.path.join(mysomaticdir,'varScan')
				create_dir(myvarscandir)
				##  varscan
				varscan_cmd,order1 = aSomatic.varscan()
				add_items(orders,order1)
				new_jobs['_'.join(['varscan',eachsample])] = \
					{'name' : '_'.join(['varscan',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['varscan'],
					'cmd' : 'sh '+os.path.join(myvarscandir,'_'.join(['varscan',eachsample])+'.sh')}
				safe_open(os.path.join(myvarscandir,'_'.join(['varscan',eachsample])+'.sh'),'w').write(varscan_cmd)
				job_points['varscan'].append('_'.join(['varscan',eachsample]))
				##  somaticvarscanannovar
				varscanann_cmd,order1 = aSomatic.somaticvarscanannovar()
				add_items(orders,order1)
				new_jobs['_'.join(['somaticvarscanannovar',eachsample])] = \
					{'name' : '_'.join(['somaticvarscanannovar',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['somaticvarscanannovar'],
					'cmd' : 'sh '+os.path.join(myvarscandir,'_'.join(['somaticvarscanannovar',eachsample])+'.sh')}
				safe_open(os.path.join(myvarscandir,'_'.join(['somaticvarscanannovar',eachsample])+'.sh'),'w').write(varscanann_cmd)
				job_points['somaticvarscanannovar'].append('_'.join(['somaticvarscanannovar',eachsample]))
			elif analysis[6] == 6.3:

			## muTect and Strelka
				sosnp_method = 'muTect'
				soindel_method = 'Strelka'
				mymutectdir = os.path.join(mysomaticdir,'muTect')
				mystrelkadir = os.path.join(mysomaticdir,'Strelka')
				create_dir(mymutectdir)
				create_dir(mystrelkadir)
				##  muTect
				if not argv['socallbychrom']:
					mutect_cmd,order1 = aSomatic.somatic_muTect(FFPE=argv['ffpe'])
					add_items(orders,order1)
					new_jobs['_'.join(['somatic_muTect',eachsample])] = \
						{'name' : '_'.join(['somatic_muTect',eachsample]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s' % queue_list,
						'memory' : memory['somatic_muTect'],
						'cmd' : 'sh '+os.path.join(mymutectdir,'_'.join(['somatic_muTect',eachsample])+'.sh')}
					safe_open(os.path.join(mymutectdir,'_'.join(['somatic_muTect',eachsample])+'.sh'),'w').write(mutect_cmd)
					job_points['somatic_muTect'].append('_'.join(['somatic_muTect',eachsample]))
				if argv['socallbychrom']:
					mutect_chr1to5_cmd,order1 = aSomatic.somatic_muTect_chr1to5(FFPE=argv['ffpe'])
					add_items(orders,order1)
					new_jobs['_'.join(['somatic_muTect_chr1to5',eachsample])] = \
						{'name' : '_'.join(['somatic_muTect_chr1to5',eachsample]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s' % queue_list,
						'memory' : memory['somatic_muTect_chr1to5'],
						'cmd' : 'sh '+os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr1to5',eachsample])+'.sh')}
					safe_open(os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr1to5',eachsample])+'.sh'),'w').write(mutect_chr1to5_cmd)
					job_points['somatic_muTect_chr1to5'].append('_'.join(['somatic_muTect_chr1to5',eachsample]))
					mutect_chr6to12_cmd,order1 = aSomatic.somatic_muTect_chr6to12(FFPE=argv['ffpe'])
					add_items(orders,order1)
					new_jobs['_'.join(['somatic_muTect_chr6to12',eachsample])] = \
						{'name' : '_'.join(['somatic_muTect_chr6to12',eachsample]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s' % queue_list,
						'memory' : memory['somatic_muTect_chr6to12'],
						'cmd' : 'sh '+os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr6to12',eachsample])+'.sh')}
					safe_open(os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr6to12',eachsample])+'.sh'),'w').write(mutect_chr6to12_cmd)
					job_points['somatic_muTect_chr6to12'].append('_'.join(['somatic_muTect_chr6to12',eachsample]))
					mutect_chr13toMT_cmd,order1 = aSomatic.somatic_muTect_chr13toMT(FFPE=argv['ffpe'])
					add_items(orders,order1)
					new_jobs['_'.join(['somatic_muTect_chr13toMT',eachsample])] = \
						{'name' : '_'.join(['somatic_muTect_chr13toMT',eachsample]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s' % queue_list,
						'memory' : memory['somatic_muTect_chr13toMT'],
						'cmd' : 'sh '+os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr13toMT',eachsample])+'.sh')}
					safe_open(os.path.join(mymutectdir,'_'.join(['somatic_muTect_chr13toMT',eachsample])+'.sh'),'w').write(mutect_chr13toMT_cmd)
					job_points['somatic_muTect_chr13toMT'].append('_'.join(['somatic_muTect_chr13toMT',eachsample]))
					mutect_annot_cmd,order1 = aSomatic.somatic_muTect_annot(FFPE=argv['ffpe'])
					add_items(orders,order1)
					new_jobs['_'.join(['somatic_muTect_annot',eachsample])] = \
						{'name' : '_'.join(['somatic_muTect_annot',eachsample]),
						'status' : 'waiting',
						'sched' : '-V -cwd %s' % queue_list,
						'memory' : memory['somatic_muTect_annot'],
						'cmd' : 'sh '+os.path.join(mymutectdir,'_'.join(['somatic_muTect_annot',eachsample])+'.sh')}
					safe_open(os.path.join(mymutectdir,'_'.join(['somatic_muTect_annot',eachsample])+'.sh'),'w').write(mutect_annot_cmd)
					job_points['somatic_muTect_annot'].append('_'.join(['somatic_muTect_annot',eachsample]))
				##  Strelka
				strelka_cmd,order1 = aSomatic.somatic_Strelka()
				add_items(orders,order1)
				new_jobs['_'.join(['somatic_Strelka',eachsample])] = \
					{'name' : '_'.join(['somatic_Strelka',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['somatic_Strelka'],
					'cmd' : 'sh '+os.path.join(mystrelkadir,'_'.join(['somatic_Strelka',eachsample])+'.sh')}
				safe_open(os.path.join(mystrelkadir,'_'.join(['somatic_Strelka',eachsample])+'.sh'),'w').write(strelka_cmd)
				job_points['somatic_Strelka'].append('_'.join(['somatic_Strelka',eachsample]))

		### somatic SV
		if set([1,2,7]).issubset(includes):
			sosv_method = 'breakdancer'
			if analysis[7] == 7.1:
			# breakdancer
				mybreakdancerdir = os.path.join(mysomaticdir,'breakdancer')
				create_dir(mybreakdancerdir)
				# sobam2cfg
				bam2cfg_cmd,order1 = aSomatic.sobam2cfg()
				add_items(orders,order1)
				new_jobs['_'.join(['sobam2cfg',eachsample])] = \
					{'name' : '_'.join(['sobam2cfg',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['sobam2cfg'],
					'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['sobam2cfg',eachsample])+'.sh')}
				safe_open(os.path.join(mybreakdancerdir,'_'.join(['sobam2cfg',eachsample])+'.sh'),'w').write(bam2cfg_cmd)
				job_points['sobam2cfg'].append('_'.join(['sobam2cfg',eachsample]))
				# sobreakdancerSV
				breakdancer_cmd,order1 = aSomatic.sobreakdancerSV()
				add_items(orders,order1)
				new_jobs['_'.join(['sobreakdancerSV',eachsample])] = \
					{'name' : '_'.join(['sobreakdancerSV',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['sobreakdancerSV'],
					'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['sobreakdancerSV',eachsample])+'.sh')}
				safe_open(os.path.join(mybreakdancerdir,'_'.join(['sobreakdancerSV',eachsample])+'.sh'),'w').write(breakdancer_cmd)
				job_points['sobreakdancerSV'].append('_'.join(['sobreakdancerSV',eachsample]))
				# bam2cfg_cmd
				annovar_cmd,order1 = aSomatic.sobreakdancerSvAnno()
				add_items(orders,order1)
				new_jobs['_'.join(['sobreakdancerSvAnno',eachsample])] = \
					{'name' : '_'.join(['sobreakdancerSvAnno',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['sobreakdancerSvAnno'],
					'cmd' : 'sh '+os.path.join(mybreakdancerdir,'_'.join(['sobreakdancerSvAnno',eachsample])+'.sh')}
				safe_open(os.path.join(mybreakdancerdir,'_'.join(['sobreakdancerSvAnno',eachsample])+'.sh'),'w').write(annovar_cmd)
				job_points['sobreakdancerSvAnno'].append('_'.join(['sobreakdancerSvAnno',eachsample]))
			elif analysis[7] == 7.2:
			## crest
				sosv_method = 'crest'
				mycrestdir = os.path.join(mysomaticdir,'crest')
				create_dir(mycrestdir)
				# sclip diff
				diff_cmd,order1 = aSomatic.diff_sclip()
				add_items(orders,order1)
				new_jobs['_'.join(['diff_sclip',eachsample])] = \
					{'name' : '_'.join(['diff_sclip',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['diff_sclip'],
					'cmd' : 'sh '+os.path.join(mycrestdir,'_'.join(['diff_sclip',eachsample])+'.sh')}
				safe_open(os.path.join(mycrestdir,'_'.join(['diff_sclip',eachsample])+'.sh'),'w').write(diff_cmd)
				job_points['diff_sclip'].append('_'.join(['diff_sclip',eachsample]))
				# crest
				crest_cmd,order1 = aSomatic.crest_somaticSV()
				add_items(orders,order1)
				new_jobs['_'.join(['crest_somaticSV',eachsample])] = \
					{'name' : '_'.join(['crest_somaticSV',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=8' % queue_list,
					'memory' : memory['crest_somaticSV'],
					'cmd' : 'sh '+os.path.join(mycrestdir,'_'.join(['crest_somaticSV',eachsample])+'.sh')}
				safe_open(os.path.join(mycrestdir,'_'.join(['crest_somaticSV',eachsample])+'.sh'),'w').write(crest_cmd)
				job_points['crest_somaticSV'].append('_'.join(['crest_somaticSV',eachsample]))
				# annotation
				ann_cmd,order1 = aSomatic.crest_somaticSVann()
				add_items(orders,order1)
				new_jobs['_'.join(['crest_somaticSVann',eachsample])] = \
					{'name' : '_'.join(['crest_somaticSVann',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['crest_somaticSVann'],
					'cmd' : 'sh '+os.path.join(mycrestdir,'_'.join(['crest_somaticSVann',eachsample])+'.sh')}
				safe_open(os.path.join(mycrestdir,'_'.join(['crest_somaticSVann',eachsample])+'.sh'),'w').write(ann_cmd)
			elif analysis[7] == 7.3:
			##delly
				sosv_method = 'delly'
				mydellydir = os.path.join(mysomaticdir,'delly')
				create_dir(mydellydir)
				#so_sv_calling
				delly_sv_del_cmd,order1 = aSomatic.delly_somatic_sv_del()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_del',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_del',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=2' % queue_list,
					'memory' : memory['delly_somatic_sv_del'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_del',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_del',eachsample])+'.sh'),'w').write(delly_sv_del_cmd)
				job_points['delly_somatic_sv_del'].append('_'.join(['delly_somatic_sv_del',eachsample]))
				delly_sv_dup_cmd,order1 = aSomatic.delly_somatic_sv_dup()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_dup',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_dup',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=2' % queue_list,
					'memory' : memory['delly_somatic_sv_dup'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_dup',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_dup',eachsample])+'.sh'),'w').write(delly_sv_dup_cmd)
				job_points['delly_somatic_sv_dup'].append('_'.join(['delly_somatic_sv_dup',eachsample]))
				delly_sv_ins_cmd,order1 = aSomatic.delly_somatic_sv_ins()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_ins',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_ins',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=2' % queue_list,
					'memory' : memory['delly_somatic_sv_ins'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_ins',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_ins',eachsample])+'.sh'),'w').write(delly_sv_ins_cmd)
				job_points['delly_somatic_sv_ins'].append('_'.join(['delly_somatic_sv_ins',eachsample]))
				delly_sv_inv_cmd,order1 = aSomatic.delly_somatic_sv_inv()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_inv',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_inv',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=2' % queue_list,
					'memory' : memory['delly_somatic_sv_inv'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_inv',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_inv',eachsample])+'.sh'),'w').write(delly_sv_inv_cmd)
				job_points['delly_somatic_sv_inv'].append('_'.join(['delly_somatic_sv_inv',eachsample]))
				delly_sv_tra_cmd,order1 = aSomatic.delly_somatic_sv_tra()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_tra',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_tra',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=2' % queue_list,
					'memory' : memory['delly_somatic_sv_tra'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_tra',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_tra',eachsample])+'.sh'),'w').write(delly_sv_tra_cmd)
				job_points['delly_somatic_sv_tra'].append('_'.join(['delly_somatic_sv_tra',eachsample]))
				#annotation
				delly_somatic_sv_annot_cmd,order1 = aSomatic.delly_somatic_sv_annot()
				add_items(orders,order1)
				new_jobs['_'.join(['delly_somatic_sv_annot',eachsample])] = \
					{'name' : '_'.join(['delly_somatic_sv_annot',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['delly_somatic_sv_annot'],
					'cmd' : 'sh '+os.path.join(mydellydir,'_'.join(['delly_somatic_sv_annot',eachsample])+'.sh')}
				safe_open(os.path.join(mydellydir,'_'.join(['delly_somatic_sv_annot',eachsample])+'.sh'),'w').write(delly_somatic_sv_annot_cmd)
				job_points['delly_somatic_sv_annot'].append('_'.join(['delly_somatic_sv_annot',eachsample]))				
			else:
			## lumpy
				sosv_method = 'lumpy'
				mylumpydir = os.path.join(mysomaticdir,'lumpy')
				create_dir(mylumpydir)
				# lumpy
				lumpy_cmd,order1 = aSomatic.lumpy_somaticSV()
				add_items(orders,order1)
				new_jobs['_'.join(['lumpy_somaticSV',eachsample])] = \
					{'name' : '_'.join(['lumpy_somaticSV',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['lumpy_somaticSV'],
					'cmd' : 'sh '+os.path.join(mylumpydir,'_'.join(['lumpy_somaticSV',eachsample])+'.sh')}
				safe_open(os.path.join(mylumpydir,'_'.join(['lumpy_somaticSV',eachsample])+'.sh'),'w').write(lumpy_cmd)
				job_points['lumpy_somaticSV'].append('_'.join(['lumpy_somaticSV',eachsample]))
				# annotation
				ann_cmd,order1 = aSomatic.lumpy_somaticSVann()
				add_items(orders,order1)
				new_jobs['_'.join(['lumpy_somaticSVann',eachsample])] = \
					{'name' : '_'.join(['lumpy_somaticSVann',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['lumpy_somaticSVann'],
					'cmd' : 'sh '+os.path.join(mylumpydir,'_'.join(['lumpy_somaticSVann',eachsample])+'.sh')}
				safe_open(os.path.join(mylumpydir,'_'.join(['lumpy_somaticSVann',eachsample])+'.sh'),'w').write(ann_cmd)
				job_points['lumpy_somaticSVann'].append('_'.join(['lumpy_somaticSVann',eachsample]))

		### somatic CNV
		if set([1,2,8]).issubset(includes):
			socnv_method = 'freec'
			if analysis[8] == 8.1:
				myfreecdir = os.path.join(mysomaticdir,'freec')
				create_dir(myfreecdir)
				freec_cmd,order1 = aSomatic.freec_somaticcnv(seq)
				add_items(orders,order1)
				new_jobs['_'.join(['freec_somaticcnv',eachsample])] = \
					{'name' : '_'.join(['freec_somaticcnv',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['freec_somaticcnv'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freec_somaticcnv',eachsample])+'.sh')}
				safe_open(os.path.join(myfreecdir,'_'.join(['freec_somaticcnv',eachsample])+'.sh'),'w').write(freec_cmd)
				job_points['freec_somaticcnv'].append('_'.join(['freec_somaticcnv',eachsample]))
				# annotation
				#freecann_cmd,order1 = aSomatic.freec_somaticCNVann()
				#add_items(orders,order1)
				#new_jobs['_'.join(['freec_somaticCNVann',eachsample])] = \
				#	{'name' : '_'.join(['freec_somaticCNVann',eachsample]),
				#	'status' : 'waiting',
				#	'sched' : '-V -cwd %s' % queue_list,
				#	'memory' : memory['freec_somaticCNVann'],
				#	'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freec_somaticCNVann',eachsample])+'.sh')}
				#safe_open(os.path.join(myfreecdir,'_'.join(['freec_somaticCNVann',eachsample])+'.sh'),'w').write(freecann_cmd)
				#job_points['freec_somaticCNVann'].append('_'.join(['freec_somaticCNVann',eachsample]))
			elif analysis[8] == 8.4:
				myfreecdir = os.path.join(mysomaticdir,'freec_loh')
				create_dir(myfreecdir)
				freec_cmd,order1 = aSomatic.freec_somaticcnv_loh(seq)
				add_items(orders,order1)
				new_jobs['_'.join(['freec_somaticcnv_loh',eachsample])] = \
					{'name' : '_'.join(['freec_somaticcnv_loh',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['freec_somaticcnv_loh'],
					'cmd' : 'sh '+os.path.join(myfreecdir,'_'.join(['freec_somaticcnv_loh',eachsample])+'.sh')}
				safe_open(os.path.join(myfreecdir,'_'.join(['freec_somaticcnv_loh',eachsample])+'.sh'),'w').write(freec_cmd)
				job_points['freec_somaticcnv_loh'].append('_'.join(['freec_somaticcnv_loh',eachsample]))
			elif analysis[8] == 8.2:
			# ExomeCNV
				socnv_method = 'ExomeCNV'
				myexomecnvdir = os.path.join(mysomaticdir,'ExomeCNV')
				create_dir(myexomecnvdir)
				# exomeCNV_Ndepth
				ndepth_cmd,order1 = aSomatic.exomeCNV_Ndepth()
				add_items(orders,order1)
				new_jobs['_'.join(['exomeCNV_Ndepth',eachsample])] = \
					{'name' : '_'.join(['exomeCNV_Ndepth',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['exomeCNV_Ndepth'],
					'cmd' : 'sh '+os.path.join(myexomecnvdir,'_'.join(['exomeCNV_Ndepth',eachsample])+'.sh')}
				safe_open(os.path.join(myexomecnvdir,'_'.join(['exomeCNV_Ndepth',eachsample])+'.sh'),'w').write(ndepth_cmd)
				job_points['exomeCNV_Ndepth'].append('_'.join(['exomeCNV_Ndepth',eachsample]))
				# exomeCNV_Tdepth
				ndepth_cmd,order1 = aSomatic.exomeCNV_Tdepth()
				add_items(orders,order1)
				new_jobs['_'.join(['exomeCNV_Tdepth',eachsample])] = \
					{'name' : '_'.join(['exomeCNV_Tdepth',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['exomeCNV_Tdepth'],
					'cmd' : 'sh '+os.path.join(myexomecnvdir,'_'.join(['exomeCNV_Tdepth',eachsample])+'.sh')}
				safe_open(os.path.join(myexomecnvdir,'_'.join(['exomeCNV_Tdepth',eachsample])+'.sh'),'w').write(ndepth_cmd)
				job_points['exomeCNV_Tdepth'].append('_'.join(['exomeCNV_Tdepth',eachsample]))
				# exomeCNV_CNV
				exomecnv_cmd,order1 = aSomatic.exomeCNV_CNV()
				add_items(orders,order1)
				new_jobs['_'.join(['exomeCNV_CNV',eachsample])] = \
					{'name' : '_'.join(['exomeCNV_CNV',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['exomeCNV_CNV'],
					'cmd' : 'sh '+os.path.join(myexomecnvdir,'_'.join(['exomeCNV_CNV',eachsample])+'.sh')}
				safe_open(os.path.join(myexomecnvdir,'_'.join(['exomeCNV_CNV',eachsample])+'.sh'),'w').write(exomecnv_cmd)
				job_points['exomeCNV_CNV'].append('_'.join(['exomeCNV_CNV',eachsample]))
			else:
				myvarscandir = os.path.join(mysomaticdir,'varScanCNV')
				create_dir(myvarscandir)
				varscan_cmd,order1 = aSomatic.varscan_somaticCNV()
				add_items(orders,order1)
				new_jobs['_'.join(['varscan_somaticCNV',eachsample])] = \
					{'name' : '_'.join(['varscan_somaticCNV',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['varscan_somaticCNV'],
					'cmd' : 'sh '+os.path.join(myvarscandir,'_'.join(['varscan_somaticCNV',eachsample])+'.sh')}
				safe_open(os.path.join(myvarscandir,'_'.join(['varscan_somaticCNV',eachsample])+'.sh'),'w').write(varscan_cmd)
				job_points['varscan_somaticCNV'].append('_'.join(['varscan_somaticCNV',eachsample]))
				# annotation
				varscanann_cmd,order1 = aSomatic.varscan_somaticCNVann()
				add_items(orders,order1)
				new_jobs['_'.join(['varscan_somaticCNVann',eachsample])] = \
					{'name' : '_'.join(['varscan_somaticCNVann',eachsample]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s' % queue_list,
					'memory' : memory['varscan_somaticCNVann'],
					'cmd' : 'sh '+os.path.join(myvarscandir,'_'.join(['varscan_somaticCNVann',eachsample])+'.sh')}
				safe_open(os.path.join(myvarscandir,'_'.join(['varscan_somaticCNVann',eachsample])+'.sh'),'w').write(varscanann_cmd)
				job_points['varscan_somaticCNVann'].append('_'.join(['varscan_somaticCNVann',eachsample]))
	## job status
	if lastjobs:
		for eachsample in list_in_qc:
			### somatic SNP/INDEL
			if set([1,2,6]).issubset(includes):
				if analysis[6] == 6.1:
				# samtools
					mpileupjob = '_'.join(['somaticsamtoolsMpileup',eachsample])
					if mpileupjob in lastjobs and mpileupjob not in new_jobs:
						add_jobs[mpileupjob] = lastjobs[mpileupjob]['txt']
						somatic_samlist.append(eachsample)  # somatic sample list
					samtoolsjob = '_'.join(['somaticsamtoolscalling',eachsample])
					if samtoolsjob in lastjobs and samtoolsjob not in new_jobs:
						add_jobs[samtoolsjob] = lastjobs[samtoolsjob]['txt']
					annovarjob = '_'.join(['somaticsamtoolsannovar',eachsample])
					if annovarjob in lastjobs and annovarjob not in new_jobs:
						add_jobs[annovarjob] = lastjobs[annovarjob]['txt']
				elif analysis[6] == 6.2:
				# varScan
					varscanjob = '_'.join(['varscan',eachsample])
					if varscanjob in lastjobs and varscanjob not in new_jobs:
						add_jobs[varscanjob] = lastjobs[varscanjob]['txt']
						somatic_samlist.append(eachsample)  # somatic sample list
					varscanannjob = '_'.join(['somaticvarscanannovar',eachsample])
					if varscanannjob in lastjobs and varscanannjob not in new_jobs:
						add_jobs[varscanannjob] = lastjobs[varscanannjob]['txt']
				elif analysis[6] == 6.3:
				# muTect and Strelka
					if argv['socallbychrom']:
						mutectjob_1to5 = '_'.join(['somatic_muTect_chr1to5',eachsample])
						if mutectjob_1to5 in lastjobs and mutectjob_1to5 not in new_jobs:
							add_jobs[mutectjob_1to5] = lastjobs[mutectjob_1to5]['txt']
							somatic_samlist.append(eachsample)  # somatic sample list
						mutectjob_6to12 = '_'.join(['somatic_muTect_chr6to12',eachsample])
						if mutectjob_6to12 in lastjobs and mutectjob_6to12 not in new_jobs:
							add_jobs[mutectjob_6to12] = lastjobs[mutectjob_6to12]['txt']
						mutectjob_13toMT = '_'.join(['somatic_muTect_chr13toMT',eachsample])
						if mutectjob_13toMT in lastjobs and mutectjob_13toMT not in new_jobs:
							add_jobs[mutectjob_13toMT] = lastjobs[mutectjob_13toMT]['txt']
						somatic_muTect_annot_job = '_'.join(['somatic_muTect_annot',eachsample])
						if somatic_muTect_annot_job in lastjobs and somatic_muTect_annot_job not in new_jobs:
							add_jobs[somatic_muTect_annot_job] = lastjobs[somatic_muTect_annot_job]['txt']
					else:
						mutectjob = '_'.join(['somatic_muTect',eachsample])
						if mutectjob in lastjobs and mutectjob not in new_jobs:
							add_jobs[mutectjob] = lastjobs[mutectjob]['txt']
							somatic_samlist.append(eachsample)  # somatic sample list
					strelkajob = '_'.join(['somatic_Strelka',eachsample])
					if strelkajob in lastjobs and strelkajob not in new_jobs:
						add_jobs[strelkajob] = lastjobs[strelkajob]['txt']
			### somatic SV
			if set([1,2,7]).issubset(includes):
				if analysis[7] == 7.1:
					cfgjob = '_'.join(['sobam2cfg',eachsample])
					if cfgjob in lastjobs and cfgjob not in new_jobs:
						add_jobs[cfgjob] = lastjobs[cfgjob]['txt']
					breakdancerjob = '_'.join(['sobreakdancerSV',eachsample])
					if breakdancerjob in lastjobs and breakdancerjob not in new_jobs:
						add_jobs[breakdancerjob] = lastjobs[breakdancerjob]['txt']
					annovarjob = '_'.join(['sobreakdancerSvAnno',eachsample])
					if annovarjob in lastjobs and annovarjob not in new_jobs:
						add_jobs[annovarjob] = lastjobs[annovarjob]['txt']
				elif analysis[7] == 7.2:
					diffjob = '_'.join(['diff_sclip',eachsample])
					if diffjob in lastjobs and diffjob not in new_jobs:
						add_jobs[diffjob] = lastjobs[diffjob]['txt']
					crestjob = '_'.join(['crest_somaticSV',eachsample])
					if crestjob in lastjobs and crestjob not in new_jobs:
						add_jobs[crestjob] = lastjobs[crestjob]['txt']
					crestannjob = '_'.join(['crest_somaticSVann',eachsample])
					if crestannjob in lastjobs and crestannjob not in new_jobs:
						add_jobs[crestannjob] = lastjobs[crestannjob]['txt']
				elif analysis[7] == 7.3:
					delly_somatic_sv_del_job = '_'.join(['delly_somatic_sv_del',eachsample])
					if delly_somatic_sv_del_job in lastjobs and delly_somatic_sv_del_job not in new_jobs:
						add_jobs[delly_somatic_sv_del_job] = lastjobs[delly_somatic_sv_del_job]['txt']
					delly_somatic_sv_dup_job = '_'.join(['delly_somatic_sv_dup',eachsample])
					if delly_somatic_sv_dup_job in lastjobs and delly_somatic_sv_dup_job not in new_jobs:
						add_jobs[delly_somatic_sv_dup_job] = lastjobs[delly_somatic_sv_dup_job]['txt']
					delly_somatic_sv_ins_job = '_'.join(['delly_somatic_sv_ins',eachsample])
					if delly_somatic_sv_ins_job in lastjobs and delly_somatic_sv_ins_job not in new_jobs:
						add_jobs[delly_somatic_sv_ins_job] = lastjobs[delly_somatic_sv_ins_job]['txt']
					delly_somatic_sv_inv_job = '_'.join(['delly_somatic_sv_inv',eachsample])
					if delly_somatic_sv_inv_job in lastjobs and delly_somatic_sv_inv_job not in new_jobs:
						add_jobs[delly_somatic_sv_inv_job] = lastjobs[delly_somatic_sv_inv_job]['txt']
					delly_somatic_sv_tra_job = '_'.join(['delly_somatic_sv_tra',eachsample])
					if delly_somatic_sv_tra_job in lastjobs and delly_somatic_sv_tra_job not in new_jobs:
						add_jobs[delly_somatic_sv_tra_job] = lastjobs[delly_somatic_sv_tra_job]['txt']
					delly_somatic_sv_annot_job = '_'.join(['delly_somatic_sv_annot',eachsample])
					if delly_somatic_sv_annot_job in lastjobs and delly_somatic_sv_annot_job not in new_jobs:
						add_jobs[delly_somatic_sv_annot_job] = lastjobs[delly_somatic_sv_annot_job]['txt']						
				else:
					lumpyjob = '_'.join(['lumpy_somaticSV',eachsample])
					if lumpyjob in lastjobs and lumpyjob not in new_jobs:
						add_jobs[lumpyjob] = lastjobs[lumpyjob]['txt']
					lumpyannjob = '_'.join(['lumpy_somaticSVann',eachsample])
					if lumpyannjob in lastjobs and lumpyannjob not in new_jobs:
						add_jobs[lumpyannjob] = lastjobs[lumpyannjob]['txt']
			
			### somatic CNV
			if set([1,2,8]).issubset(includes):
				if analysis[8] == 8.1:
					freecjob = '_'.join(['freec_somaticcnv',eachsample])
					if freecjob in lastjobs and freecjob not in new_jobs:
						add_jobs[freecjob] = lastjobs[freecjob]['txt']
					#freecannjob = '_'.join(['freec_somaticCNVann',eachsample])
					#if freecannjob in lastjobs and freecannjob not in new_jobs:
					#	add_jobs[freecannjob] = lastjobs[freecannjob]['txt']
				elif analysis[8] == 8.4:
					freecjob = '_'.join(['freec_somaticcnv_loh',eachsample])
					if freecjob in lastjobs and freecjob not in new_jobs:
						add_jobs[freecjob] = lastjobs[freecjob]['txt']
				elif analysis[8] == 8.2:
					ndepthjob = '_'.join(['exomeCNV_Ndepth',eachsample])
					if ndepthjob in lastjobs and ndepthjob not in new_jobs:
						add_jobs[ndepthjob] = lastjobs[ndepthjob]['txt']
					tdepthjob = '_'.join(['exomeCNV_Tdepth',eachsample])
					if tdepthjob in lastjobs and tdepthjob not in new_jobs:
						add_jobs[tdepthjob] = lastjobs[tdepthjob]['txt']
					exomecnvjob = '_'.join(['exomeCNV_CNV',eachsample])
					if exomecnvjob in lastjobs and exomecnvjob not in new_jobs:
						add_jobs[exomecnvjob] = lastjobs[exomecnvjob]['txt']
				else:
					varscanjob = '_'.join(['varscan_somaticCNV',eachsample])
					if varscanjob in lastjobs and varscanjob not in new_jobs:
						add_jobs[varscanjob] = lastjobs[varscanjob]['txt']
					varscanannjob = '_'.join(['varscan_somaticCNVann',eachsample])
					if varscanannjob in lastjobs and varscanannjob not in new_jobs:
						add_jobs[varscanannjob] = lastjobs[varscanannjob]['txt']
					
	### spetrum
	print "Spetrum ...\n"
	if set([1,2,6]).issubset(includes):
		create_dir(spectrumdir) 
#		if genome_version =='mm9':
#			spectrum_cmd,order1 = aSomatic.spectrum_m(somatic_samlist,spectrumdir,analydir,sosnp_method,soindel_method,sample_order)
#		else:
		spectrum_cmd,order1 = aSomatic.spectrum(somatic_samlist,spectrumdir,analydir,sosnp_method,soindel_method,sample_order)
		add_items(orders,order1)
		new_jobs['spectrum'] = \
			{'name' : 'spectrum',
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : memory['spectrum'],
			'cmd' : 'sh '+os.path.join(spectrumdir,'spectrum.sh')}
		safe_open(os.path.join(spectrumdir,'spectrum.sh'),'w').write(spectrum_cmd)
		job_points['spectrum'] = ['spectrum']

create_dir(reportdir)

exclude_jobs = ['primary_report','result']
####byebye
with safe_open(os.path.join(analydir, 'byebye.sh'),'w') as shell:
    shell.write('sh '+'/PUBLIC/source/HW/CANCER/HW_v2.1/byebye.sh '+analydir)
shell.close
## primary report
print "Report and Results ...\n"
#if set([1,2]).issubset(includes): #20180428 bysun
if set([1,2,3]).issubset(includes) or set([1,2,4]).issubset(includes) or set([1,2,5]).issubset(includes) or set([1,2,6]).issubset(includes) or set([1,2,7]).issubset(includes) or set([1,2,8]).issubset(includes):
	primaryreportdir = os.path.join(reportdir,'primary')
	create_dir(primaryreportdir)
	new_jobs['primary_report'] = {'name' : 'primary_report',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : memory['primary_report'],
		'cmd' : 'sh '+os.path.join(primaryreportdir,'primary_report.sh')}
	shell = safe_open(os.path.join(primaryreportdir,'primary_report.sh'),'w')
	if TumorGermline:
		shell.writelines('python %s --TumorGermline --project %s --projdir %s --seqsty %s --analy_array %s --qclist %s --outdir %s --reptype primary --samplelist %s --ref %s --PCRFree %s --jobname %s' \
		% (report_script, os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),primaryreportdir,os.path.abspath(argv['infile']),ref,pcrfree,newjob))
	else:
		shell.writelines('python %s --project %s --projdir %s --seqsty %s --analy_array %s --qclist %s --outdir %s --reptype primary --samplelist %s --ref %s --PCRFree %s --jobname %s' \
		% (report_script, os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),primaryreportdir,os.path.abspath(argv['infile']),ref,pcrfree,newjob))
	shell.writelines('\n')
	shell.close()
	## order
	for eachjob in new_jobs:
		if eachjob in exclude_jobs:
			continue
		order1 = 'order primary_report after %s' % eachjob
		add_items(orders,order1)
	for eachjob in add_jobs:
		if eachjob in exclude_jobs:
			continue
		order1 = 'order primary_report after %s' % eachjob
		add_items(orders,order1)
        order1 = 'order result after primary_report'   #20180423
        add_items(orders,order1)   

## PROJECT INFO backup
print "Project information BACKUP ...\n"
backupdir = backup_dir['WES']
if 'WGS' in seqstrag:
	backupdir = backup_dir['WGS']
if not somatic:
	backupdir = backup_dir['NonSomatic']
script = '''
python %s \\
	--pwd %s \\
	--qclist_info %s \\
	--seqsty %s --date %s \\
	--odir %s
''' % (backup_script,analydir,os.path.join(analydir,'qc_list'),seqstrag,time.strftime('%m.%d',time.localtime()),backupdir)
open(os.path.join(analydir,'backup_project.sh'),'w').write(script)

## Result
create_dir(resultdir)
new_jobs['result'] = {'name' : 'result',
	'status' : 'waiting',
	'sched' : '-V -cwd %s' % queue_list,
	'memory' : memory['result'],
	'cmd' : 'sh '+os.path.join(resultdir,'result_'+newjob+'.sh')}
shell = safe_open(os.path.join(resultdir,'result_'+newjob+'.sh'),'w')
create_dir(os.path.join(resultdir,newjob))
create_dir(os.path.join(resultdir,newjob,'release'))
create_dir(os.path.join(resultdir,newjob,'release_noclean'))
if argv['adv_array_cancer']:  #by sun 20180330
       assert not os.system('rm -rf '+os.path.join(resultdir,newjob)+' '+os.path.join(resultdir,'result_'+newjob+'.sh')) 
if argv['adv_array_disease']:
       assert not os.system('rm -rf '+os.path.join(resultdir,newjob)+' '+os.path.join(resultdir,'result_'+newjob+'.sh')) 

if TumorGermline:
	shell.writelines('cd %s\npython %s --TumorGermline --projdir %s --type %s --analy_array %s --qclist %s --outdir %s --readmedir %s  --ref %s' \
	% (os.path.join(resultdir,newjob),release_script,analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),os.path.join(resultdir,newjob),READMEdir,ref))
else:
	shell.writelines('cd %s\npython %s --projdir %s --type %s --analy_array %s --qclist %s --outdir %s --readmedir %s --ref %s' \
	% (os.path.join(resultdir,newjob),release_script,analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),os.path.join(resultdir,newjob),READMEdir,ref))
shell.writelines(' && \\\n')
shell.writelines('tar -chzvf %s/Primary_analysis_result.tar.gz Primary_analysis_result && \\\n' % (os.path.join(resultdir,newjob)))
shell.writelines('md5sum Primary_analysis_result.tar.gz >>MD5.txt && \\\n')
shell.writelines('cp %s %s && \\\n' % (data_relase_READMEFILE,os.path.join(resultdir,newjob)))
shell.writelines('ln -sf %s/{Bam,RawData,CleanData,MD5.txt,README.txt,Primary_analysis_result.tar.gz} release && \\\n' % os.path.join(resultdir,newjob))
shell.writelines('ln -sf %s/{Bam,RawData,Primary_analysis_result.tar.gz} release_noclean && \\\n' % os.path.join(resultdir,newjob))
if set([1,2,3]).issubset(includes) or set([1,2,4]).issubset(includes) or set([1,2,5]).issubset(includes) or set([1,2,6]).issubset(includes) or set([1,2,7]).issubset(includes) or set([1,2,8]).issubset(includes):
        shell.writelines('md5sum *primary_*_report.tar.gz >>MD5.txt && \\\n')
        shell.writelines('ln -sf %s/*primary_*_report.tar.gz release && \\\n' % os.path.join(resultdir,newjob))
        shell.writelines('ln -sf %s/*primary_*_report.tar.gz release_noclean && \\\n' % os.path.join(resultdir,newjob))
#shell.writelines('ln -sf %s/{Bam,RawData,CleanData,MD5.txt,README.txt,Primary_analysis_result.tar.gz,*primary_*_report.tar.gz} release && \\\n' % os.path.join(resultdir,newjob))
#shell.writelines('ln -sf %s/{Bam,RawData,Primary_analysis_result.tar.gz,*primary_*_report.tar.gz} release_noclean && \\\n' % os.path.join(resultdir,newjob))
shell.writelines('cp %s/MD5.txt release_noclean && \\\n' % os.path.join(resultdir,newjob) + 'sed -i \'/  CleanData/d\' release_noclean/MD5.txt && \\\n' + 'cp /PUBLIC/database/HW/CANCER/HW_v2.0/Module/Release/HW_v2.2/README_noclean.txt release_noclean/README.txt')
#if not 'TR' in seqstrag:
#	shell.writelines(script)
shell.close()

## order
for eachjob in new_jobs:
	if eachjob in exclude_jobs:
		continue
	order1 = 'order result after %s' % eachjob
	add_items(orders,order1)
for eachjob in add_jobs:
	if eachjob in exclude_jobs:
		continue
	order1 = 'order result after %s ' % eachjob
	add_items(orders,order1)


## qc
if set([1]).issubset(includes):
	qcreportdir = os.path.join(reportdir,'qc')
	create_dir(qcreportdir)
	for eachjob in new_jobs:
		if eachjob.startswith('qc_'):
			#order1 = 'order %s before qc_report' % eachjob
			order1 = 'order qc_report after %s' % eachjob
			add_items(orders,order1)
		if eachjob.startswith('pollution_'):
			order1 = 'order qc_report after %s' % eachjob
			add_items(orders,order1)
	for eachjob in add_jobs:
		if eachjob.startswith('qc_'):
			#order1 = 'order %s before qc_report' % eachjob
			order1 = 'order qc_report after %s' % eachjob
			add_items(orders,order1)
		if eachjob.startswith('pollution_'):
			order1 = 'order qc_report after %s' % eachjob
			add_items(orders,order1)
	new_jobs['qc_report'] = {'name' : 'qc_report',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : memory['qc_report'],
		'cmd' : 'sh '+os.path.join(qcreportdir,'qc_report.sh')}
	shell = safe_open(os.path.join(qcreportdir,'qc_report.sh'),'w')
	shell.writelines('python %s --project %s --projdir %s --seqsty %s --analy_array %s --qclist %s --outdir %s --reptype qc --samplelist %s --ref %s --PCRFree %s' % \
	(report_script,os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),qcreportdir,os.path.abspath(argv['infile']),ref,pcrfree))
	if argv['pollution']:
		shell.writelines(' --pollution\n')
#	if argv['english']:
#		shell.writelines(' --english')
	shell.writelines('\n')
	shell.close()

## mapping
if set([1,2]).issubset(includes):
	mappingreportdir = os.path.join(reportdir,'mapping')
	create_dir(mappingreportdir)
	for eachjob in new_jobs:
		if eachjob.startswith('combine_'):
			#order1 = 'order %s before mapping_report' % eachjob
			order1 = 'order mapping_report after %s' % eachjob
			add_items(orders,order1)
	for eachjob in add_jobs:
		if eachjob.startswith('combine_'):
			#order1 = 'order %s before mapping_report' % eachjob
			order1 = 'order mapping_report after %s' % eachjob
			add_items(orders,order1)
	new_jobs['mapping_report'] = {'name' : 'mapping_report',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : memory['mapping_report'],
		'cmd' : 'sh '+os.path.join(mappingreportdir,'mapping_report.sh')}
	shell = safe_open(os.path.join(mappingreportdir,'mapping_report.sh'),'w')
	shell.writelines('python %s --project %s --projdir %s --seqsty %s --analy_array %s --qclist %s --outdir %s --reptype mapping --samplelist %s --ref %s --PCRFree %s' % \
	(report_script,os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),mappingreportdir,os.path.abspath(argv['infile']),ref,pcrfree))
	if argv['pollution']:
		shell.writelines(' --pollution\n')
#	if argv['english']:
#		shell.writelines(' --english')
	shell.writelines('\n')
	shell.close()
	## result after mapping_report
	add_items(orders,'order result after mapping_report')


## get recursively downstream jobs
def get_succeeds (job):
	recursiveSucceeds = [job]
	succeeds = []
	if order_relation.has_key(job):
		succeeds = order_relation[job]
	if len(succeeds) >0:
		for each in succeeds:
			recursiveSucceeds.extend(get_succeeds(each))
	return recursiveSucceeds
###  for startpoint
if startpoint:
	print "Startpoint ..."
	print '   ... from %s\n' % ', '.join(startpoints)
	order_relation = {}  
	##  relationship; jobA before jobB  ==>  order_relation[jobA] = jobB;
	for one in orders:
		## order ln_EC1001_N_DHE00213_C4D7FACXX_L6 before qc_EC1001_N_DHE00213_C4D7FACXX_L6
		if (not re.search('before',one)) and (not re.search('after',one)): continue
		flag,jobA,relation,jobB = one.strip().split()
		if relation == 'before':
			if jobA not in order_relation:
				order_relation[jobA] = []
			order_relation[jobA].append(jobB)
		else:
			if jobB not in order_relation:
				order_relation[jobB] = []
			order_relation[jobB].append(jobA)
	## get the succeed jobs in startpoints list
	point_succeeds = set()
	for eachstart in startpoints:
		if eachstart not in job_points:
			print 'POINT '+eachstart+' not in you analysis, Please check your script!\nAvailable POINTS are:\n'+'  \n'.join(job_points.keys())+'\n'
			sys.exit(1)
		for each in job_points[eachstart]:
			tmp_succeeds = get_succeeds(each)
			point_succeeds |= set(tmp_succeeds)
	## get undone jobs in add_jobs dependent
	undone_depend_jobs = set()
	for eachjob in add_jobs:
		if not lastjobs[eachjob]['status'] == 'done':
			undone_depend = get_succeeds(eachjob)
			undone_depend_jobs |= undone_depend
	for eachjob in undone_depend_jobs:
		if eachjob in lastjobs:
			## undone depended jobs must be not 'done'
			assert not lastjobs[eachjob]['status'] == 'done'
	## change job status
	for eachjob in new_jobs:
		if (eachjob not in point_succeeds) and (eachjob not in undone_depend_jobs):
			new_jobs[eachjob]['status'] = 'done'
	## finalbam startpoint and mpileup, when mpileup run in mapping step
	if 'finalbam' in startpoints and 'mpileup' not in startpoints and raw_mpileup:
		if 'mpileup' not in job_points:
			print 'Startpoints finalbam wrong!'
			sys.exit(1)
		for each in job_points['mpileup']:
			new_jobs[each]['status'] = 'done'
		for each in job_points['finalbam']:
			new_jobs[each]['status'] = 'done'
	
## jobfile
if not argv['adv_array_cancer']:
	jobfile = safe_open(os.path.join(analydir,newjob),'w')
	jobfile.write('log_dir %s\n' % logdir)
	for eachjob in new_jobs:
		sjmcmd = '\n  '.join(['  name %s' % new_jobs[eachjob]['name'],
			'memory %s' % new_jobs[eachjob]['memory'],
			'status %s' % new_jobs[eachjob]['status'],
			'sched_options %s' % new_jobs[eachjob]['sched'],
			'cmd_begin',
			'  %s' % new_jobs[eachjob]['cmd'],
			'cmd_end'])
		jobfile.write('job_begin\n'+sjmcmd+'\njob_end\n')
	for eachjob in add_jobs:
		jobfile.write(add_jobs[eachjob]+'\n')
	## orders
	jobfile.write('\n'.join(orders)+'\n')
	# lastorders
	for each in lastorders:
		flag,jobA,relation,jobB = each.strip().split()
		if (jobA in new_jobs or jobA in add_jobs) and (jobB in new_jobs or jobB in add_jobs):
			if each not in orders:
				jobfile.write(each+'\n')
	jobfile.close()

## remove mpileup files
if set([5,8]).issubset(includes):
#	sorted_bams = [os.path.join(mapdir,eachsample,eachsample+'.merge.bam') for eachsample in list_in_qc]
#	sorted_bais = [os.path.join(mapdir,eachsample,eachsample+'.merge.bam.bai') for eachsample in list_in_qc]
#	rm_cmd = '\n'.join(['rm -rf '+each for each in sorted_bams])+'\n'
#	rm_cmd += '\n'.join(['rm -rf '+each for each in sorted_bais])+'\n'
	print "remove_mpileupgz.sh generated..."
	pilupgz = [os.path.join(mapdir,eachsample,eachsample+'.final.mpileup.gz') for eachsample in list_in_qc]
	rm_cmd = '\n'.join(['rm -rf '+each for each in pilupgz])+'\n'
	open(os.path.join(analydir,'remove_mpileupgz.sh'),'w').write(rm_cmd)

################################################################################################
##########################        ADVANCE ANALYSIS               ###############################
################################################################################################
#give a value to analysis steps
advance_cancer = {'1'   :  'predispose_filter',
		'2'   :  'signature_spectrum',
		'3'   :  'drivergenes_filter',
		'4'   :  'smg',
		'5'   :  'gistic',
		'6'   :  'fusion_gene',
		'7'   :  'absolutes',
		'8.1' :  'pyclones',
		'8.2' :  'sciclones',
		'9'   :  'phylip_evolution',
		'11'  :  'mrt_music',
		'12'  :  'oncodriveclust',
		'15'  :  'noncoding_filter',
                '17'  :  'mutationshow'}

adv_relation_cancer = {'1'  : [1,2,3],
		'2'  : [1,2,6],
		'3'  : [1,2,6],
		'4'  : [1,2,6],
		'5'  : [1,2,8],
		'6'  : [1,2,7],
		'7'  : [1,2,6,8],
		'8.1': [1,2,6,8],
		'8.2': [1,2,6,8],
		'9'  : [1,2,6],
		'11' : [1,2,6],
		'12' : [1,2,6],
		'15' : [1,2,6],
                '17' : [1,2,6]}
#Note: advanced analysis '11' and '12' depends on adv '4'

adv_relation_disease ={6:[1,2,3],
               6.1:[1,2,3],
               7:[1,2,3,6],
               7.1:[1,2,3,6],
               7.2:[1,2,3,6],
               8:[1,2],
               8.1:[1,2],
               8.2:[1,2],
               8.3:[1,2,3],
               8.4:[1,2,3],
               8.5:[1,2,4],
               8.6:[1,2,5],
               9:[1,2],
               9.1:[1,2],
               10:[1,2],
               10.1:[1,2,3]}  #by sun

adv_array_cancer = []
adv_symbol = []
if argv['adv_array_cancer']:
	adv_array_cancer = argv['adv_array_cancer'].strip().split(',')
	adv_symbol = [symbol_code[each] for each in analy_array]

adv_startpoint = []
if argv['adv_array_cancer'] and argv['startpoint']:
	adv_startpoints = [str(x) for x in argv['startpoint'].strip().split(',')]

will_quit = False
for i in adv_array_cancer:
	s = set(adv_relation_cancer[i])
	if not s.issubset(set(analys)):
		print 'you need do %s before %f' %(', '.join(adv_relation_cancer[i]), i)
		will_quit = True
	else:pass
if will_quit:
	sys.exit(0)

adv_array_disease = [] #by sun
if argv['adv_array_disease']:
        adv_array_disease = argv['adv_array_disease'].strip().split(',')
        adv_symbol = [symbol_code[each] for each in analy_array]

### jobstatus file
lastorders2 = []
lastjobs2 = {}
if jobstat and argv['adv_array_cancer']:
	logdir,lastjobs2,lastorders2 = parse_jobstatus(jobstat)

prepareDir = os.path.join(advdir_cancer,'00.prepare')
predisposeDir = os.path.join(advdir_cancer,'PredisposeGenes')
signatureDir = os.path.join(advdir_cancer,'Signatures')
driverGeneDir = os.path.join(advdir_cancer,'DriverGenes')
#driverDir = os.path.join(advdir_cancer,'DriverGenes')
musicDir = os.path.join(advdir_cancer,'Smg')
gisticDir = os.path.join(advdir_cancer,'Gistic')
fusionDir = os.path.join(advdir_cancer,'FusionGenes')
absoluteDir = os.path.join(advdir_cancer,'Absolute')
cloneDir = os.path.join(advdir_cancer,'Clone')
pycloneDir = os.path.join(cloneDir,'pyclone')
scicloneDir = os.path.join(cloneDir,'sciclone')
expandsDir = os.path.join(cloneDir,'expands')
evolutionDir = os.path.join(advdir_cancer,'Evolution')
circosDir = os.path.join(advdir_cancer,'Circos')
mrtDir = os.path.join(advdir_cancer,'Mrt')
oncodriveDir = os.path.join(advdir_cancer,'OncodriveClust')
targetDir = os.path.join(advdir_cancer,'DrugTarget')
resistenceDir = os.path.join(advdir_cancer,'DrugResistance')
noncodingDir = os.path.join(advdir_cancer,'Noncoding')
bubbletreeDir = os.path.join(advdir_cancer,'Bubbletree')
mutshowDir = os.path.join(advdir_cancer,'Mut_Onshow')

new_jobs2 = {}
orders2 = []
add_jobs2 = {}  ## for job status
job_points2 = {} ## for startpoints


Tsamples = [each for each in somatic_samlist if each in sample_order]
T2pid = dict([(each,sam2pid[each][0]) for each in somatic_samlist])
T2Nsample = dict([(each,patientInfo[sam2pid[each][0]]['N']) for each in somatic_samlist])
## __init__(self,Tsamples,T2pid,T2Nsample,advbin,advdir_cancer,mapDir,mutDir,svDir,somaticDir,mutM,svM,cnvM,sosnpM,soindelM,sosvM,socnvM,seqstrag,genome_info)
thisAdv = AdvAnalysis(Tsamples, Nsamples, T2pid, T2Nsample, softwares['advbin'], advdir_cancer, mapdir, mutdir, svdir, somaticdir, mut_method, sv_method, cnv_method, sosnp_method, soindel_method, sosv_method, socnv_method, seqstrag, genome_files[argv['genome']])

### prepare dir
if argv['adv_array_cancer']:
	create_dir(advdir_cancer)
	create_dir(prepareDir)
	combine_maf_cmd = thisAdv.combine_maf()
	new_jobs2['combine_maf'] = \
		{'name' : 'combine_maf',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['combine_maf']),
		 'memory' : memory['combine_maf'],
		 'cmd' : 'sh '+os.path.join(prepareDir,'combine_maf.sh')}
	safe_open(os.path.join(prepareDir,'combine_maf.sh'),'w').write(combine_maf_cmd)
	job_points2['combine_maf'] = ['combine_maf']
	if lastjobs2:
		if 'combine_maf' not in new_jobs2 and 'combine_maf' in lastjobs2:
			add_jobs2['combine_maf'] = lastjobs2['combine_maf']['txt']

### predispose
if set(['1']).issubset(adv_array_cancer):
	print 'predisposing genes ...'
	create_dir(predisposeDir)
	## predispose gene
	predispose_filter_cmd = thisAdv.predispose_filter()
	new_jobs2['predispose_filter'] = \
		{'name' : 'predispose_filter',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['predispose_filter']),
		 'memory' : memory['predispose_filter'],
		 'cmd' : 'sh '+os.path.join(predisposeDir,'predispose_filter.sh')}
	safe_open(os.path.join(predisposeDir,'predispose_filter.sh'),'w').write(predispose_filter_cmd)
	job_points2['predispose_filter'] = ['predispose_filter']
	if lastjobs2:
		if 'predispose_filter' not in new_jobs2 and 'predispose_filter' in lastjobs2:
			add_jobs2['predispose_filter'] = lastjobs2['predispose_filter']['txt']
## Mutation signature and spectrum
if set(['2']).issubset(adv_array_cancer):
	print 'Mutational signature and spectrum ...'
	create_dir(signatureDir)
	signature_spectrum_cmd,order1 = thisAdv.signature_spectrum()
	add_items(orders2,order1)
	new_jobs2['signature_spectrum'] = \
		{'name' : 'signature_spectrum',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['signature_spectrum']),
		 'memory' : memory['signature_spectrum'],
		 'cmd' : 'sh '+os.path.join(signatureDir,'signature_spectrum.sh')}
	safe_open(os.path.join(signatureDir,'signature_spectrum.sh'),'w').write(signature_spectrum_cmd)
	job_points2['signature_spectrum'] = ['signature_spectrum']
	if lastjobs2:
		if 'signature_spectrum' not in new_jobs2 and 'signature_spectrum' in lastjobs2:
			add_jobs2['signature_spectrum'] = lastjobs2['signature_spectrum']['txt']
	
## Cancer Gene
if set(['3']).issubset(adv_array_cancer):
	print 'known driver genes ...'
	## oncogene, tsg
	create_dir(driverGeneDir)
	drivergenes_filter_cmd = thisAdv.drivergenes_filter()
	new_jobs2['drivergenes_filter'] = \
		{'name' : 'drivergenes_filter',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['drivergenes_filter']),
		 'memory' : memory['drivergenes_filter'],
		 'cmd' : 'sh '+os.path.join(driverGeneDir,'drivergenes_filter.sh')}
	safe_open(os.path.join(driverGeneDir,'drivergenes_filter.sh'),'w').write(drivergenes_filter_cmd)
	job_points2['drivergenes_filter'] = ['drivergenes_filter']
	if lastjobs2:
		if 'drivergenes_filter' not in new_jobs2 and 'drivergenes_filter' in lastjobs2:
			add_jobs2['drivergenes_filter'] = lastjobs2['drivergenes_filter']['txt']
#	create_dir(driverDir)
#	driverGene_cmd,order1 = thisAdv.driverGene()
#	add_items(orders2,order1)
#	new_jobs2['driverGene'] = \
#		{'name' : 'driverGene',
#		 'status' : 'waiting',
#		 'sched' : '-V -cwd -q human.q -q all.q -l p=%s' % threads['driverGene'],
#		 'memory' : memory['driverGene'],
#		 'cmd' : 'sh '+os.path.join(driverDir,'driverGene.sh')}
#	safe_open(os.path.join(driverDir,'driverGene.sh'),'w').write(driverGene_cmd)
#	job_points2['driverGene'] = ['driverGene']
#	if lastjobs2:
#		if 'driverGene' not in new_jobs2 and 'driverGene' in lastjobs2:
#			add_jobs2['driverGene'] = lastjobs2['driverGene']['txt']
	
## SMG
if set(['4']).issubset(adv_array_cancer):
	print 'SMG ...'
	create_dir(musicDir)
	smg_cmd,order1 = thisAdv.smg()
	add_items(orders2,order1)
	new_jobs2['smg'] = \
		{'name' : 'smg',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['smg']),
		 'memory' : memory['smg'],
		 'cmd' : 'sh '+os.path.join(musicDir,'smg.sh')}
	safe_open(os.path.join(musicDir,'smg.sh'),'w').write(smg_cmd)
	job_points2['smg'] = ['smg']
	if lastjobs2:
		if 'smg' not in new_jobs2 and 'smg' in lastjobs2:
			add_jobs2['smg'] = lastjobs2['smg']['txt']
	
## gistic
if set(['5']).issubset(adv_array_cancer):
	print 'gistic ...'
	create_dir(gisticDir)
	gistic_cmd,order1 = thisAdv.gistic()
	add_items(orders2,order1)
	new_jobs2['gistic'] = \
		{'name' : 'gistic',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['gistic']),
		 'memory' : memory['gistic'],
		 'cmd' : 'sh '+os.path.join(gisticDir,'gistic.sh')}
	safe_open(os.path.join(gisticDir,'gistic.sh'),'w').write(gistic_cmd)
	job_points2['gistic'] = ['gistic']
	if lastjobs2:
		if 'gistic' not in new_jobs2 and 'gistic' in lastjobs2:
			add_jobs2['gistic'] = lastjobs2['gistic']['txt']
	
## fusion gene
if set(['6']).issubset(adv_array_cancer):
	print 'fusion gene ...'
	create_dir(fusionDir)
	fusion_gene_cmd,order1 = thisAdv.fusion_gene()
	add_items(orders2,order1)
	new_jobs2['fusion_gene'] = \
		{'name' : 'fusion_gene',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['fusion_gene']),
		 'memory' : memory['fusion_gene'],
		 'cmd' : 'sh '+os.path.join(fusionDir,'fusion_gene.sh')}
	safe_open(os.path.join(fusionDir,'fusion_gene.sh'),'w').write(fusion_gene_cmd)
	job_points2['fusion_gene'] = ['fusion_gene']
	if lastjobs2:
		if 'fusion_gene' not in new_jobs2 and 'fusion_gene' in lastjobs2:
			add_jobs2['fusion_gene'] = lastjobs2['fusion_gene']['txt']
	
## absolute
if set(['7']).issubset(adv_array_cancer):
	print 'absolute ...'
	create_dir(absoluteDir)
	maf_cnv = ['%s\t%s\t%s\n' %(each,thisAdv.contruct_somaf(each,sosnp_method,'snv'),thisAdv.contruct_socnv(each,socnv_method)) for each in Tsamples]
	safe_open(os.path.join(prepareDir,'absolute.maf_cnv.list'),'w').write(''.join(maf_cnv))
	safe_open(os.path.join(absoluteDir,'maf_cnv.list'),'w').write(''.join(maf_cnv))
	absolute_cmd,order1 = thisAdv.absolutes()
	add_items(orders2,order1)
	new_jobs2['absolutes'] = \
		{'name' : 'absolutes',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['absolutes']),
		 'memory' : memory['absolutes'],
		 'cmd' : 'sh '+os.path.join(absoluteDir,'absolutes.sh')}
	safe_open(os.path.join(absoluteDir,'absolutes.sh'),'w').write(absolute_cmd)
	job_points2['absolutes'] = ['absolutes']
	if lastjobs2:
		if 'absolutes' not in new_jobs2 and 'absolutes' in lastjobs2:
			add_jobs2['absolutes'] = lastjobs2['absolutes']['txt']

## clone == pyclone
if set(['8.1']).issubset(adv_array_cancer):
	print 'clone pyclone ...'
	create_dir(cloneDir)
	create_dir(pycloneDir)
	maf_cnv = ['%s\t%s\t%s\n' %(each,thisAdv.contruct_somaf(each,sosnp_method,'snv'),thisAdv.contruct_socnv(each,socnv_method)) for each in Tsamples]
	safe_open(os.path.join(prepareDir,'pyclone.maf_cnv.list'),'w').write(''.join(maf_cnv))
	safe_open(os.path.join(pycloneDir,'maf_cnv.list'),'w').write(''.join(maf_cnv))
	pyclones_cmd,order1 = thisAdv.pyclones()
	add_items(orders2,order1)
	new_jobs2['pyclones'] = \
		{'name' : 'pyclones',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['pyclones']),
		 'memory' : memory['pyclones'],
		 'cmd' : 'sh '+os.path.join(pycloneDir,'pyclones.sh')}
	safe_open(os.path.join(pycloneDir,'pyclones.sh'),'w').write(pyclones_cmd)
	job_points2['pyclones'] = ['pyclones']
	if lastjobs2:
		if 'pyclones' not in new_jobs2 and 'pyclones' in lastjobs2:
			add_jobs2['pyclones'] = lastjobs2['pyclones']['txt']
	
	pyclone_connect_cmd,order2 = thisAdv.pyclone_connect()
	add_items(orders2,order2)
	new_jobs2['pyclone_connect'] = \
		{
			'name' : 'pyclone_connect',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['pyclones']),
			'memory' : memory['pyclone_connect'],
			'cmd' : 'sh '+os.path.join(pycloneDir,'connect/%s'%(pyclone_connect_cmd))
		}
	job_points2['pyclones'].append('pyclone_connect')
	if lastjobs2:
		if 'pyclone_connect' not in new_jobs2 and 'pyclone_connect' in lastjobs2:
			add_jobs2['pyclone_connect'] = lastjobs2['pyclone_connect']['txt']
	
## clone == sciclone
if set(['8.2']).issubset(adv_array_cancer):
	print 'clone sciclone ...'
	create_dir(cloneDir)
	create_dir(scicloneDir)
	maf_cnv = ['%s\t%s\t%s\n' %(each,thisAdv.contruct_somaf(each,sosnp_method,'snv'),thisAdv.contruct_socnv(each,socnv_method)) for each in Tsamples]
	safe_open(os.path.join(prepareDir,'sciclone.maf_cnv.list'),'w').write(''.join(maf_cnv))
	safe_open(os.path.join(scicloneDir,'maf_cnv.list'),'w').write(''.join(maf_cnv))
	sciclones_cmd,order1 = thisAdv.sciclones()
	add_items(orders2,order1)
	new_jobs2['sciclones'] = \
		{'name' : 'sciclones',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['sciclones']),
		 'memory' : memory['sciclones'],
		 'cmd' : 'sh '+os.path.join(scicloneDir,'sciclones.sh')}
	safe_open(os.path.join(scicloneDir,'sciclones.sh'),'w').write(sciclones_cmd)
	job_points2['sciclones'] = ['sciclones']
	if lastjobs2:
		if 'sciclones' not in new_jobs2 and 'sciclones' in lastjobs2:
			add_jobs2['sciclones'] = lastjobs2['sciclones']['txt']

## clone == expands
if set(['8.3']).issubset(adv_array_cancer):
	print 'clone expands ...'
	create_dir(cloneDir)
	create_dir(expandsDir)
	maf_cnv = ['%s\t%s\t%s\n' %(each,thisAdv.contruct_somaf(each,sosnp_method,'snv'),thisAdv.contruct_socnv(each,socnv_method)) for each in Tsamples]
	safe_open(os.path.join(prepareDir,'expands.maf_cnv.list'),'w').write(''.join(maf_cnv))
	safe_open(os.path.join(expandsDir,'maf_cnv.list'),'w').write(''.join(maf_cnv))
	expands_cmd,order1 = thisAdv.expands()
	add_items(orders2,order1)
	new_jobs2['expands'] = \
		{'name' : 'expands',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['expands']),
		 'memory' : memory['expands'],
		 'cmd' : 'sh '+os.path.join(expandsDir,'expands.sh')}
	safe_open(os.path.join(expandsDir,'expands.sh'),'w').write(expands_cmd)
	job_points2['expands'] = ['expands']
	if lastjobs2:
		if 'expands' not in new_jobs2 and 'expands' in lastjobs2:
			add_jobs2['expands'] = lastjobs2['expands']['txt']

patients = patientInfo.keys()

## evolution
if set(['9']).issubset(adv_array_cancer):
	print 'evolution ...'
	create_dir(evolutionDir)
	for eachP in patientInfo.keys():
		tumors_tmp = list(set(patientInfo[eachP]['T']) & set(Tsamples))
		if len(tumors_tmp) <= 1:
			continue
		evo_tmpdir = os.path.join(evolutionDir,eachP)
		create_dir(evo_tmpdir)
		phylip_evolution_cmd,order1 = thisAdv.phylip_evolution(eachP,tumors_tmp,func='coding')
		add_items(orders2,order1)
		phylip_evolution_name = '_'.join(['phylip_evolution',eachP])
		new_jobs2[phylip_evolution_name] = \
			{'name' : phylip_evolution_name,
			 'status' : 'waiting',
			 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['phylip_evolution']),
			 'memory' : memory['phylip_evolution'],
			 'cmd' : 'sh '+os.path.join(evo_tmpdir,'phylip_evolution_'+eachP+'.sh')}
		safe_open(os.path.join(evo_tmpdir,'phylip_evolution_'+eachP+'.sh'),'w').write(phylip_evolution_cmd)
		job_points2[phylip_evolution_name] = ['phylip_evolution']

		if lastjobs2:
			if phylip_evolution_name not in new_jobs2 and phylip_evolution_name in lastjobs2:
				add_jobs2[phylip_evolution_name] = lastjobs2[phylip_evolution_name]['txt']
		
	
## circos
if set(['10']).issubset(adv_array_cancer):
	print 'circos ...'
	do_sv = 'somatic_sv' in include_symbol
	create_dir(circosDir)
	for eachsample in Tsamples:
		circos_tmpdir = os.path.join(circosDir,eachsample)
		create_dir(circos_tmpdir)
		circos_cmd,order1 = thisAdv.circos(eachsample,do_sv)
		add_items(orders2,order1)
		circos_name = 'circos_%s'%eachsample
		new_jobs2[circos_name] = \
			{'name' : circos_name,
			 'status' : 'waiting',
			 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['circos']),
			 'memory' : memory['circos'],
			 'cmd' : 'sh '+os.path.join(circos_tmpdir,'circos_'+eachsample+'.sh')}
		safe_open(os.path.join(circos_tmpdir,'circos_'+eachsample+'.sh'),'w').write(circos_cmd)
		job_points2[circos_name] = ['circos']
		if lastjobs2:
			if circos_name not in new_jobs2 and circos_name in lastjobs2:
				add_jobs2[circos_name] = lastjobs2[circos_name]['txt']

## mrt_music
if set(['11']).issubset(adv_array_cancer):
	print 'mrt music ... '
	create_dir(mrtDir)
	mrt_cmd,order1 = thisAdv.mrt_music()
	add_items(orders2,order1)
	new_jobs2['mrt_music'] = \
		{'name' : 'mrt_music',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['mrt_music']),
		 'memory' : memory['mrt_music'],
		 'cmd' : 'sh '+os.path.join(mrtDir,'mrt_music.sh')}
	safe_open(os.path.join(mrtDir,'mrt_music.sh'),'w').write(mrt_cmd)
	job_points2['mrt_music'] = ['mrt_music']
	if lastjobs2:
		if 'mrt_music' not in new_jobs2 and 'mrt_music' in lastjobs2:
			add_jobs2['mrt_music'] = lastjobs2['mrt_music']['txt']

## oncodriveclust
if set(['12']).issubset(adv_array_cancer):
	print 'oncodriveclust ... '
	create_dir(oncodriveDir)
	oncodrive_cmd,order1 = thisAdv.oncodriveclust()
	add_items(orders2,order1)
	new_jobs2['oncodriveclust'] = \
		{'name' : 'oncodriveclust',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['oncodriveclust']),
		 'memory' : memory['oncodriveclust'],
		 'cmd' : 'sh '+os.path.join(oncodriveDir,'oncodriveclust.sh')}
	safe_open(os.path.join(oncodriveDir,'oncodriveclust.sh'),'w').write(oncodrive_cmd)
	job_points2['oncodriveclust'] = ['oncodriveclust']
	if lastjobs2:
		if 'oncodriveclust' not in new_jobs2 and 'oncodriveclust' in lastjobs2:
			add_jobs2['oncodriveclust'] = lastjobs2['oncodriveclust']['txt']

## drug_target
if set(['13']).issubset(adv_array_cancer):
	print 'drug target ...'
	create_dir(targetDir)
	target_cmd,order1 = thisAdv.drug_target()
	add_items(orders2,order1)
	new_jobs2['drug_target'] = \
		{'name' : 'drug_target',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['drug_target']),
		 'memory' : memory['drug_target'],
		 'cmd' : 'sh '+os.path.join(targetDir,'drug_target.sh')}
	safe_open(os.path.join(targetDir,'drug_target.sh'),'w').write(target_cmd)
	job_points2['drug_target'] = ['drug_target']
	if lastjobs2:
		if 'drug_target' not in new_jobs2 and 'drug_target' in lastjobs2:
			add_jobs2['drug_target'] = lastjobs2['drug_target']['txt']

## drug_resistence
if set(['14']).issubset(adv_array_cancer):
	print 'drug resistence ...'
	create_dir(resistenceDir)
	resistence_cmd,order1 = thisAdv.drug_resistance()
	add_items(orders2,order1)
	new_jobs2['drug_resistance'] = \
		{'name' : 'drug_resistance',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['drug_resistance']),
		 'memory' : memory['drug_resistance'],
		 'cmd' : 'sh '+os.path.join(resistenceDir,'drug_resistance.sh')}
	safe_open(os.path.join(resistenceDir,'drug_resistance.sh'),'w').write(resistence_cmd)
	job_points2['drug_resistance'] = ['drug_resistance']
	if lastjobs2:
		if 'drug_resistance' not in new_jobs2 and 'drug_resistance' in lastjobs2:
			add_jobs2['drug_resistance'] = lastjobs2['drug_resistance']['txt']

## noncoding_filter
if set(['15']).issubset(adv_array_cancer):
	print 'noncoding filter...'
	create_dir(noncodingDir)
	noncoding_cmd,order1 = thisAdv.noncoding_filter()
	add_items(orders2,order1)
	new_jobs2['noncoding_filter'] = \
		{'name' : 'noncoding_filter',
		 'status' : 'waiting',
		 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['noncoding_filter']),
		 'memory' : memory['noncoding_filter'],
		 'cmd' : 'sh '+os.path.join(noncodingDir,'noncoding_filter.sh')}
	safe_open(os.path.join(noncodingDir,'noncoding_filter.sh'),'w').write(noncoding_cmd)
	job_points2['noncoding_filter'] = ['noncoding_filter']
	if lastjobs2:
		if 'noncoding_filter' not in new_jobs2 and 'noncoding_filter' in lastjobs2:
			add_jobs2['noncoding_filter'] = lastjobs2['noncoding_filter']['txt']

##mutation show   #by sun
if set(['17']).issubset(adv_array_cancer):
        print 'mutation show...'
        create_dir(mutshowDir)
        if set(['4']).issubset(adv_array_cancer) and not set(['12']).issubset(adv_array_cancer):
#                create_dir(os.path.join(mutshowDir,'Smg_Mutationshow'))
                x = 's'
        if set(['12']).issubset(adv_array_cancer) and not set(['4']).issubset(adv_array_cancer):
#                create_dir(os.path.join(mutshowDir,'Oncodrive_Mutationshow'))
                x = 'o'
        if set(['4']).issubset(adv_array_cancer) and set(['12']).issubset(adv_array_cancer):
#                create_dir(os.path.join(mutshowDir,'Oncodrive_Mutationshow'))
#                create_dir(os.path.join(mutshowDir,'Smg_Mutationshow'))
                x='os'
        else:
                print "before mutationshow,smg or oncodriveclust should be added"
        mutshow_cmd,order1 = thisAdv.mutshow(x)
        add_items(orders2,order1)
        new_jobs2['mutationshow'] = \
                {'name' : 'mutationshow',
                'status' : 'waiting',
                'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['mutationshow']),
                'memory' : memory['mutationshow'],
                'cmd' : 'sh '+os.path.join(mutshowDir,'mutshow.sh')}
        safe_open(os.path.join(mutshowDir,'mutshow.sh'),'w').write(mutshow_cmd)
        job_points2['mutationshow'] = ['mutationshow']
        if lastjobs2:
                if 'mutationshow' not in new_jobs2 and 'mutationshow' in lastjobs2:
                        add_jobs2['mutationshow'] = lastjobs2['mutationshow']['txt']


## bubbletree
if set(['16']).issubset(adv_array_cancer):
	print 'Bubbletree...'
	for each in Tsamples:
		create_dir(bubbletreeDir)
		bubbletree_cmd,order1 = thisAdv.bubbletree(each)
		add_items(orders2,order1)
		bubbletree_name = '_'.join(['bubbletree',each])
		new_jobs2[bubbletree_name] = \
			{'name' : bubbletree_name,
			 'status' : 'waiting',
			 'sched' : '-V -cwd %s -l p=%s' % (queue_list,threads['bubbletree']),
			 'memory' : memory['bubbletree'],
			 'cmd' : 'sh '+os.path.join(bubbletreeDir,bubbletree_name+'.sh')}
		safe_open(os.path.join(bubbletreeDir,bubbletree_name+'.sh'),'w').write(bubbletree_cmd)
		job_points2[bubbletree_name] = ['bubbletree']
		if lastjobs2:
			if bubbletree_name not in new_jobs2 and bubbletree_name in lastjobs2:
				add_jobs2[bubbletree_name] = lastjobs2[bubbletree_name]['txt']

########################################
######disease advancer#######  #by sun #
########################################
if set(['8.1']).issubset(adv_array_disease) or set(['8.2']).issubset(adv_array_disease) or set(['8.3']).issubset(adv_array_disease) or set(['8.4']).issubset(adv_array_disease) or set(['8.5']).issubset(adv_array_disease) or set(['8.6']).issubset(adv_array_disease):
        if dnmf == '':
                exit('Error: Do denovo mutation detect need family config, please see the help.')

if set(['7.1']).issubset(adv_array_disease) or set(['7.2']).issubset(adv_array_disease) or set(['10.1']).issubset(adv_array_disease):
        if samp_info == 'Null':
                exit('Error: Model Filter need sample_info file to provide family information.')
        else:
                FinFo = {}
                count = 0
                for line in safe_open(samp_info,'r'):
                        count += 1
                        if count == 1:
                                firstline = line.strip().split('\t')
                                if '#FamilyID' not in firstline or 'SampleID' not in firstline or 'Normal/Patient' not in firstline:
                                        exit('Wrong sample info format: #FamilyID, SampleID and Normal/Patient must in the firt row.')
                                else:
                                        pfid = firstline.index('#FamilyID')
                                        psid = firstline.index('SampleID')
                        else:
                                uline = line.strip().split('\t')
                                fid = uline[pfid]
                                sid = uline[psid]
                                if fid not in FinFo:
                                        FinFo[fid] = []
                                if sid not in FinFo[fid]:
                                        FinFo[fid].append(sid)
                print FinFo

#merge vcf and filterDB
pythonDIR = '/PUBLIC/software/public/Python-2.7.6/bin'
if set(['6.1']).issubset(adv_array_disease):
        print 'FileterDB...'
        if not os.path.exists(advdir_disease):
                os.makedirs(advdir_disease)
        AdvanceAna = Advance(SamID,mutdir,advdir_disease,times,snvsoft,ref,seqstrag,moduledir, pythonDIR)
        FiltDir = os.path.join(advdir_disease,'Merged_vcf')
        if not os.path.exists(FiltDir):
                os.mkdir(FiltDir)
        FiltDB_cmd = AdvanceAna.merge_filter()
        with safe_open(os.path.join(FiltDir,'merged_vcf.sh'),'w') as shell:
                shell.writelines(FiltDB_cmd+'\n')
        new_jobs2['merged_vcf'] = \
                {'name' : 'merged_vcf',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % queue_list,
                'memory' : memory['merged_vcf'],
                'cmd' : 'sh '+os.path.join(FiltDir,'merged_vcf.sh')}
        job_points2['merged_vcf'] = ['merged_vcf']
        if lastjobs2:
                if 'merged_vcf' not in new_jobs2 and 'merged_vcf' in lastjobs2:
                        add_jobs2['merged_vcf'] = lastjobs2['merged_vcf']['txt']        

        if set(['7.1']).issubset(adv_array_disease) or set(['7.2']).issubset(adv_array_disease):
                print 'Model Fileter...'
                MFdir = os.path.join(advdir_disease,'ModelF')
                if not os.path.exists(MFdir):
                        os.mkdir(MFdir)
                dbFSNP = os.path.join(advdir_disease,'Merged_vcf','Filter','SNP','snp.merged.1000g.func.syn.deleterious.xls')
                dbFindel = os.path.join(advdir_disease,'Merged_vcf','Filter','INDEL','indel.merged.1000g.func.syn.deleterious.xls')
                dbFsnpindel = os.path.join(advdir_disease,'Merged_vcf','Filter','SNP','snp.indel.merged.1000g.func.syn.deleterious.xls')
                if set(['7.1']).issubset(adv_array_disease):
                        #snp
                        modelsnpD,order1 = AdvanceAna.MFilter(samp_info,'D','snp',dbFSNP,MFdir)
                        with safe_open(MFdir+'/ModelF_D_snp'+'.sh','w') as shell:
                                shell.write(modelsnpD)
                        new_jobs2['ModelF_D_snp'] = \
                                {'name' : 'ModelF_D_snp',
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['ModelF'],
                                'cmd' : 'sh '+os.path.join(MFdir,'ModelF_D_snp.sh')}
                        add_items(orders2,order1)
                        job_points2['ModelF_D_snp'] = ['ModelF_D_snp']
                        if lastjobs2:
                                if 'ModelF_D_snp' not in new_jobs2 and 'ModelF_D_snp' in lastjobs2:
                                        add_jobs2['ModelF_D_snp'] = lastjobs2['ModelF_D_snp']['txt'] 
                        #indel
                        modelindelD,order1 = AdvanceAna.MFilter(samp_info,'D','indel',dbFindel,MFdir)
                        with safe_open(MFdir+'/ModelF_D_indel'+'.sh','w') as shell:
                                shell.write(modelindelD)
                        new_jobs2['ModelF_D_indel'] = \
                                {'name' : 'ModelF_D_indel',
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['ModelF'],
                                'cmd' : 'sh '+os.path.join(MFdir,'ModelF_D_indel.sh')}
                        add_items(orders2,order1)
                        job_points2['ModelF_D_indel'] = ['ModelF_D_indel']
                        if lastjobs2:
                                if 'ModelF_D_indel' not in new_jobs2 and 'ModelF_D_indel' in lastjobs2:
                                        add_jobs2['ModelF_D_indel'] = lastjobs2['ModelF_D_indel']['txt']   

                if set(['7.2']).issubset(adv_array_disease):
                        #snp
                        modelsnpR,order1 = AdvanceAna.MFilter(samp_info,'R','snp',dbFSNP,MFdir)
                        with safe_open(MFdir+'/ModelF_R_snp'+'.sh','w') as shell:
                            shell.write(modelsnpR)
                        new_jobs2['ModelF_R_snp'] = \
                                {'name' : 'ModelF_R_snp',
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['ModelF'],
                                'cmd' : 'sh '+os.path.join(MFdir,'ModelF_R_snp.sh')}
                        add_items(orders2,order1)
                        job_points2['ModelF_R_snp'] = ['ModelF_R_snp']
                        if lastjobs2:
                                if 'ModelF_R_indel' not in new_jobs2 and 'ModelF_R_indel' in lastjobs2:
                                        add_jobs2['ModelF_R_indel'] = lastjobs2['ModelF_R_indel']['txt']
                        #ModelF_C
                        modelsnpC,order1 = AdvanceAna.MFilter(samp_info,'C','snpindel',dbFsnpindel,MFdir)
                        with safe_open(MFdir+'/ModelF_C_snpindel'+'.sh','w') as shell:
                            catcmd = 'cat %s %s |sort -n -k2 -k3 | awk \'{if(NR!=1)print}\'> %s' %(dbFSNP,dbFindel,dbFsnpindel)
                            rmcat = 'rm -f %s' %dbFsnpindel
                            shell.write(catcmd+' && \\\n'+ modelsnpC +'\n'+rmcat)
                        new_jobs2['ModelF_C_snpindel'] =  \
                                {'name' : 'ModelF_C_snpindel',
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['ModelF'],
                                'cmd' : 'sh '+os.path.join(MFdir,'ModelF_C_snpindel.sh')}
                        add_items(orders2,order1)
                        job_points2['ModelF_C_snpindel'] = ['ModelF_C_snpindel']
                        if lastjobs2:
                                if 'ModelF_C_snpindel' not in new_jobs2 and 'ModelF_C_snpindel' in lastjobs2:
                                        add_jobs2['ModelF_C_snpindel'] = lastjobs2['ModelF_C_snpindel']['txt']
                        #indel
                        modelindelR,order1 = AdvanceAna.MFilter(samp_info,'R','indel',dbFindel,MFdir)
                        with safe_open(MFdir+'/ModelF_R_indel'+'.sh','w') as shell:
                            shell.write(modelindelR)
                        new_jobs2['ModelF_R_indel'] = \
                                {'name' : 'ModelF_R_indel',
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['ModelF'],
                                'cmd' : 'sh '+os.path.join(MFdir,'ModelF_R_indel.sh')}
                        add_items(orders2,order1)
                        job_points2['ModelF_R_indel'] = ['ModelF_R_indel']
                        if lastjobs2:
                                if 'ModelF_R_indel' not in new_jobs2 and 'ModelF_R_indel' in lastjobs2:
                                        add_jobs2['ModelF_R_indel'] = lastjobs2['ModelF_R_indel']['txt']

#Denovos_samtools
if set(['8.1']).issubset(adv_array_disease):
        print 'Denove Samtools...'
        if not os.path.exists(advdir_disease):
                os.makedirs(advdir_disease)
        sdsoft = 'samtools'
        SDir = os.path.join(advdir_disease,'DenovoSam')
        if not os.path.exists(SDir):
                os.mkdir(SDir)
        for k in dnminfo.keys():
                for j in dnminfo[k].keys():
                        samSDir = os.path.join(SDir,k,j)
                        if not os.path.exists(samSDir):
                                os.makedirs(samSDir)
                        if dnminfo[k][j][2] == 'F' or dnminfo[k][j][2] == 'Female':
                                c_sex = str(2)
                        elif dnminfo[k][j][2] == 'M' or dnminfo[k][j][2] == 'Male':
                                c_sex = str(1)
                        else:
                                c_sex = str(dnminfo[k][j][2])
                        with safe_open(os.path.join(samSDir,j+'.ped'),'w') as ped:
                                ped.write('\n'.join(['\t'.join([k,dnminfo[k][j][0],'0','0','1','1']),'\t'.join([k,dnminfo[k][j][1],'0','0','2','1']),'\t'.join([k,j,dnminfo[k][j][0],dnminfo[k][j][1],c_sex,'2'])])+'\n')
                        pedf = os.path.join(samSDir,j+'.ped')
                        with safe_open(os.path.join(samSDir,'mkdenovo_s_'+j+'.sh'),'w') as shell:
                                cbam = os.path.join(mapdir,j,j+'.final.bam')  #child
                                fbam = os.path.join(mapdir,dnminfo[k][j][0],dnminfo[k][j][0]+'.final.bam') #father
                                mbam = os.path.join(mapdir,dnminfo[k][j][1],dnminfo[k][j][1]+'.final.bam') #mather
                                shell.write(' \\\n\t'.join(['python %s/Varition/DNM/Human_reseq_dnm_v7_1.py' % (moduledir),  #new_anno
                                        '--ped %s' %pedf,
                                        '--bams %s,%s,%s' %(cbam,fbam,mbam),
                                        '--soft %s' %sdsoft,
                                        '--moduledir %s' %moduledir,
                                        '--database_dir %s' %database_dir,
                                        '--familyID %s' %j,
                                        '--ref %s' % ref,
                                        '--outdir %s' %samSDir])+'\n')
                        new_jobs2['_'.join(['mkdenovo','s_'+j])] =\
                                {'name' : '_'.join(['mkdenovo','s_'+j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['mkdenovo'],
                                'cmd' : 'sh '+os.path.join(samSDir,'mkdenovo_s_'+j+'.sh')}
                        order1 = 'order mkdenovo_s_%s after merged_vcf' % j
                        add_items(orders2,order1)

                        if c_sex == '1':
                                dechr = chrlist[:24]
                        else:
                                dechr = chrlist[:23]
                        for ch in dechr:
                                new_jobs2['_'.join(['samtoolsMpileup_chr',ch,j])] = \
                                {'name' : '_'.join(['samtoolsMpileup_chr',ch,j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['samtoolsMpileup_chr'],
                                'cmd' : 'sh '+os.path.join(samSDir,'samtoolsMpileup_chr_'+ch+'_'+j+'.sh')}
                                order1 = 'order samtoolsMpileup_chr_%s_%s after mkdenovo_s_%s' %(ch,j,j)
                                add_items(orders2,order1)
                                order2 = 'order catchr_%s after samtoolsMpileup_chr_%s_%s' %(j,ch,j)
                                add_items(orders2,order2)

                        order3 = 'order %s_DNM_%s after catchr_%s' %(sdsoft,j,j)
                        add_items(orders2,order3)
                        new_jobs2['_'.join(['catchr',j])] = \
                                {'name' : '_'.join(['catchr',j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['catchr'],
                                'cmd' : 'sh '+os.path.join(samSDir,'catchr_'+j+'.sh')}
                        new_jobs2['_'.join([sdsoft,'DNM',j])] =\
                                {'name' : '_'.join([sdsoft,'DNM',j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['samtools_DNM'],
                                'cmd' : 'sh '+os.path.join(samSDir,'samtools_DNM_'+j+'.sh')}
                        new_jobs2['_'.join([sdsoft,'DNM_snpanno',j])] = \
                                {'name' : '_'.join([sdsoft,'DNM_snpanno',j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['samtools_DNM_snpanno'],
                                'cmd' : 'sh '+os.path.join(samSDir,'samtools_DNM_snpanno_'+j+'.sh')}
                        order4 = 'order %s_DNM_snpanno_%s after %s_DNM_%s' %(sdsoft,j,sdsoft,j)
                        add_items(orders2,order4)
                        new_jobs2['_'.join([sdsoft,'DNM_indelanno',j])] = \
                                {'name' : '_'.join([sdsoft,'DNM_indelanno',j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['samtools_DNM_indelanno'],
                                'cmd' : 'sh '+os.path.join(samSDir,'samtools_DNM_indelanno_'+j+'.sh')}
                        order4 = 'order %s_DNM_indelanno_%s after %s_DNM_%s' %(sdsoft,j,sdsoft,j)
                        add_items(orders2,order4)
    
#Denovos_denovogear
if set(['8.2']).issubset(adv_array_disease):
        print 'Denovo Denovogear...'
        if not os.path.exists(advdir_disease):
                os.makedirs(advdir_disease)
        DDir = os.path.join(advdir_disease,'DenovoGear')
        gdsoft = 'denovogear'
        if not os.path.exists(DDir):
                os.mkdir(DDir)
        for k in dnminfo.keys():
                samGDir = os.path.join(DDir,k)
                if not os.path.exists(samGDir):
                        os.mkdir(samGDir)
                if dnminfo[k][3] == 'F' or dnminfo[k][3] == 'Female':
                        c_sex = str(2)
                elif dnminfo[k][3] == 'M' or dnminfo[k][3] == 'Male':
                        c_sex = str(1)
                else:
                        c_sex = str(dnminfo[k][3])
                with safe_open(os.path.join(samGDir,k+'.ped'),'w') as ped:
                        ped.write('\n'.join(['\t'.join([k,dnminfo[k][1],'0','0','1','1']),'\t'.join([k,dnminfo[k][2],'0','0','2','1']),'\t'.join([k,dnminfo[k][0],dnminfo[k][1],dnminfo[k][2],c_sex,'2'])])+'\n')
                pedf = os.path.join(samGDir,k+'.ped')
                with safe_open(os.path.join(samGDir,'mkdenovo_g_'+k+'.sh'),'w') as shell:
                        cbam = os.path.join(mapdir,dnminfo[k][0],dnminfo[k][0]+'.final.bam')
                        fbam = os.path.join(mapdir,dnminfo[k][1],dnminfo[k][1]+'.final.bam')
                        mbam = os.path.join(mapdir,dnminfo[k][2],dnminfo[k][2]+'.final.bam')
                        shell.write(' \\\n\t'.join(['python %s/Varition/DNM/Human_reseq_dnm_v7.py' % (moduledir),  #new_anno
                                '--ped %s' %pedf,
                                '--bams %s,%s,%s' %(cbam,fbam,mbam),
                                '--soft %s' %gdsoft,
                                '--moduledir %s' %moduledir,
                                '--database_dir %s' %database_dir,
                                '--outdir %s' %samGDir])+'\n')
                new_jobs2['_'.join(['mkdenovo','g_'+k])] = \
                                {'name' : '_'.join(['mkdenovo','g_'+k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['mkdenovo'],
                                'cmd' : 'sh '+os.path.join(samGDir,'mkdenovo_s_'+j+'.sh')}
                order1 = 'order mkdenovo_g_%s after merged_vcf' % k
                add_items(orders2,order1)

                new_jobs2['_'.join([gdsoft,'DNM',k])] = \
                                {'name' : '_'.join([gdsoft,'DNM',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['denovogear_DNM'],
                                'cmd' : 'sh '+os.path.join(samGDir,'denovogear_DNM_'+j+'.sh')}
                order4 = 'order %s_DNM_%s after mkdenovo_g_%s' %(gdsoft,k,k)
                add_items(orders2,order4)
                new_jobs2['_'.join([gdsoft,'DNM_snpanno',k])] = \
                                {'name' : '_'.join([gdsoft,'DNM_snpanno',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['denovogear_DNM_snpanno'],
                                'cmd' : 'sh '+os.path.join(samGDir,'denovogear_DNM_snpanno_'+j+'.sh')}
                order4 = 'order %s_DNM_snpanno_%s after %s_DNM_%s' %(gdsoft,k,gdsoft,k)
                add_items(orders2,order4)
                new_jobs2['_'.join([gdsoft,'DNM_indelanno',k])] =  \
                                {'name' : '_'.join([gdsoft,'DNM_indelanno',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['denovogear_DNM_indelanno'],
                                'cmd' : 'sh '+os.path.join(samGDir,'denovogear_DNM_indelanno_'+j+'.sh')}
                order4 = 'order %s_DNM_indelanno_%s after %s_DNM_%s' %(gdsoft,k,gdsoft,k)
                add_items(orders2,order4)
#DenovoF 
if set(['8.3']).issubset(adv_array_disease):
        if not os.path.exists(advdir_disease):
                os.makedirs(advdir_disease)        
        print 'DenovoF...'
        FDir = os.path.join(advdir_disease,'DenovoF')
        if not os.path.exists(FDir):
                os.mkdir(FDir)
        os.system('python %s/Varition/DNM/Denovo_v7_1.py --projdir %s --dnmf %s --soft %s --out %s --moduledir %s --ref %s' %(moduledir, analydir, dnmf, snvsoft, FDir, moduledir, ref))
        for k in dnminfo.keys():
                for j in dnminfo[k].keys():
                        DFSNPdir = os.path.join(FDir,k,j,'SNP')
                        DFInDelDir = os.path.join(FDir,k,j,'INDEL')
                        new_jobs2['_'.join(['Denovo','Fsnp_'+j])] = \
                                {'name' : '_'.join(['Denovo','Fsnp_'+j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['Denovo'],
                                'cmd' : 'sh '+os.path.join(DFSNPdir,'Denovo_Fsnp_'+j+'.sh')}

                        new_jobs2['_'.join(['Denovo','Findel_'+j])] = \
                                {'name' : '_'.join(['Denovo','Findel_'+j]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['Denovo'],
                                'cmd' : 'sh '+os.path.join(DFInDelDir,'Denovo_Findel_'+j+'.sh')}
                        order1 = 'order Denovo_Fsnp_%s after merged_vcf' % j
                        add_items(orders2,order1)
                        order2 = 'order Denovo_Findel_%s after merged_vcf' % j
                        add_items(orders2,order2)
                
#DenovoR
if set(['8.4']).issubset(adv_array_disease):
        print 'DenovoR...'
        if not os.path.exists(advdir_disease):
                os.makedirs(advdir_disease)
        RDir = os.path.join(advdir_disease,'DenovoR')
        if not os.path.exists(RDir):
                os.mkdir(RDir)
        os.system('python %s/Varition/DNM/Denovo_v7.py --projdir %s --dnmf %s --soft %s --fromRaw Y --out %s --moduledir %s' %(moduledir, analydir, dnmf, snvsoft , RDir, moduledir))
        for k in dnminfo.keys():
                DRSNPdir = os.path.join(RDir,k,'SNP')
                DRInDelDir = os.path.join(RDir,k,'INDEL')
                new_jobs2['_'.join(['Denovo','Rsnp_'+k])] = \
                                {'name' : '_'.join(['Denovo','Rsnp_'+k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['Denovo'],
                                'cmd' : 'sh '+os.path.join(DRSNPdir,'Denovo_Rsnp_'+k+'.sh')}
                new_jobs2['_'.join(['Denovo','Rindel_'+k])] = \
                                {'name' : '_'.join(['Denovo','Rindel_'+k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['Denovo'],
                                'cmd' : 'sh '+os.path.join(DRInDelDir,'Denovo_Rindel_'+k+'.sh')}

                order1 = 'order Denovo_Rsnp_%s after annotatVcf_snp_%s_%s' %(k,snvsoft,dnminfo[k][0])
                add_items(orders2,order1)
                order2 = 'order Denovo_Rsnp_%s after annotatVcf_snp_%s_%s' %(k,snvsoft,dnminfo[k][1])
                add_items(orders2,order2)
                order3 = 'order Denovo_Rsnp_%s after annotatVcf_snp_%s_%s' %(k,snvsoft,dnminfo[k][2])
                add_items(orders2,order3)
                order4 = 'order Denovo_Rindel_%s after annotatVcf_indel_%s_%s' %(k,snvsoft,dnminfo[k][0])
                add_items(orders2,order4)
                order5 = 'order Denovo_Rindel_%s after annotatVcf_indel_%s_%s' %(k,snvsoft,dnminfo[k][1])
                add_items(orders2,order5)
                order6 = 'order Denovo_Rindel_%s after annotatVcf_indel_%s_%s' %(k,snvsoft,dnminfo[k][2])
                add_items(orders2,order6)
#DenovoSV
#if set(['8.5']).issubset(adv_array_disease):
#        if not os.path.exists(advdir_disease):
#                os.makedirs(advdir_disease)
#        ##order added when write every sample's sv shell
#        DSVDir = os.path.join(advdir_disease,'DenovoSV')
#        if not os.path.exists(DSVDir):
#                os.makedirs(DSVDir)
#        dscmd = 'python %s/Varition/DNM/denovoSV_seeker_v0.4.py --proj %s --DNM %s --outdir %s --soft %s' %(moduledir, analydir, dnmf, advdir_disease, svsoft)
#        DSVname = DSVDir+'/'+'DenovoSV.sh'
#        with open(DSVname,'w') as shell:
#                shell.write(dscmd)
#        new_jobs['DenovoSV'] = write_jobt(DSVDir,'DenovoSV')
#        new_jobs_name['hehe'].append('DenovoSV')

#DenovoCNV
#if set(['8.6']).issubset(adv_array_disease):
#        if not os.path.exists(advdir_disease):
#                os.makedirs(advdir_disease)
#        ##order added when write every sample's cnv shell
#        DCNVDir = os.path.join(advdir_disease,'DenovoCNV')
#        if not os.path.exists(DCNVDir):
#                os.makedirs(DCNVDir)
#        dccmd = 'python %s/Varition/DNM/denovoSV_seeker_v0.4.py --proj %s --DNM %s --outdir %s --soft %s' %(moduledir, analydir, dnmf, advdir_disease, cnvsoft)
#        DCNVname = DCNVDir+'/'+'DenovoCNV.sh'
#        with open(DCNVname,'w') as shell:
#                shell.write(dccmd)
#        new_jobs['DenovoCNV'] = write_jobt(DCNVDir,'DenovoCNV')
#        new_jobs_name['hehe'].append('DenovoCNV')

#Linkage
if set(['9.1']).issubset(adv_array_disease):
        MLinDir = os.path.join(advdir_disease,'MerLinkage')
        MLdic = {}
        for line in safe_open(mlin,'r'):
                uline = line.strip().split('\t')
                if line.startswith('#'):continue
                if uline[0] not in MLdic:
                        MLdic[uline[0]] = uline[1].strip().split(',')
                else:
                        exit('There a familyID named %s already.Do not include same familyID more than once.' %uline[0])
        if not os.path.exists(MLinDir):
                os.mkdir(MLinDir)
        os.system('python %s/Linkage/Linkage_for_Pep_v2.py --con %s --out %s --proj %s --database_dir %s --moduledir %s' %(moduledir, mlin,MLinDir,analydir, database_dir, moduledir))
        for k in MLdic.keys():
                LSdir = os.path.join(MLinDir,'Linkage_'+k+'_me','shell')
                new_jobs2['_'.join(['merlin',k])] = \
                                {'name' : '_'.join(['merlin',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['merlin'],
                                'cmd' : 'sh '+os.path.join(LSdir,'merlin_'+k+'.sh')}
                new_jobs2['_'.join(['merlin2R',k])] = \
                                {'name' : '_'.join(['merlin2R',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['merlin2R'],
                                'cmd' : 'sh '+os.path.join(LSdir,'merlin2R_'+k+'.sh')}
                new_jobs2['_'.join(['merlin2Excel',k])] = \
                                {'name' : '_'.join(['merlin2Excel',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['merlin2Excel'],
                                'cmd' : 'sh '+os.path.join(LSdir,'merlin2Excel_'+k+'.sh')}
                new_jobs2['_'.join(['linkdatagen',k])] = \
                                {'name' : '_'.join(['linkdatagen',k]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['linkdatagen'],
                                'cmd' : 'sh '+os.path.join(LSdir,'linkdatagen_'+k+'.sh')}
                for s in range(0,len(MLdic[k])):
                        new_jobs2['_'.join(['samtoolsMpileup',k+'_'+MLdic[k][s]])] = \
                                {'name' : '_'.join(['samtoolsMpileup',k+'_'+MLdic[k][s]]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s' % queue_list,
                                'memory' : memory['samtoolsMpileup'],
                                'cmd' : 'sh '+os.path.join(LSdir,'samtoolsMpileup_'+k+'_'+MLdic[k][s]+'.sh')}
                        ordern = 'order samtoolsMpileup_%s_%s before linkdatagen_%s' %(k,MLdic[k][s],k)
                        add_items(orders2,ordern)
                order1 = 'order linkdatagen_%s before merlin_%s' %(k,k)
                add_items(orders2,order1)
                order1 = 'order merlin_%s before merlin2R_%s' %(k,k)
                add_items(orders2,order1)
                order1 = 'order merlin2R_%s before merlin2Excel_%s' %(k,k)
                add_items(orders2,order1)

#ROH
if set(['10.1']).issubset(adv_array_disease):
        ROHDir = os.path.join(advdir_disease,'ROH')
        if not os.path.exists(ROHDir):
                os.mkdir(ROHDir)
        os.system('python %s/ROH/ROH_analysis_module_v2.2.py --cm %s --software %s --info %s --pwd %s --odir %s --seqstrag %s --ref %s' %(softwares['advbin'],snvsoft,ROH_soft,samp_info,analydir,ROHDir,seqstrag, ref))

        sam_info={}
        fam_info={}
        for i in safe_open(samp_info,'r'):
                if i.startswith("#"):
                        head=i.lower().strip("#").split()
                else:
                        ii=i.strip().split()
                        sam_info[ii[head.index('sampleid')]]={}
                        if ii[head.index('familyid')]==".":
                                family=ii[head.index('sampleid')]
                        else:
                                family=ii[head.index('familyid')]
                        
                        if fam_info.has_key(family):
                                fam_info[family].append(ii[head.index('sampleid')])
                        else:
                                fam_info[family]=[];
                                fam_info[family].append(ii[head.index('sampleid')])
                        
                        sam_info[ii[head.index('sampleid')]]["family"]=family
                        sam_info[ii[head.index('sampleid')]]["sex"]=ii[head.index('sex')]
                        sam_info[ii[head.index('sampleid')]]["affect"]=ii[head.index('normal/patient')]

        if ROH_soft == "PLINK":
                for sam in sam_info.keys():
                        sampdir = "/".join([ROHDir,sam_info[sam]["family"],sam])
                        new_jobs2['_'.join([sam,'plink_roh'])] = \
                                        {'name' : '_'.join([sam,'plink_roh']),
                                        'status' : 'waiting',
                                        'sched' : '-V -cwd %s' % queue_list,
                                        'memory' : memory['ROH'],
                                        'cmd' : 'sh '+os.path.join(sampdir,sam+'_plink_roh'+'.sh')}
                        order1 = "order roh.statistics after %s" %(sam+"_plink_roh")
                        add_items(orders2,order1)

        elif ROH_soft == "H3M2":
                for sam in sam_info.keys():
                        sampdir = "/".join([ROHDir,sam_info[sam]["family"],sam])
                        new_jobs2['_'.join([sam,'h3m2_roh'])] = \
                                        {'name' : '_'.join([sam,'h3m2_roh']),
                                        'status' : 'waiting',
                                        'sched' : '-V -cwd %s' % queue_list,
                                        'memory' : memory['ROH'],
                                        'cmd' : 'sh '+os.path.join(sampdir,sam+'_h3m2_roh'+'.sh')}
                        order1 = "order roh.statistics after %s" %(sam+"_h3m2_roh")
                        add_items(orders2,order1)

        for fam in fam_info.keys():
                if len(fam_info[fam])>1:
                        combine = os.path.join(ROHDir,fam,"CombinedAnalysis")
                        new_jobs2['.'.join([fam,'roh.analysis'])] = \
                                        {'name' : '.'.join([fam,'roh.analysis']),
                                        'status' : 'waiting',
                                        'sched' : '-V -cwd %s' % queue_list,
                                        'memory' : memory['ROH'],
                                        'cmd' : 'sh '+os.path.join(combine,fam+'.roh.analysis'+'.sh')}
                        order1 = "order roh.statistics after %s" %(fam+".roh.analysis")
                        add_items(orders2,order1)
                        
                        for sam in fam_info[fam]:
                                order2="order %s after %s" % (fam+".roh.analysis",sam+"_"+ROH_soft.lower()+"_roh")
                                add_items(orders2,order2)

        new_jobs2['roh.statistics'] = \
                        {'name' : 'roh.statistics',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s' % queue_list,
                        'memory' : memory['ROH'],
                        'cmd' : 'sh '+os.path.join(ROHDir,'roh.statistics.sh')}

#####################################	
if argv['adv_array_cancer']:
	## Report
	create_dir(adv_reportdir)
	for eachjob in new_jobs2:
		order1 = 'order adv_report after %s' % eachjob
		add_items(orders2,order1)
	for eachjob in add_jobs2:
		order1 = 'order adv_report after %s' % eachjob
		add_items(orders2,order1)
	new_jobs2['adv_report'] = {'name' : 'adv_report',
        	'status' : 'waiting',
		'sched' : '-V -cwd %s' % (queue_list),
        	'memory' : memory['adv_report'],
	        'cmd' : 'sh '+os.path.join(adv_reportdir,'advance_report.sh')}
	shell = safe_open(os.path.join(adv_reportdir,'advance_report.sh'),'w')
	if argv['order']:
		shell.writelines('python %s --project %s --projdir %s --seqsty %s --adv_array_cancer %s --outdir %s --qclist %s --order %s\n' % \
		(advreport_script,os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['adv_array_cancer'],adv_reportdir,os.path.join(analydir,'qc_list'),argv['order']))
	else:
		shell.writelines('python %s --project %s --projdir %s --seqsty %s --adv_array_cancer %s --outdir %s --qclist %s\n' % \
		(advreport_script,os.path.join(analydir,'pn.txt'),analydir,seqstrag,argv['adv_array_cancer'],adv_reportdir,os.path.join(analydir,'qc_list')))


	## Result
	create_dir(adv_resultdir)
	new_jobs2['adv_result'] = {'name' : 'adv_result',
	        'status' : 'waiting',
		'sched' : '-V -cwd %s' % (queue_list),
	        'memory' : memory['adv_result'],
        	'cmd' : 'sh '+os.path.join(adv_resultdir,'cancer_adv_result.sh')}
	shell = safe_open(os.path.join(adv_resultdir,'cancer_adv_result.sh'),'w')
	create_dir(os.path.join(adv_resultdir,'cancer_adv_results'))
	shell.writelines('cd %s/cancer_adv_results &&\\\npython %s --projdir %s --type %s --analy_array %s --qclist %s --outdir %s/cancer_adv_results --readmedir %s --ref %s && \\\n' % \
(adv_resultdir,release_script,analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),adv_resultdir,READMEdir,ref))
	shell.writelines('tar -chzvf %s/cancer_adv_results/Primary_analysis_result.tar.gz Primary_analysis_result && \\\n' % (resultdir))
        shell.writelines('md5sum Primary_analysis_result.tar.gz >>MD5.txt && \\\n')
        shell.writelines('cp /PUBLIC/database/HW/CANCER/HW_v2.0/Module/Release/HW_v2.2/README_with_Advance.txt %s/cancer_adv_results/README.txt && \\\n' %(adv_resultdir)) 
        shell.writelines('python %s --projdir %s --type %s --analy_array %s --adv_array_cancer %s --qclist %s --outdir %s/cancer_adv_results --readmedir %s --ref %s && \\\n' % \
(release_script,analydir,seqstrag,argv['analy_array'],argv['adv_array_cancer'],os.path.join(analydir,'qc_list'),adv_resultdir,READMEdir,ref))
        shell.writelines('tar -chzvf %s/cancer_adv_results/Advanced_analysis_result.tar.gz Advanced_analysis_result && \\\n' % (adv_resultdir))
	shell.writelines('md5sum Advanced_analysis_result.tar.gz >>MD5.txt && \\\n')
        shell.writelines('mkdir %s/cancer_adv_results/release %s/cancer_adv_results/release_noclean && \\\n' % (adv_resultdir,adv_resultdir))
        shell.writelines('ln -sf %s/*.job/*_primary_*_report.tar.gz %s/cancer_adv_results/. && \\\n' % (resultdir,resultdir))
        shell.writelines('md5sum *_primary_*_report.tar.gz >>MD5.txt && \\\n')
        shell.writelines('md5sum *_AdvanceReport_*.tar.gz >>MD5.txt && \\\n')
	shell.writelines('ln -sf %s/cancer_adv_results/{Bam,RawData,CleanData,MD5.txt,README.txt,Primary_analysis_result.tar.gz,Advanced_analysis_result.tar.gz,*_primary_*_report.tar.gz,*_AdvanceReport_*.tar.gz} release && \\\n' % adv_resultdir)
        shell.writelines('ln -sf %s/cancer_adv_results/{Bam,RawData,Primary_analysis_result.tar.gz,Advanced_analysis_result.tar.gz,*_primary_*_report.tar.gz,*_AdvanceReport_*.tar.gz} release_noclean && \\\n' % adv_resultdir)
        shell.writelines('cp %s/cancer_adv_results/MD5.txt release_noclean && \\\n' % adv_resultdir)
        shell.writelines('sed -i "/  CleanData/d" release_noclean/MD5.txt && \\\n')
        shell.writelines('cp /PUBLIC/database/HW/CANCER/HW_v2.0/Module/Release/HW_v2.2/README_noclean_with_Advance.txt release_noclean/README.txt ')
	shell.close()
	exclude_jobs2 = ['adv_result','adv_report']
	## order
	for eachjob in new_jobs2:
		if eachjob in exclude_jobs2:
			continue
		order1 = 'order adv_result after %s' % eachjob
		add_items(orders2,order1)
	for eachjob in add_jobs2:
		if eachjob in exclude_jobs:
			continue
		order1 = 'order adv_result after %s ' % eachjob
		add_items(orders2,order1)
        order1 = 'order adv_result after adv_report '
        add_items(orders2,order1)

if argv['adv_array_disease']:
        ## Report
        create_dir(adv_reportdir)
        for eachjob in new_jobs2:
                order1 = 'order adv_report after %s' % eachjob
                add_items(orders2,order1)
        for eachjob in add_jobs2:
                order1 = 'order adv_report after %s' % eachjob
                add_items(orders2,order1)
        new_jobs2['adv_report'] = {'name' : 'adv_report',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % (queue_list),
                'memory' : memory['adv_report'],
                'cmd' : 'sh '+os.path.join(adv_reportdir,'advance_report.sh')}
        shell = safe_open(os.path.join(adv_reportdir,'advance_report.sh'),'w')
        shell.writelines('python %s --project %s --projdir %s  --sample %s --repsty advance --moduledir %s --english Y --analy_array %s --odir %s \n' % \
        (advreport_disease,os.path.join(analydir,'pn.txt'),analydir,samp_info,moduledir,argv['adv_array_disease'],os.path.join(adv_reportdir,'advance_report')))

 
        ## Result
        create_dir(adv_resultdir)
        new_jobs2['adv_result'] = {'name' : 'adv_result',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % (queue_list),
                'memory' : memory['adv_result'],
                'cmd' : 'sh '+os.path.join(adv_resultdir,'disease_adv_result.sh')}
        shell = safe_open(os.path.join(adv_resultdir,'disease_adv_result.sh'),'w')
        create_dir(os.path.join(adv_resultdir,'disease_adv_results'))
        shell.writelines('cd %s/disease_adv_results &&\\\npython %s --projdir %s --type %s --analy_array %s --qclist %s --outdir %s/disease_adv_results --readmedir %s --ref %s && \\\n' % \
(adv_resultdir,release_script,analydir,seqstrag,argv['analy_array'],os.path.join(analydir,'qc_list'),adv_resultdir,READMEdir,ref))
        shell.writelines('tar -chzvf %s/disease_adv_results/Primary_analysis_result.tar.gz Primary_analysis_result && \\\n' % (resultdir))
        shell.writelines('md5sum Primary_analysis_result.tar.gz >>MD5.txt && \\\n')
        shell.writelines('python %s --projdir %s --analy_array 3.2,%s --odir %s/disease_adv_results --moduledir %s --ref %s --cleandata N --ER Y && \\\n' % \
(disease_adv_release_script,analydir,argv['adv_array_disease'],adv_resultdir,moduledir,ref))
        shell.writelines('tar -chzvf %s/disease_adv_results/Advance.tar.gz Advance && \\\n' % (adv_resultdir))
        shell.writelines('md5sum Advance.tar.gz >>MD5.txt && \\\n')
        shell.writelines('ln -sf %s/*.job/*_primary_*_report.tar.gz %s/disease_adv_results/. && \\\n' % (resultdir,adv_resultdir))
        shell.writelines('md5sum *_primary_*_report.tar.gz >>MD5.txt && \\\n')
        shell.writelines('md5sum *.advance_report.zip >>MD5.txt && \\\n')
        shell.writelines('ln -sf %s/disease_adv_results/{Bam,RawData,CleanData,MD5.txt,Primary_analysis_result.tar.gz,Advance.tar.gz,*_primary_*_report.tar.gz,*.advance_report.zip} release && \\\n' % adv_resultdir)
        shell.writelines('ln -sf %s/disease_adv_results/{Bam,RawData,Primary_analysis_result.tar.gz,Advance.tar.gz,*_primary_*_report.tar.gz,*.advance_report.zip} release_noclean && \\\n' % adv_resultdir)
        shell.writelines('cp %s/disease_adv_results/MD5.txt release_noclean && \\\n' % adv_resultdir)
        shell.writelines('sed -i "/  CleanData/d" release_noclean/MD5.txt && \\\n')
        shell.writelines('cp /PUBLIC/source/HW/Disease/software/Data_Release_Advance/Readme_release.txt release/README.txt && \\\n')
        shell.writelines('cp /PUBLIC/source/HW/Disease/software/Data_Release_Advance/Readme_noclean.txt release_noclean/README.txt')
        shell.close()
        exclude_jobs2 = ['adv_result','adv_report']
        ## order
        for eachjob in new_jobs2:
                if eachjob in exclude_jobs2:
                        continue
                order1 = 'order adv_result after %s' % eachjob
                add_items(orders2,order1)
        for eachjob in add_jobs2:
                if eachjob in exclude_jobs:
                        continue
                order1 = 'order adv_result after %s ' % eachjob
                add_items(orders2,order1)
        order1 = 'order adv_result after adv_report '
        add_items(orders2,order1)


if adv_startpoint:
	print "Startpoint ..."
	print " ... from %s\n" % ', '.join(adv_startpoints)
	order_relation = {}
	for one in orders2:
		## order ln_EC1001_N_DHE00213_C4D7FACXX_L6 before qc_EC1001_N_DHE00213_C4D7FACXX_L6
		flag,jobA,relation,jobB = one.strip().split()
		if relation == 'before':
			if jobA not in order_relation:
				order_relation[jobA] = []
			order_relation[jobA].append(jobB)
		else:
			if jobB not in order_relation:
				order_relation[jobB] = []
			order_relation[jobB].append(jobA)
	## get the succeed jobs in startpoints list
	point_succeeds = set()
	for eachstart in adv_startpoints:
		if eachstart not in job_points2:
			print 'POINT '+eachstart+' not in you analysis, Please check your script!\nAvailable POINTS are:\n'+'  \n'.join(job_points2.keys())+'\n'
			sys.exit(1)
		for each in job_points2[eachstart]:
			tmp_succeeds = get_succeeds(each)
			point_succeeds |= set(tmp_succeeds)
	## get undone jobs in add_jobs dependent
	undone_depend_jobs = set()
	for eachjob in add_jobs2:
		if not lastjobs2[eachjob]['status'] == 'done':
			undone_depend = get_succeeds(eachjob)
			undone_depend_jobs |= undone_depend
	for eachjob in undone_depend_jobs:
		if eachjob in lastjobs2:
			## undone depended jobs must be not 'done'
			assert not lastjobs2[eachjob]['status'] == 'done'
	## change job status
	for eachjob in new_jobs2:
		if (eachjob not in point_succeeds) and (eachjob not in undone_depend_jobs):
			new_jobs2[eachjob]['status'] = 'done'

if argv['adv_array_cancer'] or argv['adv_array_disease']:
	## jobfile
	jobfile = safe_open(os.path.join(analydir,newjob),'w')
	jobfile.write('log_dir %s\n' % logdir)
	for eachjob in new_jobs2:
		sjmcmd = '\n  '.join(['  name %s' % new_jobs2[eachjob]['name'],
			'memory %s' % new_jobs2[eachjob]['memory'],
			'status %s' % new_jobs2[eachjob]['status'],
			'sched_options %s' % new_jobs2[eachjob]['sched'],
			'cmd_begin',
			'  %s' % new_jobs2[eachjob]['cmd'],
			'cmd_end'])
		jobfile.write('job_begin\n'+sjmcmd+'\njob_end\n')
	for eachjob in add_jobs2:
		jobfile.write(add_jobs2[eachjob]+'\n')

	## orders
	jobfile.write('\n'.join(orders2)+'\n')

	# lastorders
	for each in lastorders2:
		flag,jobA,relation,jobB = each.strip().split()
		if (jobA in new_jobs2 or jobA in add_jobs2) and (jobB in new_jobs2 or jobB in add_jobs2):
			if each not in orders2:
				jobfile.write(each+'\n')
	jobfile.close()
        
## Help
help = '''
\033[1;32;40m
1. Please carefully check the paramters and sample_list file. This is VERY important.
2. Submit the job file with "sjm" command: sjm myjobfile.
3. Prepare the "pn.txt" file, with project no., project name and contract no. seperated by tab.
4. If your job is not correctly finished, Please carefully check the log files.
5. After project finished, execute the shell "backup_project.sh" to backup your project information.
6. Donot forget to execute the shell "remove_mpileupgz.sh".
\033[0m
'''

print "DONE!"
print help
known = os.path.join(analydir,'KNOWN_RULES')
if not os.path.exists(known):
	open(known,'w').write(help)
