#!/usr/bin/python
#coding=utf-8
'''
sunxiuqiang 
20180905
'''
import sys
import os
import argparse
import ConfigParser
import gzip
import HTSeq
import itertools
##################################################################################
#parse the arguments
parser = argparse.ArgumentParser(description="Proportional removal of rRNA")
parser.add_argument('--leftseq',help='The left seq of sample,for example:/BJPROJ/HWUS/project/seq/sunxiuqiang/C101HW18020629.RRC.PE150.20180528.P101HW18020629-01-01/RawData/H3_KO_1_1.fq.gz ',required=True)
parser.add_argument('--rightseq',help='The right seq of sample, for example:/BJPROJ/HWUS/project/seq/sunxiuqiang/C101HW18020629.RRC.PE150.20180528.P101HW18020629-01-01/RawData/H3_KO_1_2.fq.gz',required=True)
parser.add_argument('--rate',help='输入要剩余多少比例的rRNA,请输入数字')
parser.add_argument('--outdir',help="project analysis dir",required=True)
parser.add_argument('--sample',help="输入样品名")
argv = vars(parser.parse_args())


leftseq  = argv['leftseq'].strip()
rightseq = argv['rightseq'].strip() 
rate = argv['rate'].strip() 
outdir = argv['outdir'].strip()
sample = argv['sample'].strip()

os.system('if [ ! -d '+ outdir +' ];then\nmkdir -p '+ outdir +'\nfi')
#if not os.path.exists(outdir):
#    os.makedirs(outdir)
os.system("/PUBLIC/software/RNA/bowtie2-2.2.3/bowtie2 -p 10 -x /PUBLIC/database/RNA/rRNA/Eukaryote_rRNA_NR  --1"+leftseq+" --2"+rightseq+" --un-conc-gz "+outdir+"/"+sample+"_unmapped.gz --al-conc-gz "+ outdir+"/"+sample+"_mapped.gz >"+outdir+"/"+sample+"_fq.bowtie 2>"+outdir+"/"+sample+"_rRNA.bowtie.stat")
os.system("tail -n 1 "+ outdir+"/"+sample+"_rRNA.bowtie.stat | awk '{print $1}' | awk -F '%' '{print $1}'> "+outdir+"/"+sample+"_rate.txt")
os.system("head  -n 1 "+ outdir+"/"+sample+"_rRNA.bowtie.stat | awk -F ' ' '{print $1}'> "+outdir+"/"+sample+"_reads_num.txt")
rate_mapped = open(outdir+"/"+sample+"_rate.txt")
r = rate_mapped.readline().strip()
rate_map = float(r)
reads_num = open(outdir+"/"+sample+"_reads_num.txt")
reads_number= reads_num.readline().strip()                     #原始fq序列总的reads数

file = gzip.open(outdir+"/"+sample+"_mapped.1.gz", 'rb')
file2 = gzip.open(outdir+"/"+sample+"_mapped.2.gz", 'rb') 
f = file.readlines()
f2 = file2.readlines()
num1 = 4* (rate_map-float(rate))/100*float(reads_number)                      #要删除的reads数目
dic ={}
dic2 = {}
if float(rate_map) <= 15:
    print "rRNA低于15%，未达到去除rRNA的阈值！"
else:
    a = 0
    for i in range(a,int(num1),4):
        id = f[i].strip().replace("@","")
	dic[id]= ''
	id2 = f2[i].strip().replace("@","")
	dic2[id2]= ''
file.close()
file2.close()

file7 = open(outdir+"/"+sample+"_1.after.fq",'w')                   
file8 = open(outdir+"/"+sample+"_2.after.fq",'w')
R1 = HTSeq.FastqReader(leftseq)
R2 = HTSeq.FastqReader(rightseq)

for r1,r2 in itertools.izip(R1,R2):
	if r1.name not in dic:
		r1.write_to_fastq_file(file7)
	if r2.name not in dic2:
		r2.write_to_fastq_file(file8)
	else:
		continue

#os.system("gzip "+outdir+"/"+sample+"_1.after.fq")
#os.system("gzip "+outdir+"/"+sample+"_2.after.fq")
#os.system("rm "+outdir+"/"+sample+"_mapped.1.gz")
#os.system("rm "+outdir+"/"+sample+"_mapped.2.gz")
#os.system("rm "+outdir+"/"+sample+"_mapped.1_1.gz")
#os.system("rm "+outdir+"/"+sample+"_mapped.2_1.gz")
#os.system("rm "+outdir+"/"+sample+"_unmapped.1.gz")
#os.system("rm "+outdir+"/"+sample+"_unmapped.2.gz")
#os.system("rm "+outdir+"/"+sample+"_fq.bowtie")
#os.system("rm "+outdir+"/"+sample+"_rate.txt")
