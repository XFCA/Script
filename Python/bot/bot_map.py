#!/usr/bin/python
#coding=utf-8
"""
Spyder Editor

This is a temporary script file.
"""
import sys
import os

def sample_read(file_name = 'sample_path'):
    sample_list = []
    f_in = open(file_name)
    for line in f_in:
        sample = line.strip('\n').split('\t')
        sample_list.append(sample)
    f_in.close()
    return sample_list

def shell_write(sample_list, pwd, rate = '15'):
    for sample in sample_list:
        print pwd
        print sample
        f_out = open(pwd + '/log/map_hm_'+ sample[0]+'.sh','w')
        f_out.write('rf=/RLNAS01/HW/project/seq/wangjinhao/data/b37/human_g1k_v37_decoy\n')
        f_out.write('fc='+sample[1]+'\n')
        f_out.write('pwd_path='+pwd+'\n')
        f_out.write('novo='+sample[0]+'\n')
        f_out.write("if [ ! -d ${pwd_path}'/'${novo} ];then\nmkdir -p ${pwd_path}'/'${novo}\nfi\n")
        f_out.write('/PUBLIC/software/RNA/bowtie2-2.2.3/bowtie2 -p 10 -x ${rf}  --1 ${fc}/${novo}/${novo}_L1_1.fq.gz,${fc}/${novo}/${novo}_L2_1.fq.gz --2 ${fc}/${novo}/${novo}_L1_2.fq.gz,${fc}/${novo}/${novo}_L2_2.fq.gz -S ${pwd_path}/${novo}/${novo}.sam 2>${pwd_path}/${novo}/${novo}.stat\n')

        f_out.write('echo 全部完成。\n')
        f_out.close()
        os.system('cd '+ pwd +'/log;qsub -cwd -l vf=10G '+ 'map_hm_'+ sample[0]+'.sh')

pwd = os.popen('pwd').read().strip('\n')
os.system('if [ ! -d "log" ];then\nmkdir log\nfi')
sample_list = sample_read()
print pwd
print sample_list
try:
    rate = sys.argv[1]
    shell_write(sample_list, pwd, rate)
except:
    shell_write(sample_list, pwd)

