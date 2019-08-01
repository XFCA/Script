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
        s0_list = sample[0].split('_')
        sample.append('_'.join(s0_list[:-1]))
        sample.append(sample[1].split('/')[-1])
        sample_list.append(sample)
    f_in.close()
    return sample_list

def shell_write(sample_list, pwd, rate = '15'):
    for sample in sample_list:
        f_out = open(pwd + '/log/rm_'+ sample[0]+'.sh','w')
        f_out.write('sample_lane='+ sample[0] + '\n')
        f_out.write('fc_path='+ sample[1] + '\n')
        f_out.write('sample='+ sample[2] + '\n')
        f_out.write('fc='+ sample[3] + '\n')
        f_out.write('root_path='+ pwd + '\n')
        f_out.write('rate='+ rate + '\n')
        f_out.write('python /HWPROJ2/HW/wangjinhao/software/script/rRNA_remove_lastv_all.py \\\n')
        f_out.write('\t--sample ${sample_lane} \\\n')
        f_out.write('\t--leftseq ${fc_path}/${sample}/${sample_lane}_1.fq.gz \\\n')
        f_out.write('\t--rightseq ${fc_path}/${sample}/${sample_lane}_2.fq.gz \\\n')
        f_out.write('\t--outdir ${root_path}/${fc}/${sample} \\\n')
        f_out.write('\t--rate ${rate}\n')
        f_out.write('echo 去rRNA完成，正在压缩。。\n')
        f_out.write('cd ${root_path}/${fc}/${sample}\n')
        f_out.write('rename .after.fq .fq ${sample_lane}*.after.fq\n')
        f_out.write('gzip ${sample_lane}*.fq\n')
        f_out.write('echo 压缩完成，正在创建软连接。。\n')
        f_out.write('ln -s ${fc_path}/${sample}/${sample_lane}*.adap* ./\n')
        f_out.write('echo 全部完成。\n')
        f_out.close()
        os.system('cd '+ pwd +'/log;qsub -cwd -l vf=10G '+ 'rm_'+ sample[0]+'.sh')

pwd = os.popen('pwd').read().strip('\n')
os.system('if [ ! -d "log" ];then\nmkdir log\nfi')
sample_list = sample_read()
try:
    rate = sys.argv[1]
    shell_write(sample_list, pwd, rate)
except:
    shell_write(sample_list, pwd)
