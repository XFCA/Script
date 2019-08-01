# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 18:15:43 2019

@author: jinhaow
"""
import os

class bowtie2_map():
    def __init__(self,pwd):
        self.file_name = 'sample_path'
        self.sample_list = self.sample_read()
        self.pwd = pwd
        self.mem = '10G'
        self.rf = '/RLNAS01/HW/project/seq/wangjinhao/data/b37/human_g1k_v37_decoy'
    
    def sample_read(self):
        """
        sample_path: smaple lib lane path
        sample_list: {sample:{lib lf rg}}
        
        """
        sample_list = {}
        f_in = open(self.file_name)
        for line in f_in:
            sample = line.strip('\n').split('\t')
            if sample[0] in sample_list.keys():
                sample_list[sample[0]]['lft'] += ','+sample[3]+'/'+sample[1]+'/'+sample[1]+'_L'+sample[2]+'_1.fq.gz'
                sample_list[sample[0]]['rgt'] += ','+sample[3]+'/'+sample[1]+'/'+sample[1]+'_L'+sample[2]+'_2.fq.gz'
            else:
                sample_list[sample[0]] = {}
                sample_list[sample[0]]['lib'] = sample[1]
                sample_list[sample[0]]['lft'] = sample[3]+'/'+sample[1]+'/'+sample[1]+'_L'+sample[2]+'_1.fq.gz'
                sample_list[sample[0]]['rgt'] = sample[3]+'/'+sample[1]+'/'+sample[1]+'_L'+sample[2]+'_2.fq.gz'
        f_in.close()
        return sample_list
    
    def shell_write(self):
        for key,value in zip(self.sample_list.keys(),self.sample_list.values()):
            print self.pwd
            print key
            f_out = open(self.pwd + '/log/map_hm_'+ value['lib']+'.sh','w')
            f_out.write('rf='+self.rf+'\n')
            f_out.write('lft='+value['lft']+'\n')
            f_out.write('rgt='+value['rgt']+'\n')
            f_out.write('pwd_path='+self.pwd+'\n')
            f_out.write('novo='+value['lib']+'\n')
            f_out.write('sample='+key+'\n')
            f_out.write("if [ ! -d ${pwd_path}'/'${novo} ];then\nmkdir -p ${pwd_path}'/'${novo}\nfi\n")
            f_out.write('/PUBLIC/software/RNA/bowtie2-2.2.3/bowtie2 -p 10 -x ${rf}  --1 ${lft} --2 ${rgt} -S ${pwd_path}/${novo}/${novo}.sam 2>${pwd_path}/${novo}/${novo}.stat\n')
            f_out.write('echo 全部完成。\n')
            f_out.close()
            #os.system('cd '+ self.pwd +'/log;qsub -cwd -l vf='+self.mem+' '+ 'map_hm_'+ value['lib']+'.sh')

    def xls_make_wrte(self):
        mk_xls = r"""

import os
class xls_stat():
    def __init__(self,pwd):
        self.file_name = 'sample_path'
        self.sample_list = self.sample_read()
        self.pwd = pwd
        
    def sample_read(self):
        '''
        sample_path: smaple lib lane path
        sample_list: {sample{lib}}
        
        '''
        sample_list = {}
        f_in = open(self.file_name)
        for line in f_in:
            sample = line.strip('\n').split('\t')
            if sample[0] not in sample_list.keys():
                sample_list[sample[0]]={}
                sample_list[sample[0]]['lib'] = sample[1]
        
        f_in.close()
        return sample_list
        

    def get_list(self,lib):
        
        f_out = open(self.pwd+'/'+lib+'/'+lib+'.stat')
        lines = f_out.readlines()
        total = lines[0].split(' ')[0]
        per = lines[-1].split(' ')[0]
        return [total,per]
    
    def xls_make(self):
        '''
        xls: sample lib all mapped per
        '''
        f_out = open(self.pwd+'/mapping_stat.xls','w')
        f_out.write('sample\tlibID\treads\tmapped\tper\n')
        for key,value in zip(self.sample_list.keys(),self.sample_list.values()):
            sample = self.get_list(value['lib'])
            mapped = int(round(float(sample[0])*float(sample[1][:-1])/100))
            f_out.write(key +'\t' + value['lib']+'\t' + sample[0]+'\t' +str(mapped)+'\t'+ sample[1]+'\n')
        f_out.close()
        
        
pwd = os.popen('pwd').read().strip('\n')
new_xls = xls_stat(pwd)
new_xls.xls_make()
        
        """

        f_out = open('mk_xls.py','w')
        f_out.write(mk_xls)
        f_out.close()
        
               
pwd = os.popen('pwd').read().strip('\n')
os.system('if [ ! -d "log" ];then\nmkdir log\nfi')
new_map = bowtie2_map(pwd)
print pwd
print new_map.sample_list
new_map.shell_write()
new_map.xls_make_wrte()
