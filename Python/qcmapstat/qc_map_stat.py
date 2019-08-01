#!/HWPROJ2/HW/wangjinhao/software/anaconda2/envs/reseq/bin/python
from collections import defaultdict
import argparse
import os,sys,fnmatch
import time
class Qc_Map():
    def __init__(self,pjdir,pjno,type):
        self.pj_type = type
        self.pjdir = pjdir
        self.pjno = pjno
        self.file_time = ''
        self.qc_list = self.get_qc_list()
        self.pj_info = self.get_pjinfo()
    
    def get_qc_list(self):
        '''
        HKG33DSXX_L3  patient  sample  lib  novo  index  path  type  sex  path_L
        '''
        '''
        HKG33DSXX_L3  sample  lib  novo  index  PE150  path  1
        '''
        qc_list = {}
        if self.pj_type == 'qc' :
            mapfile_all_file = '{}/mapfile_all'.format(self.pjdir)
            self.file_time = self.get_time(mapfile_all_file)
            f_in = open(mapfile_all_file)
            for line in f_in:
                line_list = line.strip().split('\t')
                if line_list[0][0] != '#':
                    try:
                        runid = line_list[3].split('_')[-1][1:10]
                    except:
                        continue
                    if line_list[5][:4] == 'lane':
                        line_list[5] = line_list[5][4:]
                    key = '_'.join([line_list[1],line_list[2],runid,'L{}'.format(line_list[5])])
                    qc_list[key] = ['_'.join([runid,'L{}'.format(line_list[5])]),line_list[1],line_list[1],line_list[2],line_list[2]]
            f_in.close()
        else:
            qc_list_file = '{}/qc_list'.format(self.pjdir)
            self.file_time = self.get_time(qc_list_file)
            f_in = open(qc_list_file)
            for line in f_in:
                line_list = line.strip().split('\t')
                if self.pj_type == 'animal':
                    line_list = line_list[:2] + line_list[1:]
                key = '{}_{}_{}'.format(line_list[2],line_list[4],line_list[0])
                #if not qc_list.has_key(key):
                qc_list[key] = line_list
            f_in.close()
        return qc_list

    def get_pjinfo(self):
        try:
            f_in = open('{}/pn.txt'.format(self.pjdir))
            line_list = f_in.readline().strip().split('\t')
            line_list = [line_list[0],line_list[1],'-'.join(line_list[2].split('-')[:2])]
            pj_info = [self.file_time] + line_list
        except:
            pj_info = [self.file_time,'-','-',self.pjno]
        return pj_info

    def get_stat(self):
        my_stat = {}
        for s_n_r_l in self.qc_list.keys():
            sample = self.qc_list[s_n_r_l][2]
            if self.pj_type == 'qc':
                novo = '-'
            else:
                novo = self.qc_list[s_n_r_l][4]
            qc_stat_list = self.get_qc_stat(self.qc_list[s_n_r_l])
            map_stat_list = self.get_map_stat(self.qc_list[s_n_r_l])
            #print self.pjdir
            #print self.pj_info
            #print s_n_r_l
            #print self.qc_list[s_n_r_l]
            #map_stat_list = self.get_map_stat('{0}/Alnstat/{1}/{1}_mapping_coverage.txt'.format(self.pjdir,sample))
            info = [self.pjdir,self.pj_info[0],s_n_r_l,self.pj_info[3],sample,novo,self.qc_list[s_n_r_l][3],self.qc_list[s_n_r_l][0].split('_')[0],self.qc_list[s_n_r_l][0].split('_')[1][1:]]
            my_stat[s_n_r_l] = info + qc_stat_list + map_stat_list

        print '#Pjdir\tTime\tQC stat\tProject No.\tSample\tNovo ID\tLibrary ID\tRun\tLane\tRaw reads\tError Rate\tQ20\tQ30\tMapped\tDup(%)\tAverage depth'
        for key in sorted(my_stat.keys()):
            #print my_stat[key]
            print '\t'.join(my_stat[key])

    def get_qc_stat(self,value):
        stat_dic = defaultdict(int)
        statfile_human = '{0}/QC/{1}/{1}_{2}_{3}.stat'.format(self.pjdir,value[2],value[4],value[0])
        statfile_qc = '{0}/QC/{1}/{1}.stat'.format(self.pjdir,value[2])
        statfile_animal = '{0}/QC/{1}/{1}_{2}.stat'.format(self.pjdir,value[2],value[0])
        try:
            if self.pj_type == 'human' or self.pj_type == 'qc':
                try:
                    f_in = open(statfile_human)
                except:
                    f_in = open(statfile_qc)
            elif self.pj_type == 'animal':
                f_in = open(statfile_animal)

            for line in f_in:
                line_list = line.strip().split('\t')
                if line_list[0].startswith('Number of Reads'):
                    stat_dic['Raw_reads'] += int(line_list[1])
                elif line_list[0].startswith('Error of fq'):
                    stat_dic['Error_Rate'] += float(line_list[2].strip('%'))/2
                elif line_list[0].startswith('Q20 of fq'):
                    stat_dic['Q20'] += float(line_list[2].strip('%'))/2
                elif line_list[0].startswith('Q30 of fq'):
                    stat_dic['Q30'] += float(line_list[2].strip('%'))/2
                else:
                    continue
            f_in.close()
        except:
            stat_dic={'Raw_reads':'-','Error_Rate':'-','Q20':'-','Q30':'-'}
        
        stat_list = [str(stat_dic['Raw_reads']),str(stat_dic['Error_Rate']),str(stat_dic['Q20']),str(stat_dic['Q30'])]
        return stat_list

    def get_map_stat(self,value):
        stat_dic = defaultdict(int)
        if self.pj_type == 'human':
            try:
                stat_file = '{0}/Alnstat/{1}/{1}_mapping_coverage.txt'.format(self.pjdir,value[2])
                f_in = open(stat_file)
                for line in f_in:
                    line_list = line.strip().split('\t')
                    if line_list[0].startswith('Average_sequencing_depth'):
                        stat_dic['Average_depth'] = line_list[1]
                    elif line_list[0].startswith('Duplicate'):
                        start = line_list[1].index('(') + 1
                        end = line_list[1].index('%')
                        stat_dic['Dup'] = line_list[1][start:end]
                    elif line_list[0].startswith('Mapped'):
                        start = line_list[1].index('(') + 1
                        end = line_list[1].index('%')
                        stat_dic['Mapped'] = line_list[1][start:end]
                    else:
                        continue
                f_in.close()
            except:
                stat_dic = {'Mapped':'-','Dup':'-','Average_depth':'-'}
            
        elif self.pj_type == 'animal':
            try:
                depthstat_file = '{0}/Alnstat/{1}/depthWithN/summary.txt'.format(self.pjdir,value[2])
                f_in_depth = open(depthstat_file)
                for line in f_in_depth:
                    line_list = line.strip().split('\t')
                    if line_list[0].startswith('Average_sequencing_depth'):
                        stat_dic['Average_depth'] = line_list[1]
                    else:
                        continue
                f_in_depth.close()
            except:
                stat_dic['Average_depth']='-'
            
            try:
                mapped_file = '{0}/Alnstat/{1}/{1}.alninfo'.format(self.pjdir,value[2])
                f_in_map = open(mapped_file)
                for line in f_in_map:
                    line_list = line.strip().split('\t')
                    if line_list[0].startswith('mapping rate'):
                        stat_dic['Mapped'] = line_list[1].strip('%')
                    else:
                        continue
                f_in_map.close()
            except:
                stat_dic['Mapped'] = '-'

            try:
                for root,dirs,files in os.walk('{}/log'.format(self.pjdir),topdown=False):
                    for file in files:
                        if fnmatch.fnmatch(file, 'filter_rmdup_'+value[2]+'_e*'):
                            dupstat_file = os.path.join(root,file)
                            f_in_dup = open(dupstat_file)
                            dupstat_file_lines = f_in_dup.readlines()
                            dup = str(float(dupstat_file_lines[-1].split(' ')[5])*100)
                            stat_dic['Dup'] = dup
                            f_in_dup.close()
                        else:
                            continue
            except:
                stat_dic['Dup'] = '-'
        
        elif self.pj_type == 'qc':
            stat_dic = {'Mapped':'-','Dup':'-','Average_depth':'-'}

        stat_list = [stat_dic['Mapped'],stat_dic['Dup'],stat_dic['Average_depth']]
        return stat_list
    
    def get_time(self,file_path):
        file_time = os.popen('stat {} |tail -1'.format(file_path)).read().strip()[8:27]
        return file_time

class Start():
    def __init__(self,pwd,argv):
        self.pwd = pwd
        self.argv = argv
    def start(self):
        if self.argv[1] and (self.argv[1] == 'animal'):
            try:
                new_stat = Qc_Map(self.pwd,self.argv[2],'animal')
                new_stat.get_stat()
            except:
                if os.path.exists('pn.txt'):
                    new_stat = Qc_Map(self.pwd,self.argv[2],'human')
                    new_stat.get_stat()
                else:
                    print '#cannot make'
        elif self.argv[1] and (self.argv[1] == 'human'):
            if not os.path.exists('pn.txt'):
                try:
                    new_stat = Qc_Map(self.pwd,self.argv[2],'animal')
                    new_stat.get_stat()
                except:
                    new_stat = Qc_Map(self.pwd,self.argv[2],'human')
                    new_stat.get_stat()
            else:
                new_stat = Qc_Map(self.pwd,self.argv[2],'human')
                new_stat.get_stat()
        elif self.argv[1] and (self.argv[1] == 'qc'):
            new_stat = Qc_Map(self.pwd,self.argv[2],'qc')
            new_stat.get_stat()

if __name__ == "__main__":
    pwd = os.popen('pwd').read().strip('\n')
    starter = Start(pwd,sys.argv)
    starter.start()
