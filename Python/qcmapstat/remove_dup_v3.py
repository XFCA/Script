import time,math,random,re
import sys
from collections import defaultdict
import os
import argparse

class sjm_job():
    def __init__(self,job_dir):
        self.pre = 'rmdup'
        self.job_dir = job_dir
        self.log_dir = job_dir+'/log'
        self.job_sh = []
        self.job_order = []

    def add_sh(self,name,memory,cmd_shell):
        my_job = 'job_begin\n'
        my_job += '  name {}\n'.format(name)
        my_job += '  memory {}\n'.format(memory)
        my_job += '  status waiting\n'
        my_job += '  sched_options -V -cwd -q hwus1.q,hwus2.q,joyce.q,novo.q,all.q\n'
        my_job += '  cmd_begin\n'
        my_job += '    sh {}\n'.format(cmd_shell)
        my_job += '  cmd_end\n'
        my_job += 'job_end\n'
        self.job_sh.append(my_job)

    def add_order(self,name,after_name):
        self.job_order.append('order {} after {}\n'.format(name,after_name))

    def job_writer(self):
        if not os.path.exists(self.log_dir):
            os.mkdir(self.log_dir)
        job_name = '{}/{}_{}.job'.format(self.job_dir,self.pre,time.strftime('%Y%m%d%H%M%S'))
        f_out = open(job_name,'w')
        f_out.write('log_dir {}\n'.format(self.log_dir))
        for job in self.job_sh:
            f_out.write(job)
        for order in self.job_order:
            f_out.write(order)
        f_out.close()

class MyList():
    def __init__(self,root_dir,work_dir):
        self.root_dir = root_dir
        self.work_dir = work_dir
        self.max = '0.15'
        self.min = '0.10'
        self.threshold = 0.2
        self.sl_header = ['lane','patient','sample','lib','novo','index','path','type','sex']
        self.myheader = []
        self.sl = {}
        self.slac_dup = {}
        self.sltmp = {}

    def sample_list_read(self,mylist):
        header = self.sl_header
        f_in = open(mylist)
        self.myheader = f_in.readline().strip().split('\t')[:9]
        for line in f_in :
            line_list = line.strip().split('\t')
            try:
                runid = line_list[6].strip().split('_')[-1][1:10]#190129_A00564_0075_AHGY2NDSXX-new
                #print runid
            except:
                runid = '-'
            key = '{}_{}_{}_L{}'.format(line_list[2],line_list[4],runid,line_list[0])
            
            if line[0] == '#':
               self.sltmp[key] = {header[x]:line_list[x] for x in range(len(header))}
               continue
            self.sltmp[key] = {header[x]:line_list[x] for x in range(len(header))}
            self.sl[key] = {header[x]:line_list[x] for x in range(len(header))}
            self.sltmp[key],self.sl[key]['run'] = runid,runid
            if len(line_list) >= 10:
                self.slac_dup[key] = float(line_list[9])
            else:
                self.slac_dup[key] = '-'
        f_in.close()

    def dup_mk(self,sample_list,min_data=0):
        dup_stat = {}
        max1,min1,min_data,threshold = (float(self.max),float(self.min),int(min_data)*math.pow(10,9)/300,(self.threshold)*math.pow(10,9)/300)
        print 'Min Reads: {}\tThreshold: {}'.format(int(min_data),int(threshold))
        self.sample_list_read(sample_list)
        for key,sl_value in zip(self.sl.keys(),self.sl.values()):
            sample = sl_value['sample']
            qc_stat_file = '{0}/QC/{1}/{1}_{2}_{3}_L{4}.stat'.format(self.work_dir,sl_value['sample'],sl_value['novo'],sl_value['run'],sl_value['lane'])
            map_stat_file = '{0}/Alnstat/{1}/{1}_mapping_coverage.txt'.format(self.work_dir,sl_value['sample'])
            raw,clean_reads = self.get_qc(qc_stat_file)
            dup = self.get_map(map_stat_file)
            self.sl[key]['raw'] = raw
            self.sl[key]['clean_r'] = clean_reads
            #print dup_reads
            if dup_stat.has_key(sample):
                dup_stat[sample]['clean_r'] += clean_reads
                dup_stat[sample]['raw'] += raw
            else:
                dup_stat[sample] = {'raw':raw,'clean_r':clean_reads,'old_d':dup}#,'raw_ndr':raw - clean_reads*dup,
                                    #'clean_ndr':clean_reads*(1-dup),'old_d':dup,'dup_dr':clean_reads*dup}
            dup_stat[sample]['raw_ndr'] = dup_stat[sample]['raw'] - dup_stat[sample]['clean_r']*dup_stat[sample]['old_d']
            dup_stat[sample]['clean_ndr'] = dup_stat[sample]['clean_r']*(1-dup_stat[sample]['old_d'])
            dup_stat[sample]['dup_dr'] = dup_stat[sample]['clean_r']*dup_stat[sample]['old_d']
            dup_stat[sample]['raw_other'] = dup_stat[sample]['raw'] - dup_stat[sample]['clean_r']
        
        print 'Id\tF\tDup\tThisReads\tSample'
        for key in self.sl.keys():
            sample = self.sl[key]['sample']
            if self.slac_dup[key] != '-':
                if dup_stat[sample]['old_d'] <= self.slac_dup[key]:
                    self.sl[key]['new_d'] = dup_stat[sample]['old_d']
                    self.sl[key]['f'] = '-'
                    print 'a\t{}\t{:.4f}\t{}\t{}'.format('-',self.sl[key]['new_d'],int(self.sl[key]['raw']),key)
                else:
                    self.sl[key]['new_d'] = self.slac_dup[key]
                    self.sl[key]['f'] = (dup_stat[sample]['old_d']-self.slac_dup[key])/(dup_stat[sample]['old_d']*(1-self.slac_dup[key]))
                    this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                    print 'b\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                    
            else:
                #mindr = dup_stat[sample]['clean_r']*float(self.min)
                #maxdr = dup_stat[sample]['clean_r']*float(self.max)
                min_clean_reads = min_data - dup_stat[sample]['raw_other'] + threshold
                mind1,maxd1 = float(self.min),float(self.max)
                if dup_stat[sample]['clean_r'] <= min_clean_reads:
                    #print '{}\t{}'.format(dup_stat[sample]['clean_r'],min_clean_reads)
                    self.sl[key]['new_d'] = dup_stat[sample]['old_d']
                    self.sl[key]['f'] = '-'
                    print 'c\t{}\t{:.4f}\t{}\t{}'.format('-',self.sl[key]['new_d'],int(self.sl[key]['raw']),key)
                elif dup_stat[sample]['old_d'] <= maxd1:
                    self.sl[key]['new_d'] = dup_stat[sample]['old_d']
                    self.sl[key]['f'] = '-'
                    print 'd\t{}\t{:.4f}\t{}\t{}'.format('-',self.sl[key]['new_d'],int(self.sl[key]['raw']),key)
                elif dup_stat[sample]['clean_ndr'] >= min_clean_reads:
                    new_d = float(random.randint(int(mind1*math.pow(10,4)),int(maxd1*math.pow(10,4))))/math.pow(10,4)
                    self.sl[key]['new_d'] = new_d
                    self.sl[key]['f'] = (dup_stat[sample]['old_d']-new_d)/(dup_stat[sample]['old_d']*(1-new_d))
                    this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                    print 'e\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                else:
                    #min data dup
                    mind2 = float(dup_stat[sample]['clean_r']*(dup_stat[sample]['old_d']-1)+min_clean_reads)/min_clean_reads
                    if mind2 <= mind1:
                        #print 'mind2 {}\tmind1 {}\t old_dup {} \tmaxd1 {}'.format(mind2,mind1,dup_stat[sample]['old_d'],maxd1)
                        new_d = float(random.randint(int(mind1*math.pow(10,4)),int(maxd1*math.pow(10,4))))/math.pow(10,4)
                        self.sl[key]['new_d'] = new_d
                        self.sl[key]['f'] = (dup_stat[sample]['old_d']-new_d)/(dup_stat[sample]['old_d']*(1-new_d))
                        this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                        print 'f\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                    else:
                        self.sl[key]['new_d'] = mind2
                        self.sl[key]['f'] = (dup_stat[sample]['old_d']-mind2)/(dup_stat[sample]['old_d']*(1-mind2))
                        this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                        print 'g\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                '''
                mindr = dup_stat[sample]['clean_ndr']*float(self.min)/(1-float(self.min))
                maxdr = dup_stat[sample]['clean_ndr']*float(self.max)/(1-float(self.max))
                min_clean_nodup = min_clean_reads*(1-float(self.max))
                max_clean_nodup = min_clean_reads*(1-float(self.min))
                if dup_stat[sample]['raw'] <=  min_data + self.power*math.pow(10,9)/300:
                    self.sl[key]['new_d'] = dup_stat[sample]['old_d']
                    self.sl[key]['f'] = '-'
                    print 'c\t{}\t{:.4f}\t{}\t{}'.format('-',self.sl[key]['new_d'],'-\t',key)
                elif dup_stat[sample]['dup_dr'] <= mindr + self.power*math.pow(10,9)/300:
                    self.sl[key]['new_d'] = dup_stat[sample]['old_d']
                    self.sl[key]['f'] = '-'
                    print 'd\t{}\t{:.4f}\t{}\t{}'.format('-',self.sl[key]['new_d'],'-\t',key)
                elif dup_stat[sample]['dup_dr'] <= maxdr:
                    new_dup_dr = random.randint(int(mindr),int(dup_stat[sample]['dup_dr']))
                    self.sl[key]['new_d'] = float(new_dup_dr)/(new_dup_dr+dup_stat[sample]['clean_ndr'])
                    self.sl[key]['f'] = float(dup_stat[sample]['dup_dr'] - new_dup_dr)/dup_stat[sample]['dup_dr']
                    this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                    print 'e\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                elif dup_stat[sample]['dup_dr'] > maxdr:
                    new_dup_dr = random.randint(int(mindr),int(maxdr))
                    self.sl[key]['new_d'] = float(new_dup_dr)/(new_dup_dr+dup_stat[sample]['clean_ndr'])
                    self.sl[key]['f'] = float(dup_stat[sample]['dup_dr'] - new_dup_dr)/dup_stat[sample]['dup_dr']
                    this_reads = self.sl[key]['raw'] - dup_stat[sample]['old_d']*self.sl[key]['f']*self.sl[key]['clean_r']
                    print 'f\t{:.4f}\t{:.4f}\t{}\t{}'.format(self.sl[key]['f'],self.sl[key]['new_d'],int(this_reads),key)
                '''
                    
        self.sldupwriter()
    

    def get_qc(self,qc_stat_file):
        try:
            f_in = open(qc_stat_file)
            #print qc_stat_file
            for line in f_in:
                line_list = line.strip().split('\t')
                #print line_list
                if line_list[0].startswith('Number of Reads'):
                    #print line_list
                    raw,clean = int(line_list[1]),int(line_list[2])
                else:
                    #print raw,clean
                    continue
            f_in.close()
        except:
            raw,clean = False,False
        return (raw,clean)

    def get_map(self,map_stat_file):
        try:
            f_in = open(map_stat_file)
            for line in f_in:
                line_list = line.strip().split('\t')
                if line_list[0].startswith('Duplicate'):
                    mygroup = re.search('\((.*)%\)',line_list[1])
                    
                    dup = float(mygroup.group(1))/100
            f_in.close()
        except:
            dup=False
        return dup

    def sldupwriter(self):
        f_out = open('{}/sample_list_dup'.format(self.root_dir),'w')
        f_out.write('\t'.join(self.myheader)+'\tf\tdup\n')
        for key,tmp_value in zip(self.sltmp.keys(),self.sltmp.values()):
            if self.sl.has_key(key):
                sl_value = self.sl[key]
                line = '\t'.join([str(sl_value[x]) for x in self.sl_header] + [str(sl_value['f']),str(sl_value['new_d'])])
            else:
                line = '\t'.join([tmp_value[x] for x in self.sl_header])
            f_out.write(line + '\n')
        f_out.close()

    def relister(self,sl_dup):
        new_path = self.root_dir
        f_in = open(sl_dup)
        f_out = open('{}/sample_list_new'.format(self.root_dir),'w')
        for line in f_in:
            line_list = line.split('\t')
            if line[0] == '#':
                f_out.write(line)
            elif (len(line_list)>=10) and (line_list[9] == '-'):
                f_out.write(line)
            else:
                line_list[6] = os.path.join(self.root_dir,line_list[6].strip().split('/')[-1])
                f_out.write('\t'.join(line_list))
        f_in.close()
        f_out.close()
        
class remove_dup():
    def __init__(self,root_dir):
        self.root_dir = root_dir
        self.sample_list = ''
        self.work_dir = self.mk_work_dir(root_dir)
        self.lister = MyList(root_dir,self.work_dir)
        self.new_job = sjm_job(root_dir)
    
    def sample_list_read(self,sample_list_dup='sample_list_dup'):
        '''
        sample_list{fc:{f:....,}}
        '''
        
        sample_list = {}
        f_in = open(self.root_dir+'/'+sample_list_dup)
        for line in f_in:
            line_list = line.strip().split('\t')
            
            if line_list[0][0] == '#':
                continue
            if (len(line_list) >= 10) and (line_list[9] == '-'): 
                continue
            fc = line_list[6].split('/')[-1]
            dup_f = line_list[9].strip()
            if sample_list.has_key(fc):
                if sample_list[fc].has_key(dup_f):
                    sample_list[fc][dup_f].append(line_list)
                else:
                    sample_list[fc][dup_f] = [line_list]
            else:
                sample_list[fc] = {dup_f:[line_list]}
        f_in.close()
        
        for key,value in zip(sample_list.keys(),sample_list.values()):
            if not os.path.exists(self.root_dir+'/'+key):
                os.mkdir(self.root_dir+'/'+key)
                os.mkdir(self.root_dir+'/'+key+'/tmp')
            for key2,value2 in zip(value.keys(),value.values()):
                #print value2
                r_dir = self.root_dir+'/'+key+'/tmp/rmdup_'+key2[2:]
                if not os.path.exists(r_dir):
                    os.mkdir(r_dir)
                f_out = open('{}/sample_list_{}'.format(r_dir,key2[2:]),'w')
                for sample in value2:
                    f_out.write('{}\n'.format('\t'.join(sample)))
                f_out.close()
        return sample_list
    
    def mk_work_dir(self,root_dir):
        dir_list = root_dir.split('/')
        return '/'.join(dir_list[:-1])

    def mk_sh(self,sample_list_dup):
        self.sample_list = self.sample_list_read(sample_list_dup)
        for key,value in zip(self.sample_list.keys(),self.sample_list.values()):
            for key2,value2 in zip(value.keys(),value.values()):
                key2_dir = '{}/{}/tmp/rmdup_{}'.format(self.root_dir,key,key2[2:])
                f_out = open('{}/mk_redup_{}.sh'.format(key2_dir,key2[2:]),'w')
                f_out.write('python /ifs/TJPROJ3/HWUS/USER/wangjinhao/software/process/mk_rmdup.py \\\n')
                #f_out.write('python /PUBLIC/source/HW/CANCER/HW_v1.0/script/mk_rmdup.py \\\n')
                f_out.write('-li {}/sample_list_{} \\\n'.format(key2_dir,key2[2:]))
                f_out.write('-p {} \\\n'.format(self.work_dir))
                f_out.write('-f {} \\\n'.format(key2))
                f_out.write('-o {}\n'.format(key2_dir))
                f_out.close()

                self.new_job.add_sh('mk_redup_{}_{}'.format(key2[2:],key),'3G','{}/mk_redup_{}.sh'.format(key2_dir,key2[2:]))
                for sample in value2:
                    self.new_job.add_sh('{}_L{}_1_{}'.format(sample[2],sample[0],key),'3G','{}/RemoveDup/{}/{}ReDup_L{}_1.sh'.format(key2_dir,sample[3],sample[2],sample[0]))
                    self.new_job.add_sh('{}_L{}_2_{}'.format(sample[2],sample[0],key),'3G','{}/RemoveDup/{}/{}ReDup_L{}_2.sh'.format(key2_dir,sample[3],sample[2],sample[0]))
                    self.new_job.add_order('{}_L{}_1_{}'.format(sample[2],sample[0],key),'mk_redup_{}_{}'.format(key2[2:],key))
                    self.new_job.add_order('{}_L{}_2_{}'.format(sample[2],sample[0],key),'mk_redup_{}_{}'.format(key2[2:],key))
        
    def job_mk(self):
        self.new_job.job_writer()
    
    def link_sh(self):
        for key,value in zip(self.sample_list.keys(),self.sample_list.values()):
            f_out = open('{}/{}/mk_link.sh'.format(self.root_dir,key),'w')
            for value2 in value.values():
                for sample in value2:
                    fc_dir = '{}/{}/'.format(self.root_dir,key)
                    lib_dir = '{}/{}/tmp/rmdup_{}/RemoveDup/{}'.format(self.root_dir,key,sample[9][2:],sample[3])
                    f_out.write('ln -s {} {}\n'.format(lib_dir,fc_dir))
                    self.new_job.add_order('ln_{}'.format(key),'{}_L{}_1_{}'.format(sample[2],sample[0],key))
                    self.new_job.add_order('ln_{}'.format(key),'{}_L{}_2_{}'.format(sample[2],sample[0],key))

            self.new_job.add_sh('ln_{}'.format(key),'1G','{}/{}/mk_link.sh'.format(self.root_dir,key))        
            f_out.close()

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="remove duplicate script")
    subpar = argparser.add_subparsers(help='son command')
    sub_start = subpar.add_parser('start', help='make start script and new sample_list file')
    sub_start.set_defaults(cmd='start')
    sub_start.add_argument('-s','--sample_list_dup',help=r'the sample_list_dup you need to input. default:"sample_list_dup"',default='sample_list_dup')

    sub_mkdup = subpar.add_parser('mkdup', help=r'Calculate the value of "f" and generate "sample_list_dup"(keep data size to make dup at 0.10~0.15)')
    sub_mkdup.set_defaults(cmd='mkdup')
    sub_mkdup.add_argument('-d','--data_size',default='0',help=r'Minimum amount of data. default:"0"')
    sub_mkdup.add_argument('-s','--sample_list',help=r'the sample_list you need to input. default:"sample_list"',default='sample_list')

    sub_relist = subpar.add_parser('relist', help=r'make new sample_list file')
    sub_relist.set_defaults(cmd='relist')
    sub_relist.add_argument('-s','--sample_list_dup',help=r'the sample_list_dup you need to input. default:"sample_list_dup"',default='sample_list_dup')
    #mkdup2argv = sub_mkdup2.parse_args()
    #print mkdup2argv.sample_list
    myargv = argparser.parse_args()
    print myargv

    pwd = os.popen('pwd').read().strip('\n')
    new_rmdup = remove_dup(pwd)
    if myargv.cmd == 'start':
        new_rmdup.lister.relister(myargv.sample_list_dup)
        new_rmdup.mk_sh(myargv.sample_list_dup)
        new_rmdup.link_sh()
        new_rmdup.job_mk()
    elif myargv.cmd == 'mkdup':
        new_rmdup.lister.dup_mk(myargv.sample_list,myargv.data_size)
    elif myargv.cmd == 'relist':
        new_rmdup.lister.relister(myargv.sample_list_dup) 
