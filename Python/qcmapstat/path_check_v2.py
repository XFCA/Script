#!/HWPROJ2/HW/wangjinhao/software/anaconda2/envs/reseq/bin/python
import os,fnmatch,subprocess,time,argparse
from collections import defaultdict
import sys
import re
import qc_map_stat_v2
class Path_Search():
    def __init__(self,pjno=''):
        self.rootdir = ''
        self.pjtemplate = "[PX]202SC.{8}-Z?\d{2}"#P202SC19060005-01-01 X202SC19031796-Z01-F003
        self.pjdir = ['cancer','disease','seq','reseq']
        self.pj_analysis = ['cancer','disease','reseq']
        #self.pj_analysis = ['seq']
        self.qc_map_stat_py = "/HWPROJ2/HW/wangjinhao/software/script/qcmapstat/qc_map_stat_v2.py"
        self.pjfound = set()
        self.pjnofound = set()
        self.pjno = pjno

    def pathsearch(self,rootdir):
        self.rootdir = rootdir
        if self.pjno != '':
            self.pjno = self.get_pj(self.pjno)
            self.pjpath = self.get_path()
        else:
            self.pjpath = self.get_path_all(self.pj_analysis)

        f_out = open('pjstat_{}.xls'.format(time.strftime('%Y%m%d%H%M%S')),'w')
        mydir = os.popen('pwd').read().strip()
        #pjstat_file = '{}/pjstat_{}.xls'.format(mydir,time.strftime('%Y%m%d%H%M%S'))
        mystart = qc_map_stat_v2.Start()
        mystart.out_file = f_out
        for mypj in sorted(self.pjpath.keys()):
            print '{} {}'.format(self.pjpath[mypj][0],self.pjpath[mypj][1])
            mystart.pwd = self.pjpath[mypj][0]
            mystart.argv={"type":self.pjpath[mypj][1],"pjno":self.pjpath[mypj][2]}
            mystart.start()
            mystart.myprint()
            #os.system('cd {} && python {} -t {} -p {} >> {}'.format(self.qc_map_stat_py,self.pjpath[mypj][0],self.pjpath[mypj][1],self.pjpath[mypj][2],pjstat_file))
            #subprocess.Popen('cd {} && python /HWPROJ2/HW/wangjinhao/software/script/qc_map_stat.py {} {}'.format(self.pjpath[mypj][0],self.pjpath[mypj][1],self.pjpath[mypj][2]),shell=True,stdout=f_out)
            #os.system('cd {} && QM.sh {} {}'.format(self.pjpath[mypj][0],self.pjpath[mypj][1],self.pjpath[mypj][2]))
            print '{} OVER!'.format(mypj)
        self.set_print()
        f_out.close()
        print "ALL OVER!"

    def get_pj(self,pjno_file):
        pj_set = set()
        f_in = open(pjno_file)
        for line in f_in:
            pjno = line.strip().split('\t')[0]
            #print pjno
            pj_set.add(pjno)
        f_in.close()
        return sorted(pj_set)

    def get_path(self):
        seq_list = self.get_ls('{}/seq'.format(self.rootdir))
        self.pjdir += ['seq/{}'.format(name) for name in seq_list]
        self.pjdir += ['seq/{}/auto_seq'.format(name) for name in seq_list]
        reseq_list = self.get_ls('{}/reseq'.format(self.rootdir))
        self.pjdir += ['reseq/{}'.format(name) for name in reseq_list]
        pj_path = defaultdict(str)
        for mytype in self.pjdir:
            mypj_path = os.path.join(self.rootdir,mytype)
            #print mypj_path
            new_pjtype = mypj_path.split('/')[4]
            #print new_pjtype
            if new_pjtype == 'cancer' or new_pjtype == 'disease':
                flag = 'human'
            elif new_pjtype == 'seq':
                flag = 'qc'
            elif new_pjtype == 'reseq':
                flag = 'animal'

            pj_list = self.get_ls(mypj_path)
            #print pj_list
            for pj in pj_list:
                #print pj
                if flag == 'human':
                    try:
                        f_in = open(os.path.join(mypj_path,pj,'pn.txt'))
                        pn_list = f_in.readline().strip().split('\t')
                        for mypj in self.pjno:
                            if fnmatch.fnmatch(pn_list[2],'*{}*'.format(mypj)) or fnmatch.fnmatch(pn_list[0],'*{}*'.format(mypj)):
                                pj_path[pj] = [os.path.join(mypj_path,pj),flag,mypj]
                        f_in.close()
                    except:
                        for mypj in self.pjno:
                            if fnmatch.fnmatch(pj,'*{}*'.format(mypj)):
                                pj_path[pj] = [os.path.join(mypj_path,pj),flag,mypj]
                else:
                    for mypj in self.pjno:
                        #print mypj
                        if fnmatch.fnmatch(pj,'*{}*'.format(mypj)):
                            #print mypj
                            pj_path[pj] = [os.path.join(mypj_path,pj),flag,mypj]
        #for x in self.pjdir:
            #print x
        all_set = set(self.pjno)
        self.pjfound = set([pj[2] for pj in pj_path.values()])
        self.pjnofound = all_set - self.pjfound
        return pj_path

    def get_path_all(self,dir_type):
        pj_path = []
        for mytype in dir_type:
            type_path = os.path.join(self.rootdir,mytype)
            if mytype == 'cancer' or  mytype == 'disease':
                flag = 'human'
                pj_list_1 = self.get_ls(type_path)
                human_path_1 = self.pj_path_filter(type_path,pj_list_1,flag)
                pj_path += human_path_1
                for mypath in [os.path.join(type_path,x) for x in pj_list_1]:
                    pj_list_2 = self.get_ls(mypath)
                    human_path_2= self.pj_path_filter(mypath,pj_list_2,flag)
                    pj_path += human_path_2

            elif mytype == 'reseq':
                flag = 'animal'
                #print 'kk'
                bi_name_dirs = self.get_ls(type_path)
                for bi_name_path in [os.path.join(type_path,x) for x in bi_name_dirs]:
                    pj_list_1 = self.get_ls(bi_name_path)
                    animal_1 = self.pj_path_filter(bi_name_path,pj_list_1,flag)
                    pj_path += animal_1
                    for mypath in [os.path.join(bi_name_path,x) for x in pj_list_1]:
                        pj_list_2 = self.get_ls(mypath)
                        animal_2= self.pj_path_filter(mypath,pj_list_2,flag)
                        pj_path += animal_2

            elif mytype == 'seq':
                flag = 'qc'
                bi_name_dirs = self.get_ls(type_path)
                for bi_name_path in [os.path.join(type_path,x) for x in bi_name_dirs]:
                    bi_auto_path = os.path.join(bi_name_path,"auto_seq")
                    if os.path.exists(bi_auto_path):
                        bi_auto_pj_list = self.get_ls(bi_auto_path)
                        bi_auto_pj_qc_1 = self.pj_path_filter(bi_auto_path,bi_auto_pj_list,flag)
                        pj_path += bi_auto_pj_qc_1
                        if bi_auto_pj_qc_1:
                            for auto_pj in bi_auto_pj_qc_1:
                                auto_pj_path_list = self.get_ls(auto_pj[0])
                                bi_auto_pj_qc_2 = self.pj_path_filter(auto_pj[0],auto_pj_path_list,flag)
                                pj_path += bi_auto_pj_qc_2

                    pj_list_1 = self.get_ls(bi_name_path)
                    qc_1 = self.pj_path_filter(bi_name_path,pj_list_1,flag)
                    pj_path += qc_1
                    if pj_list_1:
                        for mypath in [os.path.join(bi_name_path,x) for x in pj_list_1]:
                            pj_list_2 = self.get_ls(mypath)
                            qc_2= self.pj_path_filter(mypath,pj_list_2,flag)
                            pj_path += qc_2

        pj_path_dic = {x[0]:x for x in pj_path}
        return pj_path_dic


    def get_ls(self,mypath):
        my_list = os.popen('ls {}'.format(mypath)).read().strip().split('\n')
        if my_list[0] == mypath:
            my_list = []
        return my_list

    def pj_path_filter(self,rootdir,pj_list,flag):
        filter_pj_path = []
        if pj_list == []:
            return []

        if flag == 'human':
            for pj in pj_list:
                #print rootdir
                #print pj
                pj_path = os.path.join(rootdir,pj)
                if os.path.exists('{}/qc_list'.format(pj_path)):
                    try:
                        f_in = open('{}/pn.txt'.format(pj_path))
                        pj_no = re.search(self.template,f_in.readline().strip().split('\t')[2])
                        if pj_no :
                            pj_no = pj_no.group()
                        else:
                            pj_no = '-'
                        f_in.close()
                    except:
                        pj_no = '-'
                    filter_pj_path.append((pj_path,flag,pj_no))
                else:
                    continue
        elif flag == 'animal':
            for pj in pj_list:
                pj_path = os.path.join(rootdir,pj)
                #print 'a'
                if os.path.exists('{}/qc_list'.format(pj_path)):
                    #print 'b'
                    pj_no_raw = re.search(self.pjtemplate,pj)
                    if pj_no_raw :
                        #print 'c'
                        pj_no = pj_no_raw.group()
                        filter_pj_path.append((pj_path,flag,pj_no))
                        #print filter_pj_path[-1]
                    else:
                        pj_no_raw = re.search(self.pjtemplate,rootdir.split('/')[-1])
                        if pj_no_raw :
                            pj_no = pj_no_raw.group()
                            filter_pj_path.append((pj_path,flag,pj_no))
                        else:
                            filter_pj_path.append((pj_path,flag,'-'))
                else:
                    continue
        elif flag == 'qc':
            for pj in pj_list:
                pj_path = os.path.join(rootdir,pj)
                if os.path.exists('{}/mapfile_all'.format(pj_path)) and os.path.exists('{}/QC'.format(pj_path)):
                    pj_no_raw = re.search(self.pjtemplate,pj)
                    if pj_no_raw :
                        pj_no = pj_no_raw.group()
                        filter_pj_path.append((pj_path,flag,pj_no))
                    else:
                        pj_no_raw = re.search(self.pjtemplate,rootdir.split('/')[-1])
                        if pj_no_raw :
                            pj_no = pj_no_raw.group()
                            filter_pj_path.append((pj_path,flag,pj_no))
                        else:
                            filter_pj_path.append((pj_path,flag,'-'))
                else:
                    continue
        return filter_pj_path

    def set_print(self):
        print '\npj be found:'
        for pj in sorted(self.pjfound):
            print pj
        print '\npj not be found:'
        for pj in sorted(self.pjnofound):
            print pj

    def filter(self,filter_file):
        libs = set()
        libs_dic = {}
        qc_stat = set()
        lib_in = set()
        lib_notin = set()
        f_in = open(filter_file)
        for line in f_in:
            if line[0] == '#':
                continue
            line_list = line.strip().split('\t')
            libs_dic[line_list[2]] = line_list
            libs.add(line_list[2])
        f_in.close()
        f_in_stat = open('/HWPROJ2/HW/wangjinhao/PJSTAT/pjstat.xls')
        f_out_stat = open('/HWPROJ2/HW/wangjinhao/PJSTAT/pjstat_filter.xls','w')
        f_out_stat.write('Project No.\tNovo ID\tLibrary ID\tPC\tProductivity %\tLib QC Classification\tSample QC Classification\tRaw reads\tError Rate\tQ20\tQ30\tDup(%)\tAverage depth\n')
        for line in f_in_stat:
            if line[0] == '#':
                continue
            elif line.split('\t')[6] in libs:
                if line.split('\t')[2] in qc_stat:
                    continue
                else:
                    qc_stat.add(line.split('\t')[2])
                    lib_in.add(line.split('\t')[6])
                    f_out_stat.write('\t'.join(libs_dic[line.split('\t')[6]] + line.split('\t')[9:]))
            else:
                continue
        lib_notin = libs - lib_in
        for lib in sorted(lib_notin):
            f_out_stat.write('\t'.join(libs_dic[lib] + ['-','-','-','-','-','-\n']))
        f_in_stat.close()
        f_out_stat.close()

        print '\nlib be found:'
        for lib in sorted(lib_in):
            print lib
        print '\nlib not be found:'
        for lib in sorted(lib_notin):
            print lib
          
    def filter_v2(self,pjstat_file):
        pj_filter_dic = {}
        pj_filter_list = []
        raw_name_list = pjstat_file.split('_')
        f_in = open(pjstat_file.strip())
        f_out = open('{}_filter_{}'.format('_'.join(raw_name_list[:-1]),raw_name_list[-1]),'w')
        for line in f_in:
            if line[0] != '#':
                #print line
                line_list = line.strip().split('\t')
                primary = line_list[2]
                pj_time = line_list[1]
                new_line = '\t'.join(line_list[1:])+'\n'
                if pj_filter_dic.has_key(primary):
                    if pj_filter_dic[primary].split('\t')[0] > pj_time:
                        pj_filter_dic[primary] = new_line
                    else:
                        continue
                else:
                    pj_filter_dic[primary] = new_line

        f_out.write('#Time\tQC stat\tProject No.\tSample\tNovo ID\tLibrary ID\tRun\tLane\tRaw reads\tError Rate\tQ20\tQ30\tMapped\tDup(%)\tAverage depth\n')
        for value in reversed(sorted(pj_filter_dic.values())):
            f_out.write(value)
        f_in.close()
        f_out.close()

    def filter_v3(self,wbi_file,wobi_file):
        wbi_dic = {}
        wobi_file_list = wobi_file.split("_")
        wobi_file_new = '_'.join(wobi_file_list[:-1]+["novo",wobi_file_list[-1]])
        f_wbi = open(wbi_file)
        for line in f_wbi:
            if line[0] == "#":
                continue
            line_list = line.strip().split("\t")
            wbi_dic['_'.join([line_list[5],line_list[6],line_list[7]])] = line_list[4]
        f_wbi.close()

        f_wobi = open(wobi_file)
        f_wobi_new = open(wobi_file_new,'w')
        f_wobi_new.write(f_wobi.readline())
        for line in f_wobi:
            line_list = line.split('\t')
            key = '_'.join([line_list[5],line_list[6],line_list[7]])
            if wbi_dic.has_key(key):
                #print key
                f_wobi_new.write('\t'.join(line_list[:4]+[wbi_dic[key]]+line_list[5:]))
                #print '\t'.join(line_list[:4]+[wbi_dic[key]]+line_list[5:])
            else:
                f_wobi_new.write(line)
        f_wobi.close()
        f_wobi_new.close()



if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="qc mapping stat")
    subpar = argparser.add_subparsers(help='son command')
    sub_stat = subpar.add_parser('stat', help='stat command')
    sub_stat.set_defaults(cmd='stat')
    sub_stat.add_argument('-d','--pj_dir',help='project path',default='/ifs/TJPROJ3/HWUS/project')
    sub_stat.add_argument('-p','--pj_table',help="project table")
    sub_filter = subpar.add_parser('filter', help='filter command')
    sub_filter.set_defaults(cmd='filter')
    sub_filter.add_argument('-s','--stat_table',help='the file need to be filter')
    sub_filter = subpar.add_parser('getid', help='wobi file get novo id command')
    sub_filter.set_defaults(cmd='getnovoid')
    sub_filter.add_argument('-b','--wbi',help='wbi file')
    sub_filter.add_argument('-o','--wobi',help='wobi file')
    myargv = argparser.parse_args()

    if myargv.cmd == 'stat':
        if myargv.pj_table :
            new_pathcheck = Path_Search(myargv.pj_table)
        else:
            new_pathcheck = Path_Search()
        new_pathcheck.pathsearch(myargv.pj_dir)
    elif myargv.cmd == 'filter':
        new_pathcheck = Path_Search()
        new_pathcheck.filter_v2(myargv.stat_table)
    elif myargv.cmd == 'getnovoid':
        new_pathcheck = Path_Search()
        new_pathcheck.filter_v3(myargv.wbi,myargv.wobi)
