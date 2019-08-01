#!/HWPROJ2/HW/wangjinhao/software/anaconda2/envs/reseq/bin/python
import os,fnmatch,subprocess
from collections import defaultdict
import sys
class Path_Search():
    def __init__(self,rootdir,pjno):
        self.rootdir = rootdir
        self.pjdir = ['cancer','disease','seq','reseq']
        self.pjfound = set()
        self.pjnofound = set()
        self.pjno = self.get_pj(pjno)
        self.pjpath = self.get_path()

    def pathsearch(self):
        f_out = open('pjstat.xls','w')
        for mypj in sorted(self.pjpath.keys()):
            print self.pjpath[mypj][0]
            subprocess.Popen('cd {} && python /HWPROJ2/HW/wangjinhao/software/script/qc_map_stat.py {} {}'.format(self.pjpath[mypj][0],self.pjpath[mypj][1],self.pjpath[mypj][2]),shell=True,stdout=f_out)
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
        seq_list = os.popen('ls {}/seq'.format(self.rootdir)).read().strip().split('\n')
        self.pjdir += ['seq/{}'.format(name) for name in seq_list]
        self.pjdir += ['seq/{}/auto_seq'.format(name) for name in seq_list]
        reseq_list = os.popen('ls {}/reseq'.format(self.rootdir)).read().strip().split('\n')
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

            pj_list = os.popen('ls {}'.format(mypj_path)).read().strip().split('\n')
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
          

if __name__ == "__main__":
    new_pathcheck = Path_Search('/HWPROJ2/HW/project',sys.argv[1])
    new_pathcheck.pathsearch()
    new_pathcheck.filter('/HWPROJ2/HW/wangjinhao/PJSTAT/pjno.xls')
