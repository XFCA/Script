#!/usr/bin/python
import os,argparse
import time
import pysam
import re
import gzip
def safe_open(file_name,mode):
    try:
        if not file_name.endswith('.gz'):
            return open(file_name,mode)
        else:
            import gzip
            return gzip.open(file_name,mode)
    except IOError:
        print file_name + ' do not exist!'

def run_lane(path,lib,lane):
    f_in = gzip.open('{0}/{1}/{1}_{2}_1.fq.gz'.format(path,lib,lane))
    run_name = ':'.join(f_in.readline()[1:].split(':')[:3])
    f_in.close()
    return run_name

parser = argparse.ArgumentParser(description="Remove part of duplication,and get new rawdata.\nContact:yuhuan@novogene.com,hanqingmei@novogene.com")
parser.add_argument('-li','--samli',help="Sample list you want to remove duplication.Required.",required=True)
parser.add_argument('-p','--projdir',help="The project dir.Required.",required=True)
parser.add_argument('-o','--outpath',help="Out put path after duplication removed,default = ./",required=False)
parser.add_argument('-f','--freq',help="Frequency of dup you want remove.Float,default=0.5.",default='0.5')
parser.add_argument('-pv','--pepversion',help="Pepline version of your project after v4.2 or not.Required.",choices=['Y','N'],default='Y')
argv = parser.parse_args()
samli = os.path.abspath(argv.samli)
proj = os.path.abspath(argv.projdir)
pv = argv.pepversion
freq = float(argv.freq)

#get output path
if argv.outpath:
    outp = os.path.abspath(argv.outpath)
else:
    outp = os.getcwd()

Dpath = outp+'/RemoveDup'
if not os.path.exists(Dpath):
    os.makedirs(Dpath) 
       
if pv == 'Y':
    lpos,spos,bpos,ppos = 0,2,3,6
else:
    lpos,spos,bpos,ppos = 0,2,3,5

BamLi = []
LaneLi = [] 
for line in safe_open(samli,'r'):
    if line.startswith('#'):
        title = line.strip().split('\t')
        if 'Ori_Lane' in title:
            lpos = title.index('Ori_Lane')
        elif 'lane' in title:
            lpos = title.index('lane')
        if 'SampleID' in title:
            spos = title.index('SampleID')
        if 'Path' in title:
            ppos = title.index('Path')
        if 'LibID' in title:
            bpos = title.index('LibID')
        continue
    uline = line.split('\t')
    #get sampleID,path,lane,libraryID 
    samID = uline[spos].strip()
    path = uline[ppos].strip()
    lane = 'L'+uline[lpos][-1]
    lib = uline[bpos].strip()
    fq1 = path+'/'+lib+'/'+lib+'_'+lane+'_'+'1.fq.gz'
    fq2 = path+'/'+lib+'/'+lib+'_'+lane+'_'+'2.fq.gz'
    ada1 = path+'/'+lib+'/'+lib+'_'+lane+'_'+'1.adapter.list.gz'
    ada2 = path+'/'+lib+'/'+lib+'_'+lane+'_'+'2.adapter.list.gz'
    ofq1 = Dpath+'/'+lib+'/'+lib+'_'+lane+'_'+'1.fq'
    ofq2 = Dpath+'/'+lib+'/'+lib+'_'+lane+'_'+'2.fq'
    if os.path.exists(os.path.join(proj,'Mapping','%s.%s' %(samID,samID))):
        bam = os.path.join(proj,'Mapping','%s.%s' %(samID,samID),'%s.final.bam' %samID)
    elif os.path.exists(os.path.join(proj,'Mapping','%s' %samID)):
        bam = os.path.join(proj,'Mapping','%s' %samID,'%s.final.bam' %samID)
    else:
        print 'Bam file not exists: %s or %s neither exists.' %(os.path.join(proj,'Mapping','%s.%s' %(samID,samID),'%s.final.bam' %samID),os.path.join(proj,'Mapping','%s' %samID,'%s.final.bam' %samID))
    print bam 
    SDpath = Dpath+'/'+lib
    if not os.path.exists(SDpath):
        os.mkdir(SDpath)
    if bam not in BamLi:
        dupf=SDpath+'/%s.DupListOri.xls' %samID
        #print "Out ori list: %s" %dupf
        dupIdLi = []
        BamLi.append(bam)
        bamf=pysam.Samfile(bam,'rb')
        dupfo = safe_open(dupf,'w')
        for read in bamf.fetch():
            #skip the reads is not primary.Get dup only from primary reads.
            if read.flag & 0x900 != 0:continue
            if read.is_duplicate:
                dup_id=read.qname
                dupfo.write(dup_id+'\n')
        bamf.close()
        dupfo.close()
    run_id = run_lane(path,lib,lane)
    udupli = SDpath+'/'+'%s.DupLi_%s.xls' %(samID,lane)
    os.system('sed -n \'/%s:%s:/\'p %s > %s' %(run_id,uline[lpos][-1],dupf,udupli))
    tmpdup = open(udupli,'r')
    tmpdupli = []
    for tn in tmpdup:
        tmpdupli.append(tn.strip())
    FidupLi = sorted(list(set(tmpdupli)))
    lenthdup = int(len(FidupLi))
    print "Length:",lenthdup
    print "Freq:",freq
    uselen = int(lenthdup*freq)
    print "Use length:",uselen
    fdupLi = FidupLi[0:uselen]
    flist = SDpath+'/'+'%s.FiDupLi_%s.xls' %(samID,lane)
    wflist = open(flist,'w')
    wflist.write('\n'.join(fdupLi))
    Fshell = open(SDpath+'/'+samID+'ReDup_'+lane+'_1.sh','w')
    Sshell = open(SDpath+'/'+samID+'ReDup_'+lane+'_2.sh','w')
    Fshell.write('echo "Begin:" `date` && \\\npython /PUBLIC/source/HW/CANCER/HW_v1.0/script/rm_dup.py -id %s -i %s -o %s' %(flist,fq1,ofq1))
    Fshell.write(' && \\\ngzip -f %s && \\\necho "End:" `date`' %ofq1)
    Fshell.write(' && \\\nln -s %s %s' %(ada1,SDpath))
    Sshell.write('echo "Begin:" `date` && \\\npython /PUBLIC/source/HW/CANCER/HW_v1.0/script/rm_dup.py -id %s -i %s -o %s' %(flist,fq2,ofq2))
    Sshell.write(' && \\\ngzip -f %s && \\\necho "End:" `date`' %ofq2)
    Sshell.write(' && \\\nln -s %s %s' %(ada2,SDpath))
    os.system('rm -f %s' %udupli)
    dupfbak = SDpath+'/%s.DupListOri_bak.xls' %samID
    if not os.path.exists(dupfbak):
        os.system('cp %s %s' %(dupf,dupfbak))
    #qusbf = open(SDpath+'/'+'qsub_'+samID+'.sh','a')
    #qusbf.write('qsub -S /bin/bash -q all.q -q hw.q -cwd -V -l vf=3G %s\nqsub -S /bin/bash -q disease.q -q all.q -q cancer.q -cwd -V -l vf=3G %s' %(SDpath+'/'+samID+'ReDup_'+lane+'_1.sh',SDpath+'/'+samID+'ReDup_'+lane+'_2.sh\n'))
    Fshell.close()
    Sshell.close()
