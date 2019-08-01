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
        
