'''
Created on 2019年6月24日

@author: 18482
'''
import re
class Ssr_call():
    def __init__(self,file):
        self.ssr_file = file
        self.ssrs_dic = {}
        self.seq_dic = {}
        
    def call(self):
        self.get_seq()
        for key in self.seq_dic.keys():
            self.ssrs_dic[key] = self.get_ssr(key,self.seq_dic[key])
            
        
    
    def get_seq(self):
        f_in = open(self.ssr_file)
        for line in f_in:
            if line.startswith('>'):
                key = line.strip()[1:]
                self.seq_dic[key] = ''
            else:
                self.seq_dic[key] += line.strip()
        f_in.close()
        
    def get_ssr(self,seq):
        ssrs_list = []
        x = 0
        while(x < len(seq)):
            y = 6
            while y >0 :
                patten = "("+seq[x:(x+y)]+"){2,}"
                com = re.match(patten,seq[x:])
                if com:
                    ssrs_list.append([com.group(),com.group(1),com.end()+x])
                    x += com.end
                    print ssrs_list[-1]
                    break
                else:
                    y -= 1
        
        return ssrs_list
    
if __name__ == "__main__":
    filename = ""
    new_ssr = Ssr_call(filename)
    new_ssr.call()