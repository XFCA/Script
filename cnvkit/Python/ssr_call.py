import re
class Ssr_call():
    def __init__(self,file):
        self.ssr_file = file
        self.ssrs_dic = {}
        self.seq_dic = {}
        
    def call(self):
        self.get_seq()
        for key in self.seq_dic.keys():
            self.ssrs_dic[key] = self.get_ssr(self.seq_dic[key])
        #self.ssrs_dic['Cluster-1484.26669'] = self.get_ssr(self.seq_dic['Cluster-1484.26669'])
        self.myprint()
    def myprint(self):
        f_out = open('ssrs.xls','w')
        f_out.write('ID\tstart\tend\tseq\tunit\tunit len\tunint number\n')
        for key in  self.ssrs_dic.keys():
            value = self.ssrs_dic[key]
            for ssr in value:
                end = ssr[2]+len(ssr[0])-1
                unit_len = len(ssr[1])
                unit_no = len(ssr[0])/unit_len
                f_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,ssr[2],end,ssr[0],ssr[1],unit_len,unit_no))
        f_out.close()
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
        while x < len(seq):
            y = 1
            while y<=6 :
                patten = r"(%s){2,}" %(seq[x:(x+y)])
                #print patten
                com = re.match(patten,seq[x:])
                if (com) and (len(com.group())>=6):
                    ssrs_list.append([com.group(),com.group(1),x+1])
                    x += com.end()
                    #print y
                    y = 1
                    #print ssrs_list[-1]
                    break
                elif y == 6:
                    y += 1
                    x += 1
                else:
                    y +=1
        
        return ssrs_list
    
if __name__ == "__main__":
    filename = "/ifs/TJPROJ3/HWUS/USER/zhangjie/test/KAOHE/19Q2/transcript.fasta"
    new_ssr = Ssr_call(filename)
    new_ssr.call()
