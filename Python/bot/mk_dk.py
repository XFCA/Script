# -*- coding: utf-8 -*-
"""
Created on Tue Feb 05 14:51:53 2019

@author: jinhaow
"""

class value_base():
    def __init__(self):
        self.snp = {}#{'exon':num,...}
        self.indel = {}
        self.cnv = {}

class sample_value():
    def __init__(self):
        self.base = value_base()
        self.snp = {}#{'exon':[sample1,sample2,...],...}
        self.indel = {}
        self.cnv = {}
    
    
    def read_html(self):
        """
        读取html文件，提取表格内容
        """
        
        pass
    
    def mk_xls(self):
        """
        生成xls文件
        """
        self.read_html()
        f_out = open('for_dk.xls','w')
        for key1,value1 in zip(self.base.snp.keys(),self.base.snp.values()):
            f_out.write(key1 + value1 + '\n')
            for key2,value2 in zip(self.snp[key1].keys(),self.snp[key1].values()):
                f_out.write(key2 + value2 + '\n')
        for key1,value1 in zip(self.base.indel.keys(),self.base.indel.values()):
            f_out.write(key1 + value1 + '\n')
            for key2,value2 in zip(self.indel[key1].keys(),self.indel[key1].values()):
                f_out.write(key2 + value2 + '\n')
        for key1,value1 in zip(self.base.cnv.keys(),self.base.cnv.values()):
            f_out.write(key1 + value1 + '\n')
            for key2,value2 in zip(self.cnv[key1].keys(),self.cnv[key1].values()):
                f_out.write(key2 + value2 + '\n')
        f_out.close()
        