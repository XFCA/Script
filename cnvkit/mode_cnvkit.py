import os,system,argparse,fnmatch
class Cnvkiter:
    def __init__(self):
        self.script = ""
    
    def batch_mk(self,argv):
        script_batch = "echo \"Cnvkit Somatic CNV batch Start: \" `date` && \\\n"
        script_batch += "cnvkit.py batch \\\n"
        script_batch += "\t{} \\\n".format(argv.t_bam)
        script_batch += "\t-n {} \\\n".format(argv.n_bam)
        script_batch += "\t-p {} \\\n".format(argv.process)
        script_batch += "\t--annotate {} \\\n".format(argv.refFlat)
        script_batch += "\t-f {} \\\n".format(argv.gene)
        script_batch += "\t-g {} \\\n".format(argv.access)
        if argv.n_sex == "M":
            script_batch += "\t-y \\\n"
        if argv.type == "WES":
            script_batch += "\t-t {} \\\n".format(argv.target)
            script_batch += "\t-m hybrid \\\n"
        elif argv.type == "WGS":
            script_batch += "\t-m wgs \\\n"
        elif argv.type == "TS":
            script_batch += "\t-t {} \\\n".format(argv.target)
            script_batch += "\t-m amplicon \\\n"
            
        script_batch += "\t-d {} \\\n".format(argv.outdir)
        script_batch +=  "--diagram --acatter && \\\n"
        script_batch += "echo \"Cnvkit Somatic CNV batch End: \" `date` && \\\n"
        
        return script_batch
        
    def call_mk(self,argv):
        script_call = "echo \"Cnvkit Somatic CNV call Start: \" `date` && \\\n"
        script_call += "cnvkit.py call \\\n"
        script_call += "\t{}/{}.cns \\\n".format(argv.outdir)
        script_call += "\t-v {} \\\n".format(argv.vcf)
        script_call += "\t-m {} \\\n".format(argv.methd)
        script_call += "\t-i {} \\\n".format(argv.tumor)
        script_call += "\t-n {} \\\n".format(argv.normal)
        if argv.t_sex == "M":
            script_call += "\t-x y \\\n"
        elif argv.t_sex == "F":
            script_call += "\t-x f \\\n"
        if argv.n_sex == "M":
            script_call += "\t-y \\\n"
        script_call += "\t-o {}/{}.call.cns && \\\n".format(argv.outdir,argv.tumor)
        script_call += "echo \"Cnvkit Somatic CNV call End: \" `date` && \\\n"
        
        return script_call
    
    def to_vcf_mk(self,argv):
        script_vcf = "echo \"Cnvkit Somatic CNV to vcf Start: \" `date` && \\\n"
        script_vcf += "cnvkit.py export vcf \\\n"
        script_vcf += "\t{}/{}.call.cns \\\n".format(argv.outdir,argv.tumor)
        script_vcf += "\t--cnr {}/{}.cnr \\\n".format(argv.outdir,argv.tumor)
        script_vcf += "\t-i {} \\\n".format(argv.tumor)
        if argv.t_sex == "M":
            script_vcf += "\t-x y \\\n"
        elif argv.t_sex =="F":
            script_vcf += "\t-x x \\\n"
        if argv.n_sex == "M":
            script_vcf += "\t-y \\\n"
        script_vcf += "\t-o {}/{}.cnvkit.vcf && \\\n".format(argv.outdir,argv.tumor)
        script_vcf += "gzip {}.cnvkit.vcf && \\\n".format(argv.tumor)
        script_call += "echo \"Cnvkit Somatic CNV to vcf End: \" `date` && \\\n"
        return script_vcf
    
    def main(self,myargv):
        os.mkdir(myargv.outdir)
        self.script += "cd {} \\\n".format(myargv.outdir)
        self.script += self.batch_mk(myargv)
        self.script += self.call_mk(myargv)
        self.script += self.to_vcf_mk(myargv)
        f_out = open("{}/{}.cnvkit.sh".format(myargv.outdir,myargv.tumor),"w")
        fout.write(self.script)
        f_out.close()
    
    def filter(self,myargv):
        pass
    
class Myparser():
    def __init__(self):
        myargv=self.argv_maker()
        self.start(myargv)

    def argv_maker(self):
        argparser = argparse.ArgumentParser(description="cnvkit mode")
        subpar = argparser.add_subparsers(help='command')
        sub_make = subpar.add_parser('make', help='make command')
        sub_make.set_defaults(cmd='make')
        sub_make.add_argument('--type',help='analysis type',defalut="WES")
        sub_make.add_argument('--gene',help="reference gene fasta")
        sub_make.add_argument('--target',help='a bed file for sequence region',defalut=False)
        sub_make.add_argument('--access',help='gene access bed')
        sub_make.add_argument('--refFlat',help='refFlat.txt')
        sub_make.add_argument('--tumor',help='tumor sample')
        sub_make.add_argument('--t_bam',help='tumor sample bam')
        sub_make.add_argument('--t_sex',help='tumor sample sex')
        sub_make.add_argument('--normal',help='normal sample')
        sub_make.add_argument('--n_bam',help='normal sample bam')
        sub_make.add_argument('--n_sex',help='normal sample sex')
        sub_make.add_argument('--vcf',help='mutect output')
        sub_make.add_argument('-d','--outdir',help='cnvkit analysis dir')
        sub_make.add_argument('-p','--process',help='process number',default='8')
        sub_filter = subpar.add_parser('filter', help='filter command')
        sub_filter.set_defaults(cmd='filter')
        myargv = argparser.parse_args()
        
    def sart(self,myargv):
        cnvkit_new = Cnvkiter()
        if myargv.cmd == "make":
            cnvkit_new.main(myargv)
        elif myargv.cmd == "filter":
            cnvkit_new.filter(myargv)
    
if __name__ == "__main__":
    argv_new = Myparser()
    argv_new.sart(argv_new.myargv)
    