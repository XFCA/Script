cnvkit.py batch \
.....final.bam \
-n ......final.bam \
--targets /PUBLIC/source/HW/Disease/database/Cap/Exome_bed/Agilent/SureSelect.Human.All.Exon.V6.r2/agilent_region.bed \
--annotate /HWPROJ2/HW/wangjinhao/data/cnvkit/refFlat.b37.txt \
--fasta /HWPROJ2/HW/wangjinhao/data/cnvkit/human_g1k_v37_decoy.fasta \
--access /HWPROJ2/HW/wangjinhao/data/cnvkit/access.b37.bed \
--output-reference my_reference.cnn --output-dir ....../results/
--diagram --scatter &&\
