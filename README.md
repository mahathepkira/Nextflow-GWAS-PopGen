# Nextflow-GWAS-PopGen

## For vcf file
### GWAS
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --analyze GWAS \
     --output result
```
### PopGen
```bash
nextflow run -profile gb main.nf --vcfgzFile data/inputonebase.vcf.gz --traitsFile /nbt_main/share/pachyderm/gwasrice/datauserupload/sampleinonebase/inputtrat2.txt \
     --QCTools BCFTools \
     --method MLM \
     --PruneLDTools PLINK \
     --analyze PopGen \
     --PopGenTools phylo,faststructure,ipcaps \
     --output result
```
