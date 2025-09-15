# Nextflow-GWAS-PopGen
nextflow run -profile gb main.nf --vcfgzFile data/cucumber_core_SNP.vcf.gz --traitsFile data/cucumber_core_pheno_mock.txt \
    --method MLM \
    --analyze GWAS \
    --QCTools BCFTools \
    --output resultVCFcucumber2\
