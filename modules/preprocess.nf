process QualityControlByPLINK {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${prefix}.clean.map"), file("${prefix}.clean.ped")

  script:
  prefix=map.baseName
  """
  plink --file $prefix \
    --no-fid \
    --no-parents \
    --no-sex \
    --no-pheno \
    --geno ${params.geno} \
    --recode \
    --allow-extra-chr \
    --out ${prefix}.clean \
    --maf ${params.maf} \
  """
}

process QualityControlByBCFTools {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${prefix}.clean.vcf.gz")

  script:
  prefix=vcfgz.baseName
  """
  bcftools view -i 'MAF >= ${params.geno} && F_MISSING <= ${params.maf}' \
    -Oz -o ${prefix}.clean.vcf.gz ${vcfgz}
  """
}

process PruneLDByPLINK {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${prefix}.final.bed"), file("${prefix}.final.bim"), file("${prefix}.final.fam")

  script:
  prefix=map.baseName
  """
  plink --file ${prefix} --indep-pairwise 100 5 0.2 --allow-extra-chr --out ${prefix}.pruneLD
  plink --file ${prefix} --extract ${prefix}.pruneLD.prune.in --make-bed --allow-extra-chr --out ${prefix}.final
  """
}

process PruneLDByBCFTools {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${prefix}.final.vcf.gz")

  script:
  prefix=vcfgz.baseName
  """
  bcftools +prune ${vcfgz} -m 0.2 -w 100 -N maxAF -Oz -o ${prefix}.final.vcf.gz
  """
}

