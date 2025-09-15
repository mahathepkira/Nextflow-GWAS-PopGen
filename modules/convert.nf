include './nbt/utils'

process ConvertHapMapToMapPed {

  tag { key }

  input:
  tuple key, file(hmp)

  output:
  tuple key, file("${key}.mod.plk.map"), file("${key}.mod.plk.ped")

  script:
  """
  run_pipeline.pl -Xmx${task.memory.toGiga()}G -h \
    $hmp \
    -sortPositions \
    ${params.homozygous ? '-homozygous': ''} \
    -export ${key} \
    -exportType Plink
  perl -e ' \$sep=" "; while(<>) { s/\t/\$sep/g; print \$_; } warn "Changed tab to \$sep on \$. lines\n" ' ${key}.plk.ped \
    | cut -d' ' -f2,7- > ${key}.mod.plk.ped
  mv ${key}.plk.map ${key}.mod.plk.map
  """
}

process ConvertBedBimFamToMapPed {

  tag { key }

  input:
  tuple key, file(bed), file(bim), file(fam)

  output:
  tuple key, file("${key}.map"), file("${key}.ped")

  script:
  prefix=bed.baseName
  """
  plink --bfile ${prefix} --recode  --allow-extra-chr --out ${key}
  """
}

process ConvertMapPedToHapMap {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${key}.hmp.txt")

  script:
  """
  run_pipeline.pl -Xmx${task.memory.toGiga()}G \
    -plink \
    -map $map \
    -ped $ped \
    -sortPositions \
    -export ${key} \
    -exportType Hapmap
  """
}

process ConvertMapPedToVCF {

  tag { key }

  input:
  tuple key, file(map), file(ped)

  output:
  tuple key, file("${key}.vcf.gz")

  script:
  prefix=map.baseName
  """
  plink --file ${prefix} --recode vcf --allow-extra-chr --out ${key}
  bgzip ${key}.vcf
  """
}

process ConvertHapMapToVCF {

  tag { key }

  input:
  tuple key, file(hmp)

  output:
  tuple key, file("${key}.vcf.gz")

  script:
  """
  run_pipeline.pl -Xmx16G -h \
    $hmp \
    -sortPositions \
    ${params.homozygous ? '-homozygous' : ''} \
    -export ${key} \
    -exportType VCF
  bgzip ${key}.vcf
  """
}

process ConvertBedBimFamToVCF {

  tag { key }

  input:
  tuple key, file(bed), file(bim), file(fam)

  output:
  tuple key, file("${key}.vcf.gz")

  script:
  prefix = bed.baseName
  """
  plink --bfile ${prefix} --recode vcf --allow-extra-chr --out ${key}
  bgzip ${key}.vcf
  """
}

process ConvertVCFToHapMap {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${key}.hmp.txt")

  script:
  """
  run_pipeline.pl -Xmx${task.memory.toGiga()}G \
   -vcf ${vcfgz} \
   -export ${key} \
   -exportType Hapmap
  """
}

process ConvertVCFToBedBimFam {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${key}.bed"), file("${key}.bim"), file("${key}.fam")

  script:
  """
  plink --vcf ${vcfgz} --allow-extra-chr --make-bed --out ${key}
  """
}

process ConvertVCFToMapPed {

  tag { key }

  input:
  tuple key, file(vcfgz)

  output:
  tuple key, file("${key}.map"), file("${key}.ped")

  script:
  """
  plink --vcf ${vcfgz} --allow-extra-chr --make-bed --out ${key}
  plink --bfile ${key} --allow-extra-chr --recode --out ${key}
  """
}

