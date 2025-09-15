include './nbt/utils'

process VCFstats {

  tag { "${key}" }

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple key file(vcf)

  output:
  file("*.frq")
  file("*.lmiss")
  file("*.TsTv.summary")
  file("*.summary")
  script:
  prefix=vcf[0].baseName

  """
  bash /nbt_main/home/lattapol/nextflow-Callvariants/bin/quality.sh ${vcf}
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_AF_his.py ${key}.frq
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_lmiss_his.py ${key}.lmiss  
  """
}

