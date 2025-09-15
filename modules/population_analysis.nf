include './nbt/utils'

process IPCAPs {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  errorStrategy "ignore"

  tag { key }

  input:
  tuple key, file(bed), file(bim), file(fam), traits, file(traitTxt)

  output:
  tuple key,
    file("node*.txt"),
    file("cluster_output/images/*"),
    file("cluster_output/RData/*"),
    file("cluster_output/*.html"),
    file("cluster_output/groups.txt")

  script:
  """
  cat $traitTxt | awk 'NR > 1 { if (\$2 == "NA") print \$1" "\$1; else print \$1" "\$2; }' > sampleLabel.txt
  IPCAPS $bed sampleLabel.txt || true
  for rdataFile in \$(find . -name "node*.RData"); do
    extractIPCAPSnode $fam \$rdataFile
  done
  """
}

process Admixture {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { "$key:$k" }

  input:
  tuple key, file(bed), file(bim), file(fam)
  each k

  output:
  tuple k, key, file(bed), file(bim), file(fam), file("*.Q"), file("${key}.${k}.cv.errors.csv")

  script:
  """
  admixture $bed $k --cv=10 -j${task.cpus}
  perl -n -e'/^CV error \\(K=(\\d+)\\): ([\\d\\.]*)\$/ && print \$1,",",\$2' .command.log > ${key}.${k}.cv.errors.csv
  """
}

process SummarizeAdmixture {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { key }

  input:
  tuple optimalK, key, file(bed), file(bim), file(fam), file(minCVq), minCVerror
  file(minCVerrors)

  output:
  tuple optimalK, key, file("cv.errors.csv"), file("optimalK.csv")

  script:
  """
  echo "K,CrossValidationError" > cv.errors.csv
  for eachCVErrorByK in \$(find . -name "*.cv.errors.csv" | sort -k3,3g -t '.'); do
    cat \$eachCVErrorByK >> cv.errors.csv
    echo "" >> cv.errors.csv
  done
  echo "$optimalK,$minCVerror" >> optimalK.csv
  """
}


process GatherAdmixtureResults {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple file(cvQs), file(cvErrors)

  output:
  tuple file(cvQs), file(cvErrors)

  script:
  """
  echo "Just Gathering results to same s3 path..."
  """
}

process StructureByFastStructure {

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  errorStrategy "ignore"

  tag { "$key:$k" }

  input:
  tuple key, file(bed), file(bim), file(fam)
  each k

  output:
  tuple k, key, file("*.meanQ")

  script:
  prefix=bed.baseName

  """
  structure.py -K $k --input=$prefix --output=$key

  """
}

process Phylogenetic {

  errorStrategy 'ignore'

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

  tag { key }

  input:
  tuple key, file(vcf)

  output:
  tuple key, file("*.min4.phy.treefile")

  script:
  prefix = vcf.simpleName.replaceFirst(/\.vcf$/, "").replaceFirst(/\.gz$/, "")
  """
  vcf2phylip.py -i ${vcf} -o ${prefix}.min4.phy

  iqtree2 -s ${prefix}.min4.phy -B 1000 -T Auto

  """
}
