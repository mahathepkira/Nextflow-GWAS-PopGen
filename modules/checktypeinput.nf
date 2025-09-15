
process CheckTypeHapmap {

  //publishDir "${outputPrefixPath(params, task)}"
  //publishDir "${s3OutputPrefixPath(params, task)}"

  tag { key }

  input:
  tuple key, file(hmp)
  //tuple key, file(traits)

  output:
  //file "up${key}.traits.txt"
  tuple key, file("up${key}.hmp.txt")

  script:
  """
  length=\$(sed -n '3p' $hmp | cut -f15 | tr -d "\n" | wc -m)
  if [[ \$((length)) -ge 2 ]];
  then
    run_pipeline.pl -Xmx${task.memory.toGiga()}G -h \
    $hmp -sortPositions \
    -export up${hmp} \
    -exportType Hapmap
  else
    mv ${hmp} "up${key}.hmp.txt"
    
  fi
  """
}
//mv "${key}.txt" "up${key}.traits.txt"