
include './nbt/utils'

process WEKAConvertFromCSV {

   tag { "${prefix}"}

   input:
   file inputCSV

   output:
   file "${prefix}.arff"

   script:
   prefix=inputCSV.baseName
   """
   java -Xmx${task.memory.toGiga()}G -cp \$EBROOTWEKA/weka.jar \\
      weka.core.converters.CSVLoader \\
      ${inputCSV} > ${prefix}.arff
   """
}

process WekaFilterNumericAttributesToNominal {

   tag { "${prefix}"}

   input:
   file arff

   output:
   file "${prefix}.nominal.arff"

   script:
   prefix=arff.baseName
   """
   java -Xmx${task.memory.toGiga()}G -cp \$EBROOTWEKA/weka.jar \\
      weka.filters.unsupervised.attribute.NumericToNominal \\
      -i ${arff}  \\
      -o ${prefix}.nominal.arff \\
      -R first-last
   """
}

process WekaInfoGainAttributeEval {

   publishDir "${outputPrefixPath(params, task)}"
   publishDir "${s3OutputPrefixPath(params, task)}"

   tag { "${prefix}"}

   input:
   file nominalArff

   output:
   tuple file(nominalArff), file("${prefix}.orderSNP")

   script:
   prefix=nominalArff.baseName
   """
   java -Xmx${task.memory.toGiga()}G -cp \$EBROOTWEKA/weka.jar \\
      weka.attributeSelection.InfoGainAttributeEval \\
      -s \"weka.attributeSelection.Ranker -T -1.7976931348623157E308 -N -1\" \\
      -i ${nominalArff} | grep 'Selected attributes:' > ${prefix}.orderSNP
   """
}

process ParseFilterAndBayesClassifyByTopN {

   tag { "${prefix}:${n}"}

   input:
   tuple file(nominalArff), file(orderedSNP)
   each n

   output:
   tuple file("result.${n}.txt"), file("${prefix}.filtered.${n}.arff")

   script:
   mem=task.memory.toGiga()
   prefix=nominalArff.baseName
   """
   ParseFilter_TopN_Infogain.py $orderedSNP $n > train_filter_attribute_${n}.txt

   java -Xmx${mem}G -cp \$EBROOTWEKA/weka.jar \\
      weka.filters.unsupervised.attribute.Remove \\
      -V \\
      -R \$(cat train_filter_attribute_${n}.txt) \\
      -c last \\
      -i ${nominalArff} \\
      -o ${prefix}.filtered.${n}.arff
   
   java -Xmx${mem}G -cp \$EBROOTWEKA/weka.jar \\
      weka.classifiers.bayes.HNB \\
      -t ${prefix}.filtered.${n}.arff \\
      -c last \\
      -i > ${prefix}.result.HNB.${n}.txt

   CheckResults.py ${prefix}.result.HNB.${n}.txt \\
      | sed 's/${prefix}.result.HNB.//;s/.txt /\\t/g;s/ %//g' > result.${n}.txt
   """
}

process PickTopSNP {

   publishDir "${outputPrefixPath(params, task)}"
   publishDir "${s3OutputPrefixPath(params, task)}"

   input:
   tuple file(resultsByN), file(filteredArffs)

   output:
   tuple file("sorted.list.txt"), file("top.snp.list.txt"), file(filteredArffs)

   script:
   prefix=(filteredArffs.first() =~ /(.*).\d+.arff/)[0][1]
   """
   cat ${resultsByN} | sort -k1,1n > sorted.list.txt
   topSNPIdx=\$(cat sorted.list.txt | sort -k2,2rn | head -n1 | cut -f1)
   grep -w "@attribute" ${prefix}.\${topSNPIdx}.arff | cut -d' ' -f2 | sed '/class/d' > top.snp.list.txt
   """
}