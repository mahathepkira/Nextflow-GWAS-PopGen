
include "./nbt/utils"

process GWASAnalysis {

	tag { "$key:$method" }

  publishDir "${outputPrefixPath(params, task)}"
  publishDir "${s3OutputPrefixPath(params, task)}"

	input:
	tuple key, traits, file(traitTxt), file(hmp)
	val method

	output:
	tuple key, file("*.csv"), file("*.pdf")

	script:
	cmds="""
  cat $traitTxt | cut -f1,3- > traits.trimmed.txt
	"""
	if (method == "FASTLMM")
		cmds+"""
		GAPIT BICmodelSelection traits.trimmed.txt $hmp
		for trait in ${traits.join(" ")}; do
			optimalPCsVal=\$(optimalPCs GAPIT.MLM.\${trait}.BIC.Model.Selection.Results.csv)
			if (\$((optimalPCsVal != "0")));
			then
				cut -d "," -f 1-\$((optimalPCsVal+1)) GAPIT.PCA.csv > GAPIT.PCA.\${optimalPCsVal}.\${trait}.csv
				GAPIT $method traits.trimmed.txt $hmp GAPIT.PCA.\${optimalPCsVal}.\${trait}.csv
			else
				GAPIT $method traits.trimmed.txt $hmp
			fi
		done
		"""
	else if (method != "FASTLMM")
		cmds+"""
		GAPIT $method traits.trimmed.txt $hmp
		"""
	//else if (method == "cMLM")
	//	cmds+"""
	//	GAPIT $method traits.trimmed.txt $hmp
	//	"""
	//else if (method == "GLM")
	//	cmds+"""
	//	GAPIT $method traits.trimmed.txt $hmp
	//	"""
	else 
		throw new Exception("Unknown method...")
}



// 1.	GWAS analysis
// ## Convert input file to ped/map format
// qsub runPlinkchangeformat.sh >plink --bfile cleanoutricein88sam_plkex1 --recode --out forgapitcleanoutricein88sam_plkex1 –noweb

// ## Convert ped/map to hapmap format
// qsub conplinkModi1test.sh >tassel-5-standalone/run_pipeline.pl -Xmx32g -plink -ped forgapitcleanoutricein88sam_plkex1.ped -map forgapitcleanoutricein88sam_plkex1.map -sortPositions -export gwasforgapit88samafterqc -exportType Hapmap


// ข้อมูลที่ใช้ใน script ด้านล่าง มีส่วนที่ได้จากการ upload เอามาใช้ คือ rice_traits.txt จากหน้าเว็บ
// ##CMLM in GAPIT –Basic Scenario
// >qsub runGaptest.sh
// ##FAST-LMM
// ## Run model selection using BIC (Bayesian information criterion)
// > qsub runGaptest1.sh
// 	Output = > GAPIT.MLM.HeaderofNamephenotype.BIC.Model.Selection.Results.csv
// ตัวอย่างผลที่ได้ ทำการเลือกจำนวน PC จากการดู BIC ที่มีค่ามากที่สุด ในที่นี้คือ 2 
// 
// แต่ในตัวอย่างนี้ ผลไม่ค่อยดี เลยใช้เป็น number of PCs/Covariates = 3
// > cut –d’,’ –f1-4 GAPIT.PCA.csv > GAPIT.PCA_data1NY2_2.csv
// ##Run FastLMM
// >qsub runGaptest2.sh
// 
