nextflow.preview.dsl=2
/*
================================================================================
=                           Sinonkt Style I N I T                              =
================================================================================
*/

include './modules/nbt/utils'

if (params.exportKeySchema) exit 0, printKeySchema()
if (params.exportValueSchema) exit 0, printValueSchema()

params.MAINTAINERS = [
  'Krittin Phornsiricharoenphant (oatkrittin@gmail.com)',
  'Alisa Wilantho (alisa.wil@biotec.or.th)',
  'Sujiraporn Pakchuen (sujiraporn.pak@biotec.or.th)'
]

def schema = readAvroSchema("${workflow.projectDir}/schemas/value.avsc")
__params = getDefaultThenResolveParams(schema, params)

include './modules/nbt/log' params(__params)
include helpMessage from './modules/nbt/help' params(__params)
include './modules/convert' params(__params)
include './modules/quality' params(__params)
include './modules/preprocess' params(__params)
include VCFstats as VCFstats_GWAS from './modules/VCFstats' params(__params)
include VCFstats as VCFstats_PopGen from './modules/VCFstats' params(__params)
include './modules/population_analysis' params(__params)
include './modules/gwas_analysis' params(__params)
include Admixture as AdmixtureIter1 from './modules/population_analysis' params(__params)
include Admixture as AdmixtureIter2 from './modules/population_analysis' params(__params)
include Admixture as AdmixtureIter3 from './modules/population_analysis' params(__params)


if (__params.vcfgzFile && file(__params.vcfgzFile).exists()) {
    params.inputs = "VCF"
} else if (__params.hmpFile && file(__params.hmpFile).exists()) {
    params.inputs = "HapMap"
} else if (__params.bedbimfamFile && file(__params.bedbimfamFile).exists()) {
    params.inputs = "BedBimFam"
} else {
    error "errors"
}


__params.PopGenTools = (__params.PopGenTools instanceof List) ? __params.PopGenTools : __params.PopGenTools.tokenize(',')



if (params.version) exit 0, workflowVersionMessage()
if (params.help) exit 0, helpMessage(schema)

/*
================================================================================
=                   Sinonkt Style Workflows definitions                        =
================================================================================
*/


def pickCVerror = { 
  init = it.init()
  cv_error = it.last().text.split("\\,").last()
  return [*init, cv_error]
}

workflow GWAS {
  get:
    cleanedHapMap
    traits 
    modes
  main:
    merged = traits.merge(cleanedHapMap ) { left, right -> [ right.first(), *left.tail(), *right.tail() ]}
    gwasResults = GWASAnalysis(merged, modes)
  emit:
    gwasResults
}


workflow Phylogenetic_step {
  get:
    vcfgzCleanedPruned
  main:
    PhyloResults = Phylogenetic(vcfgzCleanedPruned)
  emit:
    PhyloResults 
}

workflow Admixture_Only {
  get:
    cleanedPrunedBedBimFams
  main:
    def ks = [(1..11), *(12..20).collate(3)]
    iter0 = Admixture(cleanedPrunedBedBimFams, ks[0])
    ksIter1 = iter0.map(pickCVerror).min { it.last() }.map { it.first() == ks[0].last() ? ks[1]: [] }.flatten()
    iter1 = AdmixtureIter1(cleanedPrunedBedBimFams, ksIter1)
    ksIter2 = iter0.mix(iter1).map(pickCVerror).min { it.last() }.map { it.first() == ks[1].last() ? ks[2]: [] }.flatten()
    iter2 = AdmixtureIter2(cleanedPrunedBedBimFams, ksIter2)
    ksIter3 = iter0.mix(iter1, iter2).map(pickCVerror).min { it.last() }.map { it.first() == ks[2].last() ? ks[3]: [] }.flatten()
    iter3 = AdmixtureIter3(cleanedPrunedBedBimFams, ksIter3)
    allIters = iter0.mix(iter1, iter2, iter3)
    minError = allIters.map(pickCVerror).min { it.last() }
    SummarizeAdmixture(minError, allIters.map { it.last() }.toList())
    best_admixture = GatherAdmixtureResults(groupTupleWithOutKey(allIters.map { [it.init().last(), it.last() ] }))
  emit:
    best_admixture
}


workflow Faststructure {
  get:
    cleanedPrunedBedBimFams
  main:
    def ks = [(1..11), *(12..20).collate(3)]
    iter0 = Admixture(cleanedPrunedBedBimFams, ks[0])
    ksIter1 = iter0.map(pickCVerror).min { it.last() }.map { it.first() == ks[0].last() ? ks[1]: [] }.flatten()
    iter1 = AdmixtureIter1(cleanedPrunedBedBimFams, ksIter1)
    ksIter2 = iter0.mix(iter1).map(pickCVerror).min { it.last() }.map { it.first() == ks[1].last() ? ks[2]: [] }.flatten()
    iter2 = AdmixtureIter2(cleanedPrunedBedBimFams, ksIter2)
    ksIter3 = iter0.mix(iter1, iter2).map(pickCVerror).min { it.last() }.map { it.first() == ks[2].last() ? ks[3]: [] }.flatten()
    iter3 = AdmixtureIter3(cleanedPrunedBedBimFams, ksIter3)
    allIters = iter0.mix(iter1, iter2, iter3)
    minError = allIters.map(pickCVerror).min { it.last() }
    SummarizeAdmixture(minError, allIters.map { it.last() }.toList())
    best_admixture = GatherAdmixtureResults(groupTupleWithOutKey(allIters.map { [it.init().last(), it.last() ] }))
    allKs = Channel.from(ks[0]).mix(ksIter1, ksIter2, ksIter3)
    FastStructureResults = StructureByFastStructure(cleanedPrunedBedBimFams, allKs)
  emit:
    best_admixture
    FastStructureResults
}

workflow IPCAPS {
  get:
    cleanedPrunedBedBimFams
    traits
  main:
    merged = cleanedPrunedBedBimFams.merge(traits) { left, right -> [ *left, *right.tail() ] }
    ipcaps = IPCAPs(merged)
  emit:
    ipcaps 
}


workflow {
  println("====================")
  println(__params)
  println("====================")
   
  if (params.inputs == "HapMap" ) {
    hmps = Channel.fromPath(__params.hmpFile).map{ [it.name.split("\\.").first(), it] }
    hmps.view()
    traits = Channel.fromPath(__params.traitsFile)
      .map{
      traitsHeader = file(it).text.readLines().first().split("\\t")
      [ it.name.split("\\.").first(), traitsHeader.takeRight(traitsHeader.size()-2), it ] }
    traits.view()
	
    if (__params.QCTools == "PLINK" && __params.analyze == "GWAS"){
      MapPed = ConvertHapMapToMapPed(hmps)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      cleanedHapMap  = ConvertMapPedToHapMap(cleanedMapPed)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      vcfgzCleaned = ConvertMapPedToVCF(cleanedMapPed)
      VCFstats_GWAS(vcfgzCleaned)
    }
    else if (__params.QCTools == "BCFTools"  && __params.analyze == "GWAS"){
      vcfgz = ConvertHapMapToVCF(hmps)
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      cleanedHapMap = ConvertVCFToHapMap(vcfgzCleaned)
      gwasResults = GWAS(cleanedHapMap , traits, Channel.from(__params.method))
      VCFstats_GWAS(vcfgzCleaned)
    }
    else if (__params.QCTools == "PLINK" && __params.analyze == "PopGen"){
      MapPed = ConvertHapMapToMapPed(hmps)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      vcfgzCleaned = ConvertMapPedToVCF(cleanedMapPed)
      if (__params.PruneLDTools == "PLINK"){
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}		
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
    else if (__params.QCTools == "BCFTools" && __params.analyze == "PopGen"){
      vcfgz = ConvertHapMapToVCF(hmps)
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      if (__params.PruneLDTools == "PLINK"){
        cleanedMapPed = ConvertVCFToMapPed(vcfgzCleaned)
	cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
	vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
	vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
    else if (__params.QCTools == "PLINK"  && __params.analyze == "BOTH"){
      MapPed = ConvertHapMapToMapPed(hmps)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      cleanedHapMap = ConvertMapPedToHapMap(cleanedMapPed)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      vcfgzCleaned = ConvertMapPedToVCF(cleanedMapPed)
      VCFstats_GWAS(vcfgzCleaned)
      if (__params.PruneLDTools == "PLINK"){
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
	vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }	
    else if (__params.QCTools == "BCFTools" && __params.analyze == "BOTH"){
      vcfgz = ConvertHapMapToVCF(hmps)
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      cleanedHapMap = ConvertVCFToHapMap(vcfgzCleaned)
      gwasResults = GWAS(cleanedHapMap , traits, Channel.from(__params.method))
      VCFstats_GWAS(vcfgzCleaned)
      if (__params.PruneLDTools == "PLINK"){
        cleanedMapPed = ConvertVCFToMapPed(vcfgzCleaned)
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
	vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
  }
  else if (params.inputs == "VCF" ) {
    vcfgz = Channel.fromPath(__params.vcfgzFile).map{ [it.name.split("\\.").first(), it] }
    vcfgz.view()
    traits = Channel.fromPath(__params.traitsFile)
      .map{
      traitsHeader = file(it).text.readLines().first().split("\\t")
      [ it.name.split("\\.").first(), traitsHeader.takeRight(traitsHeader.size()-2), it ] }
    traits.view()
    if (__params.QCTools == "PLINK" && __params.analyze == "GWAS"){
      MapPed = ConvertVCFToMapPed(vcfgz)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      cleanedHapMap = ConvertMapPedToHapMap(cleanedMapPed)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      vcfgzCleaned = ConvertMapPedToVCF(cleanedMapPed)
      VCFstats_GWAS(vcfgzCleaned)
    }
    else if (__params.QCTools == "BCFTools"  && __params.analyze == "GWAS"){
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      cleanedHapMap = ConvertVCFToHapMap(vcfgzCleaned)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      VCFstats_GWAS(vcfgzCleaned)
    }
    else if (__params.QCTools == "PLINK" && __params.analyze == "PopGen"){
      MapPed = ConvertVCFToMapPed(vcfgz)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      vcfgzCleaned = ConvertMapPedToVCF(cleanedBedBimFams)
      if (__params.PruneLDTools == "PLINK"){
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
    else if (__params.QCTools == "BCFTools" && __params.analyze == "PopGen"){
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      if (__params.PruneLDTools == "PLINK"){
        cleanedMapPed = ConvertVCFToMapPed(vcfgzCleaned)
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
         ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
    else if (__params.QCTools == "PLINK" && __params.analyze == "BOTH"){
      MapPed = ConvertVCFToMapPed(vcfgz)
      cleanedMapPed = QualityControlByPLINK(MapPed)
      cleanedHapMap = ConvertMapPedToHapMap(cleanedMapPed)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      vcfgzCleaned = ConvertHapMapToVCF(cleanedHapMap)
      VCFstats_GWAS(vcfgzCleaned)
      if (__params.PruneLDTools == "PLINK"){
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFtoBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
    else if (__params.QCTools == "BCFTools" && __params.analyze == "BOTH"){
      vcfgzCleaned = QualityControlByBCFTools(vcfgz)
      cleanedHapMap = ConvertVCFToHapMap(vcfgzCleaned)
      gwasResults = GWAS(cleanedHapMap, traits, Channel.from(__params.method))
      VCFstats_GWAS(vcfgzCleaned)
      if (__params.PruneLDTools == "PLINK"){
        cleanedMapPed = ConvertVCFToMapPed(vcfgzCleaned)
        cleanedPrunedBedBimFams = PruneLDByPLINK(cleanedMapPed)
        vcfgzCleanedPruned = ConvertBedBimFamToVCF(cleanedPrunedBedBimFams)
      }
      else if (__params.PruneLDTools == "BCFTools"){
        vcfgzCleanedPruned = PruneLDByBCFTools(vcfgzCleaned)
        cleanedPrunedBedBimFams = ConvertVCFToBedBimFam(vcfgzCleanedPruned)
      }
      VCFstats_PopGen(vcfgzCleanedPruned)
      if (__params.PopGenTools.contains("phylo")){
        PhyloRes = Phylogenetic_step(vcfgzCleanedPruned)}
      if (__params.PopGenTools.contains("admixture")){
        admixtureRes = Admixture_Only(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("faststructure")){
        admixture_and_fastsRes = Faststructure(cleanedPrunedBedBimFams)}
      if (__params.PopGenTools.contains("ipcaps")){
        ipcapsRes = IPCAPS(cleanedPrunedBedBimFams,traits)}
    }
  }
}

workflow.onComplete { handleCompleteMessage() }
workflow.onError { handleErrorMessage() }
