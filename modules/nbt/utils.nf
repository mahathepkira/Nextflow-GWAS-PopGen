import groovy.json.JsonSlurper
/*
================================================================================
=                     Start Sinonkt Style Utils                                =
================================================================================
*/


def pathResolver(workflowHomePath, paramsPath) {
  return { path -> 
    if (path instanceof ArrayList) return path
    if (path instanceof Integer) return path
    if (path instanceof BigDecimal) return path
    if (path instanceof Float) return path

    path
      .replace("env://", workflowHomePath + "/")
      .replace("launchDir://", "${workflow.launchDir}/")
      .replace("params://", paramsPath ? "${paramsPath}/" : "${workflow.launchDir}/")
  }
}

def isPathType(field) { 
  return ["path", "dir"].contains(field.__type)
}

def isArrayType(field) { 
  return ["array"].contains(field.__type)
}

def paramResolver(workflowHomePath, paramsPath) {
  def resolvePath = pathResolver(workflowHomePath, paramsPath)
  return { field, value -> 
    if (value && isPathType(field)) {
      return resolvePath(value) 
    } else if (value && isArrayType(field)) {
      return value.split(",")
    } else {
      value
    }
  }
}

def getDefaultThenResolveParams(schema, params) {
  def workflowHomePath = schema.__WORKFLOW_HOME_PATH_BY_ENVIRONMENT[workflow.profile]
  def resolveParam = paramResolver(workflowHomePath, params.paramsPath)
  def isEnumType =  { fieldType -> (fieldType instanceof Map && fieldType.type == "enum") }
  def isListEnumType = { fieldType -> (fieldType instanceof List && fieldType.find(isEnumType)) }
  def hasEnumType = { field -> (isEnumType(field.type) || isListEnumType(field.type)) }
  def getSymbolsFromEnumType = { field ->
    if (isEnumType(field.type)) {
      return field.type.symbols
    } else if (isListEnumType(field.type)) {
      def enumType = field.type.find(isEnumType) 
      if (enumType) {
        return enumType.symbols
      }
    }
    throw new Exception("Unknown case of enumType to be resolve... __need_resolve param should have at least one enum type...")
  }

  return schema.fields.inject([:]) { result, field ->
    def value = params[field.name] ? params[field.name]: field.default
    if (field.__need_resolve && hasEnumType(field) && value) {
      def symbolsToResolve = getSymbolsFromEnumType(field)
      if (symbolsToResolve.contains(value)) {
        def keys = field.__need_resolve.keySet() as String[];
        def resolved = keys.inject([:]) { acc, key -> 
          def mapKey = field.__need_resolve[key]
          def resovledDefaultValue = field[mapKey][value]
          if (resovledDefaultValue instanceof Map) {
            def nestedKeys = resovledDefaultValue.keySet() as String[];
            def nestedResolved = nestedKeys.inject([:]) { nestedAcc, nestedKey ->
              def nestedResolvedDefaultValue = resovledDefaultValue[nestedKey]
              return [*:nestedAcc, (nestedKey): resolveParam(field, nestedResolvedDefaultValue)]
            }
            return [*:acc, (key): nestedResolved]
          }
          return [ *:acc, (key): resolveParam(field, resovledDefaultValue)]
        }
        return [ *:result, (field.name): resolveParam(field, value), *:resolved]
      }
    }
    return [ *:result, (field.name): resolveParam(field, value)]
  }
}

def loadJSON(path) { (new JsonSlurper()).parse(new File(path)) }
def readAvroSchema(path) { loadJSON(path) }

def groupTupleWithoutCommonKey(channel, spreadTail) {
  channel.map { ['tempKey', *it] }
    .groupTuple()
    .map { it.tail() }
    .map { spreadTail ? [it.first(), it.tail().flatten()] : it }
}

def processName(task) { task.process.split(":").last() }

def workflowName(task) { task.process.split(":").first() }

def outputPrefixPath(params, task) { "${params.output}/${workflowName(task)}/${processName(task)}" }

def s3OutputPrefixPath(params, task) { "${params.s3output}/${workflowName(task)}/${processName(task)}" }

def outputShortPrefixPath(params, task) { "${params.output}/${processName(task)}" }

def s3OutputShortPrefixPath(params, task) { "${params.s3output}/${processName(task)}" }

def groupTupleWithOutKey(channel) {
  channel.reduce([]) { acc, mesg -> 
    if (acc.size() == 0) {
      def numElements = mesg.size()
      acc = (1..numElements).collect{[]}
    }
    mesg = (mesg instanceof String[]) ? Arrays.asList(mesg) : mesg
    mesg.withIndex().collect { elem, idx -> [ *acc[idx], elem ] }
  }
}
def groupTupleWithOutKeyFromList(list) {
  list.inject([]) { acc, mesg -> 
    if (acc.size() == 0) {
      def numElements = mesg.size()
      acc = (1..numElements).collect{[]}
    }
    mesg = (mesg instanceof String[]) ? Arrays.asList(mesg) : mesg
    mesg.withIndex().collect { elem, idx -> [ *acc[idx], elem ] }
  }
}
// https://stackoverflow.com/questions/13155127/deep-copy-map-in-groovy
// standard deep copy implementation
def deepcopy(orig) { orig.getClass().newInstance(orig) }

def isDuplicatedList(l) { l.clone().unique().size() != l.size() }

def tails(ls) {
  if (ls.size == 0) return [[]]
  [ ls, *tails(ls.tail()) ]
}

// Implmentation from https://wiki.haskell.org/99_questions/Solutions/26
def combinations(ns, k) {
  if (ns.size == 0) return [[]]
  if (k == 0) return [[]]
  tails(ns).collect {
    if (it.size == 0) return []
    def (x, xs) = [it.first(), it.tail()]
    combinations(xs, k-1).collect { ys -> [x, *ys] }
  }.inject([]) { acc, combs -> [*acc, *combs] }
    .findAll { it.size >= k }
}

def printKeySchema() {
  log.info "================================================================================"
  log.info "===================== Begin Key Schema Content ================================="
  log.info "================================================================================"
  printFileContent("${workflow.projectDir}/schemas/key.avsc")
  log.info "================================================================================"
  log.info "===================== End Key Schema Content ==================================="
  log.info "================================================================================"
}

def printValueSchema() {
  log.info "================================================================================"
  log.info "===================== Begin Value Schema Content ==============================="
  log.info "================================================================================"
  printFileContent("${workflow.projectDir}/schemas/value.avsc")
  log.info "================================================================================"
  log.info "===================== End Value Schema Content ================================="
  log.info "================================================================================"
}

def printFileContent(filePath) {
  file(filePath)
    .readLines()
    .each { println it }
}

def checkIsExtensions(input_dir_path, extensions) {
  return file("${input_dir_path}/*").any { filePath -> extensions.any { filePath.name.endsWith(it) } }
}
def IsFASTQ(input_dir_path) { checkIsExtensions(input_dir_path, ["fastq", "fastq.gz", "fq", "fq.gz"]) }
def IsBAM(input_dir_path) { checkIsExtensions(input_dir_path, ["bam", "bam.bai", "bai"]) }

def chunkOf(num_chunks, channel) {
  channel.map { it }
}

def compareFirstElem(a, b) { b.first() <=> a.first() }
def sortByFirstElement(channel, ASC) {
  channel.toSortedList{ a, b -> compareFirstElem(a, b) }
    .flatMap { ASC ? it.reverse(): it }
}

def Arrayify(list) { list.collect { it instanceof Collection ? it: [it] } }

def Singlify(list) {
  out = []
  for (i = 0; i < list.first().size(); i++) {
    out[i] = list.collect { it[i] }
  }
  return out
}

def countUniq(list) {
  def _reduce = { acc, val ->
    if (!acc) return [[val, 1]]
    else if (acc && acc.last()[0] != val) return [*acc, [val, 1]]
    else {
      return [*acc.init(), [val, acc.last()[1]+1]]
    }
  }
  return list.inject([], _reduce)
}

/*
================================================================================
=                     End Sinonkt Style Utils                                  =
================================================================================
*/