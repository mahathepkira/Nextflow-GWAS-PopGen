#!/usr/bin/env python
import os
import sys
import subprocess


(__, hmpPath, traitsPath, trimmedHmpOutPath, trimmedTraitsOutPath) = sys.argv
CWD = os.getcwd()

def resolvePath(path):
  if os.path.isabs(path):
    return path
  return f"{CWD}/{path}"

hmpHeaders=None
with open(hmpPath, "r") as hmp:
  hmpHeaders=hmp.readline()

qcColumnIdx = 11 # According to HmpMap format, the last column before followed by sampleIds
with open(traitsPath, "r") as traits:
  sampleIdColumns = hmpHeaders.split("\t")[qcColumnIdx:]
  sampleIdsFromHmp = list(map(lambda x: x.strip(), sampleIdColumns))
  sampleIdsFromTraits = list(map(lambda line: line.split("\t")[0], traits.readlines()[1:]))
  concordanceSampleIds = list(set(sampleIdsFromHmp).intersection(set(sampleIdsFromTraits)))
  concordanceSampleIds.sort()

  excludeFromHmp = list(filter(lambda x: not x in concordanceSampleIds, sampleIdsFromHmp))
  excludeFromHmp.sort()

  excludeFromTraits = list(filter(lambda x: not x in concordanceSampleIds, sampleIdsFromTraits))
  excludeFromTraits.sort()

  [hmpPath, traitsPath, trimmedHmpOutPath, trimmedTraitsOutPath] = list(map(resolvePath, [hmpPath, traitsPath, trimmedHmpOutPath, trimmedTraitsOutPath]))

  HmpColumnIdxes = [ str(idx+12) for idx, sampleId in enumerate(sampleIdsFromHmp) if sampleId in concordanceSampleIds]

  TraitsStartWithSampleIdsPatterns = [f'^{id}' for id in concordanceSampleIds]
  subprocess.run(f"cut -f1-11,{','.join(HmpColumnIdxes)} {hmpPath} > {trimmedHmpOutPath}", shell=True, check=True)
  subprocess.run(f"head -n 1 {traitsPath} > {trimmedTraitsOutPath}", shell=True, check=True)
  subprocess.run(f"tail -n +2 {traitsPath} | grep -E '{'|'.join(TraitsStartWithSampleIdsPatterns)}' >> {trimmedTraitsOutPath}", shell=True, check=True)

  file1 = open(f'{CWD}/exclude.from.hmp.list.txt', 'w') 
  file1.write('\n'.join(excludeFromHmp))
  file1.close() 

  file2 = open(f'{CWD}/exclude.from.traits.list.txt', 'w') 
  file2.write('\n'.join(excludeFromTraits))
  file2.close() 

  file1 = open(f'exclude.from.hmp.list.txt', 'w') 
  file1.write('\n'.join(excludeFromHmp))
  file1.close() 