fastq = open("reads_for_analysis.fastq")

sequences = []
qualityValues = []

line = fastq.readline()
while line:
    sequences.append(fastq.readline())
    fastq.readline()
    qualityValues.append(fastq.readline())
    line = fastq.readline()

###############################################################################
### 4 a) Nustatyti, koks kokybės kodavimas yra naudojamas pateiktame faile. ###
###############################################################################

maxQualityValue = 0
minQualityValue = 126
# maxQualityValue = 73
# minQualityValue = 10

for qualityValue in qualityValues:
    for char in qualityValue:
        ascii = ord(char)
        if ascii < minQualityValue:
            minQualityValue = ascii
        else:
            if ascii > maxQualityValue:
                maxQualityValue = ascii

from enum import Enum
class Score(Enum):
    SangerPhred33 = "Sanger Phred+33"
    SolexaSolexa64 = "Solexa Solexa+64"
    Illumina13Phred64 = "Illumina 1.3+ Phred+64"
    Illumina15Phred64 = "Illumina 1.5+ Phred+64"
    Illumina18Phred33 = "Illumina 1.8+ Phred+33"

score: Score
if maxQualityValue >= 74 and minQualityValue >= 59 and minQualityValue <= 63:
    score = Score.SolexaSolexa64
elif maxQualityValue >= 74 and maxQualityValue  <= 104 and minQualityValue >= 64:
    score = Score.Illumina13Phred64
elif maxQualityValue >= 74 and maxQualityValue  <= 105 and minQualityValue >= 67:
    score = Score.Illumina15Phred64
elif maxQualityValue == 74 and minQualityValue <= 58:
    score = Score.Illumina18Phred33
elif maxQualityValue <= 74 and minQualityValue <= 58:
    score = Score.SangerPhred33

resultFile = open("result.txt", "a")
resultFile.write(score.value + "\n")

#####################################################################
### 4 b) Pateikti C/G nukleotidų pasiskirstymo read’uose grafiką. ###
#####################################################################

sequenceIds = []
baseCounts = []
id = 1
for sequence in sequences:
    sequenceIds.append(id)
    id+=1
    totalLength = 0
    baseCount = 0
    for letter in sequence:
        totalLength += 1
        if letter == "G" or letter == "C":
            baseCount += 1

    baseCounts.append(baseCount / totalLength)

import matplotlib.pyplot as plt
plt.plot(sequenceIds, baseCounts)
plt.show()


##################################################################################
### 4 c) Paimti po 5 kiekvieno piko viršūnės sekų ir atlikti blast’o paieškas. ###
##################################################################################

maxValues = [0, 0, 0, 0]
maxIds = [0, 0, 0, 0]
# maxIds = [6048, 11480, 18521, 20096]

id = 1
for baseCount in baseCounts:
    minimum = min(maxValues)
    if baseCount > minimum:
        minimumId = maxValues.index(minimum)
        maxValues[minimumId] = baseCount
        maxIds[minimumId] = id
    id+=1

apexes = []
for id in maxIds:
    for i in [-2, -1, 0, 1, 2]:
        apexes.append(id + i)

from Bio.Blast import NCBIWWW
for apex in apexes:
    sequence = sequences[apex - 1]
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1, megablast=True, entrez_query='bacteria')
    res = result_handle.read()
    f = open("blast" + str(apex) + ".xml", "w")
    f.write(res)
    f.close()

import xml.etree.ElementTree as ET
for apex in apexes:
    fileName = "blast" + str(apex) + ".xml"
    tree = ET.parse(fileName)
    root = tree.getroot()
    hit = root.find("BlastOutput_iterations").find("Iteration").find("Iteration_hits").find("Hit")
    hitDef = ""
    if hit:
        hitDef = hit.find("Hit_def").text
    resultFile.write(str(apex) + ": " + hitDef + "\n")

resultFile.close()