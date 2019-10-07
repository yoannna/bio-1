#from Bio import SeqIO
#for seq_record in SeqIO.parse("plazmide.fasta", "fasta"):
#    print(seq_record.id)
#    print(repr(seq_record.seq))
#    print(len(seq_record))

plazmide = open("plazmide.fasta")
plazmide.readline()

startCodon = "ATG"
stopCodon1 = "TAA"
stopCodon2 = "TAG"
stopCodon3 = "TGG"

reading = False
tirpletCount = 0
end = False

sequence = []
seqCount = 0

while True:
    if end:
        break
    for i in range(1, 20):
        if end:
            break
        tirpletCount += 1
        triplet = plazmide.read(3)
        if triplet == "":
            end = True
            break
        if reading == False:
            if triplet == startCodon:
                reading = True
                print(triplet, end=" ")
                sequence.append(triplet)
        else:
            print(triplet, end=" ")
            sequence.append(triplet)
            if triplet == stopCodon1 or triplet == stopCodon2 or triplet == stopCodon2:
                print ("\n\n")
                if len(sequence) > 33:
                    seqCount += 1
                reading = False
    else:
        plazmide.readline()

print (seqCount)
