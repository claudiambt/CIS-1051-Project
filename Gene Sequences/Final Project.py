#Final Project
import Bio as Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo


#import cytb sequences
c1 = SeqIO.read("Anolis carolinensis cytb.fasta", "fasta")
c2 = SeqIO.read("Chelonia mydas cytb.txt", "fasta")
c3 = SeqIO.read("Drosophila melanogaster cytb.txt", "fasta")
c4 = SeqIO.read("Danaus plexippus cytb.txt", "fasta")
c5 = SeqIO.read("Taeniopygia guttata cytb.txt", "fasta")
c6 = SeqIO.read("Danio rerio cytb.txt", "fasta")
c7 = SeqIO.read("Erithacus rubecula cytb.txt", "fasta")
c8 = SeqIO.read("Rhincodon typus cytb.txt", "fasta")
c9 = SeqIO.read("Homo sapiens cytb.txt", "fasta")
c10 = SeqIO.read("Felis catus cytb.txt", "fasta")


#rename cytb sequences to common names
c1.id = "GreenAnole"
c2.id = "Turtle"
c3.id = "FruitFly"
c4.id = "Butterfly"
c5.id = "ZebraFinch"
c6.id = "ZebraFish"
c7.id = "EuropeanRobin"
c8.id = "Shark"
c9.id = "Human"
c10.id = "Cat"

#import 12S rRNA sequences to common names
r1 = SeqIO.read("Anolis carolinensis 12S.fasta", "fasta")
r2 = SeqIO.read("Chelonia mydas 12S.fasta", "fasta")
r3 = SeqIO.read("Drosophila melanogaster 12S.fasta", "fasta")
r5 = SeqIO.read("Danaus plexippus 12S.fasta", "fasta")
r4 = SeqIO.read("Taeniopygia guttata 12S.fasta", "fasta")
r6 = SeqIO.read("Danio rerio 12S.fasta", "fasta")
r7 = SeqIO.read("Erithacus rubecula 12S.fasta", "fasta")
r8 = SeqIO.read("Rhincodon typus 12S.fasta", "fasta")
r9 = SeqIO.read("Homo sapiens 12S.fasta", "fasta")
r10 = SeqIO.read("Felis catus 12S.fasta", "fasta")

#rename 12S rRNA sequences to common names
r1.id = "GreenAnole"
r2.id = "Turtle"
r3.id = "FruiFly"
r4.id = "Butterfly"
r5.id = "ZebraFinch"
r6.id = "ZebraFish"
r7.id = "EuropeanRobin"
r8.id = "Shark"
r9.id = "Human"
r10.id = "Cat"

#combined files
animalscytb = SeqIO.write((c1,c2,c3,c4,c5,c6,c7,c8,c9,c10), "animalscytb.txt", "fasta")
animalsrRNA = SeqIO.write((r1,r2,r3,r4,r5,r6,r7,r8,r9,r10), "animalsrRNA.txt", "fasta")

import subprocess #allows program to run clustalo 

notalignedcytb = "animalscytb.txt"
alignedanimalscytb = "animalsalignedcytb.fasta"
notalignedrRNA = "animalsrRNA.txt"
alignedanimalsrRNA = "animalsalignedrRNA.fasta"


clustalolocation= "/opt/homebrew/bin/clustalo" #so that program knows where to find clustalo

def alignsequence (notaligned, alignedanimals,clustalolocation):
    cmd = [
        clustalolocation, #tells the program where to look for clustalo
        "-i", notaligned, #specifies the input file
        "-o", alignedanimals, #specifies the putput file
        "--seqtype=DNA", #tells clustalo to treat the sequence as DNA
        "--force" #tells clustalo to overwrite the output file if it already exists
    ]

    subprocess.run(cmd, check=True)
    alignment = AlignIO.read(alignedanimals, "fasta")
    print(alignment)
    return alignment

alignedsequencecytb = alignsequence(notalignedcytb, alignedanimalscytb,clustalolocation)
alignedsequencerRNA = alignsequence(notalignedrRNA, alignedanimalsrRNA,clustalolocation)


from Bio.Phylo.TreeConstruction import DistanceCalculator
calculator = DistanceCalculator('blastn') #tells the program to use the identity model to calculate the genetic distance
distancematrix = calculator.get_distance(alignedsequencecytb) #determined the length of the branches
print(distancematrix)

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor(calculator)

treeinfo = constructor.build_tree(alignedsequencecytb)
treeinfo.rooted = True
print(treeinfo)

insects = ["FruitFly", "Butterfly"]
treeinfo.root_with_outgroup(insects)#forces program to build tree with insects as the ones that split first
asciitree = Phylo.draw_ascii(treeinfo)

tree = Phylo.write(treeinfo, "tree.xml", "phyloxml")








