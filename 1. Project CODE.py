#Final Project
import Bio as Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq


#import cytb sequences
c1 = SeqIO.read("Varanus komodoensis cytb.txt", "fasta")
c2 = SeqIO.read("Chelonia mydas cytb.txt", "fasta")
c3 = SeqIO.read("Drosophila melanogaster cytb.txt", "fasta")
c4 = SeqIO.read("Danaus plexippus cytb.txt", "fasta")
c5 = SeqIO.read("Taeniopygia guttata cytb.txt", "fasta")
c6 = SeqIO.read("Latimeria chalumnae cytb.txt", "fasta")
c7 = SeqIO.read("Aptenodytes patagonicus cytb.txt", "fasta")
c8 = SeqIO.read("Tachysurus fulvidraco cytb.txt", "fasta")
c9 = SeqIO.read("Homo sapiens cytb.txt", "fasta")
c10 = SeqIO.read("Felis catus cytb.txt", "fasta")


#rename cytb sequences to common names
c1.id = "KomodoDragon"
c2.id = "Turtle"
c3.id = "FruitFly"
c4.id = "Butterfly"
c5.id = "ZebraFinch"
c6.id = "Coelacanth"
c7.id = "Penguin"
c8.id = "YellowCatfish"
c9.id = "Human"
c10.id = "Cat"

#import elongation factor sequences 
ef1 = SeqIO.read("Varanus komodoensis eef1a2.txt", "fasta")
ef2 = SeqIO.read("Chelonia mydas eef1a2.txt", "fasta")
ef3 = SeqIO.read("Drosophila melanogaster eef1a2.txt", "fasta")
ef5 = SeqIO.read("Danaus plexippus eef1a2.txt", "fasta")
ef4 = SeqIO.read("Taeniopygia guttata eef1a2.txt", "fasta")
ef6 = SeqIO.read("Latimeria chalumnae eef1a2.txt", "fasta")
ef7 = SeqIO.read("Aptenodytes patagonicus eef1a2.txt", "fasta")
ef8 = SeqIO.read("Tachysurus fulvidraco eef1a2.txt", "fasta")
ef9 = SeqIO.read("Homo sapiens eef1a2.txt", "fasta")
ef10 = SeqIO.read("Felis catus eef1a2.txt", "fasta")

#rename elongation factor sequences to common names
ef1.id = "KomodoDragon"
ef2.id = "Turtle"
ef3.id = "FruitFly"
ef4.id = "Butterfly"
ef5.id = "ZebraFinch"
ef6.id = "Coelacanth"
ef7.id = "Penguin"
ef8.id = "YellowCatfish"
ef9.id = "Human"
ef10.id = "Cat"

#combined files
animalscytb = SeqIO.write((c1,c2,c3,c4,c5,c6,c7,c8,c9,c10), "animalscytb.txt", "fasta")
animalsef = SeqIO.write((ef1,ef2,ef3,ef4,ef5,ef6,ef7,ef8,ef9,ef10), "animalsef.txt", "fasta")

import subprocess #allows program to run clustalo 

notalignedcytb = "animalscytb.txt"
alignedanimalscytb = "animalsalignedcytb.fasta"
notalignedef= "animalsef.txt"
alignedanimalsef = "animalsalignedef.fasta"


clustalolocation= "/opt/homebrew/bin/clustalo" #so that program knows where to find clustalo

def alignsequence(notaligned, alignedanimals,clustalolocation):
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
alignedsequenceef = alignsequence(notalignedef, alignedanimalsef,clustalolocation)

import subprocess
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

insects = ["FruitFly", "Butterfly"]

def treeconstructor(alignedsequence, insects):
    calculator = DistanceCalculator('blastn') #tells the program to use the blastn model to calculate the genetic distance
    distancematrix = calculator.get_distance(alignedsequence) #determined the length of the branches
    #print(distancematrix)
    constructor = DistanceTreeConstructor(calculator)
    treeinfo = constructor.build_tree(alignedsequence)
    treeinfo.rooted = True
    #print(treeinfo)
    treeinfo.root_with_outgroup(insects)#forces program to build tree with insects as the ones that split first
    asciitree = Phylo.draw_ascii(treeinfo)
    return treeinfo
##    tree = Phylo.write(treeinfo, "tree.xml", "phyloxml")

treecytb = treeconstructor(alignedsequencecytb, insects)
Phylo.write(treecytb, "cytbtree.xml", "phyloxml")
treeef = treeconstructor(alignedsequenceef, insects)
Phylo.write(treeef, "eftree.xml", "phyloxml")

from Bio.Align import MultipleSeqAlignment

def concatenate(align1, align2):
    combinedSequences = [] #the new combined sequences will be store here
    
    for i in range(len(align1)): #loops pver the number of species
        rec1 = align1[i] #gets i from cytb
        rec2 = align2[i] #gets i from ef

        newSeq = rec1.seq + rec2.seq #combines sequences
        rec = rec1[:]         # copy the record so we can replace its sequence without changing the original rec1 alignment
        rec.seq = newSeq #replace the old with the concatenated version
        combinedSequences.append(rec) #add combined sequence into the list
        finalAlignment = MultipleSeqAlignment(combinedSequences)

    return finalAlignment

concatenatedAlignment = concatenate(alignedsequencecytb, alignedsequenceef)
AlignIO.write(concatenatedAlignment, "concatenatedAlignment.fasta", "fasta")
treecombined = treeconstructor(concatenatedAlignment, insects)
Phylo.write(treecombined, "treecombined.xml", "phyloxml")


alignment_file = "concatenatedAlignment.fasta"
iqtreelocation= "/opt/homebrew/bin/iqtree3"

def iqtree(iqtree, aligmnentfile):
    constraint_text = """(
        (FruitFly,Butterfly),
        (KomodoDragon,Turtle),
        ZebraFinch,
        Penguin,
        Human,
        Cat,
        Coelacanth,
        YellowCatfish
    );"""

    with open("constraint.nwk", "w") as f:
        f.write(constraint_text)

    cmd = [
        iqtreelocation,
        "-s", alignment_file,
        "-o", "FruitFly,Butterfly",
        "-g", "constraint.nwk", 
        "-m", "TESTNEW", 
        "-bb", "1000",
        "-alrt", "1000",
        "--redo" #run the calculation from start to finish, even if there is an existing output file
    ]
    subprocess.run(cmd, check=True)

    treefile = alignment_file + ".treefile"
    tree = Phylo.read(treefile, "newick")
    Phylo.draw_ascii(tree)
    return tree

finaltree = iqtree(iqtreelocation, alignment_file)
Phylo.write(finaltree, "IQtree.xml", "phyloxml")

import matplotlib
import matplotlib.pyplot as plt

def drawtree (asciitree, title):
    fig = plt.figure(figsize = (15,5), dpi = 100)
    matplotlib.rc("font", size = 10)
    matplotlib.rc("xtick", labelsize = 10)
    matplotlib.rc("ytick", labelsize = 10)
    axes = fig.add_subplot(1,1,1)
    plt.title(title)
    Phylo.draw(asciitree, axes = axes)
    return fig
   
drawtree(treecytb, "Cytochrome B (cytb)")
drawtree(treeef,"Elongation Factor 1 Alpha 2 (eef1a2)")
drawtree(treecombined,"cytb + eef1a")
drawtree(finaltree, "Final Tree")



