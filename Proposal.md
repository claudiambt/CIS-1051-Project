
- What will (likely) be the title of your project?
Cytochrome b Phylogenetic Tree

- In just a sentence or two, summarize your project:
Creation of a phylogenetic tree using Cytochrome b as the reference gene. The tree will show the evolutionary relationship between 2 species of reptiles, birds, fish, mammals, and     insects each. 
  
- In a paragraph or more, detail your project. What will your software do? What features will it have? How will it be executed?
This project involves the construction of a basic phylogenetic tree using Python's BioPython's library, starting with raw genetic data. The goal is to visually represent the evolutionary relationships between a small set of animals by comparing a short, common gene like Cytochrome b (cytb). The process includes three steps: getting sequence data for the chosen gene from NCBI GenBank, generating a Multiple Sequence Alignment of these sequences using an external tool like MUSCLE, and using BioPython's Bio.Phylo module to read the alignment, calculate the genetic between species, apply the Neighbor-joining (NJ) algorithm to build the tree structure, and output the result both as an ASCII tree and a graphical plot using Matplotlib. The final result is a functional, visually interpreted model of evolutionary history.

- In the world of software, most everything takes longer to implement than you expect. And so it's not uncommon to accomplish less in a fixed amount of time than you hope. In a sentence (or list of features), define a GOOD outcome for your final project. I.e., what WILL you accomplish no matter what?
Get the gene sequence for all of the species, read the files, combine the sequences into a single file, and change the id so the tree displays common names instead of scientific names.

- In a sentence (or list of features), define a BETTER outcome for your final project. I.e., what do you THINK you can accomplish before the final project's deadline?
Aligning the sequences and calculating the genetic distances between species

- In a sentence (or list of features), define a BEST outcome for your final project. I.e., what do you HOPE to accomplish before the final project's deadline?
  Being able to visualize the phylogenetic tree

- In a paragraph or more, outline your next steps. What new skills will you need to acquire? What topics will you need to research? If working with one of two classmates, who will do what?
Import BioPython, get gene code for each of the animals using GenBank, Read each file individually, Combine individual sequences into a single file, Sequence alignment (there are several   ways to do this), Calculate the distance between each of the different species, Build the tree using tree constructor in BioPython, Visualize the tree using matplotlib. By working on this project, I will learn how to use Biopython, which is useful for my career since I am a Genomic Medicine major and bioinformatics is becoming more and more important. I will need to research about how alignment and distance relate to each other to build the phylogenetic tree. 

