from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import *
import random
userinput1 = ''
userinput2 = ''

"""
This is a broken down function to reconstruct FASTA files from the tree structure given by the other functions
"""
def create_sankoff_files(njTree, align_dict, name_of_file):
  name_of_file = str(name_of_file)
  tree_nodes = njTree.get_terminals()

  f = open("output/" + name_of_file, "w")

  for i in tree_nodes:
      for k,v in align_dict.items():
          temp = k.split(" ")
          if temp[0] == str(i):
              f.write(">" + str(k) + "\n")
              f.write(str(v) + "\n")
  f.close()            

  align = AlignIO.read('output/' + name_of_file, "fasta")
  f = open('output/' + name_of_file, "w")
  f.write(format(align, "fasta"))
  f.close()
  
"""
This function will read in the alignment file and then make a dictionary out of the names and the sequences
"""
def ReadFastaFile(filename):
  fileObj = open("input/" + filename, 'r')
  sequences = dict()
  seqFragments = []
  id= 0
  for line in fileObj:
    if line.startswith('>'):
      if seqFragments:
        sequence = ''.join(seqFragments)
        sequences[id] = sequence
      seqFragments = []
      id = line.rstrip()[1:]
    else:
      seq = line.rstrip()
      seqFragments.append(seq)
  if seqFragments:
    sequence = ''.join(seqFragments)
    sequences[id] = sequence
  fileObj.close()
  return sequences

"""
This function returns a matrix of distances, the aligned sequences from the fasta file and the filename
"""
def get_dist_matrix(filename):
    align = AlignIO.read('input/' + filename,'fasta')
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)
    #print(distMatrix)
    return distMatrix, align, filename

"""
This function constructs the neighborhood-joining tree from the distance matrix

If the argument provided is 1 then the function will write to a file the tree strucutre returned by BioPython and then
print out a phylogeny tree

otherwise 

the function goes into the tree structure and grabs all nodes with sequences in them. It will then write those results to a 
a new file and then the format function will format it into a new fasta file called "phy_align.py" for later use. 
(I built sankoff kinda weird so bare with me)

This takes about a minute or two to run
"""
def build_phy_tree(distMatrix, alignment, filename):
  constructor = DistanceTreeConstructor()
  neighbor_Tree = constructor.nj(distMatrix)
  align_dict = ReadFastaFile(filename)

  constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()), neighbor_Tree)
  tree = constructor.build_tree(alignment)

  f = open("tree_structure.txt", "w")
  f.write(str(tree))
  f.close()
  Phylo.draw(neighbor_Tree)
  
  global userinput1
  userinput1 = input("Enter an output file name for a reconstructed Fasta file using this tree's information (Please ensure it is of type Fasta): ")
  create_sankoff_files(tree, align_dict, userinput1)

  return neighbor_Tree

"""
This function is incharge of finding a random node from the tree to reinsert in another place
"""
def prunning(njTree, alignment, filename):
  align_dict = ReadFastaFile(filename)
  clade = random.choice(njTree.get_nonterminals())
  prunned_tree = njTree.from_clade(clade)
  Phylo.draw(prunned_tree)
  njTree.root_with_outgroup(clade)
  
  constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()), njTree)
  tree = constructor.build_tree(alignment)

  global userinput2
  userinput2 = input("Enter an output file name for a reconstructed Fasta file using this tree's information (Please ensure it is of type Fasta): ")

  create_sankoff_files(tree, align_dict, userinput2)
  Phylo.draw(njTree)
     
def run_phy_file(filename):
  dist_matrix, alignment, filename = get_dist_matrix(filename)
  tree = build_phy_tree(dist_matrix, alignment, filename)
  prunned_tree = prunning(tree, alignment, filename)
  
