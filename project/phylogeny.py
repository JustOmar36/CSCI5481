from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import *
import sys

filename = "alignments.fasta"

"""
This function will read in the alignment file and then make a dictionary out of the names and the sequences
"""
def ReadFastaFile(filename):
  fileObj = open(filename, 'r')
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
    align = AlignIO.read(filename,'fasta')
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
def build_phy_tree(distMatrix, alignment, filename, problem):
  score = 0
  constructor = DistanceTreeConstructor()
  neighbor_Tree = constructor.nj(distMatrix)
  align_dict = ReadFastaFile(filename)

  constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()), neighbor_Tree)
  tree = constructor.build_tree(alignment)


  if problem == str(1):
      f = open("tree_structure.txt", "w")
      f.write(str(tree))
      f.close()
      Phylo.draw(neighbor_Tree)

  else:
      tree_nodes = tree.get_terminals()

      f = open("phy_align.fasta", "w")

      for i in tree_nodes:
          for k,v in align_dict.items():
              temp = k.split(" ")
              if temp[0] == str(i):
                  f.write(">" + str(k) + "\n")
                  f.write(str(v) + "\n")
      f.close()            

      align = AlignIO.read("phy_align.fasta", "fasta")
      f = open("phy_align.fasta", "w")
      f.write(format(align, "fasta"))
      f.close()

      score = ParsimonyScorer().get_score(tree, alignment)
      print(score)
  return score

"""
This is the function I call onto in phySankoff.py
"""
def run_phy_file(args):
    _, problem = args
    dist_matrix, alignment, filename = get_dist_matrix('alignments.fasta')
    score = build_phy_tree(dist_matrix, alignment, filename, problem)
    return score
