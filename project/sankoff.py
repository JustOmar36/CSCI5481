'''
The SankOff algorithm 
'''
from math import inf
from phylogeny import run_phy_file
p_inf = inf

scoring_matrix = [
   [0, 3, 4, 9, 8],
   [3, 0, 2, 4, 8],
   [4, 2, 0, 4, 8],
   [9, 4, 4, 0, 8],
   [8, 8, 8, 8, 8]
];

"""
This function will check the nucleotide we are at and then give index a number to go into the matrix later
"""
def get_from_matrix(x):
    index = 0
    
    if x == 'A':
        index = 0
    elif x == 'U':
        index = 1
    elif x == 'G':
        index = 2
    elif x == 'C':
        index = 3
    else:
        index = 4
        
    return index

"""
This function will check what index we are at and then give it the proper nucleotide
"""
def give_char(index):
    nuc = ''
    if index == 0:
        nuc = 'A'
    elif index == 1:
        nuc = 'U'
    elif index == 2:
        nuc = 'G'
    elif index == 3:
        nuc = 'C'
    else:
        nuc = '-'
    return nuc

"""
This function will go through and initalize the arrays where the proper index is 0 and the rest is infinity
"""
def init_array(first):
    array = []
    for i in range(5):
        array.append(p_inf)
    array[get_from_matrix(first)] = 0
    return array

"""
This function will go down the built tree structure and give the appropirate nucleotide based on the min values

If there are two min values the function will go into the matrix and pick the min distance of the two
"""
def go_down_tree(scored_arr):
    char_array = []
    for i in range(len(scored_arr)):
        min_idx1 = scored_arr[i].index(min(scored_arr[i]))
        value_min1 = scored_arr[i].pop(min_idx1)
        min_idx2 = scored_arr[i].index(min(scored_arr[i]))
        if value_min1 == scored_arr[i][min_idx2]:
            if scoring_matrix[min_idx1][min_idx1+1] < scoring_matrix[min_idx2][min_idx1+1]:
                char_array.append(give_char(min_idx1))
            elif scoring_matrix[min_idx1][min_idx1+1] > scoring_matrix[min_idx1][min_idx2+1]:
                char_array.append(give_char(min_idx2))
        else:
            char_array.append(give_char(min_idx1))   

    return char_array[0]

"""
This function will go up the tree one nucleotide at a time and return the appropriate array with appropriate min values of all
possible outcomes between 
A, U, C, G
"""
def traverse_tree(arr1, arr2):
    
    min_arr1 = []
    min_arr2 = []
    new_arr = [[],[],[],[],[]]
    new_arr2 = [[],[],[],[],[]]
    final_array = []

    for i in range(len(arr1)):
        for j in range(len(scoring_matrix[i])):
            temp = arr1[j]
            new_arr[i].append(temp + scoring_matrix[i][j])
   
    for i in new_arr:
        min_arr1.append(min(i))
    
    for i in range(len(arr2)):
        for j in range(len(scoring_matrix[i])):
            temp2 = arr2[j]
            new_arr2[i].append(temp2 + scoring_matrix[i][j])

    for i in new_arr2:
        min_arr2.append(min(i))

    for i in range(len(min_arr1)):
        final_array.append(min_arr1[i] + min_arr2[i])
       
    return final_array

"""
This function will use the traverse function and run the sankoff logic on each nucleotide
"""
def run_sankof_on_index(sequences_list):
    array_list = []
    going_up = []
    final = []
    final_result = []
    for i in range(len(sequences_list)):
        for j in range(len(sequences_list[0])):
            array_list.append(init_array(sequences_list[i][j]))

    inverse_array_list = array_list[::-1]
    for i in range(0, len(array_list), 2):
        final.append(inverse_array_list[i])
        final.append(inverse_array_list[i+1])
        result = traverse_tree(array_list[i], array_list[i+1])
        going_up.append(result)

    for i in range(len(going_up)):
        if len(going_up) >= 2:
            final_result = traverse_tree(going_up[0], going_up[1])
            final.append(final_result)
            final.append(going_up.pop(1))
            final.append(going_up.pop(0))
            going_up.append(final_result)
        else:
            final.append(going_up[i-1])
    final = final[::-1]

    return go_down_tree(final)

"""
This function will grab each nucleotide and run the "run_sankoff_on_each_index" function on each and then return the
appropriate nucleotide and create the next sequence in the tree.
"""
def get_sequences_back(seq1, seq2):
    new_seq = ""
    for i in range(len(seq1)):
        nuc_list = []
        nuc_list.append([seq1[i]])
        nuc_list.append([seq2[i]])

        for i in range(0, len(nuc_list), 2):
            new_seq += run_sankof_on_index(nuc_list)
    
    
    return new_seq

"""
This function runs all functions on the sequences
"""
def sankoff(sequences):
    new_seq = ""
    all_sequences = []
    odd_seq = ""

    if len(sequences) % 2 == 1:
        odd_seq = sequences.pop(len(sequences) -1)

    for i in sequences:
        all_sequences.append(i)


    for i in range(0, len(sequences), 2):
        i = 0
        new_seq = get_sequences_back(sequences[i], sequences[i+1])
        all_sequences.append(new_seq)
        sequences.append(new_seq)
        sequences.pop(0)
        sequences.pop(0)
        if len(sequences) == 2:
            new_seq = get_sequences_back(sequences[0], sequences[1])
            all_sequences.append(new_seq)
            sequences.append(odd_seq)

    if len(sequences) != 1:
        sankoff(sequences)
    
    return all_sequences


"""
This function will run only if we are running the toy sequences, but it will sort of print the sequences in tree structure
"""
def print_tree(list_of_sequences):
    array_degrees = []
    temp_array = list_of_sequences
    for i in range(int((len(list_of_sequences)/2))):
        array_degrees.append(i)
    
    for i in range(len(array_degrees)):
        array_degrees[i] = array_degrees[i]*2

    f = open("toy_result.txt", "w")
    f.write(str(list_of_sequences))
    f.write("\n")
    for i in array_degrees:
        if i == 0:
            f.write(str(temp_array[0:1]))
            f.write("\n")
            temp_array.pop(0)
        else:
            f.write(str(temp_array[:i]))
            f.write("\n")
            for j in range(i -1):
                temp_array.pop(0)

"""
This is the function that is called in phySankOff

If the arg is 2 then it will run the toy sequence problem,
otherwise it will run on the 38 sequences
"""
def run_sankoff(args):
    _, problem = args

    if problem == str(2):
        sequence = ""
        print("For part 2 the only file avaliable is: toy_alignments.txt (no quotation marks needed)")
        filename = input("Enter toy sequence filename: ")
        f = open(filename, "r")
        for i in f:
            print(i)
            sequence += i
        sequence = sequence.replace(" ", '')
        seq_list = sequence.split(',')
        print(seq_list)
        array = sankoff(seq_list)
        print_tree(array)
    
    else:
        score = run_phy_file(args)
        print(score)
        sequences = []
        sequence = ""

        f = open('phy_align.fasta', 'r')
        for line, i in enumerate(f):
            if line == 0:
                continue

            if i.startswith('>'):
                sequence = sequence.replace('\n', '')
                sequences.append(sequence)
                sequence = ""
            else:
                sequence += i
        f.close()
        
        array = sankoff(sequences)
        
        f = open('sankoff_result.txt', 'w')
        f.write("parsimony result = " + str(score) + "\n")
        for i in array:
            f.write(i)
            f.write("\n")
        f.close()

