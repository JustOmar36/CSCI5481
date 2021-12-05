import operator
import sys
# sampleFile ='sample1.fa'
# sepian = 'bwtHomoSapiens.fa'
# chrom = 'bwtChrYnew(25M-26M).fa'
# sr1 = 'SRR089545_1.fq'
# sr2 = 'SRR089545_2.fq'

"""This function reads in fa files"""
def read_file(filename):
    sequence = ''
    f = open('files/' + filename, 'r')
    for l, i in enumerate(f):
        if l != 0:
            sequence += i
    f.close()
    sequence = sequence.replace('\n', '')
    sequence = sequence.replace('N', '')
    return sequence

'''
This function does the BWT functionality
Appends '$' to the end of the sequence
'''
def bwt(args):
    _, filename, num =  args
    sequence = read_file(filename)
    arr = []
    new_sequences = list(sequence)
    new_sequences.append('$')
    seq = ''

    for i in range(len(new_sequences)):
        new_sequences = new_sequences[1:] + new_sequences[:1]
        seq = seq.join(new_sequences)
        arr.append(seq)
        seq = ''

    arr.sort()
    bwt_arr = []
    new_seq = ''
    for i in arr:
        bwt_arr.append(i[len(i) -1])
    new_seq = new_seq.join(bwt_arr)

    f= open('files/bwtSample1.fa', 'w')
    f.write(str(new_seq))
    f.close()
    return new_seq

"""This function creates the rank arrays for FM function"""
def give_rank(seq):
    dict_of_things = dict()
    ranks = []
    for i in seq:
        if i not in dict_of_things:
            dict_of_things[i] = 0
        ranks.append(dict_of_things[i])
        dict_of_things[i] += 1
    return ranks, dict_of_things

"""This maps counts to the letters"""
def mapping(dict_of_things):
    first = {}
    items_counts = 0
    for i, count in sorted(dict_of_things.items()):
        first[i] = (items_counts, items_counts + count)
        items_counts += count
    return first

'''This does the inverse functionality using the mapping and rank function written previously'''
def inverse_bwt(args):
    _, filename, num = args
    sequence = read_file(filename)
    ranks, dict_of_things = give_rank(sequence)
    first = mapping(dict_of_things)
    i = 0
    symbol = "$"
    while sequence[i] != '$':
        letter = sequence[i]
        symbol = letter + symbol
        i = first[letter][0] + ranks[i]
    
    f = open('files/rBwtSample2.fa', 'w')
    f.write(symbol)
    f.close()
    return symbol

'''This function reads fq files'''
def read_fq(filename):
    count = 0
    sequence = ''
    sequences = []
    with open('files/' + filename, 'r') as f:
        for l in f:
            if l.startswith('@'):
                count += 1
            elif l.startswith('+'):
                next(f)
            else:
                sequence = l
                sequence = sequence[:len(sequence)-1]
                sequence = sequence.replace('N', '')
                sequences.append(sequence)
    f.close()
    return sequences

'''
This runs the functionality/logic to find matching sequences. 
Sort of uses a selection sort algorithm idea to find the indices and counts of how many matching sequences exist
'''
def fullindexFM(sequence, match):
    ranks, items = give_rank(sequence)
    first = mapping(items)
    left, right = first[match[-1]]
    i = len(match)-2

    while i >= 0 and right > left:
        searching_for = match[i]
        j = left
        while j < right:
            if sequence[j] == searching_for:
                left = first[searching_for][0] + ranks[j]
                break
            j += 1
        if j == right:
            left = right
            break 
        right -= 1
        while sequence[right] != searching_for:
            right -= 1
        right = first[searching_for][0] + ranks[right] + 1
        i -= 1
    return right - left, left

'''
This use fullindexFM to print how many matches we have and their indices
The first array to print to file is for READ1 and second array is for READ2
'''
def run_fullIndexFM(args):
    _, chrom, filename1, filename2, num = args
    count = 0
    reads1 = read_fq(filename1)
    reads2 = read_fq(filename2)
    seq = read_file(chrom)
    read1_idx = []
    read2_idx = []
    total_reads = len(reads1) + len(reads2)
    for i in reads1:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read1_idx.append(q)
            count+=1
    for i in reads2:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read2_idx.append(q)
            count+=1
    
    f = open('files/mapping_fullindexFM.txt', 'w')
    f.write(str(count) + ' out of ' + str(total_reads) + ' short read pairs map on chrY:20M~40M: \n')
    f.write(str(read1_idx))
    f.write('\n')
    f.write(str(read2_idx))
    f.close()
    return

'''This function is built as a helper function that recudes the amount of ranks we have by reducing the amount of reads we use by 30%.'''
def third(arr, new_arr):
    if len(arr) == 0:
        return new_arr

    if len(arr) >= 3:
        for i in range(3):
            new_arr.append(arr.pop(0))
    else:
        for i in range(len(arr)):
            new_arr.append(arr.pop(0))
            
    return third(arr[7:], new_arr)

'''This function is built as a helper function that recudes the amount of ranks we have by reducing the amount of reads we use by 60%.'''
def sixth(arr, new_arr):
    if len(arr) == 0:
        return new_arr

    if len(arr) >= 6:
        for i in range(6):
            new_arr.append(arr.pop(0))
    else:
        for i in range(len(arr)):
            new_arr.append(arr.pop(0))

    return third(arr[4:], new_arr)

'''This function runs Hybrid Functionality'''
def run_HybridIndexFM(args):
    _, chrom, filename1, filename2, num = args
    count = 0
    reads1 = read_fq(filename1)
    reads2 = read_fq(filename2)
    reads1 = third(reads1, [])
    reads2 = third(reads2, [])
    seq = read_file(chrom)
    read1_idx = []
    read2_idx = []
    total_reads = len(reads1) + len(reads2)
    for i in reads1:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read1_idx.append(q)
            count+=1
    for i in reads2:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read2_idx.append(q)
            count+=1
    
    f = open('files/mapping_FMhybrid.txt', 'w')
    f.write(str(count) + ' out of ' + str(total_reads) + ' short read pairs map on chrY:20M~40M: \n')
    f.write(str(read1_idx))
    f.write('\n')
    f.write(str(read2_idx))
    f.close()
    return

'''This function runs checkpoint Functionality'''
def run_checkpointingIndexFM(args):
    _, chrom, filename1, filename2, num = args
    count = 0
    reads1 = read_fq(filename1)
    reads2 = read_fq(filename2)
    reads1 = sixth(reads1, [])
    reads2 = sixth(reads2, [])
    seq = read_file(chrom)
    read1_idx = []
    read2_idx = []
    total_reads = len(reads1) + len(reads2)
    for i in reads1:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read1_idx.append(q)
            count+=1
    for i in reads2:
        c, q = fullindexFM(seq, i)
        if c > 0:
            read2_idx.append(q)
            count+=1
    
    f = open('files/mapping_FMCheckpointing.txt', 'w')
    f.write(str(count) + ' out of ' + str(total_reads) + ' short read pairs map on chrY:20M~40M: \n')
    f.write(str(read1_idx))
    f.write('\n')
    f.write(str(read2_idx))
    f.close()
    return

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 3:
        if(args[1][-3:] == '.fa'):
            if (args[2].isnumeric()):
                if(int(args[2]) == 1):
                    print('Running BWT...')
                    bwt(args)
                elif(int(args[2]) == 2):
                    print('Running BWT-Reverse')
                    inverse_bwt(args)
                else: 
                    print('If you are only entering 2 arguements then the second should be either numbers 1 or 2')
            else:
                print('Second argument should be an integer')
        else:
            print('First argument should be a file of type .fa')
    elif len(args) == 5:
        if (args[1][-3:] == '.fa'):
            if (args[2][-3:] == '.fq'):
                if (args[3][-3:] == '.fq'):
                    if (args[4].isnumeric()):
                        if (int(args[4])) == 1:
                            print('Running Full-Index FM...')
                            run_fullIndexFM(args)
                        elif (int(args[4])) == 2:
                            print('Running Hybrid Full-Index FM...')
                            run_HybridIndexFM(args)
                        elif (int(args[4])) == 3:
                            print('Running Checkpoint Full-Index FM...')
                            run_checkpointingIndexFM(args)
                        else:
                            print('If you are entering 4 arguments then the 4th must be numbers 1, 2 or 3')
                    else:
                        print('Fourth argument must be an integer')
                else:
                    print('Third argument must be set to .fq')
            else:
                    print('second argument must be set to .fq')
        else:
            print('First argument must be set to .fa')
    else:
        print('You can only enter 2 arguments to run BWT or BWT-Reverse or 4 arguments to run FM variations.')
