
import operator
sampleFile ='sample1.fa'
sepian = 'bwtHomoSapiens.fa'

def read_file(filename):
    sequence = ''
    f = open(filename, 'r')
    for l, i in enumerate(f):
        if l != 0:
            sequence += i
    f.close()
    sequence = sequence.replace('\n', '')
    return sequence

def bwt(sequence):
    arr = []
    new_sequences = list(sequence)
    new_sequences.append('$')
    #new_sequences = ['b', 'a', 'n', 'a', 'n', 'a', "$"]
    seq = ''

    for i in range(len(new_sequences)):
        new_sequences = new_sequences[1:] + new_sequences[:1]
        seq = seq.join(new_sequences)
        arr.append(seq)
        seq = ''

    #print(arr)
    arr.sort()
    #print(arr)

    bwt_arr = []
    new_seq = ''
    for i in arr:
        bwt_arr.append(i[len(i) -1])
    new_seq = new_seq.join(bwt_arr)

    #print(new_seq)
    f= open('bwtSample1.fa', 'w')
    f.write(str(new_seq))
    f.close()
    return new_seq

# def inverse_bwt(encoded_seq):
#     encode_str = ''
#     encode = list(encoded_seq)
#     org_str = encoded_seq
#     encoded_seq = list(encoded_seq)
    
#     encode.sort()
#     count = 0
#     while(count < len(encoded_seq) -1):
#         for i in range(len(encode)):
#             encoded_seq[i] += encode[i]
#         encode = encoded_seq
#         encode.sort()
#         encoded_seq = list(org_str)
#         count = count + 1
    
#     #print(encode)

#     for i in encode:
#         if i[len(i) - 1] == '$':
#             encode_str = i
#     print(encode_str)
#     return encode_str

def give_rank(bw):
    dict_of_things = dict()
    ranks = []
    for i in bw:
        if i not in dict_of_things:
            dict_of_things[i] = 0
        ranks.append(dict_of_things[i])
        dict_of_things[i] += 1
    return ranks, dict_of_things

def mapping(dict_of_things):
    first = {}
    items_counts = 0
    for i, count in sorted(dict_of_things.items()):
        first[i] = (items_counts, items_counts + count)
        items_counts += count
    return first

def inverse_bwt(sequence):
    ranks, dict_of_things = give_rank(sequence)
    first = mapping(dict_of_things)
    i = 0
    symbol = "$"
    while sequence[i] != '$':
        letter = sequence[i]
        symbol = letter + symbol
        i = first[letter][0] + ranks[i]
    
    f = open('rBwtSample2.fa', 'w')
    f.write(symbol)
    f.close()
    return symbol

# if __name__ == "__main__":
#     args = sys.argv
#     if len(args) == 2:
#         if (args[1].isnumeric()):
#             if (int(args[1])) == 1:
#                 print("Running BWT on Sample1")
#                 bwt(read_file(sampleFile))
#             if(int(args[1])) == 2:
#                 print("Running reverseBWT on bwtHomoSapiens.fa")
#                 inverse_bwt(read_file(sapiensFile))
# =======

if __name__ == "__main__":
    # seq = read_file('sample1.fa')
    # encoded_seq = bwt(seq)
    # inverse_bwt(encoded_seq)
    homoSepain = read_file(sepian)
    inverse_bwt(homoSepain)

