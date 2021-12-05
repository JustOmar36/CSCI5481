
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

def inverse_bwt(encoded_seq):
    final = []
    forward = sorted(encoded_seq)
    tracker=[]
    counts = {"A":0,"G":0,"C":0,"T":0,"$":0}

    for i in encoded_seq:
        counts[i]+=1
        tracker.append((i, counts[i]))
    
    tracking_dict={"A":1,"G":0,"C":0,"T":0,"$":0}
    tracking_dict["c"] = counts["A"]+1
    tracking_dict["G"] = tracking_dict["C"]+counts["C"]
    tracking_dict["T"] = tracking_dict["G"]+counts["G"]

    unorder = tracker.copy()
    tracker.sort(key=operator.itemgetter(0))
    final.append(unorder[0])
    letter=unorder[0]

    for i in range(len(unorder)-1):
        index = tracking_dict[letter[0]]+ letter[1]-1
        final.append(unorder[index])
        letter=unorder[index]
    
    sequence = ""
    for i in final[:-1]:
        sequence = sequence + i[0]
    sequence = sequence[::-1]
    f = open('rBwtSample2.fa', 'w')
    f.write(sequence)
    f.close()

    return sequence
    
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

