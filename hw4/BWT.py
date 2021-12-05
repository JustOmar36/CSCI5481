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

    print(new_seq)
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
    print(encoded_seq[0])
    tracker=[]
    counter = {"A":0,"G":0,"C":0,"T":0,"$":0}
    for i in encoded_seq:
        counter[i]+=1
        tracker.append((i, counter[i]))
    
    starters={"A":1,"G":0,"C":0,"T":0,"$":0}
    starters["c"] = counter["A"]+1
    starters["G"] = starters["C"]+counter["C"]
    starters["T"] = starters["G"]+counter["G"]

    unorder = tracker.copy()
    tracker.sort(key=operator.itemgetter(0))
    final.append(unorder[0])
    letter=unorder[0]

    for i in range(len(unorder)-1):
        index = starters[letter[0]]+ letter[1]-1
        final.append(unorder[index])
        letter=unorder[index]
    
    sequence = ""

    for i in final[:-1]:
        sequence = sequence + i[0]

    sequence = sequence[::-1]

    f = open('rBwtSample2.fa', 'w')
    f.write(sequence)
    f.close()
    #print(sequence)
    return sequence

if __name__ == "__main__":
    # seq = read_file('sample1.fa')
    # encoded_seq = bwt(seq)
    # inverse_bwt(encoded_seq)
    homoSepain = read_file(sepian)
    inverse_bwt(homoSepain)

