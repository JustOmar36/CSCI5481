Omar Masri
Masri013
5571773
###########################################################################################################################

Answer for 3.4:

I do not have an exact measurment to answer this question but as an estimiate...
3.1 takes about 15 minutes to fully run and print
3.2 and 3.3 run in almost the same amount of timee, but 3.3 is slightly faster when running it on chromsome Y
Both 3.2 and 3.3 take about 5-8 minutes to run

###########################################################################################################################

How to run program:

To run BWT and BWT-Reverse the following format is recommended:

BWT: python3 source/BWT.py [filename] 1
BWT-Reverse: python3 source/BWT.py [filename] 2

In order to run both of these function you must only have 2 arguments

###########################################################################################################################

To run Full-Index-FM, Hybrid-Index-FM, Checkpoint-Index-FM the following format is recommended:

Full-Index-FM: python3 source/BWT.py [chromsome_file_name] [first_reads_file] [second_reads_file] 1
Hybrid: python3 source/BWT.py [chromsome_file_name] [first_reads_file] [second_reads_file] 2
checkpoint: python3 source/BWT.py [chromsome_file_name] [first_reads_file] [second_reads_file] 3

In order to run both of these function you must only have 4 arguments

###########################################################################################################################

Functions:

read_file(filename): 
takes in a file that is of type .fa 
returns it as a sequence

bwt(args): 
takes in the arguemnts provided by user
returns and writes an encoded sequence

give_rank(seq):
takes in an encoded BWT sequence
returns a list of ranks calculated and a dictionary of letters and ranks

mapping(dict_of_things):
takes in a dictionary that consist of all the ranks
returns a dictionary of tuples that consists of items and their counts

inverse_bwt(args):
takes in args from user
runs give_rank
runs mapping
returns and writes a decoded sequence

read_fq(filename):
takes in a file of type .fq
returns an array of all reads

fullindexFM(sequence, match):
takes in a chromosome sequence and a sequence to match with
returns the amount of matches found and the beginning index

run_fullIndexFM(args):
takes args from user
writes number of matched sequnces and their indices to a file
the first array printed is READ1 matches
the second array printed is READ2 matches

third(arr, new_arr):
this function returns an array of every 3 in 10 reads

sixth(arr, new_arr):
this function returns an array of every 6 in 10 reads

run_HybridIndexFM(args):
takes in args from user
writes number of matches and indices to file
this function uses function third for accuracy

run_checkpointingIndexFM(args):
takes in args from user
writes number of matches and indices to file
this function uses function sixth for accuracy

###########################################################################################################################