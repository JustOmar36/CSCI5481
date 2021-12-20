Omar Masri
Masri013
5571773

For this assignment I will have comments above the functions that describe functionality since there are too many of them.

How to run this program:

For this program I run all 3 problems in the file called "phySankoff.py"

to run any of the problems you input one of this three in the terminal

Problem 1:
"python3 phySankoff.py 1"

This specifies to the terminal to call phylogeny.py only. If the arg is one the script will generate a file called 
"tree_strucutre.txt" which contains a pairwise view of the tree and will draw an image of the Neighbor-Joining Tree and return the score

Problem 2:
"python3 phySankoff.py 2"

This will call sankoff.py on the toy sequence. If the arg is 2 then the script will ask you to input toy sequence file name.
I created a file called "toy_alignments.txt". Please put that filename in the terminal when asked. "toy_alignments.txt" only exists
in folder "Problem2". This will generate a txt file called "toy_result" with an array representation of the result 
and a descending representation of the array

Problem3:
"python3 phySankoff.py 3"

This will call sankoff.py on the 38 aligned sequences. There is nothing else you need to specficy here. it will result into a file called
"sankoff_result.txt" that prints out the score and then the sequences as nodes.

Notes:
I will provide a folder called "All_files" that will allow you to run everything in that folder if you do not wish to go
into each folder one at a time, but the critiera for this assignment will be in the individual folders.