from sankoff import run_sankoff
from  phylogeny import run_phy_file
import sys

if __name__ == "__main__":
    args = sys.argv
    #input
    #grabs 3 arguments
    if len(args) == 2:
        if (args[1].isnumeric()):
            if(int(args[1])) == 1:
                print("running problem 1")
                run_phy_file(args)
            elif (int(args[1])) == 2:
                print("running problem 2")
                run_sankoff(args)
            elif (int(args[1])) == 3:
                print("running problem 3")
                run_sankoff(args)
            else:
                print("Only acceptable numbers are 1, 2 or 3")
        else:
            print("Second arguement should be an integer")
    else:
        print('The number of arguments is invalid. This file required 1 argument.')