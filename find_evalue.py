#looks through rscape output to find PCCRs with product evalue less than treshhold
#writes pair (id - significant evalue) to file
#prints id with description 

from operator import mul
from functools import reduce
from sys import argv
import numpy as np

if argv[2]=="inf":
    treshold = np.inf
else:
    treshhold = float(argv[2])
pref = ["aa","ab","ac","ad"] #big rscape output was splitted into files

pccrs = 0
out = open(argv[1], "w") #for example, ./results/id_eval.tsv
for filename in ["./results/" + i + "_rscape.out" for i in pref]:
    with open(filename,"r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            header = line.rstrip()
            if header[:3] == "id=":
                pccrs += 1
                good_lines = []
                while True:
                    line = file.readline().rstrip()
                    if line == "":
                        break
                    elif line[0] == "#":
                        continue
                    elif line[0] == "*":
                        good_lines.append(line)
                evalues = [float(i.split()[4]) for i in good_lines]
                if len(evalues)!=0:
                    eval_prod = reduce(mul,evalues,1)
                    if eval_prod<treshold:
                        print("#" + header)
                        print("#" + " ".join(str(evalues)) + ": " + str(eval_prod))
                        for l in good_lines:
                            print(l)
                        out.write(header.split()[0] + "\t" + str(eval_prod) + "\n")
print(f"{pccrs} analysed in total")
out.close()