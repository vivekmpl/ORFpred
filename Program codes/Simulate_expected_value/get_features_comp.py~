# A program to identify sequence features and write to arff file

from __future__ import division
import re, itertools , Bio , csv
from Bio import SeqIO
from Bio import SeqUtils
import math

# ask the UTR lenght 

utr = 350 #int(raw_input("lenght of UTR regions"))
dtr = 100 #int(raw_input("lenght of DTR regions"))

# make list of combinations and append in it

combi = []
combi1 = []
combi3 = []
combi6 = []

# combination of nucleotides

# n
for p in itertools.product("ATGC",repeat = 1):
    combi.append("".join(p))
    combi1.append("".join(p))
# nnn
for p in itertools.product("ATGC",repeat = 3):
    combi.append("".join(p))
    combi3.append("".join(p))
# nnnnnn
for p in itertools.product("ATGC",repeat = 6):
    combi.append("".join(p))
    combi6.append("".join(p))

combi = combi #+["gc13" , "gc12"] 

# all sequence patterns required :
seq_pat1 = re.compile(r"[ATGC]{1}")
seq_pat2 = re.compile(r"[ATGC]{2}")
seq_pat3 = re.compile(r"[ATGC]{3}")
seq_pat4 = re.compile(r"[ATGC]{4}")
seq_pat5 = re.compile(r"[ATGC]{5}")
seq_pat6 = re.compile(r"[ATGC]{6}")

# defining seq features isolation procedure

def isf_predictor(seq):
    x_len = len(seq[utr+3:-(dtr+3)])
    f = []
    lis1 = []
    lis3 = []
    lis6 = []
    for j in seq_pat1.findall(seq[utr+3:-(dtr+3)]):
            lis1.append(j)
    for j in seq_pat3.findall(seq[utr+3:-(dtr+3)]):
            lis3.append(j)
    for j in seq_pat6.findall(seq[utr+3:-(dtr+3)]):
            lis6.append(j)
    for i in combi1:
            f.append(float(lis1.count(i)/(x_len)))
    for i in combi3:
            f.append(float(lis3.count(i)/(x_len/3)))
    for i in combi6:
            f.append(float(lis6.count(i)/(x_len/6)))
    return [str(i) for i in f]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print "--"*5+"get_fet_comp module ready"+"--"*5
# culminate everything
def get_fet_predictor(sequences):
	return isf_predictor(str(sequences.seq))
def give_names():
	return combi
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# end of the program .
