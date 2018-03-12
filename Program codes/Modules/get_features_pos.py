# A program to identify sequence features and write to arff file

from __future__ import division
import re, itertools , Bio , csv
from Bio import SeqIO
from Bio import SeqUtils

# ask the UTR lenght 

utr = 350 #int(raw_input("lenght of UTR regions"))
dtr = 100 #int(raw_input("lenght of DTR regions"))

# make list of combinations and append in it

combi = []
combi5prime = []
combi3prime = []

# combination of nucleotides


for i in range(1,utr+1):
    combi5prime.append("5p%s"%i)

for i in range(1,dtr+1):
    combi3prime.append("3p%s"%i)

combi = combi  +  combi5prime[::-1] + combi3prime +["Stop" , "Base after ATG"] 

# defining seq features isolation procedure

def isf_predictor(seq):
    f = []
    for i in seq[:utr] : # 5 prime
        f.append(i)
    for i in seq[-dtr:] : #3prime
        f.append(i)
    f.append(seq[-(dtr+3):-dtr]) # Stop
    f.append(seq[utr+3])# Base after ATG
    return [str(i) for i in f]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print "--"*5+"get_fet_pos module ready"+"--"*5
# culminate everything
def get_fet_predictor(sequences):
	return isf_predictor(str(sequences.seq))
def give_names():
	return combi
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# end of the program .
