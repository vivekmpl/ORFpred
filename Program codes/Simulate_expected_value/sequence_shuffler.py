# A program to reshuffle sequence in a random manner


import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# read file and randomize sequences
def shuffle_seq (file_name):
    shuffled = []
    for s in SeqIO.parse(file_name, "fasta"):
        shuf1 = (''.join(random.sample(s.seq,len(s.seq))))
        shuffled.append(SeqRecord(Seq(''.join(random.sample(shuf1,len(shuf1)))) ,\
                                  id = s.id , description = s.description))
    return shuffled
    #SeqIO.write(shuffled, "randomised_%s" % file_name,"fasta")
