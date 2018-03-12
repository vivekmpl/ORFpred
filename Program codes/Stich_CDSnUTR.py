#Program to stich UTR

import sys, getopt
from Bio import SeqIO, SeqRecord
from collections import defaultdict

def stich(cds,utr,dtr):
    resultant_seq = []
    cds = list(SeqIO.parse(cds,"fasta"))
    utr = SeqIO.to_dict(SeqIO.parse(utr,"fasta"))
    dtr = SeqIO.to_dict(SeqIO.parse(dtr,"fasta"))
    for i in cds:
        resultant_seq.append(SeqRecord.SeqRecord(seq = utr[i.id]+i.seq+dtr[i.id],\
                                                 id = i.id, description = i.description))
        return resultant_seq

def print_something(a,b):
    return "I am printing this - %s -- %s" % (a,b)

def main():
    d=defaultdict(list)
    for k, v in ((k.lstrip('-'), v) for k,v in (a.split('=') for a in sys.argv[1:])):
        d[k].append(v)
    print print_something(d["a"][0],d["b"][0])

def usage():
    print "Please give a correct input"

if __name__ == "__main__":
    main()
else:
    usage()

        

