# A brand new simple program to get ORFs from Genome

# Uses Simple regex patterns and re.findall


import re , operator
from Bio import SeqIO, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Defining patterns

_start = r'ATG'
_stop = r'(?:TAG|TGA|TAA)'
_nonstop = r'(?:[CAG][TCAG]{2}|T(?:[TC][TCAG]|[AG][TC])|TGG)'
_codon = r'(?:[TCAG]{3})'
_orf_re = re.compile('(ATG'+ _nonstop + '*)\
('+ _stop +')',flags = re.I) # Edited for finding only orfs
_lead_re = re.compile(r'ATG',flags = re.I)


#define all the extractors
def extract_orfs(record,min_len , max_len, utr , dtr):
    orfs = []
    start_pos = []
    seqt = record.seq
    for a in _lead_re.finditer(str(seqt)):
        start_pos.append(a.start())
    for i in start_pos:
        a = _orf_re.search(str(seqt),i)
        if a is not None:
            s =(seqt[int(a.start()-utr):a.end()+dtr],"%s_%d_%d" \
                %(record.id,a.start(),a.end()))
            if len(s[0]) in range(int(min_len+utr+dtr),int(max_len+utr+dtr)):
                orfs.append(">%s(+)\n%s\n" %(s[1],s[0]))
    return (orfs)

def extract_orfs_rev(record,min_len , max_len, utr , dtr):
    orfs = []
    start_pos = []
    seqt = record.seq.reverse_complement()
    for a in _lead_re.finditer(str(seqt)):
        start_pos.append(a.start())
    for i in start_pos:
        a = _orf_re.search(str(seqt),i)
        if a is not None:
            s =(seqt[int(a.start()-utr):a.end()+dtr],"%s_%d_%d"\
                %(record.id,len(seqt)-a.start(),len(seqt)-a.end()))
            if len(s[0]) in range(int(min_len+utr+dtr),int(max_len+utr+dtr)):
                orfs.append(">%s(-)\n%s\n" %(s[1],s[0]))
    return (orfs)

"""
#make a new list
newlist = []

#ask if reverse is needed
print "ORFs from reverse sequence needed ? (yes:1, no: 0)"
rev = str(input())

# get orfs after reading the sequence

for records in SeqIO.parse(file_handle,"fasta"):
    for i in extract_orfs(records):
        newlist.append(i)
    
if rev is "1":
    for records in SeqIO.parse(file_handle,"fasta"):
        for i in extract_orfs_rev(records):
            newlist.append(i)

with open(raw_input("Output handle :"),"w") as f:
    for i in newlist:
        f.write(i)
"""
