# A brand new simple program to get ORFs from Genome

# Uses Simple regex patterns and re.findall


import re , operator
from Bio import SeqIO,Seq,SeqRecord

# Defining patterns

_start = r'ATG'
_stop = r'(?:TAG|TGA|TAA)'
_nonstop = r'(?:[CAG][TCAG]{2}|T(?:[TC][TCAG]|[AG][TC])|TGG)'
_codon = r'(?:[TCAG]{3})'
_orf_re = re.compile('(ATG'+ _nonstop + '*)\
('+ _stop +')',flags = re.I) # Edited for finding only orfs
_lead_re = re.compile(r'ATG',flags = re.I)

file_handle = raw_input("Type input file directory :")  # Read in the input file

#define all the extractors
def extract_orfs(record):
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
            if len(s[0]) >= int(min_len+utr+dtr) and len(s[0]) <=int(max_len+utr+dtr):
                orfs.append(">%s(+)\n%s\n" %(s[1],s[0]))
    return (orfs)

def extract_orfs_rev(record):
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
            if len(s[0]) >= int(min_len+utr+dtr) and len(s[0]) <=int(max_len+utr+dtr):
                orfs.append(">%s(-)\n%s\n" %(s[1],s[0]))
    return (orfs)

# Ask for seqence min and max 
min_len = int(raw_input("Minimum lenght of Sequence"))
if min_len is None:
    min_len = 3
max_len = int(raw_input("Maximum length of Sequence"))
if max_len is None:
    max_len = 1000000

# ask for Utr and dtr to be concidered

utr = int(raw_input("UTR lenght of Sequence"))
if utr is None:
    utr = 3
dtr = int(raw_input("DTR length of Sequence"))
if dtr is None:
    dtr = 0

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
    
