#!/usr/bin/python
'''
@author: Lin An & Difei Wang, 2013, Georgetown Univ. 
'''
import os
# Standard dependancies
from optparse               import OptionParser
from urllib                 import urlopen
from os.path                import isfile
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqIO              import write
from Bio.SeqIO              import parse
from Bio.SeqRecord          import SeqRecord
from Bio.Alphabet           import generic_protein
from Bio.PDB.PDBParser      import PDBParser
from Bio                    import AlignIO
from subprocess             import Popen
from subprocess             import PIPE
from sys                    import stderr
from os.path                import splitext
from os                     import chdir
from sets		    import Set
import collections

## bypass certi problem
## https://thomas-cokelaer.info/blog/2016/01/python-certificate-verified-failed/
##
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# ===Constants===

parser = PDBParser()
# Catch user input options
def parseArguments ():
    parser = OptionParser ( "Usage: %prog -p <pdbid> -c <chainid> -u <uniprotid>" )
    parser.add_option ( "-p", "--pdb", dest = "PDBid", help = "The PDB id, e.g. 1hso." )
    parser.add_option ( "-c", "--chain", dest = "Chainid", help = "The chain id, e.g A." )
    parser.add_option ( "-u", "--uniprot", dest = "Uniprotid", help = "The protein id, e.g P07327" )
    (options, arguments) = parser.parse_args ()
    if ( not options.PDBid ) or ( not options.Chainid ) or ( not options.Uniprotid ):
        parser.error ( "Requires a PDBID and a chain id and a protein id. Run with --help for help." )
        quit ( 1 )
    return options

# Table collection
class Table(object):
    def __init__(self):
        self.d = collections.defaultdict(dict)
    def add(self, row, col, val):
        self.d[row][col] = val
    def get(self, row, col, default=None):
        return self.d[row].get(col, default)
    def incol(self, col, val):
        for z,w in self.d.iteritems():
            if w.get(col)==val:
                return w
    def inrow(self, row):
        return self.d[row]
    def printtable (self):
        diff  = set()
        e     = self.d.keys()
        s     = sorted(self.d[0].keys())
        s1    = [h[1:] for h in s]
        structure = parser.get_structure(pdbid,pdbid+'.pdb')
        model = structure[0]
        chain = model[chainid]
        num   = 0

### set flags for missing residues etc 
        mut_flag = 'False' 

        mflag1 = 'False'
        mflag2 = 'False'

        prior1 = 'False'
        prior2 = 'False'
        next1  = 'False'
        next2  = 'False'

        for i in e:
            m1  = []  # pdb list
            m2  = []  # uniprt list
            mp1 = []  # prior in pdb
            mp2 = []  # prior in uniprot
            for j in s:
### record alignment result
		if i >1:
		   mp1.append(self.d[i-1][j][0])
		   mp2.append(self.d[i-1][j][1:])
                m1.append(self.d[i][j][0])
                m2.append(self.d[i][j][1:])
## pdb seq & uniprot seq
            if m1[0]==m1[1]:
                flag = '.'
            else:
                if (m1[0] == '-') or (m1[1] == '-'):
                    flag = '-' # missing,  flag ='-'
                else:          
                    flag = 'N' # mutation, flag ='N'
                    mut_flag = 'True'

## count missing aa in this alignment
## prior1 record TRUE if previous not - but current is -.
            if mp1:
		if mp1[0] != "-" and m1[0] == "-":
			prior1 = 'True'
            if mp2:
                if mp2[0] != "-" and m2[0] == "-":
                        prior2 = 'True'

## prior1 = TRUE, current not eq -
	    if prior1 == 'True' and m1[0] != '-':
		next1 = 'True'
###                print m
	    if prior1 == 'True' and next1 == 'True':
		mflag1 = 'True'

## prior2 = TRUE, current not eq -
            if prior2 == 'True' and m2[0] != '-':
                next2 = 'True'
###                print m
            if prior2 == 'True'  and next2 == 'True':
                mflag2 = 'True'
                
## calculate the delta between unipt pos and pdb pos
## skip if "-" 
            if m1[0][0] == '-':
                pdbnum = '-'
                delta = '-'
            else:
                pdbnum = chain.child_list[num].id[1]
                num += 1
                delta = str(int(m2[1])-pdbnum)
            print pdbid+'\t'+uniprotid+'\t'+chainid+'\t'+'\t'.join(m1)+'\t'+str(pdbnum)+'\t'+m2[1]+'\t'+flag+'\t'+delta
            diff.add(delta)

###  set store only uniq values
        diff_num=list(diff)

        if "-" in diff: diff_num.remove("-")
        print '== '+pdbid+'_'+chainid+'\t'+uniprotid+'\t'+'delta='+','.join(diff)+'\t'+str(len(diff_num))+'\t'+'\t'+'mut='+mut_flag+'\t'+'missing='+mflag1+','+mflag2 

# Retrive sequences from pdb and uniprot
def retriving(a,b,c):
    pdbid     = a
    chainid   = b
    uniid     = c
    my_record = []
    log       = open('pdb.fasta','w')
    seqpy     = Popen(["python","pdb_seq.py",pdbid],stdout=PIPE,stderr=PIPE)
    stdout    = seqpy.communicate()[0]
    log.write(stdout)
    wait      = seqpy.wait()
    log.close()
    seqfile   = open("pdb.fasta")
    for seq_record in parse(seqfile, "fasta"):
        r = seq_record.id.split('_')
        if r[0][-1]==chainid:
            my_record.append(seq_record)
    seqfile.close()
    url = 'https://www.uniprot.org/uniprot/'+uniid+'.fasta'
    seqfile2 = urlopen(url)
    for seq_record in parse(seqfile2, "fasta"):
        r = seq_record.id.split('|')
        uniprot = r[1]
        my_record.append(seq_record)
    seqfile2.close()
    write(my_record, "test.fasta", "fasta")
    

# Perform multiple sequence alignment
def alignseq():
#   clustalw_exe = r"/home/wangdi/apps/clustalw2"
    clustalw_exe = r"/Users/wangdi/apps/anaconda3/envs/py2/bin/clustalw"
    cline = ClustalwCommandline(clustalw_exe, infile="test.fasta", outorder="input")
    assert isfile(clustalw_exe), "Clustal W executable missing"
    stdout = cline()    
    alignfile=open("test.aln")
    align = AlignIO.read(alignfile, "clustal")
    alignfile.close()
    return align

# Reformat alignment results
def createtable(align):
    NumofSeq  = len(align)
    letterpos = Table()
    for i in range(0, NumofSeq):
        k=0
        for pos,base in enumerate(align[i].seq):
            row = pos
            col = str(i)+align[i].id
            val = base+str(k+1)
            k+=1
            letterpos.add(row,col,val)
    return letterpos

# main
if __name__ == "__main__":
    arguments = parseArguments ()
    pdbid     = arguments.PDBid
    chainid   = arguments.Chainid
    uniprotid = arguments.Uniprotid
    rawdata   = retriving(pdbid,chainid,uniprotid)
    align     = alignseq()
    datatable = createtable(align)
    datatable.printtable()
    os.remove('test.aln')
    os.remove('test.dnd')
    os.remove('test.fasta')
    os.remove('pdb.fasta')

