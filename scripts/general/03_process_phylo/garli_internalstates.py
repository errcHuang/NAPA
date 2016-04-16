import itertools as it
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sequenceManipulation as seqManip
import argparse

def parseArguments():
	"""
	Builds the argument parser object
	"""
	parser = argparse.ArgumentParser(prog="garli_internalstates.py", 
		description="Translate cDNA internal states from GARLI output in protein AA sequence format", 
		usage="garli_internalstates.py -g garliFile -d cDNAfastaFile -c codonTableNum -p proteinFastaFile")
	parser.add_argument("-g", dest="garliFile",
		help="Name of Garli file containing inferred coding sequences at internal nodes", 
		required=True)
	parser.add_argument("-d", dest="cDNAfastaFile", 
		help="Where to output fasta formating of Garli internal nodes file",
		required=True)
	parser.add_argument("-c", dest="codonTableNum", 
		help="The number of the NCBI codon table to use.",
		required=True)
	parser.add_argument("-p", dest="proteinFastaFile",
		help="Where to output fasta of translated cDNA sequence",
		required=True)
	args = parser.parse_args()
	return args	

def translate_cDNA(cDNAstr, codonTableNum):
	cDNArecord = Seq(cDNAstr, IUPAC.unambiguous_dna)
	return str(cDNArecord.translate(table=codonTableNum,to_stop=True))
	
def translate_many_cDNA(cDNAfasta, codonTableNum):
	cDNAsequs = seqManip.get_FASTA_records(cDNAfasta)
	protSequs = {}
	for seqId in cDNAsequs:
		protSequs[seqId] = translate_cDNA(cDNAsequs[seqId], codonTableNum)
	return protSequs
	
def garliInternal_to_fasta(garliFile,fastaFile):
	cDNAfastaRecs = {}
	with open(garliFile, 'rb') as f:
		for key,group in it.groupby(f, lambda line: "node" in line):
        	if not key:
            	groupList = list(group)
            	cDNAstr = ""
            	for line in groupList:
            		if "site" in line: continue
            		cDNAstr += line.strip().split('\t')[-1].strip()
            internalNodeNum = key.strip().split('\t')[0]
            cDNAfastaRecs[">"+internalNodeNum] = cDNAstr
    simple_print_FASTA(cDNAfastaRecs,fastaFile)
    
def main():
	args = parseArguments()
	garliInternal_to_fasta(args.garliFile, args.cDNAfastaFile)
	protSequs = translate_many_cDNA(args.cDNAfastaFile, args.codonTableNum)
	simple_print_FASTA(protSequs, args.proteinFastaFile)