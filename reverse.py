#!/usr/bin/env python2.7
#takes in sequences in tab file format, reverse compliments and prints out sequence name \n then sequence
import sys 
class Gene:
	def __init__(self, creationid, creationseq):
		self.id = creationid
		self.seq = creationseq
		
	def reverse_complement(self):
		return(self.seq[::-1].replace("A", "1").replace("T", "2").replace("G", "3").replace("C", "4").replace("1", "T").replace("2", "A").replace("3", "C").replace("4", "G"))
		
for line in sys.stdin:
	seq_list = line.strip().split("\t")
	seq_id = seq_list[0]
	seq_nu = seq_list[1]
	new_seq = Gene(seq_id, seq_nu)
	print seq_id
	print new_seq.reverse_complement() 
	
