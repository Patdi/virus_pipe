#!/usr/bin/env python2.7
import sys
import re
f = open('ICTV_and_NCBI_Tax-2017-07-23.txt', 'r')
i = open('ICTV_corrected.txt', 'r')
virus_data = dict()
virus_ictv = dict()

for line in f:
	if virus_data.has_key(line.strip().split(";")[0]) == None:
		virus_data[line.strip().split(";")[0]] = line.strip().split(";")[1:]
	virus_data[line.lower().strip().split(";")[1]] = line.strip().split(";")[1:]

for line in i:
	virus_ictv[line.strip().split(";")[0]] = line.strip().split(";")[0:]

for line in sys.stdin:
	tax_id = line.lower().strip().split("\t")[1]
	virus_name = line.lower().strip().split("\t")[2]
	
	if re.search(';', tax_id):
		tax_id_list = tax_id.split(";")
		
		if (virus_data.has_key(tax_id_list[0])):
			virus_str = str(virus_data[tax_id_list[0]])
			print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
		
		elif (virus_data.has_key(tax_id_list[1])):
			virus_str = str(virus_data[tax_id_list[1]])
			print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
		
		elif (virus_ictv.has_key(virus_name)):
			virus_str = str(virus_ictv[virus_name])
			print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
		
		else:
			print line.strip(), "\t", "\t", "\t", "\t", "\t" 
		
	elif (virus_data.has_key(tax_id)):
		virus_str = str(virus_data[tax_id])
		print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
	
	elif (virus_data.has_key(virus_name)):
		virus_str = str(virus_data[virus_name])
		print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
	
	elif (virus_ictv.has_key(virus_name)):
		virus_str = str(virus_ictv[virus_name])
		print line.strip(), "\t" ,virus_str.replace("'", "").replace("[", "").replace("]", "").replace("  ", " ").replace(",", "\t")
	
	else:
		print line.strip(), "\t", "\t", "\t", "\t", "\t" 