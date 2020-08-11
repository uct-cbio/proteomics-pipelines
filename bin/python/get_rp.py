#!/usr/bin/env python

import Bio; from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import urllib2
import sys
import StringIO
import os

base = os.path.expanduser('~')
prot_folder = base + '/biotools/uniprot_proteomes/'
fasta_records = []

if len(sys.argv) == 1:
	accession = input('Enter UNIPROT proteome accession number: ')
	filename = accession

elif len(sys.argv) > 1:
	filename = '_'.join(sys.argv[1:])
	for arg in sys.argv[1:]:
		accession = arg
		present = False
		for prot_file in os.listdir(prot_folder):
			if  present == False:
				if prot_file == '{}.fasta'.format(accession):
					fasta = list(SeqIO.parse(prot_folder + prot_file,'fasta'))
					fasta_records.append(fasta)
					present = True
		if present == False:
			url = 'http://www.uniprot.org/uniprot/?query=proteome:{}&format=fasta'.format(accession)
			path = prot_folder
			f = urllib2.urlopen(url)
			page = StringIO.StringIO(f.read())
			f.close()
			prot = list(SeqIO.parse(page,'fasta'))
			SeqIO.write(prot,'{}/{}.fasta'.format(path, accession),'fasta')
			#SeqIO.write(embl,'{}/{}.fasta'.format(path,accession),'fasta')
			#fasta = SeqIO.read('{}/{}.fasta'.format(path,accession),'fasta')
			fasta_records.append(prot)

final = []
for i in fasta_records:
    final +=i

SeqIO.write(final, sys.stdout,'fasta')


 
