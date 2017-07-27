#!/usr/bin/env python

import Bio; from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import urllib2
import sys
import StringIO
import os

base = os.path.expanduser('~')
genome_folder = base + '/biotools/ebi_embl_files/'
fasta_records = []

if len(sys.argv) == 1:
	accession = raw_input('Enter EBI genome accession number: ')
	filename = accession
elif len(sys.argv) > 1:
	filename = '_'.join(sys.argv[1:])
	for arg in sys.argv[1:]:
		accession = arg
		present = False
		for genome_file in os.listdir(genome_folder):
			if  present == False:
				if genome_file == '{}.fasta'.format(accession):
					fasta = SeqIO.read(genome_folder + genome_file,'fasta')
					fasta_records.append(fasta)
					present = True
		if present == False:
			url = "http://www.ebi.ac.uk/ena/data/view/{}&display=txt&expanded=true".format(accession)
			path = genome_folder
			f = urllib2.urlopen(url)
			page = StringIO.StringIO(f.read())
			f.close()
			embl = SeqIO.read(page,'embl')
			SeqIO.write(embl,'{}/{}.embl'.format(path, accession),'embl')
			SeqIO.write(embl,'{}/{}.fasta'.format(path,accession),'fasta')
			fasta = SeqIO.read('{}/{}.fasta'.format(path,accession),'fasta')
			fasta_records.append(fasta)

if len(fasta_records) == 1:
	fasta_records = fasta_records[0]
SeqIO.write(fasta_records, '{}.fasta'.format(filename),'fasta')


 
