#/usr/bin/env python

import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate
from Bio.Data import CodonTable
from numpy import NaN
import Bio; from Bio import SeqIO
import pandas as pd
import os
import natsort
import subprocess
import tempfile
import collections; from collections import defaultdict
from io import StringIO

def alt_starts_recs(records, starts = ['ATG','GTG','TTG']):
    all= []
    for seq in records:
        alt_seqs = alt_tss(seq, starts = starts)
        all += alt_seqs
    return all

def Three_Frames(sequence): # All three frames on the forward strand (returns three sequences and includes possible stop codons)
    sequence_length = len(sequence)
    Frame1 = sequence[0:int(sequence_length/3)*3]
    Frame2 = sequence[1:int((sequence_length+2)/3)*3-2]
    Frame3 = sequence[2:int((sequence_length+1)/3)*3-1]
    Frames = (Frame1, Frame2, Frame3)
    return Frames

def Six_Frames_and_Reverse_Complement(sequence):# Do three frame translation on forward and reverse-complement of forward strand
    seq = Seq(sequence)
    reversed_seq = seq.reverse_complement()
    rev_seq = str(reversed_seq).rstrip()
    return Three_Frames(sequence) + Three_Frames(rev_seq), rev_seq

def Stop_List(table): # Return possible stop codons given translation table
    Codon_Table = table
    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    Stop_Codons = codon_table.stop_codons
    return Stop_Codons

def Start_List(table): # Possible starts
    Codon_Table = table
    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    Start_Codons = codon_table.start_codons
    return Start_Codons

def codon_list(table):
    codons = CodonTable.unambiguous_dna_by_id[table].forward_table.keys()
    return codons

def six_frame(Genome, table, peptide_length, codons = ['ATG','GTG','TTG'], translated = True):
    Minimum_Peptide_Length = peptide_length
    Codon_Table = table
    Start_Codons = codons

    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    Stop_Codons = codon_table.stop_codons
    #for i in Start_Codons:
    #        assert i in codon_table.start_codons                       # make sure the selected start codons are found in the selected translation table

    Six_Frame, Reverse_Complement_Genome = Six_Frames_and_Reverse_Complement(Genome)

    ORFs_by_Frame = {}

    for i in Six_Frame:
        count = 1
        Frame = Six_Frame.index(i)+1
        [a] = [Frame - 3 if Frame > 3 else Frame] # Adjust for frameshift to find index.
        b = a-1         # Shift the value of y to == index of last nucleotide (-1), and adjust for frame shift with 'a'.

        ORFs_by_Frame[Frame] = []
        starts = False
        index = 0
        ORF = []
        Start = NaN

        if Frame <= 3:
            String = 1

        elif Frame >3:
            String = -1

        while index < (len(i))-2:
            x = index
            y = index + 3
            if starts == True:
                if i[x:y] in Stop_Codons:
                    ORF.append(i[x:y])
                    temp = ''.join(ORF)
                    if Frame <= 3:
                        location_start = Start + a
                        location_end = y + b
                    elif Frame > 3:
                        location_start = len(Genome)-(y + b)+1 # Index by direct strand
                        location_end = len(Genome)-(Start+a)+1  # Index by direct strand

                    ORFs_by_Frame[Frame].append((count,Frame, String, location_start, location_end, temp))
                    count += 1
                    ORF = []
                    starts = False
                elif x == (len(i) - 3):
                    ORF.append(i[x:y])
                    if (len(Genome)-(y+b))<=0:
                        assert len(Genome)-(y+b) == 0
                    temp = ''.join(ORF)
                    if Frame <= 3:
                        location_start = Start + a
                        location_end = y + b
                    elif Frame > 3:
                        location_start = len(Genome)-(y + b)+1 # Index by direct strand
                        location_end = len(Genome)-(Start+a)+1  # Index by direct strand
                    ORFs_by_Frame[Frame].append((count,Frame, String,location_start, location_end, temp))
                    count += 1
                    ORF = []
                    starts = False
                else:
                    ORF.append(i[x:y])
            if starts == False:
                if i[x:y] in Start_Codons:
                    starts = True
                    Start = x
                    ORF.append(i[x:y])
            index += 3
    Six_Frame = []
    if translated == True:
        for i in ORFs_by_Frame:                          #lets check six frame ORF db
            for j in ORFs_by_Frame[i]:                   #check by frame
                if j[5][-3:] in Stop_Codons:             #check if strand ends in a stop codon
                    Seq_len = (len(j[5])-3)/3                  #subtract 3 if above true, and divide by three to get length
                elif j[5][-3:] not in Stop_Codons:       #if last codon is not a stop codon, then just divide by 3 to get petide length
                    Seq_len = len(j[5])/3
                if Seq_len >= Minimum_Peptide_Length:    #compare to minimum peptide length
                    translated = str(translate(Seq(j[5])))           #translate into peptide
                    if translated[-1] == '*':            #if end in a stop codon, then change first letter to M
                        prot_seq = translated[0:len(translated)-1]    #checked
                    elif translated[-1] != '*':
                        prot_seq = translated
                    if j[1] < 4:
                        frame_direction = '+'
                    elif j[1] > 3:
                        frame_direction = '-'
                    record = SeqRecord(Seq(prot_seq),id = "({}){}:{}".format(frame_direction,j[3],j[4]), description = "Six Frame Translated ORF")
                    Six_Frame.append(record)

    elif translated == False:
        for i in ORFs_by_Frame:
            for j in ORFs_by_Frame[i]:                   #check by frame
                if j[5][-3:] in Stop_Codons:             #check if strand ends in a stop codon
                    Seq_len = (len(j[5])-3)/3                  #subtract 3 if above true, and divide by three to get length
                elif j[5][-3:] not in Stop_Codons:       #if last codon is not a stop codon, then just divide by 3 to get petide length
                    Seq_len = len(j[5])/3
                if Seq_len >= Minimum_Peptide_Length:    #compare to minimum peptide length
                    seq = str(Seq(j[5]))
                    if j[1] < 4:
                        frame_direction = '+'
                    elif j[1] > 3:
                        frame_direction = '-'
                    record = SeqRecord(Seq(seq),id = "({}){}:{}".format(frame_direction,j[3],j[4]), description = "Six Frame ORF")
                    Six_Frame.append(record)
    return Six_Frame

def get_frame(seq, genome_seq, stops): #'''BUG'''translated seqrecord with coords, genomic seq record, possible stops eg ['TAA', 'TAG', 'TGA']
    '''In put a translated seq_record with description with format eg. (+)4454545:4454548 and a single genomic seq record to output the stop-to-stop nucleic acid sequence that contains that seq_record. NB input a list of possible stop codons, eg. stops = ['TAA', 'TAG', 'TGA'] (for bacterial table 11)'''
    frames = []
    coords = seq.description.split()[1].split('|')
    for i in coords:
        coord = i[3:].split(':')
        start = int(coord[0])
        end = int(coord[1])
        strand = i[:3]
        if strand == '(-)':
            while str(genome_seq.seq[end:end+3].reverse_complement()) not in stops:
                end = end + 3
            end = end + 3
            seq_ = genome_seq.seq[start-1:end].reverse_complement()
            description = '(-){}:{}'.format(start,end)
        elif strand == '(+)':
            while str(genome_seq.seq[start-4:start-1]) not in stops:
                start = start - 3
            start = start - 3
            seq_ = genome_seq.seq[start-1:end]
            description = '(+){}:{}'.format(start, end)
        record = SeqRecord( seq_, id = seq.id+' Genomic sequence option (stop to stop) at {}'.format(description))
        frames.append(record)
    return frames

def gmhmmp(model, path, sequence):   # --mode (native, combined, heuristic), --org (path), seq_record object
    temp_dir = tempfile.mkdtemp()
    temp_out = temp_dir + '/out'
    temp_faa = temp_dir +'/faa'
    temp_fnn = temp_dir + '/fnn'
    f = temp_dir+'/fasta'                 #temporary file to write the fasta record to
    # SeqIO.write(sequence, f, 'fasta')
    command = 'gmhmmp -m {}{} {} -r -a -s d -n -g 11 -A {}gmhmmp_proteins.fasta'.format(path, model, f, temp_dir)
    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    process.wait()
    fasta = list(SeqIO.parse(temp_dir +'gmhmmp_proteins.fasta', 'fasta'))
    return fasta                                                      # protein sequence of predicted genes

def gmsn(genome, outpath = '/tmp'):  #prokaryotes only, requires a seq-record object for query sequence
    if outpath[-1] =='/':
        outpath = outpath[:-1]
    temp_dir = tempfile.mkdtemp()
    f = temp_dir+'/fasta'                 #temporary file to write the fasta record to
    SeqIO.write(genome, f, 'fasta')
    command = 'cd {} && gmsn.pl --prok --faa {} && rm -rf {}/GeneMarkS_output_files && mv {} {}/GeneMarkS_output_files'.format(temp_dir, f, outpath,temp_dir, outpath)
    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    process.wait()
    fasta = list(SeqIO.parse('{}/{}/fasta.faa'.format(outpath, 'GeneMarkS_output_files'), 'fasta'))
    return fasta                                                      # protein sequence of predicted genes

def m_start_list(proteins):     #inport a list of amino acid sequences, and return all starting with M
    m = []
    for protein in proteins:
        seq = 'M'+str(protein.seq)[1:]
        record = SeqRecord(seq = Seq(seq), id = protein.id, description = protein.description)
        m.append(record)
    return m

def sf_contigs(contigs, table=11, peptide_length = 20, codons = 'All', translated = False): # contigs fasta file, output genomic six frame sequences
    if codons == 'All':
        codons = codon_list(table)
    six_frame_seqs = []
    seq_count = 1
    for i in contigs:
        cid = i.id    #get the contig id
        six = six_frame(str(i.seq), table = table, peptide_length = peptide_length, codons = codons, translated=translated)
        for j in six:
            id = cid +'|' + j.id
            sequence =  j.seq   #sequence seq
            record = SeqRecord(seq=sequence, description ='Six_frame_sequence_{}'.format(str(float(seq_count))), id =id)
            six_frame_seqs.append(record)
            seq_count += 1

    #sf_dct = defaultdict(list)       #the following lines combine repeating sequences in the database into one record mapped to different locations in the 					contigs
    #for i in six_frame_seqs:
    #	string = str(i.seq)
    #	sf_dct[string].append(i.id)
    
    #new_dct = defaultdict()
    #for i in sf_dct:
    #	new_id = '|'.join(sf_dct[i])
	#seq = Seq(i)
	#record = SeqRecord(seq=seq,id=new_id,description='')
	#new_dct[new_id] =record
    #keys = new_dct.keys()
    #keys = natsort.natsorted(keys)
    #combined_seqs = []
    #count = 1

    #for key in keys:
    #    record = new_dct[key]
    #    record.description = 'Six_Frame_ORF_{}'.format(float(count))
    #    count += 1
    #    combined_seqs.append(record)
    return six_frame_seqs


class gff3:
    '''Basic class to parse GFF3'''
    gffcols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    def __init__(self, GFF3):
        self.table = pd.read_csv(os.path.abspath(GFF3),sep='\t',comment='#')
        self.table.columns = self.gffcols

    def entryCount(self):
        return len(self.table)

class vcf:
    '''Basic class for vcf data'''
    def __init__(self, VCF):
        path = os.path.abspath(VCF)
        w = open(path).readlines()
        headers = [line.strip() for line in w if line.startswith('#CHROM')][0].split('\t')
        self.table = pd.read_csv(path,sep='\t',comment='#', header=None)
        self.table.columns = headers

    def entryCount(self):
        return len(self.table)

class genome:
    '''Basic class of genome object'''
    def __init__(self, FASTA):
        self.reflist = list(SeqIO.parse(os.path.abspath(FASTA),'fasta'))    
        self.refdict = {}
        for seq in self.reflist:
            self.refdict[seq.id] = seq
        

class variant:
    def __init__(self, VCF, FASTA, GFF3, name=None):  # vcf file, name, genome,gff3
        self.ref_genome = genome(FASTA)
        self.gff = gff3(GFF3)
        self.vcf = vcf(VCF)

        print(self.vcf.table.head())

        self.name = name
    
        self.var_genome = {} 
        self.var_mapping = {}

        for chrom in self.ref_genome.refdict:      
            
            start_pos = 0
            step = 0
            seen_pos = []
            rec = self.ref_genome.refdict[chrom]  
            rseq = str(rec.seq)

            _ = self.vcf.table
        
            new_seq_list = []
            new_seq_map  = {}
            
            subset = _[_['#CHROM']==chrom]  # get only vcf rows relevant to the fasta id
            for row in subset.iterrows():
                pos = row[1]['POS'] - 1
                assert pos not in seen_pos
                seen_pos.append(pos)
                ref = row[1]['REF']
                alt = row[1]['ALT']

                for p in range(start_pos, pos):
                    new_seq_map[p] = p + step
                    new_seq_list.append(rseq[p])

                new_seq_map[pos] = pos + step
                new_seq_list.append(alt)
                step += len(alt) - len(ref) 
                start_pos = pos + len(ref) 
            
            for p in range(start_pos, len(rseq)):
                new_seq_map[p] = p + step
                new_seq_list.append(rseq[p])

            vseq = ''.join(new_seq_list)
            self.var_genome[chrom] = vseq
            self.var_mapping[chrom]= new_seq_map

        





def get_end(cigar,start): # subtract 1 from start
    flags={'D':1,'I':0,'M':1,'S':0}
    end = start 
    num =[]
    for char in cigar:
        if char not in flags:
            num.append(char)
        elif char in flags:
            new_num = int(''.join(num))*flags[char]
            num=[]
            end += new_num
    return end

def sam_parse(sam):
    count = 0
    sequences = []
    orfs = defaultdict(list)
    seqdict = defaultdict(list)
    sam = sam.split('\n')
    if sam[0].startswith('@'):
        sam = sam[1:]
        #reference_title = sam[0].split('\tSN:')[1].split('\tLN:')[0].split('.')[0]+'_aligned.fasta'
        #length = sam[0].split('\tLN:')[1]

    for aln in sam:
        if len(aln.split('\t')) > 1:
            datum =  aln.split('\t')
            q_name = datum[0]
            sam_flag = datum[1]
            ref_name = datum[2]
            start = datum[3]
            
            mapq = datum[4]
            cigar = datum[5]
            
            end = get_end(cigar, int(start)-1)

            seq = Seq(datum[9])
            
            if sam_flag =='16':
                seq = seq.reverse_complement()
                strand = '-'
            else:
                strand = '+'
            
            header = '({}){}:{}'.format(strand, start, end)
            mapping = ref_name+'_'+start+'_'+sam_flag+'_'+cigar+'_'+mapq+'_'+header
            seqdict[q_name].append(str(seq))

            if sam_flag == '4':
                record = SeqRecord(seq=seq,id=q_name,description='Unmapped')
                sequences.append(record)
            else:
                record = SeqRecord(seq=seq,id=q_name,description=mapping)
                sequences.append(record)

                #orfs[q_name].append(record)
    #new_seqs = []
    #for i in orfs:
    #    description = []
    #    id = i
    #    assert len(list(set(seqdict[i])))==1
    #    seq = Seq(seqdict[i][0])
    #    for j in orfs[i]:
    #        description.append(j.description)
    #        new_seqs.append(j)   
    #    record = SeqRecord(seq=seq, description=' '.join(description),id=id)
    #   sequences.append(record)
    return sequences


'''
def alignment_code(query, query_start, alignment, ref_prot):
    query_len = len(query[query_start:]) # query len is the length of the translated region
    query_translated = str(translate(Seq(query[query_start:])))
    assert ref_prot.endswith('*')#to account for stop codon in algorithm, ensure that the function recieves the translated sequence from BioPython
    ref_prot_len = len(ref_prot)*3 #nucleotide length
    if query == alignment:
        if query_len == ref_prot_len:
            code=5        # when whole sequence match and translation matches reference proteome (after first amino acid)
        else:
            code=1     #whole sequence match but translation does not excactly match reference proteome   
    elif query[-query_len:]==alignment[-query_len:]:
        if query_len == ref_prot_len:
            code=6 #translated sequence nucleotide sequence match to genomeand matches reference protein
        else:
            code=2    # translated sequence nucleotide sequence match  but not to reference protein
    elif query[-ref_prot_len:]==alignment[-ref_prot_len:]:
        code=3   # only regions matching are length of ref prot nuc sequence
    else:
        code=4 # mapps to a region coding for a known protein but not exact match
    return code
'''
'''
def sam_cigar(description):
    new_desc = description.split('_')
    ref = new_desc[0]
    start = int(new_desc[1])
    flag = new_desc[2]
    if flag == '0':
        strand = '+'
    elif flag == '16':
        strand ='-'
    cigar=new_desc[3]
    datum=(ref,strand,start,cigar)
    return datum
'''
'''
def sf_map(six_frame, genome,ref_prot, starts=['ATG','GTG','TTG'], table=11):
    stops=Stop_List(table)
    count=0
    sf_count=0
    y=set()
    for i in six_frame:
        query = str(i.seq)
        if query[-3:] in stops:
            x = 0
            found = False
            prot_query = ''
            while x<= len(query) and found == False:
                if query[x:x+3] in starts:
                    prot_query = translate(Seq(query[x:]),
                            table=11,
                            cds=False)         
                    found = True
                else:
                    x+=3
            seq_id = i.id
            seq_description = i.description.split()[1:] 
            for j in seq_description:
                prots_found = []
                datum = sam_cigar(j)
                strand=datum[1]
                start =int(datum[2])-1
                cigar = datum[3]
                end = get_end(cigar,start)
                if strand =='+':
                    aligned = str(genome[start:end].seq)
                elif strand =='-':
                    aligned = str(genome[start:end].seq.reverse_complement())
                template =  Seq(aligned).translate(table=11, cds=False)
                if aligned[-3:] in stops:
                    prot_found = False
                    for k in ref_prot:
                        if prot_found == False:
                            ref_record = k
                            prot = str(k.seq)+'*'
                            query_prot = prot[1:]
                            if query_prot in template:
                                prot_start = template.find(query_prot)
                                code = alignment_code(query, 
                                        x ,  
                                        aligned,
                                        prot)
                                prot_found = True
                                if strand == '+':
                                    orf_start = start+x
                                    orf_end = end
                                    prot_id = k.id
                                    prot_start=template.find(str(k.seq)[1:])*3-3    
                                    prot_header = '{}|{}:{}'.format(prot_id,
                                            start+ prot_start +1, 
                                            orf_end)
                                    print prot_header.split('|')[-1]
                                    count += 1
                                elif strand=='-':
                                        pass
                    if prot_found == False:
                        code = '*'
                        prot_header = None
                else:
                    pass
                if prot_found == False:
                    header = '({}){}:{}|{}|{}|*'.format(strand, start+1, end, cigar,code)
                else:
                    header = '({}){}:{}|{}|{}|{}'.format(strand, start+1, end, cigar,code,prot_header)
                if code == 4 and strand == '+':
                    print 
                    print header
                    print 'Query:\n', query
                    print 'Alignment:\n',aligned
                    print str(aligned)==str(query)
                    print len(query), len(aligned)
                    print [i for i in xrange(len(query)) if query[i] != aligned[i]]
                    #print query[494], aligned[494]

                    print '\nTranslated Query:\n', prot_query
                    print 
                    print ref_record.format('fasta')
                    print '########################################'
        elif query[-3:] not in stops:
             pass
        sf_count+=1 
      
'''


'''
def annotate_gm(gm, sf):#not used, but can be used to annotate a gm prediction
    annotated = []
    count = 0
    found = []
    for i in gm:
        lst = i.description.split('\t>')[0].split('|')
        done = False
        strand = lst[3]
        start = lst[4]
        end = lst[5]
        description = '({}){}:{}'.format(strand, start, end)
        annos = []
        count += 1
        while done == False:
            for j in sf:
                if done == False and j not in found:
                    desc = j.description.split()[1].split('|')
                    for coord in desc:
                        j_strand = j.description.split()[1].split(')')[0][1:]
                        j_start  = j.description.split()[1].split(')')[1].split(':')[0]
                        j_end = j.description.split()[1].split(')')[1].split(':')[1]
                        if strand == j_strand:
                            if strand == '-':
                                if start == j_start:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
                            elif strand == '+':
                                if end == j_end:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
            done = True
        annotations = description+' '+'|'.join(annos)
        record = SeqRecord(i.seq, id = i.id, description = annotations)
        annotated.append(record)
        print count
return annotated
'''

'''
def m_start_dict(proteins):
	m_dict = defaultdict()
	for prot in proteins:
		assert not isinstance(proteins[prot], (list, tuple))
		record =proteins[prot]
		if prot[:1] !='M':
			id = record.id
			seq = 'M' + str(record.seq)[1:]
			description = record.description
			s_record = SeqRecord(id = id, seq = Seq(seq), description = description)
			if seq not in proteins:
				m_dict[seq] = s_record   #discard elements already known to exist
		elif prot[:1]=='M':
			s_record = SeqRecord(id = record.id, seq = Seq(str(record.seq)), description = record.description)
			m_dict[prot] =s_record     #add known reference proteins starting with 'M'
        for i in m_dict:
		rec = m_dict[i]
		rec.id = rec.id.split('|')[0]+ '|' + 'M_START_'+ rec.id.split('|')[1]
	return m_dict
'''
'''
def alt_ref(ref, seq):     # refs annotated by cds location, as dict, and a SeqRecord object (fasta, embl etc). Returns alternative start amino acid
    fasta = seq.format('fasta')
    genome = SeqIO.read(StringIO.StringIO(fasta),'fasta')
    M_refs = defaultdict(list)
    M_refs_id = defaultdict(list)
    alt_ref = []
    for i in ref:
        for j in ref[i]:
            if len(j.description.split(' Coordinates: ')) > 1:
                coords =  ref[i][0].description.split(' Coordinates: ')[1]
                for coord in coords.split('|'):
                    try:
                        start = int(coord.split(')')[1].split(':')[0])-1
                        end = int(coord.split(')')[1].split(':')[1])
                        strand = coord.split(')')[0][1:]
                        nucs = genome.seq[start:end]

                        if strand == '+':
                            prot =  str(translate(nucs, table = 11, cds = False))
                            assert prot[-1] == '*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)

                        elif strand == '-':
                            prot = str(translate(nucs.reverse_complement(), table= 11, cds = False))
                            assert prot[-1] =='*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)
                    except:
                        pass
    for i in M_refs_id:
        for j in M_refs_id[i]:
            j.id = j.id +' Alternative Start'
            j.description = j.description.split(' Coordinates: ')[0] +' Non-ATG start at coordinates: '+'|'.join(M_refs[i])
            j.seq = Seq(i)
            alt_ref.append(j)
            break
    return alt_ref
'''

'''
#5. GET NOVEL SEQUENCES
def get_novel(proteins):    # input a list of annotated and non annotated sequences, and output only the non-annotated ones
    novel = []
    for i in proteins:
        if len(i.description.split()) == 1:
            novel.append(i)
    return novel

def ldict(seqs):
	seqs_dict = defaultdict(list)
	for i in seqs:
		seqs_dict[str(i.seq)]
'''
'''
def annotate_gm(gm, sf):#not used, but can be used to annotate a gm prediction
    annotated = []
    count = 0
    found = []
    for i in gm:
        lst = i.description.split('\t>')[0].split('|')
        done = False
        strand = lst[3]
        start = lst[4]
        end = lst[5]
        description = '({}){}:{}'.format(strand, start, end)
        annos = []
        count += 1
        while done == False:
            for j in sf:
                if done == False and j not in found:
                    desc = j.description.split()[1].split('|')
                    for coord in desc:
                        j_strand = j.description.split()[1].split(')')[0][1:]
                        j_start  = j.description.split()[1].split(')')[1].split(':')[0]
                        j_end = j.description.split()[1].split(')')[1].split(':')[1]
                        if strand == j_strand:
                            if strand == '-':
                                if start == j_start:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
                            elif strand == '+':
                                if end == j_end:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
            done = True
        annotations = description+' '+'|'.join(annos)
        record = SeqRecord(i.seq, id = i.id, description = annotations)
        annotated.append(record)
        print count
return annotated
'''

'''
def m_start_dict(proteins):
	m_dict = defaultdict()
	for prot in proteins:
		assert not isinstance(proteins[prot], (list, tuple))
		record =proteins[prot]
		if prot[:1] !='M':
			id = record.id
			seq = 'M' + str(record.seq)[1:]
			description = record.description
			s_record = SeqRecord(id = id, seq = Seq(seq), description = description)
			if seq not in proteins:
				m_dict[seq] = s_record   #discard elements already known to exist
		elif prot[:1]=='M':
			s_record = SeqRecord(id = record.id, seq = Seq(str(record.seq)), description = record.description)
			m_dict[prot] =s_record     #add known reference proteins starting with 'M'
        for i in m_dict:
		rec = m_dict[i]
		rec.id = rec.id.split('|')[0]+ '|' + 'M_START_'+ rec.id.split('|')[1]
	return m_dict
'''
'''
def alt_ref(ref, seq):     # refs annotated by cds location, as dict, and a SeqRecord object (fasta, embl etc). Returns alternative start amino acid
    fasta = seq.format('fasta')
    genome = SeqIO.read(StringIO.StringIO(fasta),'fasta')
    M_refs = defaultdict(list)
    M_refs_id = defaultdict(list)
    alt_ref = []
    for i in ref:
        for j in ref[i]:
            if len(j.description.split(' Coordinates: ')) > 1:
                coords =  ref[i][0].description.split(' Coordinates: ')[1]
                for coord in coords.split('|'):
                    try:
                        start = int(coord.split(')')[1].split(':')[0])-1
                        end = int(coord.split(')')[1].split(':')[1])
                        strand = coord.split(')')[0][1:]
                        nucs = genome.seq[start:end]

                        if strand == '+':
                            prot =  str(translate(nucs, table = 11, cds = False))
                            assert prot[-1] == '*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)

                        elif strand == '-':
                            prot = str(translate(nucs.reverse_complement(), table= 11, cds = False))
                            assert prot[-1] =='*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)
                    except:
                        pass
    for i in M_refs_id:
        for j in M_refs_id[i]:
            j.id = j.id +' Alternative Start'
            j.description = j.description.split(' Coordinates: ')[0] +' Non-ATG start at coordinates: '+'|'.join(M_refs[i])
            j.seq = Seq(i)
            alt_ref.append(j)
            break
    return alt_ref
'''

'''
#5. GET NOVEL SEQUENCES
def get_novel(proteins):    # input a list of annotated and non annotated sequences, and output only the non-annotated ones
    novel = []
    for i in proteins:
        if len(i.description.split()) == 1:
            novel.append(i)
    return novel

def ldict(seqs):
	seqs_dict = defaultdict(list)
	for i in seqs:
		seqs_dict[str(i.seq)]
'''


