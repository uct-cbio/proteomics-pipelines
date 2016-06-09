#!/usr/bin/env python

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import tempfile
import shutil
import urllib2

def pblast(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/blast/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastp_cline = NcbiblastpCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out)
    stdout, stderr = blastp_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record

def pblast_pretty(seqstr, db, threshold = 0.04):
    rec = pblast(seqstr, db)
    output_list = []
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = [
                'Alignment E value: {}'.format(hsp.expect), 
                'length: {}'.format(alignment.length),
                'Feature: {}'.format(alignment.title)]
                output = ', '.join(output)
                output_list.append(output) 
    for i in output_list:
        return i
        break

def nblast(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/cdna_blast/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastn_cline = NcbiblastnCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out)
    stdout, stderr = blastn_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record


def tblastn(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/tblastn/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastn_cline = NcbiblastnCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out,
                                        strand = plus)
    stdout, stderr = blastn_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record

def nblast_pretty(seqstr, db, threshold = 0.04):
    rec = nblast(seqstr, db)
    output_list = []
    count =1
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = [
                'Alignment E value: {}'.format(hsp.expect), 
                'length: {}'.format(alignment.length),
                'Feature: {}'.format(alignment.title)]
                output = ', '.join(output)
                output_list.append(output)
                count +=1
    return output_list

def gi2up(gi_accession_number):
    try:
        url = 'http://www.uniprot.org/mapping/?from=P_GI&to=ACC&query={}'.format(gi_accession_number)
        p = urllib2.urlopen(url).read()
        splitter = 'xml:lang="en"><head><title>'
        ls = p.split(splitter)[1].split(' ')[0].split(':')[1]
        url2 = 'http://www.uniprot.org/uniprot/?query=yourlist:{}&sort=yourlist:{}&columns=yourlist%28{}%29,id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength%2Cexistence%2Ccomment%28PATHWAY%29%2Cgo%2Cgo%28biological%20process%29%2Cgo%28molecular%20function%29%2Cgo%28cellular%20component%29%2Cgo-id'.format(ls,ls,ls)
        x = urllib2.urlopen(url2).read()
        datum = x.split('class="addRemoveColumn mid"')[0].split('</script></td></tr></thead><tbody><tr ')[1]
        datum = datum.split('class=')
        biorec = {}
        biorec['Entry (Uniprot)'] = datum[5].split('"entryID"><a href="/uniprot/')[1].split('"')[0]
        biorec['Entry Name (Uniprot)'] = datum[6].split('>')[1].split('<')[0]
        biorec['Protein Names'] = datum[11].split('title="')[1].split('"')[0]
        biorec['Gene Names'] = ''.join(''.join('>'.join(datum[14].split('>')[1:]).split('<strong>')).split('</strong>')).split('</div>')[0].strip() 
        biorec['Organism'] = ''.join(datum[15].split('">')[2:]).split('</a><')[0]
        biorec['Length'] = datum[16].split('>')[1].split('<')[0]
        biorec['Protein Existence'] = datum[17].split('>')[1].split('<')[0] 
        biorec['Pathway'] = datum[18].split('td>')[1].split('<td')[0] 
        go = datum[19].split('<td style=""')[0].split('</a>')
        biorec['Gene Ontology (GO)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[20].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (biological process)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        #return entryID1, entryID2, protein_name, gene_names, organisms
        go = datum[21].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (molecular function)'] =  [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[22].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (cellular component)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[23].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology IDs'] = [i.split('>')[-1] for i in go if ((len(i.split('>')[-1]) > 0) and (i.split('>')[-1] != '<td ')) ] 
        return biorec
    except:
        return None

def bparse(rec, threshold = 0.04):
    output_list = []
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = ['****Alignment****\n',
                    'Sequence: {}'.format(alignment.title),
                    'Length: {}'.format(alignment.length),
                    'E value: {}'.format(hsp.expect),
                    'HSP Query:   {}'.format(hsp.query),
                    'HSP Match:   {}'.format(hsp.match),
                    'HSP Subject: {}'.format(hsp.sbjct)]
                output = '\n'.join(output)
                output_list.append(output)
    return output_list

