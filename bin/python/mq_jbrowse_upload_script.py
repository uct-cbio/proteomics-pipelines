#!/usr/bin/env python3

import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3
import subprocess
import pandas as pd

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
output = os.path.abspath(sys.argv[2])

outpath = os.path.abspath(config['outdir'])


cmd = 'rm -rf {}/jbrowse/local && jbrowse create {}/jbrowse/local'.format(outpath, outpath)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()

def upload_reference_genome(ref_path, name, outpath):
    cmd = "samtools faidx {}".format(ref_path)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()

    cmd = 'jbrowse add-assembly {} --out {}/jbrowse/local --load copy --force'.format(ref_path, outpath)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    
    with open('{}/jbrowse/local/config.json'.format(outpath)) as f:
        js = json.loads(f.read())
    for assembly in  js['assemblies']:
        if assembly['name'] == ref_path.split('.fasta')[0].split('/')[-1]:
            assembly['name'] = name
    with open('{}/jbrowse/local/config.json'.format(outpath), 'w') as w:
        w.write(json.dumps(js))

def upload_reference_gff3(ref_path,name,outpath):
    with open(ref_path) as f:
        headers=f.readlines()[:5]
    d = pd.read_csv(ref_path,sep='\t',skiprows=5,header=None)
    ref_folder ='{}/jbrowse/{}'.format(outpath, name)
    if not os.path.exists(ref_folder):
        os.mkdir(ref_folder)
    #top=d[:1]#.to_csv(sep='\t',header=None,index=None)
    #d = d[d[2] == 'CDS']
    def get_id(x):
        _ = str(x).split(';')
        for i in _:
            if i.startswith('ID='):
                return i.split('ID=')[1]
        return 'None'

    assert 9 not in d.columns
    
    d[9] = d[8].apply(get_id)
    
    parent_map = {}
    for row in d.iterrows():
        parent_map[row[1][9]] = row
    
    formatting = {'gene' :{'color': "rgba(61,180,94,1)", 'height':10, 'label':'locus_tag'}, 'CDS' :{'color': "rgba(134,218,32,1)", 'height':10, 'label':'locus_tag|gene'}, 'ncRNA' :{'color': "rgba(218,53,32,1)", 'height':10, 'label':'gene'}}
    seen_parents = []
    for feature, group in d.groupby(2):
        #group = pd.concat([top,group])
        newlist=[]
        newlist+=headers
        #newlist.append(top)
        #group=group.to_csv(sep='\t',header=None,index=None)
        for row in group.iterrows():
            meta = row[1][8]
            rows = [row]
            if 'Parent=' in meta:
                parent=meta.split('Parent=')[1].split(';')[0]
                if not parent in seen_parents:
                    parent_row = parent_map[parent]
                    rows = [parent_row,row]
                    seen_parents.append(parent)
            for row in rows:
                entry =[]
                for i in range(9):
                    entry.append(str(row[1][i]).strip())
                entry= '\t'.join(entry) +'\n'
                newlist.append(entry)


        #newlist.append(group)
        newlist=''.join(newlist)
        
        outfile ='{}/{}_{}.gff3'.format(ref_folder, name, feature)
        with open(outfile,'w') as w:
            w.write(newlist)
        height = 5
        color ="rgba(32,218,166,1)"
        Label=None
        if feature in formatting:
            Label = formatting[feature]['label']
            color = formatting[feature]['color']
            height = formatting[feature]['height']
        upload_gff3(outfile, outpath,feature, height, color, True , name, Label)

def update_config(outpath, name, filepath, height, color, showLabels,Label=None):
    with open('{}/jbrowse/local/config.json'.format(outpath)) as f:
        js = json.loads(f.read())
        displays = [
        {
          "type": "LinearBasicDisplay",
          "displayId": "{}-LinearBasicDisplay".format(filepath),
          "renderer": {
            "type": "SvgFeatureRenderer",
            "color1": "rgba(32,218,166,1)",
            "height": 5,
            "showLabels": False

          }
        },
        {
          "type": "LinearArcDisplay",
          "displayId": "{}-LinearArcDisplay".format(filepath)
        }
      ]
    for track in js['tracks']:
        if filepath == track['trackId']:
            track['name'] = name
            if not 'displays' in track:
                track['displays'] = displays
            for display in track['displays']:
                if display['type'] == "LinearBasicDisplay":
                    renderer = display['renderer']
                    renderer['height'] = height
                    renderer["color1"] =  color
                    if not Label is None:
                        Label = Label.split('|')[::-1]
                        base = "get(feature,'name') || get(feature,'id') || get(feature,'ID')"
                        for l in Label:
                            base= "get(feature, '{}')".format(l) + ' || ' + base
                        renderer['labels'] =  {"name": "jexl:{}".format(base)}
                    renderer["showLabels"] = showLabels

    with open('{}/jbrowse/local/config.json'.format(outpath), 'w') as w:
        w.write(json.dumps(js))
                

def minimap2(ref,query, ref_name, query_name, outpath):
    cmd='minimap2 -cx asm5 {} {} > {}_vs_{}.paf'.format(ref,query,query_name,ref_name)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    
    cmd='jbrowse add-track {}_vs_{}.paf --assemblyNames {},{} --out {}/jbrowse/local --load copy --force'.format(query_name, ref_name, query_name, ref_name, outpath)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()


def upload_gff3(gff3, outpath, name, height, color, showLabels, reference,Label=None):
    cmd = 'gt gff3 -sortlines -tidy -retainids {} > {}.sorted.gff3'.format(gff3, gff3)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    
    gff3_zip_path='{}.sorted.gff3.gz'.format(gff3)
    if os.path.exists(gff3_zip_path):
        os.remove(gff3_zip_path)

    cmd = 'bgzip {}.sorted.gff3'.format(gff3)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()

    cmd = 'tabix -f {}'.format(gff3_zip_path)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()

    cmd = "jbrowse add-track {}.sorted.gff3.gz --out {}/jbrowse/local --load copy --assemblyNames {}".format(gff3, outpath, reference)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    filepath = '{}.sorted.gff3'.format(gff3)
    filepath = filepath.split('/')[-1]
    update_config(outpath, name, filepath, height, color, showLabels,Label)

for reference in config['reference']:
    break
    assembly_id = config['reference'][reference]['assembly_id']
    genome_path = output + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
    genome_gff3 = output + '/ena/{}/{}.gff3'.format(assembly_id, assembly_id)
    upload_reference_genome(genome_path, reference , outpath)
    #upload_reference_gff3(genome_gff3, reference , outpath)
    color = "jexl:cast({CDS:'rgba(134,218,32,1)',biological_region:'rgba(32,218,166,1)',databank_entry:'rgba(218,32,87,1)', gene:'yellow',ncRNA:'rgba(218,53,32,1)',RNA:'rgba(32,218,166,1)',repeat_region:'rgba(14,155,120,1)',sequence_feature:'rgba(32,111,218,1)',tRNA:'rgba(218,32,205,1)',transcript:'rgba(243,4,128,1)'})[get(feature,'type')]"
    width=5
    width = "jexl:cast({CDS:10,biological_region:3,databank_entry:1, gene:5,ncRNA:3,RNA:3,repeat_region:3,sequence_feature:3,tRNA:3,transcript:7.5})[get(feature,'type')]"
    labels='protein_id|locus_tag|gene|gene_synonym'
    upload_gff3(genome_gff3, outpath,reference,width, color, True , reference, labels)
    
for query in config['reference']:
    break
    assembly_id = config['reference'][reference]['assembly_id']
    query_genome_path = output + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
    for target in config['reference']:
        if not target== query:
            assembly_id2 = config['reference'][target]['assembly_id']
            target_genome_path = output + '/ena/{}/{}.fasta'.format(assembly_id2, assembly_id2)
            minimap2(query_genome_path, target_genome_path, query, target, outpath)
    

for strain in config['strains']:
    assembly = config['strains'][strain]['assembly']
    upload_reference_genome(assembly, strain , outpath)
    feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_features.gff3".format(strain)
    color = "jexl:cast({open_reading_frame:'green',polypeptide:'purple'})[get(feature,'type')]"
    width = "jexl:cast({open_reading_frame:10,polypeptide:5})[get(feature,'type')]"

    showLabels = False
    upload_gff3(feature_gff3, outpath, strain + ' features', width, color, showLabels, strain )
    continue
    for reference in config['reference']:
        assembly_id = config['reference'][reference]['assembly_id']
        peptide_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_peptides.gff3".format(strain, assembly_id)
        upload_gff3(peptide_gff3, outpath, strain + ' peptides', 5, "rgba(32,218,166,1)", False , reference)
        
        orf_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_orfs.gff3".format(strain,assembly_id)
        upload_gff3(orf_gff3, outpath, strain + ' ORFs', 10, "goldenrod", False, reference )
        
        assembly_id = config['reference'][reference]['assembly_id']
        ref_fasta = outpath + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
        minimap2(ref_fasta, assembly, reference, strain, outpath)

# strain to strain syntenty
for strain in config['strains']:
    continue
    assembly = config['strains'][strain]['assembly']
    for tstrain in config['strains']:
        if not strain == tstrain:
            tassembly = config['strains'][tstrain]['assembly']
            minimap2(assembly, tassembly, strain, tstrain, outpath)

cmd = "jbrowse text-index --perTrack --out {}/jbrowse/local".format(outpath)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()




