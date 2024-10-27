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

groups = defaultdict(list)
design = pd.read_csv(config['design'])
seen_groups = []
for group in config['group_levels']:
    temp_design = design.drop_duplicates(['Strain', group])
    for row in temp_design.iterrows():
        strain_group = row[1][group]
        strain = row[1]['Strain']
        assert not strain_group in groups[strain]
        groups[strain].append(strain_group)

# {'Non-enzymatic N-terminus non-TSS peptide', 'Putative TSS peptide at
# enzymatic cleavage site', 'MAP-cleaved iMet', 'Enzymatic N-terminus non-TSS
# peptide', 'Non-ATG iMet', 'Non-enzymatic iMet', 'Ambiguous: previous codon
# and first codon (ATG) are start codons', 'N-terminal acetylated
# iMet','Globally non-specific'

def upload_reference_genome(ref_path, name, outpath):
    cmd = "samtools faidx {}".format(ref_path)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()

    cmd = 'jbrowse add-assembly {} -n {} --out {}/jbrowse/local --load copy --force'.format(ref_path, name,  outpath)
    print(cmd)
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

def update_config(outpath, name, filepath, strokeWidth, height, color, showLabels,Label=None, type_filter=None, outline=None):
    with open('{}/jbrowse/local/config.json'.format(outpath)) as f:
        js = json.loads(f.read())
        js['adminKey'] = config['JBROWSEadminKey']
        defaultSession={
              "name": "New Session",
              "view": {
                      "id": "synteny-view",
                      "type": "SyntenyView",
                      "tracks": []
                    }
        }
        js["defaultSession"] = defaultSession
        displays = [
        {
          "type": "LinearBasicDisplay",
          "displayId": "{}-LinearBasicDisplay".format(filepath),
          "renderer": {
            "type": "SvgFeatureRenderer",
            "color1": "rgba(32,218,166,1)",
            "height": 5,
            "showLabels": True

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
            if not type_filter is None:
                  track["filter"] = {"jexl": "get(feature, 'type') == '{}'".format(type_filter)}
            if not 'displays' in track:
                track['displays'] = displays
            for display in track['displays']:
                if display['type'] == "LinearBasicDisplay":
                    renderer = display['renderer']
                    if not outline is None:
                        renderer['outline'] = outline
                    renderer['height'] = height
                    renderer["color1"] =  color
                    renderer['strokeWidth'] = strokeWidth
                    #renderer["strokeWidth"]="function(f) { return f.get('nterm_acetylated') === 'true' ? 2 : 1; }"
                    #renderer["shape"]="function(f) { return f.get('nterm_acetylated') === 'true' ? 'rect' : 'ellipse'; }"  
                    if not Label is None:
                        Label = Label.split('|')[::-1]
                        base = "get(feature,'name') || get(feature,'id') || get(feature,'ID')"
                        for l in Label:
                            base= "get(feature, '{}')".format(l) + ' || ' + base
                        renderer['labels'] =  {"name": "jexl:{}".format(base)}
                    renderer["showLabels"] = showLabels
                    renderer["showDescriptions"] = False
        meta = {'type' : track['type'], "configuration":track['trackId']}
        js["defaultSession"]['view']['tracks'].append(meta)
    with open('{}/jbrowse/local/config.json'.format(outpath), 'w') as w:
        w.write(json.dumps(js))
                

def minimap2(ref,query, ref_name, query_name, outpath):
    cmd='minimap2 -cx asm5 {} {} > {}_vs_{}.paf'.format(ref,query,query_name,ref_name)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    
    cmd='jbrowse add-track {}_vs_{}.paf --assemblyNames {},{} --out {}/jbrowse/local --load copy'.format(query_name, ref_name, query_name, ref_name,  outpath)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()


def upload_gff3(gff3, outpath, name, strokeWidth, height, color, showLabels, reference,Label=None, type_filter=None, outline=None):
    
    gff3_name = '.'.join(gff3.split('.')[:-1])
    print(gff3_name)

    cmd = 'gt gff3 -sortlines -tidy -retainids {} > {}_temp.gff3'.format(gff3, gff3_name)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)
    assert process.returncode == 0, f"Command failed with return code {process.returncode}: {stderr}"

    #cmd='sort -k1,1 -k4,4n -k5,5n {}_temp.gff3 > {}_sorted.gff3'.format(gff3_name, gff3_name)
    cmd='sort -k1,1 -k4,4n -k5,5n {}.gff3 > {}_sorted.gff3'.format(gff3_name, gff3_name)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)
    assert process.returncode == 0, f"Command failed with return code {process.returncode}: {stderr}"
    

    gff3_zip_path='{}_sorted.gff3.bgz'.format(gff3_name)
    if os.path.exists(gff3_zip_path):
        os.remove(gff3_zip_path)

    cmd = 'bgzip -c {}_sorted.gff3 > {}_sorted.gff3.bgz'.format(gff3_name, gff3_name)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)
    assert process.returncode == 0, f"Command failed with return code {process.returncode}: {stderr}"

    cmd = 'tabix -f {}'.format(gff3_zip_path)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)
    assert process.returncode == 0, f"Command failed with return code {process.returncode}: {stderr}"

    cmd = "jbrowse add-track {}_sorted.gff3.bgz --out {}/jbrowse/local --load copy --assemblyNames {}".format(gff3_name, outpath, reference)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    stdout, stderr = process.communicate()
    print("STDOUT:", stdout)
    print("STDERR:", stderr)
    assert process.returncode == 0, f"Command failed with return code {process.returncode}: {stderr}"
    filepath = '{}_sorted.gff3'.format(gff3_name)
    filepath = filepath.split('/')[-1]
    update_config(outpath, name, filepath, strokeWidth, height, color, showLabels,Label, type_filter, outline)

for reference in config['reference']:
    assembly_id = config['reference'][reference]['assembly_id']
    genome_path = output + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
    genome_gff3 = output + '/ena/{}/{}.gff3'.format(assembly_id, assembly_id)
    upload_reference_genome(genome_path, reference , outpath)
    #upload_reference_gff3(genome_gff3, reference , outpath)
    color = "jexl:cast({CDS:'rgba(134,218,32,1)',biological_region:'rgba(32,218,166,1)',databank_entry:'rgba(218,32,87,1)', gene:'yellow',ncRNA:'rgba(218,53,32,1)',RNA:'rgba(32,218,166,1)',repeat_region:'rgba(14,155,120,1)',sequence_feature:'rgba(32,111,218,1)',tRNA:'rgba(218,32,205,1)',transcript:'rgba(243,4,128,1)'})[get(feature,'type')]"
    width=5
    width = "jexl:cast({CDS:10,biological_region:3,databank_entry:1, gene:5,ncRNA:3,RNA:3,repeat_region:3,sequence_feature:3,tRNA:3,transcript:7.5})[get(feature,'type')]"
    strokeWidth = 0
    labels='protein_id|locus_tag|gene|gene_synonym'
    upload_gff3(genome_gff3, outpath,reference, strokeWidth, width, color, True , reference, labels)
    
for query in config['reference']:
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
    width = "jexl:cast({open_reading_frame:10,polypeptide:5,domain:4, proteoform:3, proteoform_orf:3})[get(feature,'type')]"
    color = "jexl:(get(feature, 'type') == 'polypeptide' ? (get(feature, 'start_type') == 'non-enzymatic n-terminus non-tss peptide' ? 'purple' : get(feature, 'start_type') == 'putative tss peptide at enzymatic cleavage site' ? 'orange' : get(feature, 'start_type') == 'map-cleaved imet' ? 'yellow' : get(feature, 'start_type') == 'enzymatic n-terminus non-tss peptide' ? 'pink' : get(feature, 'start_type') == 'non-atg imet' ? 'brown' : get(feature, 'start_type') == 'non-enzymatic imet' ? 'blue' : get(feature, 'start_type') == 'ambiguous: previous codon and first codon (atg) are start codons' ? 'red' : get(feature, 'start_type') == 'n-terminal acetylated imet' ? 'green' : get(feature, 'start_type') == 'globally non-specific' ? 'grey' : 'black') : get(feature, 'type') == 'open_reading_frame' ? 'green' : get(feature, 'type') == 'domain' ? (get(feature, 'excluded') == 'true' ? 'red' : 'cyan' ) : get(feature, 'type') == 'proteoform' ? 'magenta' : 'black')"

    outline="jexl:(get(feature, 'nterm_acetylated') == 'true' ? 'black' : get(feature,'type') == 'proteoform_orf' ? 'magenta' : 'none')"
    showLabels=True
    
    upload_gff3(feature_gff3, outpath, strain + ' features', strokeWidth, width, color, showLabels, strain , outline=outline)
    feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_orfs.gff3".format(strain)
    upload_gff3(feature_gff3, outpath, strain + ' ORFs', strokeWidth, width, color, showLabels, strain )
    feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_peptides.gff3".format(strain)
    upload_gff3(feature_gff3, outpath, strain + ' peptides', strokeWidth, width, color, showLabels, strain , type_filter='polypeptide', outline=outline)
    
    feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_domains.gff3".format(strain)
    upload_gff3(feature_gff3, outpath, strain + ' domains', strokeWidth, width, color, showLabels, strain )
    feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_proteoforms.gff3".format(strain)
    upload_gff3(feature_gff3, outpath, strain + ' proteoforms', strokeWidth, width, color, False, strain, outline=outline )
    
    for group in groups[strain]:
        #continue # remove
        feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_domains.gff3".format(group)
        upload_gff3(feature_gff3, outpath, group + ' domains', strokeWidth, width, color, showLabels, strain )
        feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_features.gff3".format(group)
        upload_gff3(feature_gff3, outpath, group + ' features', strokeWidth, width, color, showLabels, strain, outline=outline)
        feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_peptides.gff3".format(group)
        upload_gff3(feature_gff3, outpath, group + ' peptides', strokeWidth, width, color, showLabels, strain, outline=outline )
        feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_proteoforms.gff3".format(group)
        upload_gff3(feature_gff3, outpath, group + ' proteoforms', strokeWidth, width, color, False, strain, outline=outline )

    #with open(output + '/annotations/global_group_identified.json') as f: 
    for reference in config['reference']:
        assembly_id = config['reference'][reference]['assembly_id']
        peptide_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_peptides.gff3".format(strain, assembly_id)
        upload_gff3(peptide_gff3, outpath, strain + ' peptides',strokeWidth, 5, "rgba(32,218,166,1)", True , reference)
        
        orf_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_orfs.gff3".format(strain,assembly_id)
        print(orf_gff3) 
        upload_gff3(orf_gff3, outpath, strain + ' ORFs', strokeWidth, 10, 'goldenrod', True, reference )
        
        proteoform_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_proteoforms.gff3".format(strain, assembly_id)
        upload_gff3(proteoform_gff3, outpath, strain + ' proteoforms', strokeWidth, width, color, True , reference, outline=outline)
        
        feature_gff3="{}/jbrowse/".format(outpath)+strain+"/{}_{}_features.gff3".format(strain, assembly_id)
        upload_gff3(feature_gff3, outpath, strain + ' features', strokeWidth, width, color, True , reference, outline=outline)
        
        assembly_id = config['reference'][reference]['assembly_id']
        ref_fasta = outpath + '/ena/{}/{}.fasta'.format(assembly_id, assembly_id)
        
        minimap2(ref_fasta, assembly, reference, strain, outpath)

# strain to strain syntenty
for strain in config['strains']:
    assembly = config['strains'][strain]['assembly']
    for tstrain in config['strains']:
        if not strain == tstrain:
            tassembly = config['strains'][tstrain]['assembly']
            minimap2(assembly, tassembly, strain, tstrain, outpath)

cmd = "jbrowse text-index --perTrack --out {}/jbrowse/local".format(outpath)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()




