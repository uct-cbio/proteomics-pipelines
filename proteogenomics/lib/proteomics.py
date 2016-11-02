#!/usr/bin/env python

import pandas as pd
import Bio; from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def tryptic_digest(protein , min_len = 7, max_len = 30):
    digest = []
    peptide = []
    
    for aa in range(len(protein)):
       
       current = protein[aa]
       peptide.append(current)
       
       if aa < len(protein) -1:
            nxt = protein[aa + 1]
            if ((current =='R') or (current=='K')) and (nxt != 'P'):
                digest.append(''.join(peptide))
                peptide = []

       else:
           digest.append(''.join(peptide))
           peptide = []
    
    filt_digest = []
    for i in digest:
        if (len(i) >= min_len) and (len(i) <= max_len):
            filt_digest.append(i)

    return filt_digest

def emPAI(sc, protein):
    peps = tryptic_digest(protein)
    print(peps)
    empai = sc / float(len(peps))
    return empai
    

