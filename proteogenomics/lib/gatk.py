#!/home/ptgmat003/ve/bin/python

import subprocess
import tempfile
import shutil

def FARM(vcf_path, chrom, start, end, strand, genome_fasta): # 0 based coordinates
    outdir = tempfile.mkdtemp()
    out_fasta = outdir + '/out.fasta'
    inputIntervals = outdir +'/temp.intervals'
    w = open(inputIntervals,'w')
    intervals = '{}:{}-{}'.format(chrom, start + 1, end)
    w.write(intervals);w.close()
    command = 'java -jar $HOME/GATK/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R {} -o {} -L {} -V {}'.format(genome_fasta, out_fasta, inputIntervals, vcf_path)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    out, err = process.communicate()
    print(out, err)
    shutil.rmtree(outdir)
