#!/usr/bin/env bash


# Needs Python 2.7.8

input_fasta=$1
outdir=$2
ve=$3

#module load compilers/perl-5.22.1  # Sometimes there are problems with PERL and compilers. Got  FILE::PATH error or something like that.

source $ve/bin/activate
interproscan.sh -i $input_fasta -d $outdir -goterms -iprlookup -pa


