#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.constants as tk_constants
import cellranger.utils as cr_utils
import os

__MRO__ = '''
stage CALL_SNPS(
    in  path   reference_path,
    in  bam    input,
    in  int    n_donors,
    out vcf[]  output,
    src py     "stages/snpclust/call_snps_pd",
    in  path    bed_file,
) split using (
    in  string locus,
)
'''

# split the .bed file and make chunks
def split(args):
    loci = [x.split() for x in open(args.bed_file)]
    chunks = [{'locus': locus, '__mem_gb': 8} for locus in loci]
    return {'chunks': chunks}

# define the reference 
def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
        
    chrom, start, stop = args.locus
    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\n")
    
# Run GATK4    
    gatk_args = ['gatk-launch', 'HaplotypeCaller', 
                 '-R', genome_fasta_path, 
                 '-I', args.input, 
                 '-O','output.vcf', 
                 '-L', bed_path,  
                 '--minimum-mapping-quality', '30', 
                 '--min-base-quality-score', '20', 
                 '--dont-use-soft-clipped-bases', 'true', 
                 '--add-output-vcf-command-line', 'false']
            
    subprocess.check_call(gatk_args)
        
def join(args, outs, chunk_defs, chunk_outs):
    outs.output = [chunk.output for chunk in chunk_outs]
    
