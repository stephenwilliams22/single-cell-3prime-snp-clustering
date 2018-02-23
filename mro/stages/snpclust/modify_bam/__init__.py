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
stage MODIFY_BAM(
    in  path   reference_path,
    in  bam    input,
    in  sample input, # need to add this to the mro
    src py     "stages/snpclust/modify_bam",
    in  path    bed_file,
)split using (
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

    
    first_bam = martian.make_path('first_bam.bam')
    second_bam = martian.make_path('second_bam.bam')
    
    
    #add the readgroup needed for GATK
    rg_make_args = ['gatk-launch', 'AddOrReplaceReadGroups', '-I', args.input, '-O', 
                      first_bam, '-LB', 'lib1', '-PL', 'illumina','-PU', 'unit1', '-SM', args.sample]
    subprocess.check_call(rg_make_args)
      
    #this corrects the STAR mapq annotation and takes care of the split reads
    mapq_make_args = ['gatk-launch', 'SplitNCigarReads', '-R', genome_fasta_path, '-I', first_bam, '-O', 
                      second_bam, '--skip-mapping-quality-transform', 'false', 
                      '--create-output-bam-index', 'false']

    subprocess.check_call(mapq_make_args)            
    
    #remove the first .bam
    os.remove(first_bam)
    
    #join the bams together. NEED TO FIGURE THIS OUT
    def join(args, outs, chunk_defs, chunk_outs):
        outs.output = [chunk.output for chunk in chunk_outs]
