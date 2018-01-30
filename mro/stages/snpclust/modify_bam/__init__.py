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
    src py     "stages/snpclust/modify_bam",
)
'''

def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    #dic_make_args = ['gatk-launch', 'CreateSequenceDictionary', '-R', genome_fasta_path]
    #subprocess.check_call(dic_make_args)
    
    first_bam = martian.make_path('output.RG.bam')
    second_bam = martian.make_path('output.RG.STARcor.bam')
    rg_make_args = ['gatk-launch', 'AddOrReplaceReadGroups', '-I', args.input, '-O', 
                      first_bam, '-LB', 'lib1', '-PL', 'illumina','-PU', 'unit1', '-SM', 'sample']
    subprocess.check_call(rg_make_args)
      
      
    #this corrects the STAR mapq annotation. Uses 8 threads.
    samtools_args = '''samtools view -@ 8 -h {}
      awk 'BEGIN{{OFS="\t"}} $5 == 255 {{ $5 = 60; print; next}} {{print}}' | 
      samtools view -Sb -@ 8 - > {}'''.format(first_bam, second_bam)
    
    subprocess.call(samtools_args, shell=True)            
    
    os.remove(first_bam)
      
    samtools_index_args = ['samtools', 'index','-@','8',second_bam]
    subprocess.call(samtools_index_args)
