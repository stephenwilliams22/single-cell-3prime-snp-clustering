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
    
    #add readgroups to the bams
    first_tumor_bam = martian.make_path('output.tumor.RG.bam')
    second_tumor_bam = martian.make_path('output.tumor.RG.STARcor.bam')
    rg_make_tumor_args = ['gatk-launch', 'AddOrReplaceReadGroups', '-I', args.tumor.input, '-O', 
                      first_tumor_bam, '-LB', 'lib1', '-PL', 'illumina','-PU', 'unit1', '-SM', 'sample']
    subprocess.check_call(rg_make_tumor_args)
    
    
    first_normal_bam = martian.make_path('output.normal.RG.bam')
    second_normal_bam = martian.make_path('output.normal.RG.STARcor.bam')
    rg_make_normal_args = ['gatk-launch', 'AddOrReplaceReadGroups', '-I', args.input, '-O', 
                      first_normal_bam, '-LB', 'lib1', '-PL', 'illumina','-PU', 'unit1', '-SM', 'sample']
    subprocess.check_call(rg_make_normal_args)
      
      
    #this corrects the STAR mapq annotation. Uses 8 threads.
    samtools_tumor_args = '''samtools view -@ 8 -h {} |
      awk 'BEGIN{{OFS="\t"}} $5 == 255 {{ $5 = 60; print; next}} {{print}}' | 
      samtools view -Sb -@ 8 - > {}'''.format(first_tumor_bam, second_tumor_bam)
    
    subprocess.check_call(samtools_tumor_args, shell=True)
    
    #this corrects the STAR mapq annotation. Uses 8 threads.
    samtools_normal_args = '''samtools view -@ 8 -h {} |
      awk 'BEGIN{{OFS="\t"}} $5 == 255 {{ $5 = 60; print; next}} {{print}}' | 
      samtools view -Sb -@ 8 - > {}'''.format(first_normal_bam, second_normal_bam)
    
    subprocess.check_call(samtools_normal_args, shell=True)   
    
    #clean up
    os.remove(first_tumor_bam)
    os.remove(first_normal_bam)

    #index  
    samtools_index_tumor_args = ['samtools', 'index',second_tumor_bam]
    subprocess.check_call(samtools_index_tumor_args)
    
    #index  
    samtools_index_normal_args = ['samtools', 'index',second_normal_bam]
    subprocess.check_call(samtools_index_normal_args)

