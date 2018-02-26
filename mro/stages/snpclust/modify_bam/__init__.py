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
    chunks = [{'locus': locus, '__mem_gb': 16, '__threads': 1} for locus in loci]
    return {'chunks': chunks}

# define the reference 
def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)
        
    chrom, start, stop = args.locus
    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\n")

    
    # Correct the STAR mapping from 255 to 60 and take care of split reads
    output_bam = martian.make_path('output.bam')
    star_args = ['gatk-launch', 'SplitNCigarReads',
                 '-R', genome_fasta_path,
                 '-I', args.input,
                 '-L', bed_path,
                 '-O', output_bam,
                 '--max-reads-in-memory', '50000',
                 '--skip-mapping-quality-transform', 'false',
                 #'--create-output-bam-index', 'true',
                 '--TMP_DIR', os.getcwd()]
    
    #star_args = ['java', 
    #             '-Djava.io.tmpdir=/mnt/home/stephen/yard',
    #             '-jar', '/mnt/opt/gatk/3.8/GenomeAnalysisTK.jar', 
    #             '-T', 'SplitNCigarReads', 
    #             '-R', genome_fasta_path, 
    #             '-I', args.input, 
    #             '-o', output_bam,
    #             '-L', bed_path,
    #             '-rf', 'ReassignOneMappingQuality', 
    #             '-RMQF', '255', 
    #             '-RMQT', '60', 
    #             '-U', 'ALLOW_N_CIGAR_READS']
                 
    subprocess.check_call(star_args)
    
    #sort and index the bam
    args = ['samtools', 'sort', 'output.bam']
    subprocess.check_call(args)
    #tk_bam.sort(output_bam)
    os.remove(outs.output)
    os.rename('output_sorted.bam', 'output.bam')
    tk_bam.index('output.bam')
    
    #join the bams together
def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    tk_bam.concatenate(outs.output, input_bams)
    #tk_bam.sort(outs.output)
    #os.remove(outs.output)
    #os.rename('output_sorted.bam', 'output.bam')
    tk_bam.index(outs.output)
