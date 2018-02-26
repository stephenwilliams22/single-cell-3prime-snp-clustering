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
    return {'chunks': chunks, "join": {'__mem_gb': 4, '__threads': 10}}

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
                 '--create-output-bam-index', 'false',
                 '--TMP_DIR', os.getcwd()]
                 
    subprocess.check_call(star_args)
    
    #join the bams together
def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    args_merge = ['samtools', 'merge', '-@', str(args.__threads), outs.output]
    args.extend(input_bams)
    subprocess.check_call(args_merge)
    #tk_bam.concatenate(outs.output, input_bams)
    #tk_bam.sort(outs.output)
    args_sort = ['samtools', 'sort', '-@', str(args.__threads), 'output_sorted.bam']
    subprocess.check_call(args_sort)
    os.remove(outs.output)
    os.rename('output_sorted.bam', 'output.bam')
    tk_bam.index(outs.output)
