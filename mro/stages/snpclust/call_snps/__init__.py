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

__MRO__ = '''
stage CALL_SNPS(
    in  path   reference_path,
    in  bam    input,
    in  int    n_donors,
    out vcf[]  output,
    src py     "stages/snpclust/call_snps_pd",
) split using (
    in  string locus,
)
'''

def split(args):
    in_bam = tk_bam.create_bam_infile(args.input)
    loci = tk_bam.generate_tiling_windows(in_bam, tk_constants.PARALLEL_LOCUS_SIZE)
    chunks = [{'locus': locus} for locus in loci]
    return {'chunks': chunks}

def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    chrom, start, stop = tk_io.get_locus_info(args.locus)
    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\n")
    
    dic_make_args = ['gatk-launch', 'CreateSequenceDictionary', '-R', genome_fasta_path]
    subprocess.call(dic_make_args)
    
    
    first_bam = martian.make_path('output.RG.bam')
    second_bam = martian.make_path('output.RG.STARcor.bam')
    rg_make_args = ['gatk-launch', 'AddOrReplaceReadGroups', '-I', args.input, '-O', 
                    first_bam, '-LB', 'lib1', '-PL', 'illumina',
                    '-PU', 'unit1', '-SM', 'example']
    subprocess.check_call(rg_make_args)
    
    samtools_args = '''samtools view -@ 4 -h {} | 
                       awk 'BEGIN{{OFS="\t"}} $5 == 255 {{ $5 = 60; print; next}} {{print}}' | 
                       samtools view -Sb -@ 4 - > {}'''.format(first_bam, second_bam)
    subprocess.call(samtools_args, shell=True)            
    
    samtools_index_args = ['samtools', 'index',second_bam]
    subprocess.call(samtools_index_args)
    
    gatk_args = ['gatk-launch', 'HaplotypeCaller', '-R', genome_fasta_path, '-I', second_bam, 
                 '--minimum-mapping-quality', '30', '--min-base-quality-score', '20', '-L', bed_path, '-O','output.vcf']
    
    #modify the GATK vcf to remove header issues
    sed_args = '''sed -i '/##FORMAT=<ID=PL/,/##INFO=<ID=AC/{//!d}' output.vcf'''
    subprocess.call(sed_args, shel=True)
    
    #os.remove first_bam
    
    with open(outs.output, 'w') as f:
        subprocess.call(gatk_args, stdout=f)

def join(args, outs, chunk_defs, chunk_outs):
    outs.output = [chunk.output for chunk in chunk_outs]
   
