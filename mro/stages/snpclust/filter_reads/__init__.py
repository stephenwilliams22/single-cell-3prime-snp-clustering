#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import itertools
import subprocess
import os
import tenkit.bam as tk_bam
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

__MRO__ = '''
stage FILTER_READS(
    in  bam[]    input_bams,  # array of BAMs
    in  tsv    cell_barcodes,
    in  map    align,
    out bam    output,
    src py     "stages/snpclust/filter_reads",
) split using (
    in bam input_bam,  # single BAM
)
'''

def split(args):
    return {'chunks': {'input_bam': x for x in args.input_bams}, "join": {'__mem_gb': 200, '__threads': 20}}

def main(args, outs):
    outs.coerce_strings()

    in_bam = tk_bam.create_bam_infile(args.input_bam)
    out_bam, _ = tk_bam.create_bam_outfile(outs.output, None, None, template=in_bam)
    cell_bcs = set(cr_utils.load_barcode_tsv(args.cell_barcodes))

    for (tid, pos), reads_iter in itertools.groupby(in_bam, key=cr_utils.pos_sort_key):
        dupe_keys = set()
        for read in reads_iter:
            if cr_utils.get_read_barcode(read) not in cell_bcs:
                continue

            if cr_utils.is_read_dupe_candidate(read, cr_utils.get_high_conf_mapq(args.align)):
                dupe_key = (cr_utils.si_pcr_dupe_func(read), cr_utils.get_read_umi(read))
                if dupe_key in dupe_keys:
                    continue

                dupe_keys.add(dupe_key)
                out_bam.write(read)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    #merg and index
    args_merge = ['sambamba', 'merge', '-t', str(args.__threads), 'output_merge.bam']
    #create an extended list to put at the end of args_merge
    args_merge.extend(input_bams)
    subprocess.check_call(args_merge)
    os.rename('output_merge.bam', outs.output)
    os.rename('output_merge.bam.bai', outs.output+'.bai')
    tk_bam.concatenate(outs.output, input_bams)
    tk_bam.index(outs.output)
