# WisecondorX

import logging

import numpy as np
import pysam
import sys

'''
Converts aligned reads file to numpy array by transforming
individual reads to counts per bin.
'''


def convert_reads(args):
    bins_per_chr = dict()
    
    for chr in range(0, 25):
        bins_per_chr[str(chr)] = None

    logging.info('Importing data ...')

    if args.infile.endswith(".bam"):
        reads_file = pysam.AlignmentFile(args.infile, 'rb')
    elif args.infile.endswith(".cram"):
        if args.reference is not None:
            reads_file = pysam.AlignmentFile(args.infile, 'rc', reference_filename=args.reference)
        else:
            logging.error("Cram support requires a reference file, please use the --reference argument")
            sys.exit(1)
    else:
        logging.error(
            "Unsupported input file type. Make sure your input filename has a correct extension ( bam or cram)")
        sys.exit(1)

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1
    larp2 = -1

    logging.info('Converting aligned reads ... This might take a while ...')
    
    for index, chr in enumerate(reads_file.references):

        chr_name = chr
        # remove the chr from chr1,chr2,chr3 etc...
        
        if chr_name[:3].lower() == 'chr':
            chr_name = chr_name[3:]

        chr_rename_dict = {
                "telomere":0,
                "NC_060925.1":1,
                "NC_060926.1":2,
                "NC_060927.1":3,
                "NC_060928.1":4,
                "NC_060929":5,
                "NC_060930":6,
                "NC_060931":7,
                "NC_060932":8,
                "NC_060933":9,
                "NC_060934":10,
                "NC_060935":11,
                "NC_060936":12,
                "NC_060937":13,
                "NC_060938":14,
                "NC_060939":15,
                "NC_060940":16,
                "NC_060941":17,
                "NC_060942":18,
                "NC_060943":19,
                "NC_060944":20,
                "NC_060945":21,
                "NC_060946":22,
                "NC_060947":23,
                "NC_060948":24
                }
        chr_name = chr_rename_dict[chr_name]
        if chr_name not in bins_per_chr and chr_name != 'X' and chr_name != 'Y':
            
            continue

        logging.info('Working at {}; processing {} bins'
                     .format(chr, int(reads_file.lengths[index] / float(args.binsize) + 1)))
        counts = np.zeros(int(reads_file.lengths[index] / float(args.binsize) + 1), dtype=np.int32)
        bam_chr = reads_file.fetch(chr)

        if chr_name == 'X':
            chr_name = '23'
        if chr_name == 'Y':
            chr_name = '24'
        
        for read in bam_chr:
            if read.is_paired:
                if not read.is_proper_pair:
                    reads_pairf += 1
                    continue
                if not args.normdup and larp == read.pos and larp2 == read.next_reference_start:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start
                reads_seen += 1
                larp = read.pos
            else:
                if not args.normdup and larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

        bins_per_chr[chr_name] = counts
        reads_kept += sum(counts)

    qual_info = {'mapped': reads_file.mapped,
                 'unmapped': reads_file.unmapped,
                 'no_coordinate': reads_file.nocoordinate,
                 'filter_rmdup': reads_rmdup,
                 'filter_mapq': reads_mapq,
                 'pre_retro': reads_seen,
                 'post_retro': reads_kept,
                 'pair_fail': reads_pairf}
    return bins_per_chr, qual_info
