#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong Univeristy, Leiden University
@contact: jiadong324@gmail.com
@time: 2019/11/4
'''


import sys

from optparse import OptionParser
import pysam
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import re

parser = OptionParser()
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
          "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

class Interval:
    def __init__(self, chrom, start, end, pattern, sample, interval_str):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.pattern = pattern
        self.sample = sample
        self.interval = interval_str

    def overlap(self, interval, max_dist, len_prop):
        this_size = self.end - self.start
        inter_size = interval.end - interval.start

        min_len = this_size * len_prop
        max_len = this_size * (2 - len_prop)

        # Two intervals overlap
        if min(self.end, interval.end) >= max(self.start, interval.start) and self.chrom == interval.chrom:
            # Check breakpoint distance
            if abs(self.start - interval.start) <= max_dist and abs(self.end - interval.end) <= max_dist:
                # Check if SV size matches
                if inter_size >= min_len and inter_size <= max_len:
                    return True

        return False

    def toString(self):

        out_str = "{0}\t{1}\t{2}\t".format(self.chrom, self.start, self.end)
        interval_tokens = self.interval.split(";")
        sample_tokens = self.sample.split(";")

        sample_str = ""
        for i in range(len(sample_tokens)):
            sample_str += "{0},{1};".format(sample_tokens[i], interval_tokens[i])

        out_str += sample_str[:-1] + "\t" + self.pattern

        return out_str


def mako_to_vcf(mako, out):
    '''
    Convert Mako raw output to standard VCF format
    :param mako:
    :param out:
    :param ref:
    :param sample:
    :return:
    '''

    print("Convert to VCF ...")

    names = ("chr", "start", "end", "type", "filter", "info", "pattern", "weight")

    call_list = list()

    calls = pd.read_csv(mako, header=None, sep="\t", names=names)

    for idx, row in calls.iterrows():
        info_tokens = row["info"].split(";")
        qual = info_tokens[0].split("=")[1]
        brkp_str = ""
        supp_str = ""
        CR = 0
        for i in range(1, len(info_tokens)):
            info = info_tokens[i]
            if "cr" in info:
                CR = info.split("=")[1]
            else:
                brkp_str += "BRK{0}={1},{2};".format(i, info.split(",")[0], info.split(",")[1])
                supp_str += "SA{0}={1},RP{2}={3};".format(i, info.split(",")[2].split("=")[1], i, info.split(",")[3].split("=")[1])

        supp_str += "CR={0}".format(CR)

        alt = "<SV>"
        if "," in row["type"]:
            alt = "<CSV>"
        svlen = int(row["end"]) - int(row["start"])

        vcf_info_str = "END={0};SVLEN={1};SVTYPE={2};{3};{4};PATTERN={5};WEIGHT={6}".format(row["end"], svlen, row['type'], brkp_str[:-1], supp_str, row['pattern'], row['weight'])

        this_call = (row["chr"], row["start"], "N", alt, qual, row["filter"], vcf_info_str)

        call_list.append(this_call)

    df_calls = pd.DataFrame(call_list, columns=["#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO"])

    sorted_df_calls = df_calls.sort_values(['#CHROM', 'POS'], ascending=[True, True])

    with open(out, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=Mako V1.0\n")
        vcf.write(
            '##REF=<Description="Reference sequence at variant region">' + "\n")
        vcf.write(
            '##ALT=<ID=<SV>, Description="Simple SV inferred from subgraph">' + "\n")
        vcf.write(
            '##ALT=<ID=<CSV>, Description="Complex SV inferred from subgraph">' + "\n")
        vcf.write(
            '##QUAL=<Description="SV complexity score derived from subgraph">' + "\n")
        vcf.write(
            '##FILTER=<ID=ARP_Span,Number=1,Type=String,Description="SVs supported by ARP derived edges connecting two subgraphs">' + "\n")
        vcf.write('##FILTER=<ID=ARP_Self,Number=1,Type=String,Description="SVs supported by ARP derived edges in the subgraph">' + "\n")
        vcf.write('##FILTER=<ID=Split,Number=1,Type=String,Description="SVs supported by edges derived from split alignment">' + "\n")
        vcf.write(
            '##FILTER=<ID=SVTYPE,Number=1,Type=String,Description="Inferred SV type from edge connections in the subgraph">' + "\n")
        vcf.write(
            '##INFO=<ID=BRK,Number=1,Type=String,Description="Internal breakpoints derived from subgraph">' + "\n")
        vcf.write(
            '##INFO=<ID=SA,Number=1,Type=String,Description="Number of split alignment supporting an edge connection in the subgraph">' + "\n")
        vcf.write(
            '##INFO=<ID=RP,Number=1,Type=String,Description="Number of discordant read-pairs supporting an edge connection in the subgraph">' + "\n")
        vcf.write(
            '##INFO=<ID=CR,Number=1,Type=Integer,Description="Length of cross matched sequence">' + "\n")
        vcf.write(
            '##INFO=<ID=PATTERN,Number=1,Type=String,Description="Node types of the subgraph">' + "\n")
        vcf.write(
            '##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description="Read support of each node in the subgraph">' + "\n")

        sorted_df_calls.to_csv(vcf, sep="\t", index=False)


def merge_multiple_makos(sample_files, mako_dir, out_file, max_dist, len_prop):
    '''
    Merge multiple Mako call set
    :param sample_files:
    :param mako_dir:
    :param out_file:
    :param max_dist:
    :param len_prop:
    :return:
    '''

    intervals = []
    for line in open(sample_files, "r"):
        file_name = line.strip()
        mako_file_path = mako_dir + file_name
        sample_name = file_name.split(".")[0]

        cur_sample_sv_num = 0

        for line in open(mako_file_path, "r"):
            tmp = line.strip().split("\t")
            chrom = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])

            cur_sample_sv_num += 1
            sv_info_tokens = tmp[4].split(';')
            pattern_str = ""
            for token in sv_info_tokens:
                if token.split('=')[0] == 'Pattern':
                    pattern_str = token.split('=')[1]
                    break

            interval_str = "{0},{1},{2}".format(chrom, start, end)

            cur_interval = Interval(chrom, start, end, pattern_str, sample_name, interval_str)

            intervals = add_interval(intervals, cur_interval, max_dist, len_prop)

        # print sample_name + ", " + str(cur_sample_sv_num) + " SVs processed .."
        print("Merge sample: {0} total entries: {1}".format(sample_name, len(intervals)))

    writer = open(out_file, "w")

    for interval in intervals:
        writer.write(interval.toString() + "\n")
    writer.close()


def add_interval(intervals, new_interval, max_dist, len_prop):
    new_intervals = []

    num = len(intervals)

    if num == 0:
        new_intervals.append(new_interval)
        return new_intervals

    if new_interval.end < intervals[0].start or new_interval.start > intervals[num - 1].end:

        if new_interval.end < intervals[0].start:
            new_intervals.append(new_interval)

        new_intervals.extend(intervals)

        if new_interval.start > intervals[num - 1].end:
            new_intervals.append(new_interval)

        return new_intervals

    for i in range(len(intervals)):
        ele = intervals[i]
        overlap = ele.overlap(new_interval, max_dist, len_prop)
        # Overlapped
        if not overlap:
            new_intervals.append(ele)

            # check if given interval lies between two intervals
            if i < num and new_interval.start > intervals[i].end and new_interval.end < intervals[i + 1].start:
                new_intervals.append(new_interval)

            continue

        new_start = min(ele.start, new_interval.start)
        new_pattern = ele.pattern
        new_end = max(ele.end, new_interval.end)
        new_sample = ele.sample
        new_interval_str = ele.interval

        while i < num and overlap:
            new_end = max(intervals[i].end, new_interval.end)
            new_pattern += ";" + new_interval.pattern
            new_sample += ";" + new_interval.sample
            new_interval_str += ";" + new_interval.interval
            if i == num - 1:
                overlap = False
            else:
                overlap = intervals[i + 1].overlap(new_interval, max_dist, len_prop)

            i += 1

        i -= 1

        new_intervals.append(
            Interval(new_interval.chrom, new_start, new_end, new_pattern, new_sample, new_interval_str))

    return new_intervals


def mako_filter(in_file, out_file, cxs, format):
    '''
    Filter mako raw call site with different evidence
    :param in_file:
    :param out_file:
    :return:
    '''
    writer = open(out_file, 'w')
    csvs_num = 0
    all_calls = 0
    for line in open(in_file, 'r'):
        if "#" in line:
            continue
        all_calls += 1
        tmp = line.strip().split("\t")
        cx_score = int(tmp[5].split(";")[0].split("=")[1])
        # if cx_score not in scores:
        #     scores.append(cx_score)
        if cx_score > cxs:
            csvs_num += 1
            if format == "mako":
                writer.write(line)
            elif format == "bed":
                sv_len = int(tmp[2]) - int(tmp[1])
                out_str = "{0}\t{1}\t{2}\t{3}\n".format(tmp[0], tmp[1], tmp[2], sv_len)
                writer.write(out_str)
    print("Number of calls after filtering: ", csvs_num)
    writer.close()
    # print(np.percentile(scores, 25))

def not_primary(aln):
    return aln.is_supplementary or aln.is_secondary

def classify_rps(bam, fai_file, min_mapq, min_insert, max_insert):
    seen_aln = {}
    npairs = 0
    genome_length = 0
    rp_type_dict = {}

    with open(fai_file, 'r') as f:
        for line in f:
            entries = line.strip().split("\t")
            chrom = entries[0]
            if "chr" not in chrom:
                chrom = "chr{0}".format(chrom)

            if chrom in CHROMS:
                genome_length += int(entries[1])
    print("Genome length: ", genome_length)

    for aln in bam.fetch(until_eof=True):
        if not_primary(aln) or aln.is_duplicate or aln.is_unmapped or aln.mate_is_unmapped:
            continue
        chrom = aln.reference_name
        if "chr" not in chrom:
            chrom = "chr{0}".format(chrom)

        if chrom not in CHROMS:
            continue

        if aln.qname not in seen_aln:
            seen_aln[aln.qname] = aln
            continue
        mate = seen_aln[aln.qname]
        npairs += 1
        del seen_aln[aln.qname]

        if npairs % 1000000 == 0:
            print("[bam summary] processed read-pairs: ", npairs)
        if aln.mapq < min_mapq or mate.mapq < min_mapq or aln.is_unmapped or \
           mate.is_unmapped or not_primary(aln) or not_primary(mate):
            continue

        ilen = abs(aln.reference_start - mate.reference_end)
        sig_type = ""
        if aln.is_reverse != mate.is_reverse:
            second = aln if aln.is_reverse else mate
            first = aln if second is mate else mate
            if ilen > max_insert:
                sig_type = 'ARP_LARGE'

            elif (first.reference_start > second.reference_start) or \
                 (first.reference_end > second.reference_end):
                sig_type = 'ARP_RF'

            elif ilen < min_insert:
                sig_type = 'ARP_SMALL'
        else:
            sig_type = 'ARP_RR' if aln.is_reverse else "ARP_FF"

        if sig_type == "":
            continue

        if sig_type in rp_type_dict:
            rp_type_dict[sig_type] += 1
        else:
            rp_type_dict[sig_type] = 1

    return rp_type_dict, genome_length


def mako_config(bam, fai_file, num_to_check, min_mapq, out, sample):

    required = 97
    restricted = 3484
    flag_mask = required | restricted

    read_length = 0
    read_counter = 0

    L = []

    bam_file = pysam.AlignmentFile(bam, "r")

    for read in bam_file.fetch():
        if read_counter >= num_to_check:
            break

        cigar = read.cigarstring
        if cigar == None:
            continue

        read_length = get_read_length(cigar)
        flag = read.flag
        refname = read.reference_name
        mate_refname = read.next_reference_name
        isize = read.template_length

        valid = mate_refname == refname and flag & flag_mask == required and isize >= 0

        if valid:
            read_counter += 1
            L.append(isize)

    L = np.array(L)
    L.sort()
    med, umad = unscaled_upper_mad(L)
    upper_cutoff = med + 30 * umad
    L = L[L < upper_cutoff]

    mean = int(np.mean(L))
    stdev = int(np.std(L))

    min_insert = mean - 3 * stdev
    max_insert = mean + 3 * stdev

    print("mean: {0}\tstd: {1}\nStart to classify disocrdant read-pairs".format(mean, stdev))

    rp_lambda, genome_length = classify_rps(bam_file, fai_file, min_mapq, min_insert, max_insert)

    out_str = "mean:{0}\nstdev:{1}\nreadlen:{2}\nworkDir:{3}\nbam:{4}\nname:{5}\n".format(mean, stdev, read_length, out, out + bam, sample)
    for rp_type, val in rp_lambda.items():
        out_str += "{0}:{1}\n".format(rp_type, val / genome_length)

    print("All discordant read pairs processed!")

    writer = open(out + sample + '.mako.cfg', 'w')
    writer.write(out_str)
    writer.close()


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    med = np.median(xs)
    return med, np.median(xs[xs > med] - med)


def get_read_length(cigar):
    cigarPattern = '([0-9]+[MIDNSHP])'
    cigarSearch = re.compile(cigarPattern)
    atomicCigarPattern = '([0-9]+)([MIDNSHP])'
    atomicCigarSearch = re.compile(atomicCigarPattern)

    readLen = 0
    if (cigar != '*'):
        cigarOpStrings = cigarSearch.findall(cigar)

        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)[0]
            readLen += int(cigarOpList[0])

    return readLen

def get_mako_sub(mako, out):

    ind_svs = []

    for line in open(mako, "r"):
        if "#" in line:
            continue

        tmp = line.strip().split("\t")
        chrom = tmp[0]

        sv_info = tmp[4]
        if ";;" in sv_info:
            sv_info_tokens = sv_info.split(";;")
            for info_token in sv_info_tokens:
                tmp_token = info_token.split(",")
                if len(tmp_token) == 2 and "-" in tmp_token[0]:
                    this_start = tmp_token[0].split("-")[0]
                    this_end = tmp_token[0].split("-")[1]
                    this_sv = (chrom, this_start, this_end, tmp_token[1])
                    ind_svs.append(this_sv)
                else:
                    for i in range(2, len(tmp_token)):
                        token = tmp_token[i]
                        if "-" in token:
                            this_start = token.split("-")[0]
                            this_end = token.split("-")[1]
                            this_sv = (chrom, this_start, this_end, tmp_token[0] + "," + tmp_token[1] + "," + tmp_token[i + 1])
                            ind_svs.append(this_sv)
        else:
            sv_info_tokens = sv_info.split(",")
            for i in range(2, len(sv_info_tokens)):
                token = sv_info_tokens[i]
                if "-" in token:
                    this_start = token.split("-")[0]
                    this_end = token.split("-")[1]
                    this_sv = (chrom, this_start, this_end, sv_info_tokens[0] + "," + sv_info_tokens[1] + "," + sv_info_tokens[i + 1])
                    ind_svs.append(this_sv)


    writer = open(out, "w")

    for sv in ind_svs:

        out_str = "{0}\t{1}\t{2}\t{3}\n".format(sv[0], sv[1], sv[2], sv[3])

        writer.write(out_str)

    writer.close()

script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('ParseMako.py         Last Update:2020-7-20\n')
    print('This script is used to process Mako raw callset\n')
    print('Usage:')
    print('ParseMako.py [options] <parameters>\n')
    print('Options:')
    print('config: Create config file for Mako input ')
    print('filter: filter Mako calls ')
    print('tovcf: convert Mako calls to standard VCF format')
    # print('merge: merge mutiple Mako calls')

    print("=======================================================")
else:
    option = sys.argv[1]

    if option == "config":
        parser.add_option("-b", dest='bam', help='BAM file to config')
        parser.add_option("-n", type=int, dest="num", help="Number of samples used for estimation")
        parser.add_option("-m", type=int, dest="mapq", help="Minimum mapping quality for aligned read-pairs", default=20)
        parser.add_option("-f", dest="fai", help="Index of reference file")
        parser.add_option("-w", dest='out', help='Working directory')
        parser.add_option("-s", dest="name", help="Name of the sample")

        (options, args) = parser.parse_args()
        mako_config(options.bam, options.fai, options.num, options.mapq, options.out, options.name)

    elif option == "tovcf":

        parser.add_option("-m", dest="mako")
        parser.add_option("-o", dest="out")
        (options, args) = parser.parse_args()

        if not options.mako:
            parser.error("Mako call not given")

        if not options.out:
            parser.error("VCF output not given")

        mako_to_vcf(options.mako, options.out)

    elif option == "filter":

        parser.add_option("-i", dest="input", help="Input Mako callset")
        parser.add_option("-o", dest="out", help="Output of filtered Mako callset by CXS")
        parser.add_option("-c", type=int, dest="cxs", help="CXS threshold")
        parser.add_option("-f", dest="format", help="output format (original, bed)")
        (options, args) = parser.parse_args()
        mako_filter(options.input, options.out, options.cxs, options.format)


# if __name__ == '__main__':
#     bam_file = "/Users/jiadonglin/Data/HG00733/HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram.bam"
#     fai_file = "/Users/jiadonglin/Data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
#     bam_stats = "/Users/jiadonglin/Data/HG00733/MakoV1/bam_summary.txt"
#     classify_rps(bam_file, fai_file, 20, 101, 1037, bam_stats)
