#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2019/12/5
'''
import pysam
from optparse import OptionParser


def convert_gridss(vcf_file, bed_out):
    sv_type_dict = {}
    sv_len_dict = {}
    writer = open(bed_out, "w")
    vcf_base = pysam.VariantFile(vcf_file, 'r')

    filters = ["ASSEMBLY_TOO_FEW_READ", "ASSEMBLY_TOO_SHORT", "INSUFFICIENT_SUPPORT",
               "LOW_QUAL"]

    trans = []
    validate_svs = 0
    for record in vcf_base.fetch():
        chrom = record.chrom
        start = record.pos
        end = record.stop

        is_filtered = False
        alt = record.alts[0]

        for fil in record.filter.keys():
            if fil in filters:
                is_filtered = True
                continue

        if is_filtered:
            continue

        rp = record.info["RP"]
        sr = record.info["SR"]
        if int(rp) == 0 and int(sr) == 0:
            continue
        if "[" in alt:
            bnd_info = alt[alt.find('[')+1:alt.rfind('[')].split(":")
            bnd_chrom = bnd_info[0]
            if bnd_chrom != chrom:
                continue
            bnd_pos = int(bnd_info[1])
        elif "]" in alt:
            bnd_info = alt[alt.find(']') + 1:alt.rfind(']')].split(":")
            bnd_chrom = bnd_info[0]
            if bnd_chrom != chrom:
                continue
            bnd_pos = int(bnd_info[1])

        sv_type = record.info["SVTYPE"]
        if sv_type in sv_type_dict:
            sv_type_dict[sv_type] += 1
        else:
            sv_type_dict[sv_type] = 1

        if sv_type == "BND":
            if bnd_pos > start:
                end = bnd_pos
            else:
                tmp = start
                start = bnd_pos
                end = tmp

            this_trans = "{0}-{1}".format(start, end)
            if this_trans in trans:
                continue
            else:
                trans.append(this_trans)

        sv_len = int(end) - int(start)
        if "SVLEN" in record.info.keys():
            sv_len = record.info["SVLEN"]

        if sv_len < 50:
            continue

        if sv_type in sv_type_dict:
            sv_type_dict[sv_type] += 1
        else:
            sv_type_dict[sv_type] = 1

        if sv_len in sv_len_dict:
            sv_len_dict[sv_len] += 1
        else:
            sv_len_dict[sv_len] = 1

        validate_svs += 1
        bed_out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrom, start, end, sv_len, sv_type, rp, sr)
        writer.write(bed_out)

    print("Obtained number of SVs: ", validate_svs)

def convert_tardis(vcf_file, bed_out):
    writer = open(bed_out, "w")
    vcf_base = pysam.VariantFile(vcf_file, 'r')
    sv_type_dict = {}
    sv_len_dict = {}

    validate_svs = 0
    for record in vcf_base.fetch():

        chrom = record.chrom
        start = record.pos
        end = record.stop
        sv_type = record.info["SVTYPE"]
        sv_len = end - start

        if "SVLEN" in record.info.keys():
            sv_len = abs(record.info["SVLEN"][0])

        if sv_len < 50:
            continue

        if sv_type in sv_type_dict:
            sv_type_dict[sv_type] += 1
        else:
            sv_type_dict[sv_type] = 1

        if sv_len in sv_len_dict:
            sv_len_dict[sv_len] += 1
        else:
            sv_len_dict[sv_len] = 1

        rp = record.info["RPSUP"]
        sr = record.info["SRSUP"]

        if start > end:
            tmp_value = start
            start = end
            end = tmp_value
        validate_svs += 1
        bed_out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrom, start, end, sv_len, sv_type, rp, sr)
        writer.write(bed_out)


def convert_svelter(vcf_file, bed_out, sample):
    writer = open(bed_out, "w")
    vcf_base = pysam.VariantFile(vcf_file, 'r')
    sv_type_dict = {}
    sv_len_dict = {}

    validate_svs = 0
    for record in vcf_base.fetch():

        chrom = record.chrom
        start = record.pos
        end = record.stop
        sv_type = record.info["SVTYPE"]

        if sample == 'yri':
            if float(record.info["VALIDATION"]) == -1:
                continue

        sv_len = end - start

        if "SVLEN" in record.info.keys():
            sv_len = abs(record.info["SVLEN"][0])

        if sv_len < 50:
            continue

        if sv_type in sv_type_dict:
            sv_type_dict[sv_type] += 1
        else:
            sv_type_dict[sv_type] = 1

        if sv_len in sv_len_dict:
            sv_len_dict[sv_len] += 1
        else:
            sv_len_dict[sv_len] = 1

        validate_svs += 1
        bed_out = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom, start, end, sv_len, sv_type)
        writer.write(bed_out)
    print("Obtained number of SVs: ", validate_svs)

def run():
    parser = OptionParser()
    parser.add_option("-s", dest="sample", help="Sample name (yri, skbr3)")
    parser.add_option("-t", dest="tool", help="Detection method (gridss, svelter, tardis)")
    parser.add_option("-v", dest="vcf", help="Path to VCF file")
    parser.add_option("-b", dest="bed", help="Path to output BED file")
    (options, args) = parser.parse_args()
    if options.tool == 'gridss':
        convert_gridss(options.vcf, options.bed)

    elif options.tool == 'svelter':
        convert_svelter(options.vcf, options.bed, options.sample)

    elif options.tool == 'tardis':
        convert_tardis(options.vcf, options.bed)

if __name__ == '__main__':
    run()