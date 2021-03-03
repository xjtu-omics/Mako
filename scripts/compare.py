'''
Thanks for the comparison code from Peter.A Audano.
And we use some code from truvari (https://github.com/spiralgenetics/truvari)

@author: Jiadong Lin, Xi'an Jiaotong Univeristy, Leiden University

@contact: jiadong324@gmail.com

@time: 2019/11/4
'''

import pandas as pd
import numpy as np
from collections import OrderedDict, defaultdict
import logging
import pysam
import re
import argparse
import sys


def read_from_bed(bed_file, sample_name):
    df_raw = pd.read_csv(bed_file, delimiter="\t", header=None, usecols=[0, 1, 2, 3], names=["#CHROM", "POS", "END", "SVLEN"])
    sv_list = list()
    for index, entry in df_raw.iterrows():
        id = "{0}_{1}_{2}_{3}".format(entry["#CHROM"], entry["POS"], entry["END"], index)
        this_sv = (entry["#CHROM"], entry["POS"], entry["END"], int(entry["SVLEN"]), id, sample_name)
        sv_list.append(this_sv)

    df = pd.DataFrame(sv_list, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))

    return df

def read_svs_from_source(source, sample_name):

    df_source_raw = pd.read_csv(source, delimiter="\t", header=None, usecols=[0,1,2,3,4], names=["#CHROM", "POS", "END", "TYPE", "COMPs"])

    complete_target_list = list()
    sub_targe_list = list()
    for index, entry in df_source_raw.iterrows():
        id = "{0}_{1}_{2}_{3}".format(entry["#CHROM"], entry["POS"], entry["END"], index)
        sv_size = int(entry["END"]) - int(entry["POS"])
        complete_sv = (entry["#CHROM"], entry["POS"], entry["END"], sv_size, id, sample_name)
        complete_target_list.append(complete_sv)

        sub_comp_info = entry['COMPs'].replace(">",";").replace("<","")[:-1]
        sub_comp_idx = 1
        for sub_comp in sub_comp_info.split(";"):
            sub_info = sub_comp.split(",")
            comp_type = sub_info[0]

            if "translocation" in comp_type:
                sorc_start = int(sub_info[1].split(":")[1].split("-")[0])
                sorc_end = int(sub_info[1].split(":")[1].split("-")[1])
                sv_size = sorc_end - sorc_start
                ins_dest = int(sub_info[2].split(":")[1])
                ins_ori = sub_info[3]
                if ins_ori == "reverse":
                    comp_type = "inverted duplication"
                else:
                    comp_type = "dispersed duplication"

                if ins_dest > sorc_start:
                    start = sorc_start
                    end = ins_dest
                else:
                    start = ins_dest
                    end = sorc_start
            else:
                start = int(sub_info[1])
                end = int(sub_info[2])
                sv_size = end - start

            this_comp_id = "{0}_{1}_{2}_{3}_{4}".format(comp_type, sub_comp_idx, entry["#CHROM"], entry["POS"], entry["END"])
            this_comp = (entry["#CHROM"], start, end, sv_size, this_comp_id, sample_name)
            sub_targe_list.append(this_comp)
            sub_comp_idx += 1

    df_complete_source = pd.DataFrame(complete_target_list, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))
    df_sub_source = pd.DataFrame(sub_targe_list, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))

    return df_complete_source, df_sub_source

def read_target_vcf(target, sample_name):

    vcf_base = pysam.VariantFile(target, 'r')
    target_list = list()
    index = 0
    for record in vcf_base.fetch():
        alt = record.alts[0]
        chrom = record.contig
        start = record.pos
        end = record.stop

        new_alt_1 = re.findall(r"\]([A-Za-z0-9_:]+)\]", alt)
        new_alt_2 = re.findall(r"\[([A-Za-z0-9_:]+)\[", alt)

        if new_alt_1 != []:
            if new_alt_1[0].split(":")[0] == chrom:
                end = int(new_alt_1[0].split(":")[1])

        elif new_alt_2 != []:
            if new_alt_2[0].split(":")[0] == chrom:
                end = int(new_alt_2[0].split(":")[1])

        size = abs(end - start) + 1
        if "SVLEN" in record.info:
            size = int(record.info["SVLEN"][0])

        id = "{0}_{1}_{2}_{3}".format(chrom, start, end, index)
        this_sv = (chrom, start, end, size, id, sample_name)
        target_list.append(this_sv)
        index += 1

    df_target = pd.DataFrame(target_list, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))

    return df_target

def read_mako(targe, sample_name):
    ind_svs = list()
    complete_svs = list()

    index = 0
    sub_id_set = {}
    for line in open(targe, "r"):
        if "#" in line:
            continue

        tmp = line.strip().split("\t")
        cx_score = int(tmp[4].split(";")[0].split("=")[1])
        # Filter Mako with CXS score, 2 is the default value
        if cx_score > 2:
            chrom = tmp[0]
            sv_len = int(tmp[2]) - int(tmp[1])
            id = "{0}_{1}_{2}_{3}".format(chrom, tmp[1], tmp[2], index)

            complete_event = (chrom, int(tmp[1]), int(tmp[2]), sv_len, id, sample_name)
            complete_svs.append(complete_event)

            sv_info_tokens = tmp[4].split(";")

            for i in range(0, len(sv_info_tokens)):
                if len(sv_info_tokens[i].split(",")) < 3:
                    continue
                sv_info = sv_info_tokens[i].split(",")[0]
                this_start = int(sv_info.split("-")[0])
                this_end = int(sv_info.split("-")[1])

                ind_sv_id = "{0}_{1}_{2}_{3}".format(chrom, this_start, this_end, index)
                if ind_sv_id in sub_id_set:
                    dup_time = len(sub_id_set[ind_sv_id])
                    sub_id_set[ind_sv_id].append(1)
                    ind_sv_id = "{0}_{1}_{2}_{3}_{4}".format(chrom, this_start, this_end, index, dup_time)
                else:
                    sub_id_set[ind_sv_id] = [1]

                this_sv = (chrom, this_start, this_end, this_end - this_start, ind_sv_id, sample_name)
                ind_svs.append(this_sv)

            index += 1

    df_complete_target = pd.DataFrame(complete_svs, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))
    df_sep_target = pd.DataFrame(ind_svs, columns=('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES'))

    return df_complete_target, df_sep_target

def reciprocal_overlap(begin_a, end_a, begin_b, end_b):
    """
    Get reciprocal overlap of two intervals. Intervals are expected to be half-open coordinates (length is end - start).

    :param begin_a: Begin of interval a.
    :param end_a: End of interval a.
    :param begin_b: Begin of interval b.
    :param end_b: End of interval b.

    :return: A value between 0 and 1 if intervals a and b overlap and some negative number if they do not.
    """

    overlap = min(end_a, end_b) - max(begin_a, begin_b)

    return round(min([overlap / (end_a - begin_a), overlap / (end_b - begin_b)]), 3)

def size_similarity(size_a, size_b):
    return round(min(size_a, size_b) / max(size_a, size_b), 3)

def sim_score(ropct, simpct, dist, size_a, size_b):
    norm_dist = (size_a + size_b) / (size_a + size_b + dist)
    return round((2 * norm_dist + 1 * ropct + 1 * simpct) / 3.0, 3)

def order_variant_columns(df, head_cols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')):
    """
    Rearrange columns with a set list first (in defined order of `head_cols`) and leave the remaining columns
    in the order they were found.

    :param df: Data frame.
    :param head_cols: Columns to move to the first columns. Set variant BED order by default.

    :return: Data frame with rearranged columns.
    """

    # Check head columns
    head_cols = list(head_cols)

    for col in head_cols:
        if col not in df.columns:
            raise RuntimeError('Missing column in variant file: {}'.format(col))

    # Define tail columns
    tail_cols = [col for col in df.columns if col not in head_cols]

    # Arrange with head columns first. Leave remaining columns in order
    return df.loc[:, head_cols + tail_cols]


def get_variant_id(df):
    """
    Get variant IDs using '#CHROM', 'POS', 'SVTYPE', and 'SVLEN' columns.

    :param df: Dataframe.

    :return: A Series of variant IDs for `df`.
    """

    # Set ID
    return df.apply(lambda row: '{}-{}-{}-{}'.format(
        row['#CHROM'], row['POS'] + 1, row['SVTYPE'], row['SVLEN']
    ), axis=1)

def nearest_by_svlen_overlap(df_source, df_target, max_dist, size_sim, restrict_samples=False):
    """
    For each variant in df_source, get the nearest variant in df_target. Both dataframes must contain fields
    "#CHROM", "POS", "END", and "SVLEN" (which is the absolute value, no negative SVLEN for DELs). Return a dataframe
    with each source variant ("ID") with the best match by distance and then by size overlap ("TARGET_ID") along with
    distance ("DISTANCE"), reciprocal overlap if variants intersect ("RO", 0 if no overlap), and size proportion
    ("TARGET_SIZE_PROP", target-size / source-size).

    :param df_source: Source dataframe.
    :param df_target: Target dataframe.
    :param size_sim: Length proportion of matches.
    :param restrict_samples: If `True` and both dataframes contain a `MERGE_SAMPLES` column, then restrict matches to
        only those that share samples.

    :return: A dataframe with "ID", "TARGET_ID", "DISTANCE", "RO", and "TARGET_SIZE_PROP".
    """

    # Check for expected columns
    if any(col not in df_source.columns for col in ('#CHROM', 'POS', 'END', 'SVLEN', 'ID')):
        raise RuntimeError('Source Dataframe missing at least one column in ("#CHROM", "POS", "END", "SVLEN"): {}')

    if any(col not in df_target.columns for col in ('#CHROM', 'POS', 'END', 'SVLEN', 'ID')):
        raise RuntimeError('Target Dataframe missing at least one column in ("#CHROM", "POS", "END", "SVLEN"): {}')

    # IDs must be unique
    if len(set(df_source['ID'])) != df_source.shape[0]:
        raise RuntimeError('Source Dataframe IDs are not unique')

    if len(set(df_target['ID'])) != df_target.shape[0]:
        raise RuntimeError('target Dataframe IDs are not unique')

    # Determine if variants are matched on MERGE_SAMPLES
    if restrict_samples and 'MERGE_SAMPLES' in df_source.columns and 'MERGE_SAMPLES' in df_target.columns:
        restrict_samples = True
    else:
        restrict_samples = False

    # Subset and cast to int16 (default int64 uses more memory and is not needed)
    if restrict_samples:
        subset_cols = ('#CHROM', 'POS', 'END', 'SVLEN', 'ID', 'MERGE_SAMPLES')
        col_map = {
            '#CHROM': np.object, 'POS': np.int32, 'END': np.int32, 'SVLEN': np.int32, 'ID': np.object,
            'MERGE_SAMPLES': np.object
        }

    else:
        subset_cols = ('#CHROM', 'POS', 'END', 'SVLEN', 'ID')
        col_map = {
            '#CHROM': np.object, 'POS': np.int32, 'END': np.int32, 'SVLEN': np.int32, 'ID': np.object
        }

    stats_box = OrderedDict()
    stats_box["TP-base"] = 0
    stats_box["TP-call"] = 0
    stats_box["FP"] = 0
    stats_box["FN"] = 0
    stats_box["precision"] = 0
    stats_box["recall"] = 0
    stats_box["f1"] = 0
    stats_box["base cnt"] = 0
    stats_box["call cnt"] = 0

    df_source = df_source.loc[:, subset_cols]
    df_target = df_target.loc[:, subset_cols]

    # df_source.set_index('ID', inplace=True)
    copy_df_targat = df_target.set_index('ID', inplace=False)

    matched_calls = defaultdict(bool)

    if size_sim <= np.float16(0.0) or size_sim >= np.float16(1.0):
        raise RuntimeError('Length proportion must be between 0 and 1 (exclusive): {}'.format(size_sim))

    # Nearest by #CHROM
    overlap_list = list()  # Dataframe for each chrom (to be merged into one dataframe)

    for chrom in sorted(set(df_source['#CHROM'])):

        df_source_chr = df_source.loc[df_source['#CHROM'] == chrom]
        df_target_chr = copy_df_targat.loc[copy_df_targat['#CHROM'] == str(chrom)]

        # Target has no values on this chrom, all calls in source file are missed
        if df_target_chr.shape[0] == 0:
            stats_box["FN"] += df_source_chr.shape[0]
            continue

        for index, source_row in df_source_chr.iterrows():
            stats_box["base cnt"] += 1

            pos = source_row['POS']
            end = source_row['END']

            source_id = source_row['ID']

            # print(df_source_chr.columns)
            min_len = np.int32(source_row['SVLEN'] * size_sim)
            max_len = np.int32(source_row['SVLEN'] * (2 - size_sim))

            # Subset target - Skip if no targets within range
            df_target_chr_len = df_target_chr.loc[
                (df_target_chr['SVLEN'] >= min_len) & (df_target_chr['SVLEN'] <= max_len)
            ]
            # Cannot find targes of similar size to the source
            if df_target_chr_len.shape[0] == 0:
                stats_box["FN"] += 1
                continue


            # Get distance to all target records
            distance = df_target_chr_len.apply(
                lambda row: row['POS'] - end if row['POS'] > end else (
                    pos - row['END'] if pos > row['END'] else 0
                ),
                axis=1
            )

            # Assign index and subset minimum
            min_distance = np.min(distance)

            # Only consider potential match within max distance
            if min_distance > max_dist:
                stats_box["FN"] += 1
                continue

            distance = distance.loc[distance == min_distance]

            # Multiple match continue or find the best match
            if len(distance) > 1:
                match_id = np.argmin(np.abs(df_target_chr_len.loc[distance.index, 'SVLEN'] - source_row['SVLEN']))
                continue
            else:
                # Only one best match, get ID
                match_id = distance.index[0]

            # Save base ID to matched calls
            if not matched_calls[source_row["ID"]]:
                stats_box["TP-base"] += 1
            matched_calls[source_row["ID"]] = True


            if not matched_calls[match_id]:
                stats_box["TP-call"] += 1
            matched_calls[match_id] = True

            # Save match record
            overlap_list.append((index, source_id, match_id, min_distance, len(distance)))

    do_stats_math = True
    if stats_box["TP-base"] == 0 and stats_box["FN"] == 0:
        logging.warning("No TP or FN calls in base!")
        do_stats_math = False
    elif stats_box["TP-call"] == 0 and stats_box["FP"] == 0:
        logging.warning("No TP or FP calls in base!")
        do_stats_math = False
    else:
        logging.info("Results peek: %d TP-base %d FN %.2f%% Recall", stats_box["TP-base"], stats_box["FN"],
                     100 * (float(stats_box["TP-base"]) / (stats_box["TP-base"] + stats_box["FN"])))

    for index, entry in df_target.iterrows():
        # print(entry)
        if not matched_calls[entry["ID"]]:
            stats_box['FP'] += 1

    if do_stats_math:
        stats_box["precision"] = float(stats_box["TP-call"]) / (stats_box["TP-call"] + stats_box["FP"])
        stats_box["recall"] = float(stats_box["TP-base"]) / (stats_box["TP-base"] + stats_box["FN"])

    # f-measure
    neum = stats_box["recall"] * stats_box["precision"]
    denom = stats_box["recall"] + stats_box["precision"]
    if denom != 0:
        stats_box["f1"] = 2 * (neum / denom)
    else:
        stats_box["f1"] = "NaN"

    stats_box["call cnt"] = stats_box["TP-base"] + stats_box["FP"]

    # Merge records
    df = pd.DataFrame(overlap_list, columns=('ID', 'BASE_ID', 'TARGET_ID', 'DISTANCE','TARGET_MATCHES'))

    # Annotate overlap and size proportion (RO and TARGET_SIZE_PROP).
    df['#CHROM'] = list(df_source.loc[df['ID'], '#CHROM'])
    df['POS'] = list(df_source.loc[df['ID'], 'POS'])
    df['END'] = list(df_source.loc[df['ID'], 'END'])
    df['SVLEN'] = list(df_source.loc[df['ID'], 'SVLEN'])



    df['TARGET_POS'] = list(copy_df_targat.loc[df['TARGET_ID'], 'POS'])
    df['TARGET_END'] = list(copy_df_targat.loc[df['TARGET_ID'], 'END'])
    df['TARGET_SVLEN'] = list(copy_df_targat.loc[df['TARGET_ID'], 'SVLEN'])
    df['RO'] = df.apply(lambda row: reciprocal_overlap(*row.loc[['POS', 'END', 'TARGET_POS', 'TARGET_END']]), axis=1)

    # # Set TARGET_SIZE_PROP
    # df['TARGET_SIZE_PROP'] = round(df['TARGET_SVLEN'] / df['SVLEN'], 3)

    df['TARGET_SIZE_SIM'] = df.apply(
        lambda row: size_similarity(*row.loc[['TARGET_SVLEN', 'SVLEN']]), axis=1
    )

    # df['SIM_SCORE'] = df.apply(
    #     lambda row: sim_score(*row.loc[['RO', 'TARGET_SIZE_SIM', 'DISTANCE', 'SVLEN', 'TARGET_SVLEN']]), axis=1
    # )

    # Sort
    df.sort_values(['#CHROM', 'POS'])
    df = df.loc[:, ('TARGET_ID', 'BASE_ID', 'TARGET_MATCHES', 'RO','TARGET_SIZE_SIM')]

    # Return dataframe
    return df, stats_box

def multibp_match_eval(sub_overlap_df, source, tool):
    source_csv_sub_dict = {}
    df_source_raw = pd.read_csv(source, delimiter="\t", header=None, usecols=[0, 1, 2, 3, 4],
                                names=["#CHROM", "POS", "END", "TYPE", "COMPs"])

    for index, entry in df_source_raw.iterrows():
        sv_id = "{0}_{1}_{2}".format(entry["#CHROM"], entry["POS"], entry["END"])

        sub_comp_info = entry['COMPs'].replace(">", ";").replace("<", "")[:-1]
        source_csv_sub_dict[sv_id] = len(sub_comp_info.split(";"))

    all_csvs = len(source_csv_sub_dict)
    source_csv_sub_matched = {}
    for index, entry in sub_overlap_df.iterrows():
        base_id_tokens = entry["BASE_ID"].split("_")
        sv_id = "{0}_{1}_{2}".format(base_id_tokens[2], base_id_tokens[3], base_id_tokens[4])

        if sv_id in source_csv_sub_matched:
            source_csv_sub_matched[sv_id] += 1
        else:
            source_csv_sub_matched[sv_id] = 1
    matched_csv = {"all_match":0, "one_miss_match":0}
    for csv, comp_num in source_csv_sub_matched.items():
        if comp_num >= source_csv_sub_dict[csv] - 1:
            matched_csv["one_miss_match"] += 1
        if comp_num == source_csv_sub_dict[csv]:
            matched_csv["all_match"] += 1

    for k, v in matched_csv.items():
        matched_csv[k] = v / float(all_csvs)
    return matched_csv

def compare_eval(source, target, sample, max_dist, size_sim, output, tool, compare_method):

    if compare_method == "both":
        df_complete_source, df_sub_source = read_svs_from_source(source, sample)
        if tool == 'mako':
            complete_target, sub_target = read_mako(target, sample)
            complete_overlap_df, complete_stats_dict = nearest_by_svlen_overlap(df_complete_source, complete_target, max_dist, size_sim)
            sub_overlap_df, sub_stats_dict = nearest_by_svlen_overlap(df_sub_source, sub_target, max_dist, size_sim)
            multibp_match_stats = multibp_match_eval(sub_overlap_df, source, tool)
            writer_both_results(complete_overlap_df, complete_stats_dict, sub_overlap_df, multibp_match_stats, output, tool)
        else:
            df_target = read_target_vcf(target, sample)
            complete_overlap_df, complete_stats_dict = nearest_by_svlen_overlap(df_complete_source, df_target, max_dist, size_sim)
            sub_overlap_df, sub_stats_dict = nearest_by_svlen_overlap(df_sub_source, df_target, max_dist, size_sim)
            multibp_match_stats = multibp_match_eval(sub_overlap_df, source, tool)
            writer_both_results(complete_overlap_df, complete_stats_dict, sub_overlap_df, multibp_match_stats, output, tool)

    elif compare_method == "base":
        df_complete_source = read_from_bed(source, sample)
        if tool == 'mako':
            complete_target, sub_target = read_mako(target, sample)
            complete_overlap_df, complete_stats_dict = nearest_by_svlen_overlap(df_complete_source, complete_target, max_dist, size_sim)
            writer_base_results(complete_overlap_df, complete_stats_dict, output, tool)
        else:
            df_target = read_target_vcf(target, sample)
            complete_overlap_df, complete_stats_dict = nearest_by_svlen_overlap(df_complete_source, df_target, max_dist, size_sim)
            writer_base_results(complete_overlap_df, complete_stats_dict, output, tool)

def writer_both_results(complete_overlap, complete_stats, sub_overlap, sub_stats, output, tool):
    complete_overlap.to_csv(output + tool + "_base.overlapped_sites.bed", index=False, sep="\t")
    stats_writer = open(output + tool + "_base.overlapped_sites.stats.txt", "w")
    for k, v in complete_stats.items():
        stat_str = "{0}: {1}\n".format(k, v)
        stats_writer.write(stat_str)

    sub_overlap.to_csv(output + tool + "_sub.overlapped_sites.bed", index=False, sep="\t")
    stats_writer = open(output + tool + "_sub.overlapped_sites.stats.txt", "w")

    for k, v in sub_stats.items():
        stat_str = "{0}: {1}\n".format(k, v)
        stats_writer.write(stat_str)

    stats_writer.close()

def writer_base_results(complete_overlap, complete_stats, output, tool):
    complete_overlap.to_csv(output + tool + "_base.overlapped_sites.bed", index=False, sep="\t")
    stats_writer = open(output + tool + "_base.overlapped_sites.stats.txt", "w")
    for k, v in complete_stats.items():
        stat_str = "{0}: {1}\n".format(k, v)
        stats_writer.write(stat_str)
    stats_writer.close()

def compare_bed_files(source, target, sample,  max_dist, size_sim, output, tool):
    df_source = read_from_bed(source, sample)
    df_target = read_from_bed(target, sample)
    complete_overlap_df, complete_stats_dict = nearest_by_svlen_overlap(df_source, df_target, max_dist, size_sim)

    writer_base_results(complete_overlap_df, complete_stats_dict, output, tool)

def run():
    parser = argparse.ArgumentParser(description='Evaluate overlaps between prediction and benchmarks')
    subparsers = parser.add_subparsers(title='commands', dest='command')  # two submodules


    all_eval = subparsers.add_parser('both', help='Evaluate both entire and sub events')

    all_eval.add_argument("-b", dest="base", help="Source BED file (usually groud truth set)")
    all_eval.add_argument("-t", dest="target", help="Target VCF file (VCF file from caller)")
    all_eval.add_argument("-w", dest="dir", help="Result output directory")
    all_eval.add_argument("-s", dest="sample", help="Sample name")
    all_eval.add_argument("-p", dest="prefix", help="Output file prefix (use tool name for convenience)")
    all_eval.add_argument("-d", dest="dist", type=int, default=500, help="Max distance between breakpoints")
    all_eval.add_argument("-r", dest="size", type=float, default=0.7, help="SV size similarity ratio")

    # base_eval = subparsers.add_parser('base', help='Only compare the entire events')
    #
    # base_eval.add_argument("-b", dest="base", help="Source BED file (usually groud truth set)")
    # base_eval.add_argument("-t", dest="target", help="Target VCF file (VCF file from caller)")
    # base_eval.add_argument("-w", dest="dir", help="Result output directory")
    # base_eval.add_argument("-s", dest="sample", help="Sample name")
    # base_eval.add_argument("-p", dest="prefix", help="Output file prefix (use tool name for convenience)")
    # base_eval.add_argument("-d", dest="dist", type=int, default=500, help="Max distance between breakpoints")
    # base_eval.add_argument("-r", dest="size", type=float, default=0.7, help="SV size similarity ratio")

    bed_eval = subparsers.add_parser('bed', help='Compare events in two BED files')

    bed_eval.add_argument("-a", dest="bed_one", help="BED file one")
    bed_eval.add_argument("-b", dest="bed_two", help="BED file two")
    bed_eval.add_argument("-w", dest="dir", help="Result output directory")
    bed_eval.add_argument("-s", dest="sample", help="Sample name")
    bed_eval.add_argument("-p", dest="prefix", help="Output file prefix")
    bed_eval.add_argument("-d", dest="dist", type=int, default=500, help="Max distance between breakpoints")
    bed_eval.add_argument("-r", dest="size", type=float, default=0.7, help="SV size similarity ratio")

    args = parser.parse_args()
    if args.command == "both":
        compare_eval(args.base, args.target, args.sample, args.dist, args.size, args.dir, args.prefix, args.command)
    elif args.command == "bed":
        compare_bed_files(args.bed_one, args.bed_two, args.sample, args.dist, args.size, args.dir, args.prefix)

script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('compare.py         Last Update:2019-11-03\n')
    print('Compare between base set and results from callers\n')
    print('Usage:')
    print('compare.py [commands] <parameters>\n')
    print('Commands:')
    print('bed: evaluate unique-interval match between two BED files (for real data)')
    print('both: evaluate all-breakpoint match and unique-interval match (for simulate data)')
    print("=======================================================")
else:
    run()

# For debug
# if __name__ == '__main__':
    # work_dir = '/mnt/d/mako_works/skbr3/'
    # work_dir = '/mnt/d/mako_works/trios/yri_csvs/illumina/'
    # compare_eval(work_dir + "csv_bench.bed", work_dir + "tardis.bed", 'skbr3', 500, 0.5, work_dir, 'tardis')
    # compare_eval(work_dir + 'integration/csvs_refinedBp.bed', work_dir + 'mako_debug/mako_bp_filtered.bed', 'NA19240', 500, 0.5, work_dir + 'mako_debug/', 'mako_filtered')
    # compare_bed_files(work_dir + 'integration/csvs_refinedBp.bed', work_dir + 'mako_debug/mako.bed', 'NA19240', 500, 0.7, work_dir + 'mako_debug/', 'mako')
    # compare_bed_files(work_dir + 'integration/csvs_refinedBp.bed', work_dir + 'tardis.bed', 'NA19240', 500, 0.5, work_dir, 'tardis')
    # compare_bed_files(work_dir + 'integration/csvs_refinedBp.bed', work_dir + 'gridss.bed', 'NA19240', 500, 0.5, work_dir, 'gridss')
    # compare_bed_files(work_dir + 'integration/csvs_refinedBp.bed', work_dir + 'gridss.bed', 'NA19240', 500, 0.5, work_dir,'gridss')
    # compare_eval("/mnt/d/mako_works/visor_sim/known_csv_chr1/version2/chr1.added_SVs.nested.bed", "/mnt/d/mako_works/visor_sim/known_csv_chr1/version2/p100_c30/lumpy.vcf", "sim",
    #              500, 0.7, "/mnt/d/mako_works/visor_sim/known_csv_chr1/version2/p100_c30/", "lumpy", "both")