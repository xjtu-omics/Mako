#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/1/7
'''

import sys
import argparse
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from statistics import median
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import math
from scipy.signal import find_peaks

cut_offs = [0.1, 5, 10, 50, 100, 500, 1000, 10000]

class Candidate:

    def __init__(self, sv_list):
        starts = []
        ends = []
        type_str = ''
        for sv in sv_list:
            starts.append(sv[1])
            ends.append(sv[2])
            self.chrom = sv[0]
            type_str += sv[3] + ","

        self.start = median(starts)
        self.end = median(ends)
        self.len = self.end - self.start
        self.type = type_str[:-1]

    def to_string(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\n".format(self.chrom, int(self.start), int(self.end), int(self.len), self.type)

def span_position_distance(bp1, bp2):
    center1 = (bp1[0] + bp1[1]) // 2
    center2 = (bp2[0] + bp2[1]) // 2
    position_distance = min(abs(bp1[0] - bp2[0]), abs(bp1[1] - bp2[1]), abs(center2 - center1)) / 1000
    return position_distance

def find_threshold(bed_file, chrom):

    svs = []
    for line in open(bed_file, 'r'):
        tmp = line.strip().split("\t")
        # chrom, start, end, type, len
        if tmp[0] != chrom:
            continue
        this_sv = (tmp[0], int(tmp[1]), int(tmp[2]), tmp[3], float(tmp[4]))
        svs.append(this_sv)

    sorted_svs = sorted(svs, key=lambda k: (k[1] + k[2]) // 2)

    data = np.array([[sv[1], sv[2]] for sv in sorted_svs])
    print("Caculate cluster cutoff ...")
    Z = linkage(data, method = "average", metric = span_position_distance)

    cut_clusters_num = []
    for cut in cut_offs:
        cluster_indices = list(fcluster(Z, cut, criterion='distance'))
        new_clusters = [[] for i in range(max(cluster_indices))]
        for sv_index, cluster_index in enumerate(cluster_indices):
            new_clusters[cluster_index - 1].append(sorted_svs[sv_index])

        valid_candidates = []
        valid_candidates_size = []
        for cluster in new_clusters:
            if len(cluster) > 1:
                this_candidate = Candidate(cluster)
                if this_candidate.len == 0:
                    continue
                valid_candidates.append(this_candidate)
                valid_candidates_size.append(math.log10(this_candidate.len))
        cut_clusters_num.append(len(valid_candidates))
        plt.show()

    peak_index = find_peaks(cut_clusters_num)[0][0]

    return cut_offs[peak_index], cut_clusters_num

def draw_dendrogram(bed_file, chrom, figure_out):
    svs = []
    for line in open(bed_file, 'r'):
        tmp = line.strip().split("\t")
        # chrom, start, end, type, len
        if tmp[0] != chrom:
            continue
        this_sv = (tmp[0], int(tmp[1]), int(tmp[2]), tmp[3], float(tmp[4]))
        svs.append(this_sv)

    sorted_svs = sorted(svs, key=lambda k: (k[1] + k[2]) // 2)

    data = np.array([[sv[1], sv[2]] for sv in sorted_svs])
    print("Caculate cluster cutoff ...")
    Z = linkage(data, method="average", metric=span_position_distance)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    dendrogram(Z, no_labels=True, ax=ax)
    ax.set_title(chrom + " hierarchical clustering tree view")
    ax.set_ylabel("Distance")
    # plt.show()
    plt.tight_layout()
    plt.savefig(figure_out, dpi=600)

def draw_threshold_curve(chroms, bed_file, figure_out):

    cluster_num_dict = {}
    for chrom in chroms:
        select_thresh, cluster_num_list = find_threshold(bed_file, chrom)
        cluster_num_dict[chrom] = cluster_num_list

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    x_values = np.arange(len(cut_offs))
    for chrom, value in cluster_num_dict.items():
        ax.plot(x_values, value, marker='.', label=chrom)

    ax.set_yticks(np.linspace(0,200, 5))
    ax.set_yticklabels([int(ele) for ele in np.linspace(0,200, 5)])
    ax.set_ylabel("Nubmer of clusters")
    ax.set_xticks(np.arange(len(cut_offs)))
    ax.set_xticklabels(cut_offs)
    ax.set_xlabel("Cluster merge cutoff")
    ax.legend(bbox_to_anchor=(0.98, 1), fancybox=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.show()
    plt.tight_layout()
    plt.savefig(figure_out, dpi=600)

def yri_bp_clustering(bed_file, chrom, cluster_thresh, out_dir):

    out_file = out_dir + chrom + ".sv_clusters_ct" + str(cluster_thresh) + ".bed"
    out_writer = open(out_file, 'w')
    svs = []
    for line in open(bed_file, 'r'):
        tmp = line.strip().split("\t")
        # chrom, start, end, type
        if tmp[0] != chrom:
            continue
        this_sv = (tmp[0], int(tmp[1]), int(tmp[2]), tmp[3])
        svs.append(this_sv)

    sorted_svs = sorted(svs, key=lambda k: (k[1] + k[2]) // 2)

    data = np.array([[sv[1], sv[2]] for sv in sorted_svs])
    print("Start clustering ...")
    Z = linkage(data, method = "average", metric = span_position_distance)
    cluster_indices = list(fcluster(Z, cluster_thresh, criterion='distance'))

    new_clusters = [[] for i in range(max(cluster_indices))]
    for sv_index, cluster_index in enumerate(cluster_indices):
        new_clusters[cluster_index - 1].append(sorted_svs[sv_index])

    candidates = []
    for cluster in new_clusters:
        if len(cluster) > 1:
            candidates.append(Candidate(cluster))
    num = len(candidates)
    sorted_candidates = sorted(candidates, key=lambda k:k.start)
    for cand in sorted_candidates:
        if cand.len > 1000000:
            num -= 1
            continue
        out_writer.write(cand.to_string())

    print("Candidates found: ", num)

def skbr3_bp_clustering(bed_file, chrom, cluster_thresh, out_dir):

    out_file = out_dir + chrom + ".sv_clusters_ct" + str(cluster_thresh) + ".bed"
    out_writer = open(out_file, 'w')
    svs = []
    for line in open(bed_file, 'r'):
        tmp = line.strip().split("\t")
        # chrom, start, end, type
        if tmp[0] != chrom:
            continue
        this_sv = (tmp[0], int(tmp[1]), int(tmp[2]), tmp[3])
        svs.append(this_sv)

    sorted_svs = sorted(svs, key=lambda k: (k[1] + k[2]) // 2)

    data = np.array([[sv[1], sv[2]] for sv in sorted_svs])
    print("Start clustering ...")
    Z = linkage(data, method = "average", metric = span_position_distance)
    cluster_indices = list(fcluster(Z, cluster_thresh, criterion='distance'))

    new_clusters = [[] for i in range(max(cluster_indices))]
    for sv_index, cluster_index in enumerate(cluster_indices):
        new_clusters[cluster_index - 1].append(sorted_svs[sv_index])

    candidates = []
    for cluster in new_clusters:
        if len(cluster) > 1:
            candidates.append(Candidate(cluster))
    num = len(candidates)
    sorted_candidates = sorted(candidates, key=lambda k:k.start)
    for cand in sorted_candidates:
        if cand.len > 1000000:
            num -= 1
            continue
        out_writer.write(cand.to_string())

    print("Candidates found: ", num)

def run():
    parser = argparse.ArgumentParser(description='Cluster SV breakpoint intervals with hierarchical clustering')
    subparsers = parser.add_subparsers(title='commands', dest='command')  # two submodules

    yri_parser = subparsers.add_parser('yri', help='Evaluate both entire and sub events')
    yri_parser.add_argument("-b", dest="bed", help="SV callset in BED format")
    yri_parser.add_argument("-c", dest="chrom", help="Chromosome to cluster")
    yri_parser.add_argument("-w", dest="dir", help="Working directory for results output")

    skbr3_parser = subparsers.add_parser('skbr3', help='Compare events in two BED files')
    skbr3_parser.add_argument("-b", dest="bed", help="SV callset in BED format")
    skbr3_parser.add_argument("-c", dest="chrom", help="Chromosome to cluster")
    skbr3_parser.add_argument("-w", dest="dir", help="Working directory for results output")

    args = parser.parse_args()
    if args.command == "yri":
        cluster_cutoff = find_threshold(args.bed, args.chrom)
        yri_bp_clustering(args.bed, args.chrom, cluster_cutoff, args.dir)
    if args.command == "skbr3":
        cluster_cutoff = find_threshold(args.bed, args.chrom)
        skbr3_bp_clustering(args.bed, args.chrom, cluster_cutoff, args.dir)

script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('clusterBp.py         Last Update:2019-12-03\n')
    print('Cluster SV breakpoint intervals with hierarchical clustering\n')
    print('Usage:')
    print('clusterBp.py [commands] <parameters>\n')
    print('Commands:')
    print('yri: cluster NA19240 SV callset from HGSVC')
    print('skbr3: cluster SK-BR-3 SV callset')
    print("=======================================================")
else:
    run()