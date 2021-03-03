import pysam
import re
from consensus import simplify_aligns, joint_aligns
from consensus.simplify_aligns import cigar_to_list

min_sv_size = 20


def fetch_read_seq(seq, cigar, ref_start, target_ref_start, target_ref_end):
    ops, lengths = cigar_to_list(cigar)
    read_pos = 0
    ref_pos = ref_start

    read_cut_start = -1
    read_cut_end = -1

    for i in range(len(ops)):


        op = ops[i]
        op_len = lengths[i]

        if op == "I":
            read_pos += op_len

        elif op == "D":
            if ref_pos + op_len > target_ref_start and read_cut_start == -1:
                read_cut_start = read_pos
                # print(0, op, op_len)

                if read_cut_start != -1 and read_cut_end != -1:
                    break
            if ref_pos + op_len > target_ref_end and read_cut_end == -1:
                read_cut_end = read_pos
                # print(1, op, op_len)
                if read_cut_start != -1 and read_cut_end != -1:
                    break

            ref_pos += op_len

        elif op == 'M':
            if ref_pos + op_len > target_ref_start and read_cut_start == -1:
                read_cut_start = read_pos + (target_ref_start - ref_pos)
                # print(0, op, op_len)

                if read_cut_start != -1 and read_cut_end != -1:
                    break
            if ref_pos + op_len > target_ref_end and read_cut_end == -1:
                # print(1, op, op_len)

                read_cut_end = read_pos + (target_ref_end - ref_pos)
                if read_cut_start != -1 and read_cut_end != -1:
                    break

            ref_pos += op_len
            read_pos += op_len
        elif op in ["X", "E"]:
            pass
        else:
            pass
    # print(read_cut_start, read_cut_end, ref_pos, read_pos)
    return seq[read_cut_start: read_cut_end]


def create_align(ref_id, old_align,):
    nm = old_align.get_tag('NM')


    new_align = pysam.AlignedSegment()

    new_align.reference_id = ref_id
    new_align.reference_start = old_align.reference_start

    new_align.query_name = old_align.query_name
    #

    # if not old_align.is_reverse:
    #     new_align.flag = 2048
    # else:
    #     new_align.flag = 2064

    new_align.is_supplementary = old_align.is_supplementary
    new_align.is_reverse = old_align.is_reverse

    if not new_align.is_supplementary:
        new_align.query_sequence = old_align.query_sequence

    try:
        new_align.mapping_quality = old_align.mapq
    except OverflowError:
        new_align.mapping_quality = 0

    new_align.cigarstring = str(old_align.cigarstring).replace('H', "S")
    new_align.next_reference_id = -1
    new_align.next_reference_start = -1
    new_align.template_length = 0
    # new_align.query_qualities = old_align.query_qualities
    new_align.set_tags([("NM", nm, "i")])


    return new_align

def run(reads_dict, bam_file, ref_path, spec_chrom, spec_start, spec_end):



    candidate_seqs = []
    for qname in reads_dict.keys():

        # if qname != 'cluster24_000001F':
        #     continue

        # separate pm and sa
        pm_align = None
        supp_aligns = []
        for align in reads_dict[qname]:
            if not align.is_supplementary:
                pm_align = align
            else:
                supp_aligns.append(align)

        if pm_align == None:
            continue
        # set sa's query sequence to pm align's query sequence
        for sa in supp_aligns:
            sa.query_sequence = pm_align.query_sequence

        whole_read_seq = pm_align.query_sequence

        # get simplified aligns
        simplied_aligns = simplify_aligns.run(pm_align, supp_aligns, whole_read_seq, bam_file, ref_path, min_sv_size)
        simplied_aligns_srt_by_read = sorted(simplied_aligns, key=lambda aln: (aln.read_start, aln.read_end))


        covered_aligns = []
        for align in simplied_aligns_srt_by_read:
            if spec_start < align.ref_start < spec_end or spec_start < align.ref_end < spec_end \
                    or align.ref_start < spec_start < align.ref_end or align.ref_start < spec_end < align.ref_end:
                covered_aligns.append(align)
                # print(align.to_string())

        if len(covered_aligns) == 0:
            continue
        elif len(covered_aligns) == 1:
            # print('--------------------------------', pm_align.qname, spec_start, spec_end)

            align = covered_aligns[0]
            seq = fetch_read_seq(align.simplified_seq, align.simplified_cigar, align.ref_start, spec_start, spec_end)
            candidate_seqs.append(seq)

        # more than one aligns, then we need joint them
        elif len(covered_aligns) > 1:
            # print('--------------------------------', pm_align.qname, spec_start, spec_end)

            jointed_align = joint_aligns.run(covered_aligns, bam_file.getrname(pm_align.reference_id), pm_align.qname, whole_read_seq)
            seq = fetch_read_seq(jointed_align.seq, jointed_align.cigar, jointed_align.ref_start, spec_start, spec_end)
            candidate_seqs.append(seq)
        else:
            pass

    return candidate_seqs


















    # # simplify cigar test code
    # bam_path = '/mnt/d/Data/Bams/TGS/CCS/NA19240/ngmlr/chr20-2240289-2242427.bam'
    # ref_path = '/mnt/d/Data/Bams/REF/GRCh38/GRCh38_chr1-x.fa'
    #
    # bam = pysam.AlignmentFile(bam_path)
    #
    # out_bam = pysam.AlignmentFile('/mnt/d/Data/Bams/TGS/CCS/NA19240/ngmlr/chr20-2240289-2242427.simplified.bam', "wb", header=bam.header)
    #
    # # print(bam.header)
    # for align in bam:
    #
    #     ref_name = align.reference_name
    #
    #     ref_start = align.reference_start
    #     ref_end = align.reference_end
    #     ref_seq = fetch_ref_seq(ref_path, ref_name, ref_start, ref_end)
    #     read_seq = align.query
    #
    #     cigarstring = align.cigarstring
    #
    #     print(align.qname, ref_start)
    #
    #     simplified_seq, simplified_cigar = simplify_cigar(cigarstring, read_seq, ref_seq, 50)


        # # write new seq to bam file
        # a = pysam.AlignedSegment()  # 定义一个AlignedSegment对象用于存储比对信息
        # a.query_name = align.qname
        # a.query_sequence = simplified_seq
        # a.flag = align.flag
        # a.reference_id = 19
        # a.cigarstring = simplified_cigar
        # a.reference_start = align.reference_start
        # a.mapping_quality = 20
        #
        # out_bam.write(a)
    #
