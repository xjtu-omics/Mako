import pysam
from consensus import traverse_reads
from aligner import run_aligner
from consensus.simplify_aligns import fetch_ref_seq
from consensus.traverse_reads import create_align


if __name__ == '__main__':

    bam_paths = ['/mnt/d/Data/Bams/TGS/CCS_HAP_HHU_ASSM_v12/HG00514-hap1.bam', '/mnt/d/Data/Bams/TGS/CCS_HAP_HHU_ASSM_v12/HG00514-hap2.bam']

    ref_path = '/mnt/d/Data/Bams/REF/GRCh38/GRCh38_chr1-x.fa'
    bed_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/mako/valid/HG00514/HG00514.mako_csvs.tsv'


    for hapltype in range(2):
        # hapltype = 1
        out_file = open('/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/mako/valid/HG00514/HG00514.mako_csvs.tsv.hap{0}'.format(hapltype + 1), 'w')
        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('CHROM', "OriginStart", "OriginEnd", "ValidStart", "ValidEnd", 'Note'))

        bam_path = bam_paths[hapltype]
        bam_file = pysam.AlignmentFile(bam_path)

        # chroms = ['chr1']
        chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX']
        for spec_chrom in chroms:

            bam_file = pysam.AlignmentFile(bam_path)
            all_possible_chrs = pysam.FastaFile(ref_path).references

            spec_aligns = bam_file.fetch(spec_chrom)

            reads_dict = {}
            for align in spec_aligns:
                # # no cigar, then pass this align
                if align.cigarstring == None:
                    continue
                # # unmapped or secondary or low mapping quality, then pass this align
                if align.is_unmapped or align.is_secondary or align.mapping_quality < 20:
                    continue

                # # align to a ref that not in genome reference
                align_chr = align.reference_name
                if align_chr not in all_possible_chrs:
                    print("[Warning]: '{0}' not in reference's .fa file, skip this read".format(align_chr))
                    continue

                align_ref_id = bam_file.get_tid(align.reference_name)

                # create new align
                new_align = create_align(align_ref_id, align)

                if align.qname not in reads_dict.keys():
                    reads_dict[align.qname] = [new_align]
                else:
                    reads_dict[align.qname].append(new_align)

            for line in open(bed_path).readlines():
                line_split = line.strip().split('\t')
                spec_chr = line_split[0]
                spec_start = int(line_split[1])
                spec_end = int(line_split[2])
                if spec_chr != spec_chrom:
                    continue

                # spec_chr = 'chr1'
                # spec_start = 235078305
                # spec_end = 235078561

                spec_region_str = '{0}-{1}-{2}-hap{3}'.format(spec_chr, spec_start, spec_end, hapltype + 1)
                print('------------------------------------------------', spec_region_str)

                extend_length = 500
                spec_start -= extend_length
                spec_end += extend_length


                candidate_seqs = traverse_reads.run(reads_dict, bam_file, ref_path, spec_chr, spec_start, spec_end)

                if len(candidate_seqs) == 0:
                    print('[ERROR]: No reads in that region')
                    out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, -1, -1, 'No reads in that region'))

                    continue

                # if len(candidate_seqs) != 1:
                #     print('[ERROR]: More than one candidate seq')
                #     out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, -1, -1,'More than one candidate seq'))
                #
                #     continue

                valid_flag = 0
                for read_seq in candidate_seqs:
                    if len(read_seq) == 0:
                        continue

                    valid_flag = 1
                    ref_seq = fetch_ref_seq(ref_path, spec_chr, spec_start, spec_end)

                    ref_length = len(ref_seq)
                    read_length = len(read_seq)


                    # run hash lineplot to draw dotplots
                    segments, ref_segs = run_aligner.run(ref_seq, read_seq)
                    segments_sorted_by_read = sorted(segments, key=lambda aln: (aln.x_start, aln.x_end))

                    # for i in range(len(segments_sorted_by_read)):
                    #     print(segments_sorted_by_read[i].to_string())

                    if len(segments_sorted_by_read) == 0 or len(segments_sorted_by_read) == 1:
                        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, -1, -1, 'No obvious SV'))

                        print('[ERROR]: No obvious SV')
                    elif len(segments_sorted_by_read) == 2:

                        valid_start = spec_start + segments_sorted_by_read[0].y_end
                        valid_end = spec_start + segments_sorted_by_read[1].y_start

                        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, valid_start, valid_end, 'Simple Event'))

                        print(len(segments_sorted_by_read), valid_start, valid_end)
                    else:
                        valid_start = spec_start + min([s.y_start for s in segments_sorted_by_read[1: -1]])
                        valid_end = spec_start + max([s.y_end for s in segments_sorted_by_read[1: -1]])

                        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, valid_start, valid_end, 'Valid'))

                        print(len(segments_sorted_by_read), valid_start, valid_end)

                if valid_flag == 0:
                    print('[ERROR]: All read seq is empty')
                    out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(spec_chr, spec_start, spec_end, -1, -1, 'All read seq is empty'))

                    continue