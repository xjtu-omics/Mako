
import pysam
from consensus.simplify_aligns import cigar_to_list, list_to_cigar

class JointAlign:
    def __init__(self, qname, chr, ref_start, ref_end, is_reverse, seq, cigar, align_type):
        self.read_name = qname
        self.chr = chr
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.is_reverse = is_reverse
        self.seq = seq
        self.cigar = cigar
        self.align_type = align_type

    def to_string(self):
        return '{0} {1} {2} {3} {4} {5}'.format(self.chr, self.ref_start, self.ref_end, self.is_reverse, self.cigar, self.align_type)



def annealing(gap_on_read, anneal_length, cur_align):
    cur_align.read_end = cur_align.read_end - anneal_length - 1

    # anneal length is longer than cur align's lenght, then we add nothing to seq and cigar
    if anneal_length >= len(cur_align.simplified_seq):
        return "", ""
    else:
        old_cigar = cur_align.simplified_cigar
        old_op, old_lenght = cigar_to_list(old_cigar)

        # traverse cigar from the last op
        index = len(old_op) - 1
        while anneal_length > 0 and index >= 0:
            # this op's lenght is longer than anneal length, then anneal finished, and anneal lengh is set to 0, besides, we reduce the op's length
            if old_lenght[index] > anneal_length:
                old_lenght[index] = old_lenght[index] - anneal_length
                if old_op[index] == 'M':
                    cur_align.ref_end = cur_align.ref_end - anneal_length
                anneal_length = 0
            # this op's lenght is smaller than anneal length, then we update anneal length and remove the ops from cigar
            else:
                anneal_length = anneal_length - old_lenght[index]

                if old_op[index] == 'M':
                    cur_align.ref_end = cur_align.ref_end - old_lenght[index]

                old_op = old_op[: index]
                old_lenght = old_lenght[: index]

            index -= 1

        new_cigar = list_to_cigar(old_op, old_lenght)
        cur_align.simplified_cigar = new_cigar

        return cur_align.simplified_seq[: gap_on_read], new_cigar

def run(simplied_aligns_srt_by_read, chr, qname, whole_read_seq):

    # for align in simplied_aligns_srt_by_read:
    #     print(align.to_string())

    # # merge all simplified_aligns to get final-whole-seq
    final_whole_seq = ''
    final_cigar = ''

    final_ref_start = simplied_aligns_srt_by_read[0].ref_start
    final_ref_end = simplied_aligns_srt_by_read[0].ref_end
    for i in range(len(simplied_aligns_srt_by_read) - 1):
        cur_align = simplied_aligns_srt_by_read[i]
        next_align = simplied_aligns_srt_by_read[i + 1]
        # print(next_align.simplified_seq)

        # assign merged ref start and end
        final_ref_start = min(cur_align.ref_start, final_ref_start)
        final_ref_start = min(next_align.ref_start, final_ref_start)
        final_ref_end = max(cur_align.ref_end, final_ref_end)
        final_ref_end = max(next_align.ref_end, final_ref_end)

        # add unmapped seq between two mapped aligns
        gap_on_read = next_align.read_start - cur_align.read_end

        # there is no gap or overlap on read
        if gap_on_read == 0:
            final_whole_seq += cur_align.simplified_seq
            final_cigar += cur_align.simplified_cigar

        # there is a overlap on read
        elif gap_on_read < 0:

            print( cur_align.read_end, next_align.read_start, gap_on_read)
            # we perform annealing
            anneal_length = abs(gap_on_read)

            annealed_seq, annealed_cigar = annealing(gap_on_read, anneal_length, cur_align)

            # add annealed seq and cigar to final
            final_whole_seq += annealed_seq
            final_cigar += annealed_cigar

        elif gap_on_read > 0:
            gap_seq = whole_read_seq[cur_align.read_end: next_align.read_start]
            final_whole_seq += cur_align.simplified_seq + gap_seq
            final_cigar += cur_align.simplified_cigar + '{0}X'.format(gap_on_read)

        else:
            pass
    # # add the last mapped align to list
    final_whole_seq += simplied_aligns_srt_by_read[-1].simplified_seq
    final_cigar += simplied_aligns_srt_by_read[-1].simplified_cigar

    # # add the unmapped seq located after the last mapped align
    unmapped_seq = whole_read_seq[simplied_aligns_srt_by_read[-1].read_end: ]
    if len(unmapped_seq) != 0:
        final_whole_seq += unmapped_seq
        final_cigar += '{0}X'.format(len(unmapped_seq))

    return JointAlign(qname, chr, final_ref_start, final_ref_end, False, final_whole_seq, final_cigar, 'JT')


