import re
import pysam


class SimplifiedAlign:
    def __init__(self, chr, ref_start, ref_end, read_start, read_end, is_reverse, simplified_seq, simplified_cigar, align_type):
        self.chr = chr
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.is_reverse = is_reverse
        self.simplified_seq = simplified_seq
        self.simplified_cigar = simplified_cigar
        self.align_type = align_type

    def to_string(self):
        return '{0} {1} {2} {3} {4} {5} {6} {7}'.format(self.chr, self.ref_start, self.ref_end, self.read_start, self.read_end, self.is_reverse, self.align_type, self.simplified_cigar)


def get_reverse_complement(seq):
    inv_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        bp = seq[i]
        inv_bp = ''
        if bp == 'A':
            inv_bp = 'T'
        elif bp == 'T':
            inv_bp = 'A'
        elif bp == 'C':
            inv_bp = 'G'
        elif bp == 'G':
            inv_bp = 'C'
        else:
            inv_bp = 'N'

        inv_seq += inv_bp

    return inv_seq

def fetch_ref_seq(ref_path, chr, start, end):
    ref = pysam.FastaFile(ref_path)

    ref_cutted = ref.fetch(chr, start, end)
    return ref_cutted



def list_to_cigar(ops, lengths):
    if len(ops) != len(lengths):
        # print(ops, lengths, len(ops), len(lengths))
        print('[ERROR]: length of ops and its length are not equal')
        exit()

    cigar_string = ''
    for i in range(len(ops)):
        cigar_string += "{0}{1}".format(lengths[i], ops[i])

    return cigar_string
def cigar_to_list(cigar):

    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    return ops, lengths


def simplify_cigar(cigarstring, read_seq, ref_seq, min_sv_size, ref_start):

    # I_and_Ds = []

    simplified_seq = ''
    new_cigar = ''
    ops, lengths = cigar_to_list(cigarstring)
    read_pos = 0
    ref_pos = 0

    for i in range(len(ops)):
        op = ops[i]
        op_len = lengths[i]


        if op == "I":
            if op_len >= min_sv_size:

                new_cigar += '{0}I'.format(op_len)

                simplified_seq += read_seq[read_pos: read_pos + op_len]
            else:
                pass
            read_pos += op_len

        elif op == "D":
            if op_len >= min_sv_size:

                new_cigar += '{0}D'.format(op_len)
            else:
                simplified_seq += ref_seq[ref_pos: ref_pos + op_len]
                new_cigar += '{0}M'.format(op_len)

            ref_pos += op_len

        elif op in ["M", "X", "E"]:
            simplified_seq += ref_seq[ref_pos: ref_pos + op_len]
            ref_pos += op_len
            read_pos += op_len

            new_cigar += '{0}M'.format(op_len)

        else:
            pass

    simplified_cigar = ''


    new_ops, new_lengths = cigar_to_list(new_cigar)

    previous_op = new_ops[0]
    previous_length = new_lengths[0]
    for i in range(1, len(new_ops)):
        if new_ops[i] == previous_op:
            previous_length += new_lengths[i]
        else:
            simplified_cigar += '{0}{1}'.format(previous_length, previous_op)
            previous_op = new_ops[i]
            previous_length = new_lengths[i]
    simplified_cigar += '{0}{1}'.format(previous_length, previous_op)
    return simplified_seq, simplified_cigar


def run(pm_align, supp_aligns, whole_read_seq, bam_file, ref_path, min_sv_size):
    simplied_aligns = []

    pm_chr = bam_file.getrname(pm_align.reference_id)
    pm_is_reverse = pm_align.is_reverse

    # # simplify primary's seq and cigar
    pm_read_seq = whole_read_seq[pm_align.query_alignment_start: pm_align.query_alignment_end]
    pm_ref_seq = fetch_ref_seq(ref_path, pm_chr, pm_align.reference_start, pm_align.reference_end)
    pm_simplified_read_seq, pm_simplified_cigar = simplify_cigar(pm_align.cigarstring, pm_read_seq, pm_ref_seq, min_sv_size, pm_align.reference_start)

    # add to list
    simplied_aligns.append(SimplifiedAlign(pm_align.reference_id, pm_align.reference_start, pm_align.reference_end,
                                           pm_align.query_alignment_start, pm_align.query_alignment_end, False,
                                           pm_simplified_read_seq, pm_simplified_cigar, 'PM'))

    # ref_seq = fetch_ref_seq(ref_path, pm_chr, 2232676, 2244552)
    # print('> ref')
    # print(ref_seq)
    # print('> seq')
    # print(pm_align.query_sequence[24785692: 25416958])
    # print(pm_align.reference_id, pm_align.reference_start, pm_align.reference_end, pm_align.query_alignment_start, pm_align.query_alignment_end, pm_align.is_reverse,
    #       pm_align.is_supplementary)
    # print(pm_align.cigarstring)
    # # process supplementary aligns
    for sa in supp_aligns:
        sa_chr = bam_file.getrname(sa.reference_id)

        # # focus on primary's forward, if other reads' forward is not same as primary, then adjust cords
        if sa.is_reverse != pm_is_reverse:
            sa_q_start = sa.query_length - sa.query_alignment_end
            sa_q_end = sa.query_length - sa.query_alignment_start
            sa_is_reverse = True
        else:
            sa_q_start = sa.query_alignment_start
            sa_q_end = sa.query_alignment_end
            sa_is_reverse = False

        # print(sa.reference_id, sa.reference_start, sa.reference_end, sa_q_start, sa_q_end, sa.is_reverse,
        #       sa.is_supplementary, sa.query_alignment_start, sa.query_alignment_end)

        # fetch sa's read and ref seq, if it is a reverse seq, then we get the reverse complement seq
        sa_read_seq = whole_read_seq[sa_q_start: sa_q_end]
        if sa_is_reverse == True:
            sa_read_seq = get_reverse_complement(sa_read_seq)
        sa_ref_seq = fetch_ref_seq(ref_path, sa_chr, sa.reference_start, sa.reference_end)


        # simplify sa's seq and cigar
        sa_simplified_read_seq, sa_simplified_cigar = simplify_cigar(sa.cigarstring, sa_read_seq, sa_ref_seq,
                                                                     min_sv_size, sa.reference_start)
        if sa_is_reverse == True:
            sa_simplified_read_seq = get_reverse_complement(sa_simplified_read_seq)
        # add to list
        simplied_aligns.append(SimplifiedAlign(sa.reference_id, sa.reference_start, sa.reference_end, sa_q_start, sa_q_end, sa_is_reverse, sa_simplified_read_seq, sa_simplified_cigar, 'SA'))
        # print(sa_simplified_cigar)
    return simplied_aligns