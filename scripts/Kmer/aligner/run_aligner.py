from aligner.sequence import Sequence
from aligner.segment import Segment
from aligner.hash_aligner import HashAligner



def run(ref, seq):
    k = 32
    window_size = 20
    repeat_thresh = 10

    ref = Sequence(ref)
    read = Sequence(seq)

    aligner_ref = HashAligner(k, window_size,  0, repeat_thresh)
    aligner_ref.run(ref, ref)
    ref_segs = aligner_ref.getMergeSegments()
    diff_segs = aligner_ref.getSelfDiffSegs()
    avoid_mers = aligner_ref.getAvoidKmer()

    y_hashvalues = aligner_ref.getHashValues()
    # y_hashvalues = None

    aligner_merge = HashAligner(k, window_size, 0, repeat_thresh)
    aligner_merge.run(read, ref, diff_segs, y_hashvalues, avoid_mers)
    segments_merge = aligner_merge.getMergeSegments()

    return segments_merge, ref_segs


