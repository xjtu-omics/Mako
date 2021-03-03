from aligner.segment import Segment


class HashAligner():

    def __init__(self, k, windowSize, mismatchNum, repeat_thresh):
        self.k = k
        self.windowSize = windowSize
        self.mismatchNum = mismatchNum

        self.segments = []
        self.selfDiffSegs = []
        self.compareDiffSegs = []
        self.N = 'N'
        self.G = 'G'
        self.A = 'A'
        self.T = 'T'
        self.C = 'C'
        self.avoid_kmers = []
        self.repeat_thresh = repeat_thresh
    def getSegments(self):
        return self.segments

    def getSelfDiffSegs(self):
        return self.selfDiffSegs

    def run(self, x, y, compareDiffSegs=None, y_hashvalue=None, avoid_kmers_from_ref=None):
        self.ref_length = len(y.get_bases())
        self.compareDiffSegs = compareDiffSegs
        self.y_hashvalues = y_hashvalue

        self.makePairwiseAlignment(x, y, avoid_kmers_from_ref)

    def extendKmersForward(self, xBases, yBases, matchPositions, p, i, segId):
        matchLength = self.k
        mismatch = 0

        # print(start, end)
        while mismatch <= self.mismatchNum:

            # beyond the length of x
            if matchPositions[p] + matchLength >= len(xBases) - 1:
                break
            # beyond the length of y
            if i + matchLength >= len(yBases) - 1:
                break

            xBase = xBases[matchPositions[p] + matchLength]
            yBase = yBases[i + matchLength]

            # stop extension when meet 'N'
            if xBase == self.N or yBase == self.N:
                break

            # stop extension when different
            if xBase != yBase:
                mismatch += 1

            matchLength += 1

        # when longer than windowSize
        if matchLength >= self.windowSize:

            d = Segment(matchPositions[p], i, matchLength, True, segId)

            if self.compareDiffSegs == None:
                self.segments.append(d)
                # if seq is ref(compareDiff = None), then calculate the ref's diff list
                if self.calDiffForRef(d):
                    self.selfDiffSegs.append(d)
            else:
                # all seg is diff

                if self.compareWithDiffSegs(d) is False:
                    # if d.get_xStart() == 0 and d.get_yStart() == 0:
                    #     print(d.to_string())
                    self.segments.append(d)

    def extendKmersReverse(self, reverseXbases, yBases, matchPosition, i, segId):
        matchLength = self.k
        mismatch = 0
        while mismatch <= self.mismatchNum:
            if matchPosition + matchLength >= len(reverseXbases) - 1:
                break
            if i + matchLength >= len(yBases) - 1:
                break
            xBase = reverseXbases[matchPosition + matchLength]
            yBase = yBases[i + matchLength]

            if xBase == self.N or yBase == self.N:
                break

            # stop extension when different
            if xBase != yBase:
                mismatch += 1

            matchLength += 1
        if matchLength >= self.windowSize:
            d = Segment((len(reverseXbases) - 1) - matchPosition, i, matchLength, False, segId)

            if self.compareDiffSegs == None:
                self.segments.append(d)

                # if seq is ref(compareDiff = None), then calculate the ref's diff list
                # if not (d.xStart() == 0 and d.yStart() == 0):
                #     self.selfDiffSegs.append(d)

                if self.calDiffForRef(d):
                    self.selfDiffSegs.append(d)
            else:
                if self.compareWithDiffSegs(d) is False:
                    self.segments.append(d)
                # if self.calDiff(d) or (abs(d.yEnd() - d.yStart()) < y.length() * 0.1):
                #     if self.compareWithDiffSegs(d) is False:
                #         self.segments.append(d)
                #
                # else:
                #     self.segments.append(d)

    def calHash(self, kmer):

        # hashValue = ''.join([str(bp) for bp in kmer])
        hashValue = kmer
        # hashValue = 0
        # for i in range(len(kmer)):
        #     hashValue += (int(kmer[9 - i] * math.pow(self.k, i)))
        #
        # hashValue = 0
        # power = 1
        # for d in range(9, -1, -1):
        #     if kmer[d] == self.G:
        #         hashValue += 0 * power
        #     elif kmer[d] == self.A:
        #         hashValue += 1 * power
        #     elif kmer[d] == self.T:
        #         hashValue += 2 * power
        #     elif kmer[d] == self.C:
        #         hashValue += 3 * power
        #
        #     power = power * 4
        return hashValue

    def makePairwiseAlignment(self, x, y, avoid_kmers_from_ref):

        xBases = x.get_bases()
        reverseXbases = x.get_reverse_complement()
        # store hash value, format: hashvalue: [pos1, po2, ....]
        hashedPositions = {}

        # for sequence x, make hash
        for i in range(0, len(xBases) - (self.k + 1)):
            kmer = xBases[i: i + self.k]
            hashValue = self.calHash(kmer)

            if hashValue not in hashedPositions.keys():
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(i)


        # get reverse seq, and make hash table
        for i in range(0, len(reverseXbases) - (self.k + 1)):
            kmer = reverseXbases[i: i + self.k]
            hashValue = self.calHash(kmer)

            if hashValue not in hashedPositions.keys():
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(-1 - i)

        yBases = y.get_bases()
        segId = 0

        if self.y_hashvalues == None:

            self.hashvalues = []

            for i in range(0, len(yBases) - (self.k + 1)):
                kmer = yBases[i: i + self.k]
                hashValue = self.calHash(kmer)

                self.hashvalues.append(hashValue)

                if hashValue in hashedPositions.keys():
                    matchPositions = hashedPositions[hashValue]
                    # for each possible hit position
                    if len(matchPositions) >= self.repeat_thresh:
                        self.avoid_kmers.append(hashValue)
                    else:
                        for p in range(0, len(matchPositions)):

                            # hit in the same strand
                            # Kmer extension
                            if matchPositions[p] >= 0:
                                # Already matched in previous Kmer
                                if matchPositions[p] > 0 and i > 0 and xBases[matchPositions[p] - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersForward(xBases, yBases, matchPositions, p, i, segId)
                                segId += 1
                            else:
                                matchPosition = -1 - matchPositions[p]

                                if matchPosition > 0 and i > 0 and reverseXbases[matchPosition - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersReverse(reverseXbases, yBases, matchPosition, i, segId)
                                segId += 1
        else:
            for i in range(len(self.y_hashvalues)):
                hashValue = self.y_hashvalues[i]

                if hashValue in hashedPositions.keys():

                    if hashValue not in avoid_kmers_from_ref:

                        matchPositions = hashedPositions[hashValue]

                        # for each possible hit position
                        for p in range(0, len(matchPositions)):

                            # hit in the same strand
                            # Kmer extension
                            if matchPositions[p] >= 0:
                                # Already matched in previous Kmer
                                if matchPositions[p] > 0 and i > 0 and xBases[matchPositions[p] - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersForward(xBases, yBases, matchPositions, p, i, segId)
                                segId += 1
                            else:
                                matchPosition = -1 - matchPositions[p]

                                if matchPosition > 0 and i > 0 and reverseXbases[matchPosition - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersReverse(reverseXbases, yBases, matchPosition, i, segId)
                                segId += 1


    def getMergeSegments(self):

        # # return self.segments
        # def is_merge(base_seg, target_seg):
        #     if base_seg.seg_forward != target_seg.seg_forward:
        #         return False
        #
        #     dis_on_ref = target_seg.y_start - base_seg.y_end
        #
        #     if base_seg.seg_forward == True:
        #         dis_on_read = target_seg.x_start - base_seg.x_end
        #     else:
        #         dis_on_read = base_seg.x_end - target_seg.x_start
        #
        #     if abs(dis_on_read - dis_on_ref) >= 50:
        #         return False
        #     else:
        #         return True
        #     # if abs(dis_on_ref) < 10 and abs(dis_on_read) < 10:
        #     #     return True
        #     # else:
        #     #     return False
        #
        # merged_segments = []
        # segments_sorted_by_read = sorted(self.segments, key=lambda aln: (aln.x_start, aln.x_end))
        # merged_segments.append(segments_sorted_by_read[0])
        #
        # for i in range(1, len(segments_sorted_by_read)):
        #     base_segment = merged_segments[-1]
        #     target_segment = segments_sorted_by_read[i]
        #
        #     # the two seg can merge, then adjust base's cord
        #     if is_merge(base_segment, target_segment):
        #         if base_segment.seg_forward == True:
        #             base_segment.set_x_end(target_segment.x_end)
        #             base_segment.set_y_end(target_segment.y_end)
        #         else:
        #             base_segment.set_x_start(target_segment.x_start)
        #
        #             base_segment.set_y_end(target_segment.y_end)
        #
        #
        #     else:
        #         merged_segments.append(target_segment)
        # return merged_segments

        ################################################################################
        # main_segs = []
        # for seg in self.segments:
        #     if (seg.yStart() - 0) <= 5 and seg.forward() == True:
        #         main_segs.append(seg)
        #     elif abs(seg.yEnd() - self.ref_length) <= 5 and seg.forward() == True:
        #         main_segs.append(seg)
        #
        # for main_seg in main_segs:
        #     self.segments.remove(main_seg)

        curSegNum = 1

        while curSegNum < len(self.segments):
            flag = 0
            curSeg = self.segments[curSegNum]
            for i in range(curSegNum):
                candiSeg = self.segments[i]

                if self.linearOrNot(candiSeg, curSeg):
                    # print(candiSeg.xStart(), candiSeg.xEnd(), candiSeg.yStart(), candiSeg.yEnd(), curSeg.xEnd(), curSeg.yEnd())
                    # print('merge')
                    # merge and update segments list
                    if curSeg.seg_forward == True:
                        candiSeg.set_x_end(max(curSeg.x_end, candiSeg.x_start))
                    elif curSeg.seg_forward == False:
                        candiSeg.set_x_end(min(curSeg.x_end, candiSeg.x_start))
                    candiSeg.set_y_end(max(curSeg.y_end, candiSeg.y_end))

                    candiSeg.set_length(abs(candiSeg.seg_length) + abs(curSeg.x_end - candiSeg.x_end))

                    # remove curSeg
                    self.segments.remove(curSeg)
                    flag = 1
                    break

            if flag == 0:
                curSegNum += 1

        # max_length = 0
        # for seg in self.segments:
        #     if seg.length() > max_length:
        #         max_seg = seg
        #         max_length = seg.length()

        # self.segments.extend(self.main_segs)
        after_filte_segs = []
        for seg in self.segments:
            if (seg.y_end - seg.y_start) >= 20:
                # print(seg.yEnd() - seg.yStart())
                after_filte_segs.append(seg)
        self.segments = after_filte_segs
        return self.segments

    def linearOrNot(self, i, j):

        # merge condition1: different strand, pass directly
        if i.seg_forward != j.seg_forward:
            return False

        # if i.forward() is not False:
        #     return

        # merge condition2: colinear, not then false
        DIFF = self.calDiffBetTow(i, j)
        # print(DIFF)

        if DIFF > 1.1 or DIFF < 0.9:
            return False

        # merge condition3: colinear, yes, and dis
        DIS_X = abs(i.x_end - j.x_start)
        DIS_Y = abs(i.y_end - j.y_start)
        # maxDis = (i.length() + j.length()) * 1.5
        maxDis = 50

        if DIS_X > maxDis and DIS_Y > maxDis:
            return False

        # merge condition4: merged line's k is not -1 or 1
        tmp = float(j.x_end - i.x_start)
        if tmp == 0:
            tmp = 0.0001
        k = float(j.y_end - i.y_start) / tmp

        if abs((abs(k) - 1)) > 0.2:
            return False

        return True


    def compareWithDiffSegs(self, i):

        refStart = i.y_start
        refEnd = i.y_end
        forward = i.seg_forward


        for tmpSeg in self.compareDiffSegs:
            # if tmpSeg.forward() != forward:
            #     continue

            startDis = abs(refStart - tmpSeg.y_start)
            endDis = abs(refEnd - tmpSeg.y_end)

            if (startDis <= 5 and refEnd <= tmpSeg.y_end) \
                    or (endDis <= 5 and refStart >= tmpSeg.y_start):
                return True

        return False

    def calDiffForRef(self, i):
        # diff about the start
        diff2 = float(i.x_end) / float(i.y_end)

        # diff about the center point
        centryX = float(i.x_start + i.x_end) / 2.0
        centryY = float(i.y_start + i.y_end) / 2.0
        diff3 = centryX / centryY

        if diff2 != 1 or diff3 != 1:
            return True
        else:
            return False

    def calDiff(self, i):
        diff1 = float(i.x_start + 1) / float(i.y_start + 1)
        diff2 = float(i.x_end) / float(i.y_end)

        centryX = float(i.x_start + i.x_end) / 2.0
        centryY = float(i.y_start + i.y_end) / 2.0

        diff3 = centryX / centryY

        thresh1 = 0.95
        thresh2 = 1.05

        if i.seg_forward is True and diff1 == 1.0 and diff2 == 1.0 and diff3 == 1.0:
            return True
        else:
            return False
        # if i.forward() is False and diff3 <= thresh2 and diff3 >= thresh1:
        #     return True
        #
        # if (diff1 > thresh2 or diff1 < thresh1) or (diff2 > thresh2 or diff2 < thresh1):
        #     return True
        #
        # else:
        #     return False

    def calDiffBetTow(self, i, j):
        if abs(float(i.y_start - j.y_start)) == 0:
            return 5
        return abs(float(i.x_start - j.x_start)) / abs(float(i.y_start - j.y_start))

    def getHashValues(self):
        return self.hashvalues
    def getAvoidKmer(self):
        return self.avoid_kmers