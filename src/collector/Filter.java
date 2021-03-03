package collector;

import utils.FileReader;
import htsjdk.samtools.*;
import org.apache.commons.math3.distribution.PoissonDistribution;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * User: jiadonglin
 * Date: 2020/9/8
 */

public class Filter {

    private final int readDepthContainerBuffer = 1000000;

    private final int minMapQ;
    private int maxRpDist;

    public String[] chromNameMap;
    public int[] chromLengthMap;

    private final int isizeUpper;
    private final int isizeLower;
    private final int readLen;
    private final int insertMu;
    private final int insertStd;

    public Filter(int fragMean, int fragStd, int cutStd, int readLen, int maxRpDist, int minMapQ) throws IOException {
        isizeUpper = fragMean + cutStd * fragStd;
        isizeLower = fragMean - cutStd * fragStd;


        this.readLen = readLen;
        this.insertMu = fragMean;
//        this.insertCut = cutStd;
        this.insertStd = fragStd;
        this.minMapQ = minMapQ;
        this.maxRpDist = maxRpDist;

    }

    public Set<String> runFilter(String bamFile, String fastaIndexFile, String chrom, Map<String, Float> signalSummary) throws IOException{

        int windowStart = 0;
        int windowEnd = 0;

        readFaIdxFile(fastaIndexFile);

        System.out.println("\nStart filtering background noise from BAM ...\n");
        FileReader myFileReader = new FileReader();
        final SamReader samReader = myFileReader.openBamFile(bamFile, ValidationStringency.SILENT, false);

        int length = chromLengthMap.length;
        SAMRecordIterator iterator;

        Map<String, List<DiscordantRp>> discRps;
        Map<String, SAMRecord> rpTracker;
        Set<String> filteredDiscRpNames = new HashSet<>();
        Map<String, Integer> discRpSigTypeCount = new HashMap<>();

        int totalDiscRps = 0;

        // access user specified region
        if (chrom != null){

            SAMFileHeader samFileHeader = samReader.getFileHeader();

            SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();
            SAMSequenceRecord refSequenceRecord = sequenceDictionary.getSequence(chrom);
            int refSequenceLength = refSequenceRecord.getSequenceLength();
            int nWindows = refSequenceLength / readDepthContainerBuffer;

            for (int i = 0; i < nWindows; i++){
                discRps = new HashMap<>();
                rpTracker = new HashMap<>();

                windowEnd = windowStart + readDepthContainerBuffer;
                iterator = samReader.query(chrom, windowStart, windowEnd, false);

                analysisAlignment(iterator, rpTracker, discRps, discRpSigTypeCount);
                totalDiscRps += rpTracker.size();

                int failedRps = discRpTest(discRps, signalSummary, filteredDiscRpNames);

                System.out.println("[Filter] processed region: [" + windowStart + ", " + windowEnd + "] " + "# failed rps: " + failedRps);
                windowStart = windowEnd;

                discRps.clear();
                rpTracker.clear();

            }

            // process remaining alignment in BAM
            discRps = new HashMap<>();
            rpTracker = new HashMap<>();

            iterator = samReader.query(chrom, windowStart, refSequenceLength, false);
            analysisAlignment(iterator, rpTracker, discRps, discRpSigTypeCount);
            totalDiscRps += rpTracker.size();

            int failedRps = discRpTest(discRps, signalSummary, filteredDiscRpNames);

            System.out.println("[Filter] processed region: [" + windowStart + ", " + refSequenceLength + "] " + "# failed rps: " + failedRps);

        }
        else{
            for (int i = 0;i < length; i ++){
                discRps = new HashMap<>();
                rpTracker = new HashMap<>();

                windowStart = 0;

                int refSequenceLength = chromLengthMap[i];
                String curChrom = chromNameMap[i];
                System.out.println("[Filter] processing chrom: " + curChrom + ", chrom length: " + refSequenceLength);

                int nWindows = refSequenceLength / readDepthContainerBuffer;


                for (int k = 0; k < nWindows; k++){
                    windowEnd = windowStart + readDepthContainerBuffer;
                    iterator = samReader.query(curChrom, windowStart, windowEnd, false);

                    analysisAlignment(iterator, rpTracker, discRps, discRpSigTypeCount);

                    totalDiscRps += rpTracker.size();
                    int failedRps = discRpTest(discRps, signalSummary, filteredDiscRpNames);

                    System.out.println("[Filter] processed region: [" + windowStart + ", " + windowEnd + "] " + "# failed rps: " + failedRps);

                    windowStart = windowEnd;


                    discRps.clear();
                    rpTracker.clear();
                }

                if (windowStart < refSequenceLength) {

                    iterator = samReader.query(curChrom, windowStart, refSequenceLength, false);
                    analysisAlignment(iterator, rpTracker, discRps, discRpSigTypeCount);

                    totalDiscRps += rpTracker.size();
                    int failedRps = discRpTest(discRps, signalSummary, filteredDiscRpNames);

                    System.out.println("[Filter] processed region: [" + windowStart + ", " + refSequenceLength + "] " + "# failed rps: " + failedRps);


                    rpTracker.clear();
                }

            }
        }

        System.out.println("\n===== Filter results ====== ");
        System.out.println("Total discordant rps: " + totalDiscRps);
        float failedPcrt = (float) filteredDiscRpNames.size() * 100 / totalDiscRps;
        DecimalFormat f = new DecimalFormat("#.##");
        System.out.println("Filtered discordant rps: " + f.format(failedPcrt) + "%");
        System.out.println("=============================");

        return filteredDiscRpNames;
    }

    /**
     * Analysis each BAM record through different channels
     * @param iterator
     */
    private void analysisAlignment(SAMRecordIterator iterator, Map<String, SAMRecord> rpTracker,
                                   Map<String, List<DiscordantRp>> discRps, Map<String, Integer> discRpTypeCount){
//        CigarOps corasenCigar = new CigarOps();
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();
//            int recordChrIdx = record.getReferenceIndex();
//            String recordChrName = record.getReferenceName();

            String qName = record.getReadName();
            List<CigarElement> cigarElements = record.getCigar().getCigarElements();

            if (badReads(cigarElements, record.getMappingQuality()) || record.getDuplicateReadFlag()){
                continue;
            }

            if (!rpTracker.containsKey(qName)) {
                rpTracker.put(qName, record);
                continue;
            }

            SAMRecord mate = rpTracker.get(qName);

            if (!isBadRp(record, mate) && mate.getReferenceName().equals(record.getReferenceName())) {
                processDiscRps(record, mate, discRps, discRpTypeCount);
            }



        }
        iterator.close();
    }

    /**
     * Process discordant read-pairs within current window
     * @param record
     * @param mate
     * @param discRps
     */
    private void processDiscRps(SAMRecord record, SAMRecord mate, Map<String, List<DiscordantRp>> discRps,
                                Map<String, Integer> discRpTypeCount){

        SAMRecord first = record.getReadNegativeStrandFlag() ? mate : record;
        SAMRecord second = first == mate ? record : mate;

        String firstStrand = first.getReadNegativeStrandFlag() ? "-" : "+";
        String secondStrand = second.getReadNegativeStrandFlag() ? "-" : "+";

        String pairStrand = firstStrand + secondStrand;

        if (first.getReadNegativeStrandFlag() != second.getReadNegativeStrandFlag()) {
            if (first.getInferredInsertSize() > isizeUpper) {
                DiscordantRp signal = new DiscordantRp(first, "ARP_LARGE", "DEL", pairStrand, first.getAlignmentEnd(), second.getAlignmentStart());
                addDiscSignals(signal, discRps, discRpTypeCount);

            }
            else if (first.getInferredInsertSize() < isizeLower) {
                DiscordantRp signal = new DiscordantRp(first, "ARP_SMALL", "INS", pairStrand, first.getAlignmentEnd(), second.getAlignmentStart());
                addDiscSignals(signal, discRps, discRpTypeCount);

            }
            else if(first.getAlignmentStart() > second.getAlignmentStart() || first.getAlignmentEnd() > second.getAlignmentEnd()) {
                if (first.getReadNegativeStrandFlag()){
                    DiscordantRp signal = new DiscordantRp(first, "ARP_RF", "DUP", pairStrand, second.getAlignmentStart(), first.getAlignmentEnd());
                    addDiscSignals(signal, discRps, discRpTypeCount);

                }
            }
        }
        else{
            int pos1 = record.getAlignmentStart();
            int pos2 = mate.getAlignmentStart();

            if (pos1 > pos2) {
                int tmp = pos1;
                pos1 = pos2;
                pos2 = tmp;

            }
            if (record.getReadNegativeStrandFlag()) {
                DiscordantRp signal = new DiscordantRp(first, "ARP_RR", "INV", pairStrand, pos1, pos2);
                addDiscSignals(signal, discRps, discRpTypeCount);

            }else{
                DiscordantRp signal = new DiscordantRp(first, "ARP_FF", "INV", pairStrand, pos1, pos2);
                addDiscSignals(signal, discRps, discRpTypeCount);
            }

        }
    }
    /**
     * The main function of clustering discordant read pairs within current window and performing the statistical test to the
     * randomnized read pair distribution.
     * step1. Clustering discordant read pairs by potential breakpoint types;
     * step2. Compute a null distribution based on the clusters from step1;
     * step3. Filter clusters in step1 by comparing with the null distribution.
     *
     */

    private int discRpTest(Map<String, List<DiscordantRp>> discRps, Map<String, Float> bamSummary, Set<String> qnames) {

        int failedRpCount = 0;
        for (Map.Entry<String, List<DiscordantRp>> entry : discRps.entrySet()) {

            List<DiscordantRp> signals = entry.getValue();

            String sigType = entry.getKey();
            if (sigType.equals("ARP_LARGE")) {
                List<List<DiscordantRp>> clusters = clusterDiscSigs(signals);

                for (List<DiscordantRp> cluster : clusters) {

                    Collections.sort(cluster);
                    int regionSize = clusterSpanSize(cluster);

                    PoissonDistribution pd = new PoissonDistribution(regionSize * bamSummary.get(sigType));
                    double p = 1 - pd.cumulativeProbability(cluster.size());

                    // Failed the test, should be background noise
                    if (p > 0.00001) {
                        for (DiscordantRp rp : cluster) {
                            qnames.add(rp.getqName());
                            failedRpCount += 1;

                        }
                    }

                }
            }

        }

        return failedRpCount;
    }

    /**
     * Calculate a single cluster spanned genome size
     * @param signals
     * @return
     */
    private int clusterSpanSize(List<DiscordantRp> signals){
        int leftMin = Integer.MAX_VALUE;
        int rightMax = 0;

        for (DiscordantRp signal : signals) {
            if (signal.mutPos1 < leftMin) {
                leftMin = signal.mutPos1;
            }

            if (signal.mutPos2 > rightMax) {
                rightMax = signal.mutPos2;
            }
        }

        int regionSize = Math.abs(rightMax - leftMin);
        if (regionSize == 0) {
            regionSize = signals.get(0).getReadLen();
        }

        return regionSize;

    }

    /**
     * This function is used to cluster a list of discordant read pair mutation signals under rpClusteringMax distance constrain.
     * @param signals
     * @return
     */

    private List<List<DiscordantRp>> clusterDiscSigs(List<DiscordantRp> signals){


        if (signals.size() < 3) {
            return new ArrayList<>();
        }

        List<List<DiscordantRp>> clusters = new ArrayList<>(signals.size());

        Collections.sort(signals);
        Iterator<DiscordantRp> iter = signals.iterator();
        List<DiscordantRp> cluster = new ArrayList<>(signals.size());
        cluster.add(iter.next());

        while (iter.hasNext()) {
            DiscordantRp thisSignal = iter.next();
            List<DiscordantRp> pairs = new ArrayList<>(2);
            pairs.add(cluster.get(cluster.size() - 1));
            pairs.add(thisSignal);

            if (!isCompatible(pairs, maxRpDist)) {
                clusters.add(cluster);
                cluster = new ArrayList<>();
                cluster.add(thisSignal);
            }
            else {
                cluster.add(thisSignal);
            }
        }

        return clusters;
    }

    /**
     * Test if the new add signal is compatible with the current rp cluster
     * @param signals
     * @param maxClusDist
     * @return
     */
    private boolean isCompatible(List<DiscordantRp> signals, int maxClusDist){
        int minLeft = Integer.MAX_VALUE;
        int maxLeft = 0;
        int minRight = Integer.MAX_VALUE;
        int maxRight = 0;

        for (DiscordantRp sig : signals) {
            int sigLeft = sig.mutPos1;
            int sigRight = sig.mutPos2;
            if (sigLeft < minLeft) {
                minLeft = sigLeft;
            }
            else if (sigLeft > maxLeft) {
                maxLeft = sigLeft;
            }

            if (sigRight < minRight) {
                minRight = sigRight;
            }
            else if (sigRight > maxRight) {
                maxRight = sigRight;
            }
        }

        return Math.max(maxLeft - minLeft, maxRight - minRight) <= maxClusDist && maxLeft < minRight;
    }

    /**
     * Add discordant read pair induced signals
     * @param signal
     */
    public void addDiscSignals(DiscordantRp signal, Map<String, List<DiscordantRp>> discRps, Map<String, Integer> discRpTypeCount){

        // Abnormal read-pairs
        if (signal.insertSize <= readDepthContainerBuffer){
//            String bpType = signal.bpType;
            String signalType = signal.signalType;
            List<DiscordantRp> thisTypeSignals = discRps.get(signalType);
            if (thisTypeSignals == null) {
                thisTypeSignals = new ArrayList<>();
                thisTypeSignals.add(signal);
                discRps.put(signalType, thisTypeSignals);
            }
            else{
                thisTypeSignals.add(signal);
                discRps.put(signalType, thisTypeSignals);
            }

            discRpTypeCount.put(signalType, discRpTypeCount.containsKey(signalType) ? discRpTypeCount.get(signalType) + 1 : 1);
        }

    }

    /**
     * Read reference index .fai file for chromosome length
     * @param faIdxFile
     * @throws IOException
     */
    private void readFaIdxFile(String faIdxFile) throws IOException{
        FileInputStream fin = new FileInputStream(new File(faIdxFile));
        BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
        List<Integer> chromLength = new ArrayList<>();
        List<String> chromName = new ArrayList<>();
        // Whole genome
        String thisLine;

        while ((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split("\t");
            String refSeq = tokens[0];
            // escape "M", "MT", "chrM"
            if (refSeq.contains("M") || refSeq.contains("_")){
                break;
            }
            chromLength.add(Integer.parseInt(tokens[1]));

            chromName.add(refSeq);

        }

        chromNameMap = new String[chromName.size()];
        chromLengthMap = new int[chromLength.size()];

        for (int i = 0; i < chromName.size(); i++){
            chromNameMap[i] = chromName.get(i);
            chromLengthMap[i] = chromLength.get(i);
        }

    }

    /**
     * Ignore low quality rps
     * @param first
     * @param second
     * @return
     */
    private boolean isBadRp(SAMRecord first, SAMRecord second) {
        return first.getMappingQuality() < minMapQ || second.getMappingQuality() < minMapQ || first.getReadUnmappedFlag() || second.getReadUnmappedFlag()
                || first.isSecondaryOrSupplementary() || second.isSecondaryOrSupplementary();
    }

    /**
     * Discard reads of clipped length longer than 70% of the read length and reads with low mapping quality
     * @param cigarElements
     * @return
     */
    private boolean badReads(List<CigarElement> cigarElements, int quality){
        int clippedLength = 0;
        boolean isBad = false;
        for (CigarElement element : cigarElements){
            String operation = element.getOperator().toString();
            int optLength = element.getLength();
            if (operation.equals("S") || operation.equals("H")){
                clippedLength += optLength;
            }
        }
        if (clippedLength > 0.7 * readLen){
            isBad = true;
        }
        if (quality <= minMapQ){
            isBad = true;
        }
        return isBad;
    }



}
