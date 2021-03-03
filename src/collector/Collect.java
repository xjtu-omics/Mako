/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package collector;
import htsjdk.samtools.*;

import java.io.*;
import java.util.*;

import utils.FileReader;

/**
 *
 * @author jiadonglin
 */
public class Collect {
    
    private long startTime;
    private long endTime;
    
    private BufferedWriter nodeWriter = null; // Used to write nodes.
    private BufferedWriter bndWriter = null; // Write potential BNDs to BED file.
    

    public String[] chromNameMap;
    public int[] chromLengthMap;

    private final int isizeUpper;
    private final int isizeLower;
    private final int readLen;
    private final int fragMean;
    
    private final int minMapQ;
    private int maxRpDist;
    
//    private int superitemMaxSpanRange;

    private Map<String, List<SAMRecord>> rpTracker;
    // Different channels used to parse mutational signals
    private SignalParser breakChannel;
    private SignalParser isizeLargeChannel;
    private SignalParser isizeSmallChannel;
    private SignalParser oemChannel;
    private SignalParser oriChannel;
    private SignalParser breakendChannel;

    // keep normal read per base every 1Mb length.
    private int readDepthContainerBuffer = 1000000;
    private int[] readDepthContainer = new int[readDepthContainerBuffer];
    
    private int numOfARPs;
    private int numOfRPs;
    
    private Map<String, List<int[]>> maskedRegion;
    
    public Collect(int fragMean, int fragStd, int cutStd, int readLen, int maxRpDist, int minMapQ, String fastaIndexFile) throws IOException{
        isizeUpper = fragMean + cutStd * fragStd;
        isizeLower = fragMean - cutStd * fragStd;
        readFaIdxFile(fastaIndexFile);

        this.readLen = readLen;
        this.fragMean = fragMean;
        this.minMapQ = minMapQ;
        this.maxRpDist = maxRpDist;
        if (maxRpDist == -1){
            this.maxRpDist = fragMean - 2 * readLen;
        }
//        superitemMaxSpanRange = 2 * fragMean;
                
    }
    /**
     * 
     * @param bamFile the alignment file
     * @param chrom a user specified chromosome (optional)
     * @param chromStart start of a genome region (optional)
     * @param chromEnd end of a genome region (optional)
     * @param superitemOutPath output path of the created superitems
     * @param abnormalSigOut output path of the abnormal signals (optional)
     * @throws IOException 
     */
    public void runCollect(String bamFile, Set<String> qnames, String chrom, int chromStart, int chromEnd,
                       String superitemOutPath, String abnormalSigOut, String bndOutPath) throws IOException{

        startTime = System.currentTimeMillis();
        
        // the channel constructor need a name, it can be whatever you like, just used for naming some output file.
        breakChannel = new SignalParser(maxRpDist, "break", abnormalSigOut);
        isizeLargeChannel = new SignalParser(maxRpDist, "isize_large", abnormalSigOut);
        isizeSmallChannel = new SignalParser(maxRpDist, "isize_small", abnormalSigOut);
        oemChannel = new SignalParser(maxRpDist, "oem", abnormalSigOut);
        oriChannel = new SignalParser(maxRpDist, "ori", abnormalSigOut);
        breakendChannel = new SignalParser(maxRpDist, "bnd", abnormalSigOut);

        // Create output header
        if (!superitemOutPath.isEmpty()){
            nodeWriter = new BufferedWriter(new FileWriter(superitemOutPath));
            nodeWriter.write("type\tchromIdx\tnread\tpos\tsaPos\tori\tweight\tratio\tsumMapQ\tplusRead\tminusRead\tsplitRead\tsplitMapQ\titxRead\tregion\tmateRegion\tqnames\tmConsensus\tcConsensus\n");
        }

        if (bndOutPath != null) {
            bndWriter = new BufferedWriter(new FileWriter(bndOutPath));
        }

        extractSignalsFromBAM(bamFile, qnames, chrom, chromStart, chromEnd);

        nodeWriter.close();
    }
    
    /**
     * Processing the BAM file and extract abnormal alignments
     * @param bamFile
     * @param chrom
     * @param regionStart
     * @param regionEnd
     * @throws IOException 
     */
    
    private void extractSignalsFromBAM(String bamFile, Set<String> qnames, String chrom, int regionStart, int regionEnd) throws IOException{
               
        rpTracker = new HashMap<>();
        int windowStart = 0;
        int windowEnd;
        int nodeNum = 0;

        utils.FileReader myFileReader = new FileReader();
        final SamReader samReader = myFileReader.openBamFile(bamFile, ValidationStringency.SILENT, false);

        System.out.println("\nStart extracting nodes from alignments ...\n");

        // access user specified region
        if (chrom != null){

            SAMFileHeader samFileHeader = samReader.getFileHeader();
            
            SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();
            SAMSequenceRecord refSequenceRecord = sequenceDictionary.getSequence(chrom);
            int refSequenceLength = refSequenceRecord.getSequenceLength();           
            int nWindows = refSequenceLength / readDepthContainerBuffer;
            
                  
            
            if (regionStart != 0 && regionEnd != 0){
                refSequenceLength = regionEnd - regionStart + 1;
                if (refSequenceLength <= readDepthContainerBuffer){
                    nWindows = 1;
                    readDepthContainerBuffer = refSequenceLength;
                }else{
                    nWindows = refSequenceLength/readDepthContainerBuffer;
                }
                windowStart = regionStart;
                
                
            }
            
            int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
            
            for (int i = 0; i < nWindows; i++){  
                windowEnd = windowStart + readDepthContainerBuffer;                 
                SAMRecordIterator iterator = samReader.query(chrom, windowStart, windowEnd, false);               
                
                analysisAlignment(iterator, qnames, windowStart, chrom);
                processRemainingSignals();
                int curWindowNodeCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer);
                System.out.println("[Collect] processed region: [" + windowStart + ", " + windowEnd + "] " + "#nodes: " + curWindowNodeCount);
                windowStart = windowEnd;

                nodeNum += curWindowNodeCount;
                
//                readDepthPreStepBuffer = copyFromReadDepthBuffer();
                readDepthPreStepBuffer = readDepthContainer;
                
                readDepthContainer = new int[readDepthContainerBuffer];
               
                writeAllSuperItems();
                if (bndWriter != null) {
                    breakendChannel.createBreakEndCandidates(bndWriter);
                }
            }
            // process remaining alignment in BAM
            SAMRecordIterator iterator = samReader.query(chrom, windowStart, refSequenceLength, false);
            analysisAlignment(iterator, qnames, windowStart, chrom);
            processRemainingSignals();
            int curWindowNodeCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer);

            nodeNum += curWindowNodeCount;
            System.out.println("[Collect] processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#nodes: " + curWindowNodeCount);
            
            
            writeAllSuperItems();
            if (bndWriter != null) {
                breakendChannel.createBreakEndCandidates(bndWriter);
                bndWriter.close();
            }

           
        } 
        // read whole genome
        else{            
            int length = chromLengthMap.length;
            SAMRecordIterator iterator;
            for (int i = 0;i < length; i ++){
                int refSequenceLength = chromLengthMap[i];
                String curChrom = chromNameMap[i];
                System.out.println("[Collect]Start processing chrom: " + curChrom + ", chrom length: " + refSequenceLength);
                
                int nWindows = refSequenceLength / readDepthContainerBuffer;

                int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
                for (int k = 0; k < nWindows; k++){
                    windowEnd = windowStart + readDepthContainerBuffer;                 
                    iterator = samReader.query(curChrom, windowStart, windowEnd, false);
                    
                    analysisAlignment(iterator, qnames, windowStart,curChrom);
                    processRemainingSignals();
                    int curWindowNodeCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer);

                    System.out.println("[Collect] processed region: [" + windowStart + ", " + windowEnd + "] " + "#nodes: " + curWindowNodeCount);
                    windowStart = windowEnd;

                    nodeNum += curWindowNodeCount;

                    readDepthPreStepBuffer = copyFromReadDepthBuffer();
                    readDepthContainer = new int[readDepthContainerBuffer];
                    
                    writeAllSuperItems();

                    breakendChannel.createBreakEndCandidates(bndWriter);
                }
                iterator = samReader.query(curChrom, windowStart, refSequenceLength, false);
                analysisAlignment(iterator, qnames, windowStart, curChrom);
                processRemainingSignals();
                int curWindowNodeCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer);

                nodeNum += curWindowNodeCount;
                System.out.println("[Collect] processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#nodes: " + curWindowNodeCount);
                
                writeAllSuperItems();

                windowStart = 0;

                if (bndWriter != null) {
                    breakendChannel.createBreakEndCandidates(bndWriter);
                }
                rpTracker.clear();
            }
            if (bndWriter != null) {
                bndWriter.close();
            }
        }

        endTime = System.currentTimeMillis();

        StringBuilder sb = new StringBuilder();
        sb.append("\n==============  Nodes Generation =============\n");
        sb.append("Time: " + (endTime - startTime) / 1000 + "s");
        sb.append("\nTotal nodes: " + nodeNum);
        sb.append("\nTotal number of read-pairs: " + numOfRPs);
        sb.append("\nTotal number of doiscordant read-pairs: " + numOfARPs);
        sb.append("\n===============================================");

        System.out.println(sb.toString());


    }      

    /**
     * Analysis each BAM record through different channels
     * @param iterator
     * @param windowStart
     */
    private void analysisAlignment(SAMRecordIterator iterator, Set<String> qnames, int windowStart, String curRefName){
//        CigarOps corasenCigar = new CigarOps();
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();                        
//            int recordChrIdx = record.getReferenceIndex();
//            String recordChrName = record.getReferenceName();
                       
           
            List<CigarElement> cigarElements = record.getCigar().getCigarElements();
            
            if (badReads(cigarElements, record.getMappingQuality()) || record.getDuplicateReadFlag()){
                continue;
            }
            // count the number of normal read per base
            int isGoodAlign = exactAlignment(cigarElements);
            if (isGoodAlign != -1){
                updateReadDepthArray(record.getAlignmentStart(), isGoodAlign, windowStart);
            }
            transParser(record);
            SEClippedParser(record, cigarElements, curRefName);   
            RPUnmappedParser(record, curRefName);

            // Pass read pairs that have been filtered in first step.
            if (qnames.contains(record.getReadName())) {
                continue;
            }

            RPisizeParser(record, cigarElements, curRefName);
            if (!rpTracker.containsKey(record.getReadName())){
                List<SAMRecord> records = new ArrayList<>();
                records.add(record);
                rpTracker.put(record.getReadName(), records);
            }else{
                List<SAMRecord> records = rpTracker.get(record.getReadName());
                records.add(record);
                numOfRPs += 1;
                RPoriParser(records, cigarElements, curRefName);
                rpTracker.remove(record.getReadName());
            }
        } 
        iterator.close();
    }
    
    /**
     * Copy the read depth buffer in current window and save it for next window usage
     * @return 
     */
    
    private int[] copyFromReadDepthBuffer(){
        int[] newBuffer = new int[readDepthContainerBuffer];
//        int startPosToCopy = readDepthContainerBuffer - readLen;
        for (int i = 0; i < readDepthContainerBuffer; i ++){
            int val = readDepthContainer[i];
//            newBuffer[i - startPosToCopy] = val;
            newBuffer[i] = val;
        }
        return newBuffer;
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
    /**
     * Get exact matched read
     * @param cigarElements
     * @return 
     */
    private int exactAlignment(List<CigarElement> cigarElements){
        if (cigarElements.size() == 1){
            String cigarOperation = cigarElements.get(0).getOperator().toString();
            int opLength = cigarElements.get(0).getLength();
            return cigarOperation.equals("M") ? opLength : -1;
        }
        else return -1;
    }
    
    /**
     * Keep track of read depth at each site`
     * @param pos
     * @param length
     * @param windowStart 
     */
    private void updateReadDepthArray(int pos, int length, int windowStart){
        
        for (int i = 0; i < length ; i ++){
            if ( (pos + i - windowStart) >= readDepthContainerBuffer){
                continue;
            }
            else if (pos + i < windowStart) {
                continue;
            }
            else{
                readDepthContainer[pos + i - windowStart] += 1;            
            }
            
        }
        
    }
  
    /**
     * Process soft clipped reads
     * @param record
     * @param cigarElements
     */
    
    private void SEClippedParser(SAMRecord record, List<CigarElement> cigarElements, String curRefName) {
        // For a mapped read and read of relatively high mapQ
        if (!record.getReadUnmappedFlag()){
            CigarOps cigarOper = new CigarOps(record);
            String firstOperation = cigarElements.get(0).getOperator().toString();
                        
            
            cigarOper.calQueryPosFromCigar(cigarElements, 1, record.getReadNegativeStrandFlag(),readLen);
            int qsPos = cigarOper.getqsPos();
            
            String cigarStr = cigarOper.getCigarStr();
            int mutCoord = record.getAlignmentStart();           
                                    
            if (!cigarStr.equals("M") && !cigarStr.isEmpty()){
                if (firstOperation.equals("M")){                    
                    mutCoord += qsPos;                
                }
                if (cigarOper.isCoIDread() && firstOperation.equals("S")){
                    mutCoord += qsPos;
                }
                
                String ori = record.getReadNegativeStrandFlag() ? "-" : "+";                               
                MutSignal mutSignal = new MutSignal(record, cigarStr, mutCoord, ori);
                mutSignal.setIsizeNormal(isizeUpper, isizeLower);

                breakChannel.addSignals(mutSignal, curRefName);
            }
        }
        
    }
    
     /**
     * One end unmapped read
     * @param record
     * @return 
     */
    private void RPUnmappedParser(SAMRecord record, String curRefName){
        
        // read unmapped
        if (record.getReadUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
           
            oemChannel.addSignals(mutSignal, curRefName);
            numOfARPs += 1;

        }else if (record.getMateUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
                        
            oemChannel.addSignals(mutSignal, curRefName);
            numOfARPs += 1;
        }

    }
    /**
     * Process PE of abnormal insert size
     * @param record
     * @param cigarElements 
     */
    private void RPisizeParser(SAMRecord record, List<CigarElement> cigarElements, String curRefName){
        // only process read-pair mapped on the same chrom.
       
        CigarElement leftMostCigarElement = cigarElements.get(0);
        String leftMostCigarOperator = leftMostCigarElement.getOperator().toString();
      
        int mutCoord = record.getAlignmentStart();
        if (leftMostCigarOperator.equals("M") && !record.getReadNegativeStrandFlag()){
            mutCoord += leftMostCigarElement.getLength();
        }

        int insertSize = record.getInferredInsertSize();

        String ori = record.getReadNegativeStrandFlag() ? "-" : "+";
        if (Math.abs(insertSize) >= isizeUpper){    
//                System.out.println(record.getReadName() + " isize: " +Math.abs(insertSize));

            MutSignal mutSignal = new MutSignal(record, "ARP_LARGE_INSERT", mutCoord, ori);               
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), "ARP_LARGE_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
            isizeLargeChannel.addSignals(mutSignal, curRefName);
            numOfARPs += 1;
        }
        else if (Math.abs(insertSize) <= isizeLower && insertSize != 0){

            MutSignal mutSignal = new MutSignal(record, "ARP_SMALL_INSERT", mutCoord, ori);       
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                record.getInferredInsertSize(), "ARP_SMALL_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);

            isizeSmallChannel.addSignals(mutSignal, curRefName);
            numOfARPs += 1;
            
        }
                
    }
    /**
     * Process read pairs with abnormal orientation
     * @param records
     * @param cigarElements 
     */
    private void RPoriParser(List<SAMRecord> records, List<CigarElement> cigarElements, String curRefName){

        SAMRecord leftMostRecord = records.get(0);
        SAMRecord rightMostRecord = records.get(records.size() - 1);
        // For read-pair, it should be proper paired. Its read and mate are all mapped.
        if (leftMostRecord.getReadPairedFlag() && !leftMostRecord.getReadUnmappedFlag() && !leftMostRecord.getMateUnmappedFlag()){

            int mutCoord = leftMostRecord.getAlignmentStart();
            if (leftMostRecord.getReadNegativeStrandFlag() == leftMostRecord.getMateNegativeStrandFlag()){
                String mutType;
                String ori = leftMostRecord.getReadNegativeStrandFlag() ? "-":"+";
                if (leftMostRecord.getReadNegativeStrandFlag()){
                    mutType = "ARP_RR";
                }else{
                    CigarElement leftMostCigarElement = cigarElements.get(0);
                    String leftMostCigarOperation = leftMostCigarElement.getOperator().toString();

                    if (leftMostCigarOperation.equals("M")){
                        mutCoord += leftMostCigarElement.getLength();
                    }
                    mutType = "ARP_FF";
                }
                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, ori);
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                MutSignal mateMutSignal = new MutSignal(rightMostRecord, mutType, rightMostRecord.getAlignmentEnd(), ori);
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);   

                oriChannel.addSignals(readMutSignal, curRefName);
                oriChannel.addSignals(mateMutSignal, leftMostRecord.getMateReferenceName());
                numOfARPs += 1;
            }
            else if (leftMostRecord.getReadNegativeStrandFlag() && !leftMostRecord.getMateNegativeStrandFlag()){
                String mutType = "ARP_RF";

                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, "-");
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);
                MutSignal mateMutSignal = new MutSignal(rightMostRecord, mutType, rightMostRecord.getAlignmentEnd(), "+");
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                oriChannel.addSignals(readMutSignal, curRefName);
                oriChannel.addSignals(mateMutSignal, leftMostRecord.getMateReferenceName());

                numOfARPs += 1;
            }
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
     * Collect and clustering translocation signals
     * @param record
     */

    private void transParser(SAMRecord record) {
        String contig1 = record.getReferenceName();
        String contig2 = record.getMateReferenceName();

        if (!contig1.equals(contig2) && validContig(contig1) && validContig(contig2)){
            String ori = record.getReadNegativeStrandFlag() ? "rev" : "fwd";
            MutSignal traSignal = new MutSignal(record, "tra", record.getAlignmentStart(), ori);
            breakendChannel.addSignals(traSignal, record.getReferenceName());
        }
    }

    /**
     * Only deal with cannonical chromosomes
     * @param contig
     * @return
     */
    private boolean validContig(String contig) {
        for (String chrom : chromNameMap) {
            if (contig.equals(chrom)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Create SuperItems of remaining mutational signals
     */
    private void processRemainingSignals() {
        breakChannel.processFinalSignals();
        isizeLargeChannel.processFinalSignals();
        isizeSmallChannel.processFinalSignals();
        oemChannel.processFinalSignals();
        oriChannel.processFinalSignals();
                             
    }

    /**
     * Write create SuperItems to file
     * @throws IOException
     */
    private void writeAllSuperItems() throws IOException{
        if (nodeWriter != null){
            breakChannel.writeSuperItemsInChannel(nodeWriter);
            isizeLargeChannel.writeSuperItemsInChannel(nodeWriter);
            isizeSmallChannel.writeSuperItemsInChannel(nodeWriter);
            oemChannel.writeSuperItemsInChannel(nodeWriter);
            oriChannel.writeSuperItemsInChannel(nodeWriter);

        }
        
    }
    /**
     * calculate normal read aligned at a specific position and the number of superitems that generated within this window.
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer
     * @return 
     */
    private int assignReadDepthAndCountSuperItem(int windowStart, int windowSize, int[] preReadDepthBuffer){
        breakChannel.setSuperitemWeightRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeLargeChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeSmallChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oriChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oemChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);

        
        int curWindowSuperItem = 0;
        curWindowSuperItem += breakChannel.getSuperitemCount();
        curWindowSuperItem += isizeLargeChannel.getSuperitemCount();
        curWindowSuperItem += isizeSmallChannel.getSuperitemCount();
        curWindowSuperItem += oemChannel.getSuperitemCount();
        curWindowSuperItem += oriChannel.getSuperitemCount();

        
        return curWindowSuperItem;
    }
}
