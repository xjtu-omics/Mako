/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import structures.SequenceDatabase;
import structures.Node;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import utils.SuperItemLink;
import matcher.StringMatcher;
import utils.Link;

/**
 *
 * @author jiadonglin 
 */
public class Merge {
    
//    final private double minAf;
    private int patternCount;
    private String[] chrIdxNameMap;
//    private int mapq;

    public Merge(int num, String[] idxNameMap){
        patternCount = num;
        chrIdxNameMap = idxNameMap;
    }
    /**
     * Merge patterns from whole genome
     * @param mergedPatternOut
     * @param strMatcher
     * @param database
     * @param patternCandidates
     * @param regionWriter
     * @throws IOException 
     */
    public void wgsPatternMerge(String mergedPatternOut, StringMatcher strMatcher, SequenceDatabase database,
                                List<List<PseudoSequentialPattern>> patternCandidates, ReferenceSequenceFile refSeqFile, BufferedWriter regionWriter) throws IOException{
//        System.out.println("\n[Call] Start pattern post-processing, total candidate patterns: " + patternCount);

//        BufferedWriter regionWriter = new BufferedWriter(new FileWriter(svRegionOut));       
        BufferedWriter mergedWriter;
        
        if (mergedPatternOut == null){
            mergedWriter = null;
        }else{
            mergedWriter = new BufferedWriter(new FileWriter(mergedPatternOut));
        }                        
//        strMatcher = new StringMatcher();
        EstimateBPs estimator = new EstimateBPs();
        
        int numChrs = patternCandidates.size();
        for (int i = 0; i < numChrs; i ++){
            List<PseudoSequentialPattern> allPatterns = patternCandidates.get(i);
            // For single chrom, others wont be processed.
            if (allPatterns.isEmpty()){
                continue;
            }
            Map<Integer, List<Integer>> indexMap = getPatternStartIndexMap(allPatterns, database);
            List<PseudoSequentialPattern> mergedPatterns = oneChromMerge(i, database, allPatterns, patternCandidates, indexMap);
           
            oneChrPatternLinkageAnalysis(database, estimator, mergedWriter, regionWriter, mergedPatterns);
        }
        
        regionWriter.close();
        
        if (mergedWriter != null){
            mergedWriter.close();
        }        
    }
    
     /**
     * Merge patterns at single chromosome.
     * @param chrom
     * @param patterns
     * @param patternStartIndexMap
     * @return 
     */
    private List<PseudoSequentialPattern> oneChromMerge(int chrom, SequenceDatabase database, List<PseudoSequentialPattern> patterns,
                                                        List<List<PseudoSequentialPattern>> patternCandidates, Map<Integer, List<Integer>> patternStartIndexMap){
        String chr = chrIdxNameMap[chrom];
//        System.out.println("[Call] Process: " + chr +", #raw subgraphs found: " + patternCandidates.get(chrom).size());
        List<PseudoSequentialPattern> mergedPatternCandidates = new ArrayList<>();
//        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartAndIndexMap.entrySet());
        List<Map.Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartIndexMap.entrySet());
        Collections.sort(patternIndexEntrys, new Comparator<Map.Entry<Integer, List<Integer>>>(){
            @Override
            public int compare(Map.Entry<Integer, List<Integer>> o1, Map.Entry<Integer, List<Integer>> o2){
                return o1.getKey().compareTo(o2.getKey());
            }
        
        });
        
        int entrysSize = patternIndexEntrys.size();
        Set<Integer> tracker = new HashSet<>();
        for (int i = 0; i < entrysSize - 1; i ++){
            int candidateSize = mergedPatternCandidates.size();
            Map.Entry<Integer, List<Integer>> entry = patternIndexEntrys.get(i);
            int pos = entry.getKey();
            
            List<Integer> patternIndex = entry.getValue();

            Map.Entry<Integer, List<Integer>> nextEntry = patternIndexEntrys.get(i + 1);                
            List<Integer> nextPatternIndex = nextEntry.getValue();
            int nextPos = nextEntry.getKey();
            
            PseudoSequentialPattern mergedPattern = mergePatternList(patterns, patternIndex);
            PseudoSequentialPattern nextMergedPattern = mergePatternList(patterns, nextPatternIndex);
    
//            System.out.println(entry.getKey() + ": " + patternIndex.toString() + "\t" + mergedPattern.toString(database));
//            System.out.println(nextEntry.getKey() + ": " + nextPatternIndex.toString() + "\t" + nextMergedPattern.toString(database));
            List<PseudoSuperItem> mergedSuperItems = mergedPattern.mergeTwoPattern(nextMergedPattern, database);
            if (!mergedSuperItems.isEmpty()){
                tracker.add(pos);
                tracker.add(nextPos);
                
                PseudoSequentialPattern newMergedPattern = new PseudoSequentialPattern(mergedSuperItems, database);
                
//                System.out.println("merged: " + newMergedPattern.toString(database));

                // Check if the newly merged pattern is able to Merge with last pattern saved
                if (!mergedPatternCandidates.isEmpty()){
                    
                    List<PseudoSuperItem> superitems = secondaryMerge(database, mergedPatternCandidates, newMergedPattern);
                    if (!superitems.isEmpty()){
                        PseudoSequentialPattern secondaryMergedPattern = new PseudoSequentialPattern(superitems, database);
//                        System.out.println("Added pattern: " + secondaryMergedPattern.toString(database));
                        mergedPatternCandidates.remove(candidateSize - 1);
                        mergedPatternCandidates.add(secondaryMergedPattern);
                    }else{                        
//                        System.out.println("Added pattern: " + newMergedPattern.toString(database));
                        mergedPatternCandidates.add(newMergedPattern);
                    }
                }else{
//                    System.out.println("Added pattern: " + newMergedPattern.toString(database));

                    mergedPatternCandidates.add(newMergedPattern);
                }
                               
            }else{
                if (!tracker.contains(pos)){
                    tracker.add(pos);
                    if (! mergedPatternCandidates.isEmpty()){
                        List<PseudoSuperItem> superitems = secondaryMerge(database, mergedPatternCandidates, mergedPattern);
                        if (!superitems.isEmpty()){
                            PseudoSequentialPattern secondaryMergedPattern = new PseudoSequentialPattern(superitems, database);
                            mergedPatternCandidates.remove(candidateSize - 1);
                            mergedPatternCandidates.add(secondaryMergedPattern);
                        }else{
                            mergedPatternCandidates.add(mergedPattern);
                        }
                    }else{
                        mergedPatternCandidates.add(mergedPattern);
                    }                                        
                }                                                               
            }
                                
        }

        System.out.println("[Call] Process: " + chr + "\t#subgraphs found: " + mergedPatternCandidates.size());
        return mergedPatternCandidates;
    }
    /**
     * This is used to Merge a new pattern with the last pattern in the merged pattern list
     * @param mergedPatternList
     * @param aPattern
     * @return 
     */
    private List<PseudoSuperItem> secondaryMerge(SequenceDatabase database, List<PseudoSequentialPattern> mergedPatternList, PseudoSequentialPattern aPattern){
        int candidateSize = mergedPatternList.size();
        PseudoSequentialPattern lastSPInCandidates = mergedPatternList.get(candidateSize - 1);
        List<PseudoSuperItem> superitems = lastSPInCandidates.mergeTwoPattern(aPattern, database);
        return superitems;
    }
    
     /**
     * Merged patterns of a chromosome will be processed to generate SV calls.      
     * @param regionWriter
     * @param mergedPatterns
     * @throws IOException 
     */
    private void oneChrPatternLinkageAnalysis(SequenceDatabase database, EstimateBPs estimator, BufferedWriter mergedPatternWriter, BufferedWriter regionWriter,
                                              List<PseudoSequentialPattern> mergedPatterns) throws IOException{


        SuperItemLink linkageAnalyzer = new SuperItemLink();
        Collections.sort(mergedPatterns);
        int patternNums = mergedPatterns.size();

        Map<Integer, List<Link>> linkedPatternInfo = new HashMap<>();
        // pattern link by split align
        Map<Integer, Integer> splitLinkPatternBuffer = new HashMap<>();
        Set<Integer> linkedPatternIdx = new HashSet<>();
        // pattern without any link
        Set<Integer> unLinkedPattern = new HashSet<>();

        for (int i = 0; i < patternNums; i ++){

            PseudoSequentialPattern pattern = mergedPatterns.get(i);

//            if (pattern.patternLeftMostPos == 239952717) {
//                System.out.println(pattern.toString(database, chrIdxNameMap));
//            }

            if (mergedPatternWriter != null){
                mergedPatternWriter.write(pattern.toString(database, chrIdxNameMap));
                mergedPatternWriter.newLine();
            }
            // Check split alignment pos in each pattern and get the split status.
            pattern.checkSplitLinks(database, i);
//            int splitStatus = pattern.getSplitStatus();
            boolean hasSplitAlign = pattern.hasSplitAlign();

            // Patterns do not have arp superitems
            if (!pattern.hasArpSuperItem()){
                // Pattern has split aligned reads

                if (hasSplitAlign){
                    splitLinkPatternBuffer.put(i, i);

                }
                // Pattern has no arps or split-reads are not considered for confident call
                else{
                    unLinkedPattern.add(i);
                }
            }
            else{
                // First if the pattern have matched superitems through read-pairs and estimate breakpoints
                boolean isSelfLinked = pattern.isSelfLinkedPattern(database, i);

                // Linked arp superitems in pattern
                if (isSelfLinked){

                    // A pattern may contains more than one matched superitem pairs
                    List<Link> oldVal = linkedPatternInfo.get(i);
                    linkedPatternIdx.add(i);

                    if (oldVal == null){
                        oldVal = new ArrayList<>();
                        oldVal.addAll(pattern.getArpLinksInPattern());
                        linkedPatternInfo.put(i, oldVal);
                    }else{
                        oldVal.addAll(pattern.getArpLinksInPattern());
                        linkedPatternInfo.put(i, oldVal);
                    }
                }

                boolean arpSpanned = false;

                // Search mate pattern if it exists, this should be captured in pattern growth.
                // But we check again if this pattern has a mate pattern
                Map<Integer, Integer> indexMap = new HashMap<>();
                List<PseudoSequentialPattern> removedPatternCandidates = linkageAnalyzer.minusItemAndCopyWithIndexMap(mergedPatterns, i, indexMap);


                Link betweenLink = searchMateLink(database, removedPatternCandidates, pattern, i, linkageAnalyzer, indexMap);

                // Find link between patterns
                if (betweenLink != null){
                    arpSpanned = true;
                    linkedPatternIdx.add(betweenLink.getMatePatternIdx());
                    linkedPatternIdx.add(betweenLink.getPatternIdx());
                    List<Link> thisLinks = linkedPatternInfo.get(i);
                    // This is the index in original pattern list
                    if (thisLinks == null){
                        thisLinks = new ArrayList<>();
                        thisLinks.add(betweenLink);
                        linkedPatternInfo.put(i, thisLinks);
                    }else{

                        thisLinks.add(betweenLink);
                        linkedPatternInfo.put(i, thisLinks);
                    }
                }
                // Mate pattern is found through ARPs

                // If this pattern dose not have support from abnormal read-pairs, we check its split aligns
                if (!isSelfLinked && !arpSpanned){
                    if (hasSplitAlign){
                        splitLinkPatternBuffer.put(i, i);
                    }else{
                        unLinkedPattern.add(i);
                    }
                }
            }
        }

        estimator.callSVFromLinked(linkedPatternInfo, linkedPatternIdx, splitLinkPatternBuffer, mergedPatterns, database, chrIdxNameMap, regionWriter);
//        estimator.callSVFromUnlinked(unLinkedPattern, mergedPatterns, database, chrIdxNameMap, regionWriter);

    }

//    private Link searchMaxSplitLink(SequenceDatabase database, List<PseudoSequentialPattern> sortedPatterns,
//                                    PseudoSequentialPattern sourecePattern, int sourcePatternIdx, Map<Integer, Integer> indexMap){
//        Link splitLink = null;
//        for (int i = 0; i < sortedPatterns.size(); i++ ){
//            PseudoSequentialPattern targetPattern = sortedPatterns.get(i);
//            splitLink = sourecePattern.splitAlignLink(database, targetPattern, sourcePatternIdx, indexMap.get(i));
//        }
//
//        return splitLink;
//    }

     /**
     * If a target pattern has ARP Node, we can use it to search its mate pattern.Otherwise, we need to do consensus matching of itself.
     * @param sortedPatterns
     * @param sourecePattern
     * @param linker
     * @return 
     */
      
    private Link searchMateLink(SequenceDatabase database, List<PseudoSequentialPattern> sortedPatterns,
                                PseudoSequentialPattern sourecePattern, int sourcePatternIdx, SuperItemLink linker, Map<Integer, Integer> indexMap){


        List<QueryInterval> sourecePatternMateInterval = sourecePattern.superitemMateInterval;
        int length = sortedPatterns.size();
        int startIdx = 0;
        int endIdx = length - 1;
        int mateIndex = -1;
        int noQueryInterval = sourecePatternMateInterval.size();
        
//        int linkSup = -1;
        int sourcePatternMatchSuperItemIdx = -1;
        List<Integer> matchedMatePatternIdx = new ArrayList<>();
        
        // Two values to return, one is mate index and another one is number of supported read-pairs
        
        for (int i = 0; i < noQueryInterval ; i++) {
            QueryInterval interval = sourecePatternMateInterval.get(i);

            while (startIdx <= endIdx) {
                int midIdx = startIdx + (endIdx - startIdx) / 2;

                PseudoSequentialPattern targetPattern = sortedPatterns.get(midIdx);
                List<QueryInterval> sortedIntervals = targetPattern.superitemInterval;
                // Find all linked ARPs
                List<Integer> overlapAt = hasOverlap(interval, sortedIntervals);
                if (overlapAt.size() > 0) {
                    mateIndex = midIdx;
                    sourcePatternMatchSuperItemIdx = i;
                    matchedMatePatternIdx = overlapAt;
                    break;
                }

                if (isAfterInterval(interval, sortedIntervals)) {
                    startIdx = midIdx + 1;
                }
                if (isAheadInterval(interval, sortedIntervals)) {
                    endIdx = midIdx - 1;
                } else if (overlapAt.size() == 0 && !isAfterInterval(interval, sortedIntervals) && !isAheadInterval(interval, sortedIntervals)) {
                    startIdx = midIdx + 1;
                }
            }
            if (mateIndex != -1) {

                for (Integer matchedMatePatternSuperItemIdx : matchedMatePatternIdx) {
                    Node superitemOne = sourecePattern.getSuperItemOfPatternAtPos(database, sourcePatternMatchSuperItemIdx);
                    PseudoSequentialPattern matchedSequentialPattern = sortedPatterns.get(mateIndex);
                    Node superitemTwo = matchedSequentialPattern.getSuperItemOfPatternAtPos(database, matchedMatePatternSuperItemIdx);
                    boolean isEnoughARPs = linker.twoSuperItemLinkCheck(superitemOne, superitemTwo);
                    if (isEnoughARPs) {
//                        int[] mapq = new int[]{superitemOne.getSumMapQ(), superitemTwo.getSumMapQ()};
                        double[] ratio = new double[]{superitemOne.getRatio(), superitemOne.getRatio()};
//                        String[] type = new String[]{superitemOne.getType(), superitemTwo.getType()};
                        int[] weight = new int[]{superitemOne.getWeight(), superitemTwo.getWeight()};
                        int originalIndex = indexMap.get(mateIndex);
                        String bpType = inferBpTypeFromRpLink(superitemOne, superitemTwo);


                        Link betweenLink = new Link(linker.getSupLink(), superitemOne.getPos(), superitemTwo.getPos(), sourcePatternIdx, originalIndex, "bt", bpType);
                        betweenLink.setBreakInfo(ratio, weight);
                        return betweenLink;
                    }
                }

            }

            // reset start and end for next interval match
            startIdx = 0;
            endIdx = length - 1;
        }


        return null;
    }

    /**
     * Infer breakpoint type from discordant read-pairs
     * @param nodeOne
     * @param nodeTwo
     * @return
     */

    private String inferBpTypeFromRpLink(Node nodeOne, Node nodeTwo){
        String oneSigType = nodeOne.getType();
        String twoSigType = nodeTwo.getType();

        String linkSvType = "";
        if (oneSigType.equals("ARP_SMALL_INSERT")){
            linkSvType = "INS";
        }

        else if (oneSigType.equals("ARP_LARGE_INSERT")){
            linkSvType = "DEL";
        }

        else if (oneSigType.equals("ARP_RF")){
            linkSvType = "DUP";
        }

        else if (oneSigType.equals("ARP_FF") || oneSigType.equals("ARP_RR")){
            linkSvType = "INV";
        }

        return linkSvType;
    }

    /**
     * Merge patterns in a list to a new pattern.
     * @param arpPatterns
     * @param patternIndex
     * @return 
     */
    
    private PseudoSequentialPattern mergePatternList(List<PseudoSequentialPattern> arpPatterns, List<Integer> patternIndex){
        List<PseudoSequentialPattern> patterns = new ArrayList<>();
        for (Integer idx : patternIndex){
            patterns.add(arpPatterns.get(idx));
        }
        int maxLength = 0;
        int maxLengthPatternIndex = 0;
        int patternsSize = patterns.size();
        for (int i = 0; i < patternsSize; i ++){                
            if (patterns.get(i).patternLength > maxLength){
                maxLength = patterns.get(i).patternLength;
                maxLengthPatternIndex = i;
            }
        }
        return patterns.get(maxLengthPatternIndex);
    }
    private List<Integer> hasOverlap(QueryInterval targetInterval, List<QueryInterval> intervals){
        List<Integer> overlapAtIdx = new ArrayList<>();

        for (int i = 0; i < intervals.size(); i++){
            QueryInterval interval = intervals.get(i);
            if (interval != null){
                if (reciprocalOverlap(targetInterval, interval)){
                    overlapAtIdx.add(i);

                }
            }
        }
        return overlapAtIdx;
    }
    
    private boolean reciprocalOverlap(QueryInterval a, QueryInterval b){
        int aSize = a.end - a.start;
        int bSize = b.end - b.start;
        boolean isOverlapped = false;
        
        if (b.start < a.end && b.start >= a.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }else if (a.start < b.end && a.start >= b.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }
        
        return isOverlapped;
    }
    private boolean isAheadInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAhead = false;
        QueryInterval leftMostInterval = intervals.get(0);
        if (targetInterval.end < leftMostInterval.start){
            isAhead = true;
        }
        return isAhead;        
    }
    
    private boolean isAfterInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAfter = false;
        QueryInterval lastInterval = intervals.get(intervals.size() - 1);
        if (targetInterval.start > lastInterval.end){
            isAfter = true;
        }
        return isAfter;
    }
    private Map<Integer, List<Integer>> getPatternStartIndexMap(List<PseudoSequentialPattern> arpPatterns, SequenceDatabase database){
        Map<Integer, List<Integer>> indexMap = new HashMap<>();
        int numOfPatterns = arpPatterns.size();
        for (int i = 0; i < numOfPatterns ; i++){
            int patternLeftMostPos = arpPatterns.get(i).patternLeftMostPos;
            
            List<Integer> indexList = indexMap.get(patternLeftMostPos);
            if (indexList == null) {
                indexList = new ArrayList<>();
                indexMap.put(patternLeftMostPos, indexList);
            }
            indexList.add(i);
        }
        return indexMap;
    }
}
