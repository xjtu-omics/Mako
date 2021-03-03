/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import matcher.StringIdentity;
import matcher.StringMatcher;
import structures.SequenceDatabase;
import structures.Node;
import utils.*;
import htsjdk.samtools.QueryInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;




/**
 *
 * @author jiadonglin
 */
public class PseudoSequentialPattern implements Comparable<PseudoSequentialPattern> {

    private List<PseudoSuperItem> superitems;
    private List<PseudoSuperItem> ARPSuperItems = new ArrayList<>();
    private List<PseudoSuperItem> MSSuperItems = new ArrayList<>();
    List<QueryInterval> superitemInterval;
    List<QueryInterval> superitemMateInterval;
    private List<Integer> weights;
    private List<Integer> postions;
    private List<Double> ratios;
    private List<String> oris;
//    private Map<QueryInterval, List<Integer>> indexMap = new HashMap<>();

    int ChromId = -1;
    int patternLength;
    int patternLeftMostPos; // position of the first Node in the pattern
    int patternRightMostPos; // position of the last Node in the pattern
    int patternLeftMostIntervalStart; // read position of the left most interval of the pattern

    private int[] crossSupInfo = new int[]{-1, -1, 0, 0};

    private int[] splitAlignCoords = new int[]{-1, -1};

    // Configurations of self-linked patterns
    private List<Link> arpLinksInPattern = new ArrayList<>(); // for a self linked pattern, save its linked coord for boundary estimation
    private List<Link> splitLinksInPattern = new ArrayList<>();
//    private int[] selfLinkedEstimatedBp = new int[]{-1, -1};
//    private int[] selfLinkedSuperItemMapQ; // linked ARP Node mapping quality
//    private double[] selfLinkedSuperItemAF; // linked ARP Node allele fraction
//    private String[] selfLinkedPatternBpItemType; // Breakpoint Node type


    // Configurations of pair-wise pattern linkage
    private int[] arpSpanBp = new int[]{-1, -1};
    private int[] arpSpanBpMapQ;
    private double[] arpSpanBpAF;
    private String[] arpSpanBPItem;

    private boolean arpSpanUseSplit = false;
    private boolean hasSmallInsert = false;

    @Override
    public int compareTo(PseudoSequentialPattern other) {
        return patternLeftMostIntervalStart - other.patternLeftMostIntervalStart;
    }


    public PseudoSequentialPattern(List<PseudoSuperItem> itemset, SequenceDatabase database) {

        superitems = itemset;
        patternLength = itemset.size();
        superitemInterval = new ArrayList<>();
        superitemMateInterval = new ArrayList<>();
        postions = new ArrayList<>();
        weights = new ArrayList<>();
        ratios = new ArrayList<>();
        oris = new ArrayList<>();

        patternLeftMostPos = itemset.get(0).getSuperItem(database).getPos();
        patternRightMostPos = itemset.get(itemset.size() - 1).getSuperItem(database).getPos();

        for (int i = 0; i < patternLength; i++) {
            Node superitem = superitems.get(i).getSuperItem(database);

//                System.out.println(superitem.toConciseString());
            if (superitem.getType().contains("SMALL")) {
                hasSmallInsert = true;
            }
            weights.add(superitem.getWeight());
            postions.add(superitem.getPos());
            ratios.add(superitem.getWeightRatio());
            oris.add(superitem.getReadOri());

            if (ChromId == -1) {
                ChromId = superitem.getChromIdx();
            }
            if (superitem.isARPsuperitem() && !superitem.getType().contains("OEM")) {
                ARPSuperItems.add(superitems.get(i));
                superitemInterval.add(superitem.getSuperitemRegion());
                superitemMateInterval.add(superitem.getSuperitemMateRegion());
            } else {
                MSSuperItems.add(superitems.get(i));
            }
        }
        Collections.sort(ARPSuperItems);
        Collections.sort(superitemInterval);
        if (!superitemInterval.isEmpty()) {
            patternLeftMostIntervalStart = superitemInterval.get(0).start;
        }

    }


    /**
     * Return pattern spanned genome region
     *
     * @return
     */
    public int[] getPatternSpanRegion() {
        return new int[]{patternLeftMostPos, patternRightMostPos};
    }

    public int[] getPatternSuperItemPos(SequenceDatabase database) {
        int[] pos = new int[superitems.size()];
        for (int i = 0; i < superitems.size(); i++) {
            pos[i] = superitems.get(i).getSuperItem(database).getPos();
        }
        return pos;
    }

    public String toString(SequenceDatabase database, String[] idxNameMap) {
        StringBuilder sb = new StringBuilder();

        sb.append(idxNameMap[ChromId]);
        sb.append("\t");
        sb.append(patternLeftMostPos);
        sb.append("\t");
        sb.append(patternRightMostPos);
        sb.append("\t");

        String patternStr = "";

        for (PseudoSuperItem item : superitems) {
            Node superitem = item.getSuperItem(database);
            patternStr += superitem.getType();
            patternStr += ",";
            sb.append('(');
            sb.append(superitem.toConciseString());
            sb.append(')');
        }
        patternStr = patternStr.substring(0, patternStr.length() - 1);
        sb.append(patternStr);
        return sb.toString();
    }

    public String toTypeString(SequenceDatabase database) {
        StringBuilder sb = new StringBuilder();
        for (PseudoSuperItem item : superitems) {
            Node superitem = item.getSuperItem(database);
            sb.append(superitem.getType());
            sb.append(",");
        }

        return sb.substring(0, sb.length() - 1);
    }

    public List<Node> getSuperItemsOfPattern(SequenceDatabase database) {
        List<Node> superItemsList = new ArrayList<>();
        for (PseudoSuperItem item : superitems) {
            Node superitem = item.getSuperItem(database);
            superItemsList.add(superitem);
        }
        return superItemsList;
    }


    private boolean arpLinkedPatternSplitEvidence(int[] splitCoords, PseudoSequentialPattern matePattern) {
        int linkedPatternLeftPos = patternLeftMostPos < matePattern.patternLeftMostPos ? patternLeftMostPos : matePattern.patternLeftMostPos;
        int linkedPatternRightPos = patternRightMostPos > matePattern.patternRightMostPos ? patternRightMostPos : matePattern.patternRightMostPos;
        boolean splitStatus = false;
        int[] matePatternSplitCoord = matePattern.getSplitAlignCoords();
        if ((splitCoords[0] >= linkedPatternLeftPos && splitCoords[1] <= linkedPatternRightPos) ||
                (matePatternSplitCoord[0] >= linkedPatternLeftPos && matePatternSplitCoord[1] <= linkedPatternRightPos)) {
            splitStatus = true;
        }
        return splitStatus;
    }

    public void findSplitAlignLink(SequenceDatabase database, PseudoSequentialPattern matePattern, int patternIdx, int matePatternIdx) {
        List<Integer> mateItemPos = matePattern.getPatternItemPos(database);
        int maxSA = 0;
        int splitLinkedPatternItemIdx = -1;
        int splitLinkedPatternItemPos = -1;
        int thisPatternItemIdx = -1;
        int thisPatternItemPos = -1;
        List<Link> potentialSplitLinks = new ArrayList<>();
        for (int j = 0; j < superitems.size(); j++) {
            Node item = superitems.get(j).getSuperItem(database);
            int itemSplitPos = item.getSplitAlignPos();
            int sa = item.getSplitReadCount();
            if (itemSplitPos > matePattern.getPatternRightMostPos() || itemSplitPos < matePattern.getPatternLeftMostPos()) {
                continue;
            }
            for (int i = 0; i < mateItemPos.size(); i++) {
                int matePos = mateItemPos.get(i);
                if (matePos - 100 < itemSplitPos && itemSplitPos < matePos + 100) {
                    Node mateItem = matePattern.getSuperItemFromOriginal(database, i);
                    if (sa > maxSA) {
                        splitLinkedPatternItemIdx = i;
                        splitLinkedPatternItemPos = mateItem.getPos();
                        thisPatternItemIdx = j;
                        thisPatternItemPos = item.getPos();
                        maxSA = sa;
                    }
                }
            }
        }
        if (splitLinkedPatternItemIdx != -1) {
            Node thisNode = superitems.get(thisPatternItemIdx).getSuperItem(database);
            Node linkedNode = superitems.get(splitLinkedPatternItemIdx).getSuperItem(database);

            String bpType = inferBpTypeFromNonRps(thisNode, linkedNode);
            potentialSplitLinks.add(new Link(maxSA, thisPatternItemPos, splitLinkedPatternItemPos, patternIdx, matePatternIdx, "sa", bpType));


        }
        getSplitStatus(potentialSplitLinks);

    }


    private List<Integer> getPatternItemPos(SequenceDatabase database) {
        List<Integer> itemPos = new ArrayList<>();
        for (PseudoSuperItem psItem : superitems) {
            Node item = psItem.getSuperItem(database);
            itemPos.add(item.getPos());
        }
        return itemPos;
    }

    /**
     * A pattern has links between nodes is a self-linked pattern. Match superitems in the pattern
     * Use ARP Node for boundary estimation if there dose not exit split align or clipped reads
     *
     * @param database
     * @return
     */
    public boolean isSelfLinkedPattern(SequenceDatabase database, int patternIdx) {

        SuperItemLink linker = new SuperItemLink();

        Collections.sort(ARPSuperItems, new Comparator<PseudoSuperItem>() {
            @Override
            public int compare(PseudoSuperItem o1, PseudoSuperItem o2) {
                return o1.getInterval()[0] - o2.getInterval()[0];
            }
        });
        int Arps = ARPSuperItems.size();

        List<PseudoSuperItem> searchSpace;
        int machtedSuperItemIdx = -1;

        Set<Integer> matchedItem = new HashSet<>();
        int singleSuperItemARPcount = 0;

        for (int i = 0; i < Arps; i++) {
            if (matchedItem.contains(i)) {
                continue;
            }
            PseudoSuperItem target = ARPSuperItems.get(i);
            int siArpCount = singleSuperItemARPs(target, database);
            if (siArpCount > singleSuperItemARPcount) {
                singleSuperItemARPcount = siArpCount;
            }

            Map<Integer, Integer> idxMap = new HashMap<>();
            searchSpace = minusSelf(ARPSuperItems, i, idxMap);

            int mateIndex = linker.mateSuperItemSearch(searchSpace, target, database);
            // Linked with other node
            if (mateIndex != -1) {
//                curSuperItemIdx = i;
                machtedSuperItemIdx = mateIndex;
                // Original superitem index in ARPSuperItems
                int originalIdx = idxMap.get(machtedSuperItemIdx);
                matchedItem.add(originalIdx);

                Node nodeOne = target.getSuperItem(database);
                Node nodeTwo = ARPSuperItems.get(originalIdx).getSuperItem(database);
//                String[] linkedItemType = new String[]{nodeOne.getType(), nodeTwo.getType()};
//                int[] linkedItemMapQ = new int[]{nodeOne.getSumMapQ(), nodeTwo.getSumMapQ()};
                double[] linkedItemAF = new double[]{nodeOne.getRatio(), nodeTwo.getRatio()};
                int[] linkedNodeWeight = new int[]{nodeOne.getWeight(), nodeTwo.getWeight()};
                int superItemOnePos = nodeOne.getPos();
                int superItemTwoPos = nodeTwo.getPos();
                String svType = inferBpTypeFromRpLink(nodeOne, nodeTwo);


                if (superItemTwoPos > superItemOnePos) {
                    Link thisLink = new Link(linker.getSupLink(), superItemOnePos, superItemTwoPos, patternIdx, patternIdx, "rp", svType);
                    thisLink.setBreakInfo(linkedItemAF, linkedNodeWeight);
                    arpLinksInPattern.add(thisLink);
                } else {
                    Link thisLink = new Link(linker.getSupLink(), superItemTwoPos, superItemOnePos, patternIdx, patternIdx, "rp", svType);
                    thisLink.setBreakInfo(linkedItemAF, linkedNodeWeight);
                    arpLinksInPattern.add(thisLink);
                }
            }
        }
        // check if the current pattern contains self-linked node.
        if (machtedSuperItemIdx == -1){
            for (int i = 0; i < Arps; i++) {
                PseudoSuperItem target = ARPSuperItems.get(i);
                int linkWeight = singleSuperItemARPs(target, database);
                double itemAF = target.getSuperItem(database).getRatio();
                int itemWeight = target.getSuperItem(database).getWeight();
                if (linkWeight >= 3) {
                    int leftPos = target.getInterval()[0];
                    int rightPos = target.getInterval()[1];
                    Node thisNode = target.getSuperItem(database);
                    String svType = inferBpTypeFromRpLink(thisNode, thisNode);
                    Link thisLink = new Link(linkWeight, leftPos, rightPos, patternIdx, patternIdx, "rp", svType);
                    thisLink.setBreakInfo(new double[]{itemAF, itemAF}, new int[]{itemWeight, itemWeight});
                    arpLinksInPattern.add(thisLink);
                }
            }
        }
        return arpLinksInPattern.size() > 0;
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
     * Infer breakpoint type from either split alignment or clipped signals
     * @return
     */
    private String inferBpTypeFromNonRps(Node nodeOne, Node nodeTwo) {
        String nodeOneType = nodeOne.getType();
        String nodeTwoType = nodeTwo.getType();

        int nodeOnePlusRead = nodeOne.getNumPlusRead();
        int nodeOneMinusRead = nodeOne.getNumMinusRead();

        int nodeTwoPlusRead = nodeTwo.getNumPlusRead();
        int nodeTwoMinusRead = nodeTwo.getNumMinusRead();

        String linkSvType = "";

        if ((nodeOnePlusRead > 0 && nodeOneMinusRead > 0) && (nodeTwoPlusRead > 0 && nodeTwoMinusRead > 0)) {
            if (nodeOneType.equals("MS") && nodeTwoType.equals("SM")) {
                linkSvType = "DEL";
            }
            else if (nodeOneType.equals("SM") && nodeTwoType.equals("MS")) {
                linkSvType = "DUP";
            }
        }
        if (nodeOneMinusRead * nodeOnePlusRead == 0 && nodeTwoPlusRead * nodeTwoMinusRead == 0) {
            if (nodeOnePlusRead > 0 && nodeTwoPlusRead > 0) {
                if (nodeOneType.equals("MS") && nodeTwoType.equals("SM")) {
                    linkSvType = "DEL";
                }
                else if (nodeOneType.equals("SM") && nodeTwoType.equals("MS")) {
                    linkSvType = "DUP";
                }
            }
            else {
                linkSvType = "BND";
            }
        }

        return linkSvType;
    }

    /**
     * Check if a Node contain entire abnormal read-pair
     *
     * @param ps
     * @param database
     * @return
     */

    private int singleSuperItemARPs(PseudoSuperItem ps, SequenceDatabase database) {
        Node si = ps.getSuperItem(database);
        String[] qnames = si.getQNames();
        int readPairNums = readPairCounter(qnames);
        return readPairNums;
    }

    private List<PseudoSuperItem> minusSelf(List<PseudoSuperItem> psSuperItems, int idx, Map<Integer, Integer> idxMap) {
        List<PseudoSuperItem> newSuperItems = new ArrayList<>();
        int length = psSuperItems.size();
        for (int i = 0; i < length; i++) {
            if (i != idx) {
                newSuperItems.add(psSuperItems.get(i));
                idxMap.put(newSuperItems.size() - 1, i);
            }

        }
        return newSuperItems;
    }


    public Node getSuperItemFromOriginal(SequenceDatabase database, int idx) {
        PseudoSuperItem item = superitems.get(idx);
        return item.getSuperItem(database);
    }

    public Node getSuperItemOfPatternAtPos(SequenceDatabase database, int idx) {
        PseudoSuperItem item = ARPSuperItems.get(idx);
        return item.getSuperItem(database);
    }

    public List<PseudoSuperItem> mergeTwoPattern(PseudoSequentialPattern aPattern, SequenceDatabase database) {

        List<PseudoSuperItem> mergedSuperitems = new ArrayList<>();
        List<Node> patternOneNodes = getSuperItemsOfPattern(database);
        List<Node> patternTwoNodes = aPattern.getSuperItemsOfPattern(database);

        int lengthOne = patternOneNodes.size();
        int lengthTwo = patternTwoNodes.size();

        int matchedIndexAtPatternOne = -1;
        int lastMatchedIndexAtPatternOne = lengthOne;

        Node patternTwoStartNode = patternTwoNodes.get(0);
        Node patternTwoLastNode = patternTwoNodes.get(lengthTwo - 1);

        for (int i = 0; i < lengthOne; i++) {
            Node patternOneNode = patternOneNodes.get(i);
            if (patternTwoStartNode.isEqual(patternOneNode)) {
                matchedIndexAtPatternOne = i;
            }
            if (patternTwoLastNode.isEqual(patternOneNode)) {
                lastMatchedIndexAtPatternOne = i;
            }
        }

        if (matchedIndexAtPatternOne != -1) {
            if (lastMatchedIndexAtPatternOne == lengthOne) {
                List<PseudoSuperItem> subListOfPatternOne = superitems.subList(0, matchedIndexAtPatternOne);
                mergedSuperitems.addAll(subListOfPatternOne);
                mergedSuperitems.addAll(aPattern.superitems);
            } else {
                mergedSuperitems = superitems;
            }
        }

        return mergedSuperitems;
    }

    /**
     * For clipped Super-Item, match clipped sequence with aligned sequence
     *
     * @param database
     */

    public void doCrossLink(SequenceDatabase database) {
        // info[0], info[1] start and end pos. 
        // info[2] shared string length
        // info[3] number of reads share string
        int[] info = new int[]{-1, -1, 0, 0};

        StringBuilder sb;
        StringMatcher strMatcher = new StringMatcher();
        List<String> mStrings = new ArrayList<>();
        List<String> sForwardStrings = new ArrayList<>();
        List<String> sReverseStrings = new ArrayList<>();

        int allStrLength = 0;
        for (PseudoSuperItem psItem : superitems) {
            Node node = psItem.getSuperItem(database);
            String matchedStr = node.getMachedConsensus();
            String clippedStr = node.getClippedConsensus();
            sb = new StringBuilder(clippedStr);
            String clippedRevereStr = sb.reverse().toString();

            allStrLength += matchedStr.length();
            allStrLength += clippedStr.length();
            allStrLength += clippedRevereStr.length();

            mStrings.add(matchedStr);
            sForwardStrings.add(clippedStr);
            sReverseStrings.add(clippedRevereStr);
        }

        strMatcher.strCrossMatch(mStrings, sForwardStrings, sReverseStrings);

        int numOfStrs = superitems.size();
        double avgLen = (double) allStrLength / superitems.size();
        double val = Math.log(numOfStrs * avgLen) / Math.log(4);
        int expect = (int) Math.ceil(val);

        int[] linkInfo = strMatcher.isCrossLinked();
        int infoSum = linkInfo[0] + linkInfo[1] + linkInfo[2];

        if (infoSum > 0) {
            int observedLen = strMatcher.isForwardExist() ? strMatcher.getForwardSharedStrLength() : strMatcher.getReverseSharedStrLength();

            if (observedLen > expect) {
//                System.out.println("Expect: " + expect + "\tObserve: " + observedLen); 
//                strMatcher.printSharedString(superitems, database);
//                System.out.println("\n");
                List<StringIdentity> identitys = strMatcher.isForwardExist() ? strMatcher.getForwardSharedStrIdentitys() : strMatcher.getReverseSharedStrIdentitys();
                if (validClipLink(identitys, database)) {
                    int[] coords = strMatcher.getEstimateBp(superitems, database);
                    if (coords[0] != coords[1]) {
                        info[0] = coords[0];
                        info[1] = coords[1];
                        info[2] = observedLen;
                        info[3] = infoSum;
                    }

                }

            }
        }
        crossSupInfo = info;
    }

    private boolean validClipLink(List<StringIdentity> strIdentitys, SequenceDatabase database) {
        List<String> types = new ArrayList<>();
        boolean isValid = false;
        Collections.sort(strIdentitys, new Comparator<StringIdentity>() {
            @Override
            public int compare(StringIdentity o1, StringIdentity o2) {
                return o1.getSeqId() - o2.getSeqId();
            }
        });

        for (StringIdentity id : strIdentitys) {
            PseudoSuperItem psItem = superitems.get(id.getSeqId());
            Node node = psItem.getSuperItem(database);
            types.add(node.getType());
        }
        if (types.size() > 1) {
            String firstType = types.get(0);
            String lastType = types.get(types.size() - 1);
            if ((firstType.equals("MS") && lastType.equals("SM")) || (firstType.equals("SM") && lastType.equals("MS"))) {
                isValid = true;
            }
            isValid = true;
        }
        return isValid;
    }

    /**
     * First check if the pattern contains split aligned info, we only consider split read of mapQ at least 35.
     *
     * @param database
     */
    public void checkSplitLinks(SequenceDatabase database, int patternIdx) {
        int[] pos = new int[]{-1, -1};
        List<Link> potentialSplitLinks = new ArrayList<>();
        for (int i = 0; i < superitems.size(); i++) {
            Node thisNode = superitems.get(i).getSuperItem(database);
//            String superItemType = superitem.getType();
            int primaryPos = thisNode.getPos();
            int splitAlignPos = thisNode.getSplitAlignPos();
            if (splitAlignPos != -1 && primaryPos != splitAlignPos) {
                String bpType = inferBpTypeFromNonRps(thisNode, thisNode);
                if (bpType.equals("")) {
                    bpType = inferBpTypeFromRpLink(thisNode, thisNode);
                }
                // This might indicate translocation
                if (Math.abs(primaryPos - splitAlignPos) > 1000000 || Math.abs(primaryPos - splitAlignPos) < 50) {
                    continue;
                }
                int splitReadCount = thisNode.getSplitReadCount();
                if (splitReadCount < 2) {
                    continue;
                }
                if (primaryPos < splitAlignPos) {

                    Link splitLink = new Link(splitReadCount, primaryPos, splitAlignPos, patternIdx, patternIdx, "sa", bpType);
                    potentialSplitLinks.add(splitLink);
                } else {
                    Link splitLink = new Link(splitReadCount, splitAlignPos, primaryPos, patternIdx, patternIdx, "sa", bpType);
                    potentialSplitLinks.add(splitLink);
                }

            }
        }
        getSplitStatus(potentialSplitLinks);

    }


    /**
     * pattern linked from split aligned reads
     * -3: inproper split align, split aligned pos either not in the pattern or in another pattern
     * -2: proper split align, where pos and split aligned pos are within the pattern spanned region
     * >0: split aligned pos is in another pattern
     *
     * @return
     */
    private void getSplitStatus(List<Link> links) {

//        int splitAlignLeftMatch = -1;
//        int splitAlignRightMatch = -1;

        for (Link splitLink : links) {
            int[] splitAlignPos = splitLink.getLinkedItemPos();

            // If the split aligned coords within the pattern
            if (splitAlignPos[0] >= patternLeftMostPos && splitAlignPos[1] <= patternRightMostPos) {
                splitLink.setSplitLinkStatus(-2);
                splitLinksInPattern.add(splitLink);
            }
        }
    }


//    public int[] estimateBpFromArp(SequenceDatabase database, int minQ){
//        // Use split-alignment to get the BP pos if there exist split align of a Node.
//        int[] info = new int[]{-1, -1, 0};
//        int arpsNum = ARPSuperItems.size();
//        // Used to count number of whole read-pairs within a Node.
//        int maxReadPairIdx = -1;
//        int maxReadPairNum = 1;
//        if (arpsNum > 0){
//            for (int i = 0; i < arpsNum; i++){
//                Node si = ARPSuperItems.get(i).getSuperItem(database);
//                if (si.getSumMapQ() < si.getWeight() * minQ && si.getRatio() >= 0.2){
//                    continue;
//                }
//                String[] qnames = si.getQNames();
//                int readPairNums = readPairCounter(qnames);
//                if (readPairNums > maxReadPairNum){
//                    maxReadPairNum = readPairNums;
//                    maxReadPairIdx = i;
//                }
//            }
//        }
//        if (maxReadPairIdx != -1){
//            QueryInterval interval = ARPSuperItems.get(maxReadPairIdx).getSuperItem(database).getSuperitemRegion();
//            QueryInterval mateInterval = ARPSuperItems.get(maxReadPairIdx).getSuperItem(database).getSuperitemMateRegion();
//            info[0] = interval.start;
//            info[1] = mateInterval.end;
//            info[2] = maxReadPairNum;
//        }
//
//        return info;
//    }
    
    private int readPairCounter(String[] qnames){
        int count = 0;
        Set<String> qSet = new HashSet<>();
        for (String q : qnames){
            if (qSet.contains(q)){
                count += 1;
            }else{
                qSet.add(q);
            }
        }
        return count;
    }
    
    /**
     * For patterns contain ARP superitem but do not have link info.
     * @param database
     * @param minQ
     * @return 
     */
//    public int[] unlinkedArpPatternPosEstimate(SequenceDatabase database, int minQ){
//        int[] info = new int[3];
//        // If the pattern contains 'ARP_SMALL_INSERT', return the position as left and right pos of the pattern
//        info[0] = patternLeftMostPos;
//        info[1] = patternRightMostPos;
//
//        if (hasSmallInsert){
//
//            for(PseudoSuperItem psItem : ARPSuperItems) {
//                Node item = psItem.getSuperItem(database);
//                if (item.getType().contains("SMALL")) {
//                    info[2] = item.getWeight();
//                }
//            }
//        }
//
//        if (!hasSmallInsert){
//            info = estimateBpFromArp(database, minQ);
//        }
//
//        return info;
//    }
    /**
     * OEM reads are important sign for potential SVs, but only OEM RPs are not confident enough
     * @param database
     * @return 
     */
    public int[] oemPatternPosEstimate(SequenceDatabase database){
        int[] info = new int[]{-1, -1, 0};
        for (int i = 0; i < superitems.size(); i++){
            Node si = superitems.get(i).getSuperItem(database);
            if (si.getType().equals("ARP_OEM") && si.getWeightRatio() >= 0.2){
                info[2] = si.getWeight(); 
                int[] bp = oemBpSearch(i, si.getOri(), database);
                if (bp[0]!=-1 && bp[1]!=-1){
                    info[0] = bp[0];
                    info[1] = bp[1];
                }
                
            }
        }
        return info;
    }
    
    private int[] oemBpSearch(int oemIdx, String oemOri, SequenceDatabase database){
        int[] distToIdx;        
        int[] bp; 
        
        if(oemOri.equals("+")){            
            distToIdx = new int[superitems.size() - oemIdx];
            for (int i = oemIdx; i < superitems.size(); i++){
                Node si = superitems.get(i).getSuperItem(database);
                if (si.isClipped()){                     
                    distToIdx[i - oemIdx] = si.getPos() ;
//                    lastPos = si.getPos();                                       
                }else{
                    distToIdx[i - oemIdx] = -1;
                }
            }     
            
            bp = oemBpHelper(distToIdx, oemIdx, oemOri, database);
                        
        }else{
            distToIdx = new int[oemIdx + 1];
            for (int i = oemIdx; i >= 0; i--){
                Node si = superitems.get(i).getSuperItem(database);
                if (si.isClipped()){                     
                    distToIdx[i] = si.getPos();
//                    lastPos = si.getPos();                                       
                }else{
                    distToIdx[i] = -1;
                }
            }
            bp = oemBpHelper(distToIdx, oemIdx, oemOri, database);
            
        }
                                       
        return bp;
        
    }
    
    
    private int[] oemBpHelper(int[] distArray, int oemIdx, String oemOri, SequenceDatabase database){
        int bp[] = new int[]{-1, -1};
        if (oemOri.equals("+")){
            // ARP_OEM is the last superitem in the pattern
            if (oemIdx == superitems.size() - 1){
                return bp;
            }
            for (int i = 1; i < distArray.length;i++){
                if (distArray[i] == -1){
                    continue;
                }

                int closeNum = 1;
                int sumDist = 0;
                int leftBoundIdx = -1;

                for (int j = 1; i + j < distArray.length; j++){
                    if (distArray[i + j] == -1){
                        continue;
                    }
                    sumDist += distArray[i + j] - distArray[i];
                    if (sumDist < 100){
                        closeNum += 1;
                        leftBoundIdx = i + j;
                    }
                }
                if (closeNum > 1 && leftBoundIdx != -1){
                    bp[0] = superitems.get(oemIdx).getSuperItem(database).getPos();
//                    bp[0] = superitems.get(i + oemIdx).getSuperItem(database).getPos();
                    bp[1] = superitems.get(leftBoundIdx + oemIdx).getSuperItem(database).getPos();
                    
                }
            }
        }else{
            // ARP_OEM is the first superitem in the pattern
            if (oemIdx == 0){
                return bp;
            }
    
            for (int i = oemIdx - 1; i >= 0; i--){
                if (distArray[i] == -1){
                    continue;
                }

                int closeNum = 1;
                int sumDist = 0;
                int leftBoundIdx = -1;

                for (int j = 1; i - j >= 0; j++){
                    if (distArray[i - j] == -1){
                        continue;
                    }
                    sumDist += distArray[i] - distArray[i - j];
                    if (sumDist < 100){
                        closeNum += 1;
                        leftBoundIdx = i - j;
                    }
                }
                if (closeNum > 1 && leftBoundIdx != -1){
                    bp[0] = superitems.get(leftBoundIdx).getSuperItem(database).getPos();
                    bp[1] = superitems.get(oemIdx).getSuperItem(database).getPos();                    
                }
            }
        }
        
        return bp;
    }

    public String getPatternOris() {
        StringBuilder sb = new StringBuilder();
        for (String ele : oris){
            sb.append(ele);
            sb.append(",");
        }
        String outStr = sb.toString();
        return outStr.substring(0, outStr.length() - 1);
    }
    
    public int[] getSplitAlignCoords(){
        return splitAlignCoords;
    }
    public int getPatternLeftMostPos(){
        return patternLeftMostPos;
    }
    public int getPatternRightMostPos(){
        return patternRightMostPos;
    }
    
    public boolean isCrossSup() {
        return crossSupInfo[0] > 0 && crossSupInfo[1] > 0 && crossSupInfo[2] >= 20;
    }
    
    public int[] getCrossInfo() {
        return crossSupInfo;
    }
    
    public int getCrossedLen() {
        return crossSupInfo[2];
    }
    
    public int[] getCrossedPos() {
        if (isCrossSup()) {
            return new int[]{crossSupInfo[0], crossSupInfo[1]};
        }else{
            return new int[]{-1, -1};
        }
    }

//    public List<Integer> getIndex(QueryInterval aInterval){
//        return indexMap.get(aInterval);
//    }
    public List<Integer> getWeights(){
        return weights;
    }
    
    public List<Double> getWeightRatio(){
        return ratios;
    }

    public List<String> getOris(){
        return oris;
    }


    public boolean hasArpSuperItem(){
        return !ARPSuperItems.isEmpty();
    }

//    public int getSplitSupCount(){
//        return splitReadSup;
//    }

    public boolean hasSplitAlign() {
        return splitLinksInPattern.size() > 0;
    }
    public boolean hasArpLinks() {
        return arpLinksInPattern.size() > 0;
    }

//    public int getSplitStatus() {return splitStatus;}


    public List<Link> getArpLinksInPattern() {
        return arpLinksInPattern;
    }
    public List<Link> getSplitLinksInPattern() {
        return splitLinksInPattern;
    }

    public int getSelfLinkedBpSup() {
        int supRead = 0;
        int num = 0;
        if (arpLinksInPattern.size() > 0){
            for (Link link : arpLinksInPattern) {
                num += 1;
                supRead += link.getSups();

            }
            return supRead / num;
        }
        return 0;
    }


    public int getSelfLinksCount(){
        return arpLinksInPattern.size();
    }

//    public int getSplitReadMapQ(){
//        return splitReadMapQ;
//    }
    
    public boolean arpSpanHasSplit() {
        return arpSpanUseSplit;
    }
    
    public int[] getArpSpanBp() {
        return arpSpanBp;
    }
    public int getPatternChromId(){
        return ChromId;
    }
}
