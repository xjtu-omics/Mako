/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import structures.SequenceDatabase;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

import utils.Candidate;

import java.util.Map.Entry;
import utils.Link;
import utils.LinkGraph;


/**
 *
 * @author jiadonglin
 */


public class EstimateBPs {

    public EstimateBPs(){

    }
    /**
     * Call SVs from patterns that are linked by ARP, cross or split align.
     * @param linkedPatternInfo
     * @param splitLinkPattern
     * @param mergedPatterns
     * @param database
     * @param idxNameMap
     * @param regionWriter
     * @throws IOException 
     */
    public void callSVFromLinked(Map<Integer, List<Link>> linkedPatternInfo, Set<Integer> linkedPatternNames, Map<Integer, Integer> splitLinkPattern, List<PseudoSequentialPattern> mergedPatterns,
            SequenceDatabase database, String[] idxNameMap, BufferedWriter regionWriter) throws IOException{


        // Call SVs from patterns only have split align
        for (Entry<Integer, Integer> entry : splitLinkPattern.entrySet()){
            String linkType = "Split";
            int idx = entry.getKey(); 
            PseudoSequentialPattern pattern = mergedPatterns.get(idx);
            pattern.doCrossLink(database);
            if (pattern.isCrossSup()) {
                linkType += ";Cross";
            }
            Candidate sv = new Candidate(linkType, pattern, null, null, database);
            sv.writeSV(pattern, regionWriter, idxNameMap);
                       
        }
//        System.out.println("Linked patterns: " + linkedPatternInfo.size());
        Map<Integer, List<Link>> tmp = findMaximalLinks(linkedPatternInfo, linkedPatternNames);

        // Call SVs from ARP linked patterns
        Set<Integer> linkedPatternCalled = new HashSet<>();
        for (Entry<Integer, List<Link>> entry : tmp.entrySet()){
            StringBuilder linkTypeStrBuilder = new StringBuilder();
            int patternIdx = entry.getKey();
            PseudoSequentialPattern pattern = mergedPatterns.get(patternIdx);
//            if (pattern.patternLeftMostPos == 16078651) {
//                System.out.println(pattern.toString(database, idxNameMap));
//            }
            if (linkedPatternCalled.contains(patternIdx)) {
                continue;
            }
            List<Link> matePatternLinks = entry.getValue();
            Map<Integer, PseudoSequentialPattern> matePatterns = new HashMap<>();
            List<Link> betweenLinks = new ArrayList<>();
            // This pattern has more than one arp links
            if (matePatternLinks.size() > 1){
                Set<String> linkTypes = new HashSet<>();

                for (Link aLink : matePatternLinks) {

                    if (aLink.isSelfLink()) {
                        linkTypes.add("ARP_Self");
                    }else{
                        addNonDuplicatedLinks(betweenLinks, aLink);
                        linkTypes.add("ARP_Span");
                        if (pattern.arpSpanHasSplit()) {
                            linkTypes.add("Split");
                        }
                        int mateIdx = patternIdx == aLink.getMatePatternIdx() ? aLink.getPatternIdx() : aLink.getMatePatternIdx();

                        linkedPatternCalled.add(patternIdx);
                        linkedPatternCalled.add(mateIdx);

                        PseudoSequentialPattern matePattern = mergedPatterns.get(mateIdx);

                        matePatterns.put(mateIdx, matePattern);
                    }

                }

                for (String type : linkTypes) {
                    linkTypeStrBuilder.append(type);
                    linkTypeStrBuilder.append(";");
                }
                String linkTypeStr = linkTypeStrBuilder.toString().substring(0, linkTypeStrBuilder.length() - 1);
                Candidate sv = new Candidate(linkTypeStr, pattern, matePatterns, betweenLinks, database);

                sv.writeSV(pattern, regionWriter, idxNameMap);
            }
            else{
                Link thisLink = matePatternLinks.get(0);
                int mateIdx = thisLink.getMatePatternIdx();
                if (mateIdx == patternIdx) {
                    mateIdx = thisLink.getPatternIdx();
                }
                pattern.doCrossLink(database);

                // Self linked pattern
                if (thisLink.isSelfLink()) {
                    linkTypeStrBuilder.append("ARP_Self");
                    // In pattern split support
                    if (pattern.hasSplitAlign()){
                        linkTypeStrBuilder.append(";Split");
                        // split align && cross linked
                        if (pattern.isCrossSup()){
                            linkTypeStrBuilder.append(";Cross");
                        }

                    }
                    Candidate sv = new Candidate(linkTypeStrBuilder.toString(), pattern, null, null, database);
                    linkedPatternCalled.add(patternIdx);
                    sv.writeSV(pattern, regionWriter, idxNameMap);

                }
                else{
                    linkedPatternCalled.add(patternIdx);
                    linkedPatternCalled.add(mateIdx);

                    PseudoSequentialPattern matePattern = mergedPatterns.get(mateIdx);

                    linkTypeStrBuilder.append("ARP_Span");

                    if (pattern.arpSpanHasSplit()){
                        // arp linked two pattern with split align support
                        linkTypeStrBuilder.append(";Split");
                    }
                    matePatterns.put(mateIdx, matePattern);
                    betweenLinks.add(thisLink);
                    Candidate sv = new Candidate(linkTypeStrBuilder.toString(), pattern, matePatterns, betweenLinks, database);

                    sv.writeSV(pattern, regionWriter, idxNameMap);

                }
            }


        }                
    }

//    public void callSVFromUnlinked(Set<Integer> unLinkedPattern, List<PseudoSequentialPattern> mergedPatterns, SequenceDatabase database, String[] idxNameMap,
//                                   BufferedWriter regionWriter) throws IOException{
////        int localAlignedSV = 0;
//
//        for (Integer id : unLinkedPattern){
//            String linkType = "nolinks";
//
//            PseudoSequentialPattern pattern = mergedPatterns.get(id);
//
//            int[] arpBasedEstimateInfo = pattern.unlinkedArpPatternPosEstimate(database, 20);
//
//            pattern.doCrossLink(database);
//
//            int[] oemEstimatePos = pattern.oemPatternPosEstimate(database);
//            int[] crossLinkInfo = pattern.getCrossInfo();
//
//            boolean isCrossSup = pattern.isCrossSup();
//
//            if (isCrossSup && pattern.getCrossedLen() >= 20){
//                linkType += ";Cross";
//                Candidate sv = new Candidate(crossLinkInfo[0], crossLinkInfo[1], linkType, pattern.toTypeString(database));
//                sv.writeUnlinkedSV(pattern, regionWriter, arpBasedEstimateInfo[2], idxNameMap);
//            }
//            else if (oemEstimatePos[0] > 0 || oemEstimatePos[1] > 0) {
//                linkType += ";OEM";
//                Candidate sv = new Candidate(oemEstimatePos[0], oemEstimatePos[1], linkType, pattern.toTypeString(database));
//                sv.writeUnlinkedSV(pattern, regionWriter, oemEstimatePos[2], idxNameMap);
//            }
//        }
//    }

    private void addNonDuplicatedLinks(List<Link> links, Link alink) {
        for (Link thisLink : links) {
            if (thisLink.identicalLink(alink)) {
                return;
            }
        }
        links.add(alink);
    }

    private Map<Integer, List<Link>> findMaximalLinks(Map<Integer, List<Link>> patternLinks, Set<Integer> names) {

        List<Integer> linkNames = new ArrayList<>();
        Map<Integer, Integer> nodeNameIdxMap = new HashMap<>();
        int idx = 0;
        for (Integer name : names) {
            linkNames.add(name);
            nodeNameIdxMap.put(name, idx);
            idx ++;
        }
        LinkGraph lg = new LinkGraph(linkNames);

        for (Entry<Integer, List<Link>> entry : patternLinks.entrySet()) {
            for (Link alink : entry.getValue()){
                if (alink.isSelfLink()) {
                    continue;
                }
                int patternName = alink.getPatternIdx();
                int mateName = alink.getMatePatternIdx();
                try {
                    lg.addEdge(nodeNameIdxMap.get(patternName), nodeNameIdxMap.get(mateName));
                }
                catch (Exception e) {
                    System.out.println(patternName + " " + mateName);
                }
            }
        }
        lg.connectedComponents();
        List<List<Integer>> maxComps = lg.getComps();
        Map<Integer, List<Link>> newLinks = new HashMap<>();
        for (List<Integer> comp : maxComps) {
            List<Link> thisNameLinks = new ArrayList<>();
            for (int i = 0; i < comp.size(); i++) {
                int name = comp.get(i);
                if (patternLinks.containsKey(name)) {

                    thisNameLinks.addAll(patternLinks.get(name));

                }
            }
            newLinks.put(comp.get(0), thisNameLinks);
        }

        return newLinks;
    }

}