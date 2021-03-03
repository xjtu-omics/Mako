/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import structures.SequenceDatabase;
import fspm.PseudoSequentialPattern;
import java.io.BufferedWriter;
import java.io.IOException;

import java.util.*;
import java.util.Map.Entry;

/**
 *
 * @author jiadonglin
 */
public class Candidate {

    private String linkEvi;
    private int start;
    private int end;
    private String patternStr;
    private String patternOri;
    private List<Cluster> clusters;

    private String chrom;
    private Set<String> bpTypes;


    // Involve both inside and between links, this might indicate complex SVs
    public Candidate(String linkStr, PseudoSequentialPattern pattern,
                     Map<Integer, PseudoSequentialPattern> matePatterns, List<Link> links, SequenceDatabase database) {
        linkEvi = linkStr;
        bpTypes = new HashSet<>();
        estimateBpFromLinks(pattern, matePatterns, links, database);
    }


    // Construct for unlinked patterns
//    public Candidate(int s, int e, String link, String pattern) {
//        linkEvi = link;
//        start = s;
//        end = e;
//        patternStr = pattern;
//
//    }

    /**
     * Writer detected events to output.
     * @param pattern
     * @param writer
     * @param idxNameMap
     * @throws IOException
     */
    public void writeSV(PseudoSequentialPattern pattern, BufferedWriter writer,
                        String[] idxNameMap) throws IOException{

        StringBuilder sb = new StringBuilder();
        chrom = idxNameMap[pattern.getPatternChromId()];
        String patternInfo = getPatternInfo(pattern);

        if (end - start > 50){
            sb.append(getCandidateInfo());
            sb.append("\t");
            sb.append(linkEvi);
            sb.append("\t");
            sb.append(patternInfo);
            sb.append("\t");
            sb.append(patternOri);
            writer.write(sb.toString());
            writer.newLine();
        }

    }

    /**
     * Get the coordinates, breakpoint type of a candidate derived from its corresponding pattern.
     * @return
     */
    private String getCandidateInfo() {
        StringBuilder sb = new StringBuilder();
        sb.append(chrom);
        sb.append("\t");
        sb.append(start);
        sb.append("\t");
        sb.append(end);

        if (bpTypes.size() == 0) {
            sb.append("\t");
            sb.append("BND");
        }
        else{
            StringBuilder bpSb = new StringBuilder();
            if (bpTypes.size() > 1){
                for (String bp : bpTypes) {
                    if (bp.equals("BND")) {
                        continue;
                    }
                    bpSb.append(bp);
                    bpSb.append(",");
                }
                sb.append("\t");
                sb.append(bpSb.substring(0, bpSb.length() - 1));
            }
            else{
                List<String> listBp = new ArrayList<>(bpTypes);
                sb.append("\t");
                sb.append(listBp.get(0));
            }

        }

        return sb.toString();
    }


    /**
     * Get SV pattern attributes for output.
     * @param pattern
     * @return
     */
    private String getPatternInfo(PseudoSequentialPattern pattern) {
        StringBuilder sb = new StringBuilder();
        int crossSup = pattern.getCrossedLen();
        sb.append("cxs=");
        sb.append(assignScore());
        if (crossSup != 0){
            sb.append(";cr=");
            sb.append(crossSup);
        }
        for (Cluster cluster : clusters) {
            sb.append(";");
            sb.append(cluster.toString());
        }
        sb.append("\t");
        sb.append(patternStr);

        return sb.toString();
    }


    /**
     * Estimate the breakpoint of each candidate events from discovered links
     * @param pattern
     * @param matePatterns
     * @param links
     * @param database
     */
    private void estimateBpFromLinks(PseudoSequentialPattern pattern, Map<Integer, PseudoSequentialPattern> matePatterns,
                                     List<Link> links, SequenceDatabase database) {

        List<Link> allLinks = new ArrayList<>();
        StringBuilder patternStrSb = new StringBuilder();
        patternStrSb.append(pattern.toTypeString(database));

        StringBuilder patternOriSb = new StringBuilder();
        patternOriSb.append(pattern.getPatternOris());

        // Between subgraph links dose not exist
        if (matePatterns == null || matePatterns.size() == 0) {

            if (pattern.hasSplitAlign()) {
                allLinks.addAll(pattern.getSplitLinksInPattern());
            }
            if (pattern.hasArpLinks()) {
                allLinks.addAll(pattern.getArpLinksInPattern());
            }
        }
        else{
            allLinks.addAll(links);
            if (pattern.hasArpLinks()) {
                allLinks.addAll(pattern.getArpLinksInPattern());
            }
            if (pattern.hasSplitAlign()) {
                allLinks.addAll(pattern.getSplitLinksInPattern());
            }

            for (Entry<Integer, PseudoSequentialPattern> entry : matePatterns.entrySet()) {
                PseudoSequentialPattern matePattern = entry.getValue();
                patternStrSb.append("<>");
                patternStrSb.append(matePattern.toTypeString(database));

                patternOriSb.append(",");
                patternOriSb.append(matePattern.getPatternOris());

                if(matePattern.hasSplitAlign()) {
                    allLinks.addAll(matePattern.getSplitLinksInPattern());
                }
                if(matePattern.hasArpLinks()) {
                    allLinks.addAll(matePattern.getArpLinksInPattern());
                }
            }
        }

        patternStr = patternStrSb.toString();
        patternOri = patternOriSb.toString();
        clusters = resolveLinks(allLinks);

        List<Integer> bps = new ArrayList<>();
        for (Cluster cluster : clusters) {
            bpTypes.addAll(cluster.getBpTypes());
            cluster.computeClusterStats();
            bps.add((int)cluster.getMean()[0]);
            bps.add((int)cluster.getMean()[1]);
        }
//        List<Integer> bps = getHighConfBp();
        Collections.sort(bps);

        start = bps.get(0);
        end = bps.get(bps.size() - 1);
    }

    /**
     * A help function to cluster similar links with heirachical clustering
     * @param links
     * @return
     */
    private List<Cluster> resolveLinks(List<Link> links) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < links.size(); i++) {
            Link link = links.get(i);
            Cluster cluster = new Cluster(link);
            cluster.addVector(link);
            cluster.setMean(new double[]{link.getLinkedItemPos()[0], link.getLinkedItemPos()[1]});
            clusters.add(cluster);
        }

        boolean merged = false;

        do {
            merged = mergeClosestCluster(clusters);
        } while(merged);

        return clusters;
    }

    private boolean mergeClosestCluster(List<Cluster> clusters) {
        Cluster clusterToMergeOne = null;
        Cluster clusterToMergeTwo = null;
        double minDist = Integer.MAX_VALUE;

        for (int i = 0; i < clusters.size(); i++) {
            for (int j = i + 1; j < clusters.size(); j++) {
                // calculate the distance between i and j
                double distance = calculateDistance(clusters.get(i).getMean(), clusters.get(j).getMean());
                // if the distance is less than the max distance allowed
                // and if it is the smallest distance until now
                if (distance < minDist && distance <= 100) {
                    // record this pair of clusters
                    minDist = distance;
                    clusterToMergeOne = clusters.get(i);
                    clusterToMergeTwo = clusters.get(j);
                }
            }
        }

        // if no close clusters were found, return false
        if (clusterToMergeOne == null) {
            return false;
        }

        // else, merge the two closest clusters
        for(Link link: clusterToMergeTwo.getLinks()){
            clusterToMergeOne.addVector(link);
        }
        // after mergint, we need to recompute the mean of the resulting cluster
        clusterToMergeOne.recomputeClusterMean();
        // we delete the cluster that was merged
        clusters.remove(clusterToMergeTwo);

        return true;
    }

    private double calculateDistance(double[] var1, double[] var2){
        double span1 = var1[1] - var1[0];
        double span2 = var2[1] - var2[0];

        double spanRegionDiff = Math.abs(span1 - span2);

        double center1 = (var1[0] + var1[1]) / 2;
        double center2 = (var2[0] + var2[1]) / 2;

        double curMin = Math.min(Math.abs(var1[0] - var2[0]), Math.abs(var1[1] - var2[1]));
        double positionDistDiffMin = Math.min(Math.abs(center1 - center2), curMin);

//        int typeDiff = type1.equals(type2) ? 1 : 100000;

        return (spanRegionDiff + positionDistDiffMin);
    }

    private int assignScore(){
        int numClusters = clusters.size();
        Set<String> nodeType = new HashSet<>();
        if (patternStr.contains("<>")){
            String[] patternTokens = patternStr.split("<>");
            for (String pattern : patternTokens) {
                for (String type : pattern.split(",")){
                    nodeType.add(type);
                }
            }
        }else{
            for (String type : patternStr.split(",")){
                nodeType.add(type);
            }
        }
        return numClusters * nodeType.size();
    }
}
