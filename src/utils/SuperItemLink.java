/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import fspm.PseudoSequentialPattern;
import fspm.PseudoSuperItem;
import structures.Node;
import structures.SequenceDatabase;
import htsjdk.samtools.QueryInterval;

import java.util.*;
//import java.util.HashMap;

/**
 *
 * @author jiadonglin
 */
public class SuperItemLink {
    
    int supLink;
    
    public SuperItemLink(){
        
    }
    
    public int getSupLink(){
        return supLink;
    }
    
    private List<String> getUniqueQnames(String[] qnames){
        List<String> uniqueNames = new ArrayList<>(qnames.length);
        Set<String> nameSet = new HashSet<>();
        for (String name : qnames){            
            if (!nameSet.contains(name)){
                nameSet.add(name);
                uniqueNames.add(name);
            }
        }
        return uniqueNames;
    }
    
    public boolean twoSuperItemLinkCheck(Node nodeOne, Node nodeTwo){
        boolean supported = false;

        if (!nodeOne.getType().equals(nodeTwo.getType())) {
            return false;
        }

        List<String> qnameOneList = getUniqueQnames(nodeOne.getQNames());
        List<String> qnameTwoList = getUniqueQnames(nodeTwo.getQNames());
        Set<String> uniqueQName = new HashSet<>();
               
        int supportedARPs = 0;
        for (String qname : qnameOneList){  
            uniqueQName.add(qname);
        }       
        for (String qname : qnameTwoList){              
            if (uniqueQName.contains(qname)){
                supportedARPs += 1;
                if (supportedARPs >= 1 ){
                    supported = true;                   
                }
            }
            uniqueQName.add(qname);
        }
        supLink = supportedARPs;

        return supported;
    }

    
    public int mateSuperItemSearch(List<PseudoSuperItem> superItems, PseudoSuperItem targetSuperItem,
                                   SequenceDatabase database){
        QueryInterval targetInterval = targetSuperItem.getSuperItem(database).getSuperitemMateRegion();
        int length = superItems.size();
        int mateIndex = -1;

        if (length == 1){
            QueryInterval curInterval = superItems.get(0).getSuperItem(database).getSuperitemRegion();
            if (hasOverlap(curInterval, targetInterval)){
                Node itemOne = targetSuperItem.getSuperItem(database);
                Node itemTwo = superItems.get(0).getSuperItem(database);
                boolean isEnoughARPs = twoSuperItemLinkCheck(itemOne, itemTwo);
                if (isEnoughARPs) {
                    return 0;
                }else{
                    return -1;
                }
            }else{
                return -1;
            }
        }
        List<Integer> matchedSuperItemIdx = mateSearchHelper(targetInterval, superItems);
        if (matchedSuperItemIdx.size() > 0) {
            int maxSup = 0;
            for (Integer id : matchedSuperItemIdx) {
                Node itemOne = targetSuperItem.getSuperItem(database);
                Node itemTwo = superItems.get(id).getSuperItem(database);
                boolean isEnoughARPs = twoSuperItemLinkCheck(itemOne, itemTwo);

                if (isEnoughARPs && supLink > maxSup) {
                    maxSup = supLink;
                    mateIndex = id;
                }
            }
            supLink = maxSup;

        }
        return mateIndex;
    }

    private List<Integer> mateSearchHelper(QueryInterval target, List<PseudoSuperItem> superitems) {
        List<Integer> matched = new ArrayList<>();

        IntervalTree intervalSearchTree = new IntervalTree();
        for (int i = 0; i < superitems.size(); i++) {
            PseudoSuperItem thisItem = superitems.get(i);
            int start = thisItem.getInterval()[0];
            int end = thisItem.getInterval()[1];
            intervalSearchTree.add(start, end, i);
        }
        intervalSearchTree.overlap(target.start, target.end, matched);
        return matched;
    }

    private boolean hasOverlap(QueryInterval aInterval, QueryInterval targetInterval){       
        return aInterval.overlaps(targetInterval);
    }
    private boolean isAheadInterval(QueryInterval aInterval, QueryInterval targetInterval){        
        return aInterval.start > targetInterval.end;
    }

    private boolean isAfterInterval(QueryInterval aInterval, QueryInterval targetInterval){
        return aInterval.end < targetInterval.start;
    }
    
    public List<PseudoSequentialPattern> minusItemAndCopy(List<PseudoSequentialPattern> patterns, int index){
        List<PseudoSequentialPattern> removedPatterns = new ArrayList<>();
        int length = patterns.size();
        for (int i = 0; i < length ; i ++){
            if (i != index){
                removedPatterns.add(patterns.get(i));                
            }
        }
        return removedPatterns;
    }
    public List<PseudoSequentialPattern> minusItemAndCopyWithIndexMap(List<PseudoSequentialPattern> patterns, int index, Map<Integer, Integer> indexMap){
        List<PseudoSequentialPattern> removedPatterns = new ArrayList<>();
        int length = patterns.size();
        for (int i = 0; i < length ; i ++){
            PseudoSequentialPattern pattern = patterns.get(i);
            if (i != index && pattern.hasArpSuperItem()){
                removedPatterns.add(pattern);  
                indexMap.put(removedPatterns.size() - 1, i);
            }
        }
        return removedPatterns;
    }
}
