/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matcher;

import fspm.PseudoSuperItem;
import structures.Node;
import structures.SequenceDatabase;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.HashMap;
import java.util.Map.Entry;

import utils.Candidate;

/**
 * This is a class for scoring the string matching status based on k-mer.
 * @author jiadonglin
 */
public class StringMatcher {
    
    // Cross match used variables
    Map<String, List<char[]>> stringsDB;
    sharedString forwardSharedStrings;
    sharedString reverseSharedStrings;

    // Local alignment used variables
    String maxForStrMatch;
    String maxRevStrMatch;
    
    int[] maxMatchAtRef = new int[]{-1,-1};
    int[] maxMatchAtRead = new int[]{-1,-1,-1,-1};    
    char[] baseMap;     
    
    final int protecStrLength = 10;
    final int maxMisMatch = 1;
    
    
    boolean isPlusMatch = false;
    boolean isMinusMatch = false;
    
    AlignInfo plusAlignInfo;
    AlignInfo minusAlignInfo;
    
    
    public StringMatcher(){
        baseMap = new char[256];
        baseMap[(int)'A'] = 'T';
        baseMap[(int)'T'] = 'A';
        baseMap[(int)'C'] = 'G';
        baseMap[(int)'G'] = 'C';
    }    
           
     public void doStringAlign(List<String> readStrs, List<PseudoSuperItem> superitems, int refRegionLeft, String refStr, List<Candidate> alignedInfos, SequenceDatabase database){
        for (int i = 0; i < readStrs.size(); i++){
            String curRead = readStrs.get(i);
            alignStrToRef(curRead, refStr);            
            if (maxMatchAtRef[0] != -1 && maxMatchAtRef[1] != -1){
                if (maxMatchAtRead[1] + maxMatchAtRead[3] + 2 >= curRead.length()){
                    // AlignInfo -> [startPosAtRef, endPosAtRef, supSuperItemIdx]
                    
                    int[] adjustedBpAtRead = bpAdjust(maxMatchAtRead[0], maxMatchAtRead[1],maxMatchAtRead[2],maxMatchAtRead[3], curRead.length());
//                    System.out.println(adjustedBpAtRead[0] + " " + adjustedBpAtRead[1]);
                    
                    int newLeftBpAtRef;
                    int newRightBpAtRef;
                    // left Bp coords adjust
                    if (adjustedBpAtRead[0] < maxMatchAtRead[1]){
                        newLeftBpAtRef = maxMatchAtRef[0] - (maxMatchAtRead[1] - adjustedBpAtRead[0]);
                    }else{
                        newLeftBpAtRef = maxMatchAtRef[0] + (adjustedBpAtRead[0] - maxMatchAtRead[1]);
                    }
                    
                    // right bp coords adjust
                    if (adjustedBpAtRead[1] < maxMatchAtRead[3]){
                        newRightBpAtRef = maxMatchAtRef[1] + (maxMatchAtRead[3] - adjustedBpAtRead[1]);
                    }else{
                        newRightBpAtRef = maxMatchAtRef[1] - (adjustedBpAtRead[1] - maxMatchAtRead[3]);
                    }
                    Node superitem = superitems.get(i).getSuperItem(database);
                    
                    if (newLeftBpAtRef > newRightBpAtRef){
                        continue;
                    }
                    if (superitem.getWeight() < 5 || superitem.getNumMinusRead() < 1 || superitem.getNumPlusRead() < 1){
                        continue;
                    }
                    
                            
//                    svOutInfo curSV = new svOutInfo(refRegionLeft + newLeftBpAtRef, refRegionLeft + newRightBpAtRef, 10, 
//                            superitem.getWeight(), superitem.getRatio(), superitem.getNumPlusRead(), superitem.getNumMinusRead());
//                    if (!alignedInfos.isEmpty()){
//                        boolean exist = false;
//                        for (svOutInfo preSV : alignedInfos){
//                            if (curSV.identicalSV(preSV)){
//                                exist = true;
//                            }
//                        }
//                        if (!exist){
//                            alignedInfos.add(curSV);
//                        }
//                    }else{
//                        alignedInfos.add(curSV);
//                    }
                    
                }
                
            }                                                                        
        }
                        
    }
    
    
    private String createRepeatStr(){
        return new String(new char[protecStrLength]).replace("\0", "N");
    }
    
    
    public void alignStrToRef(String readStr, String refStr){
//        int[] alignedPos = new int[]{-1,-1,-1,-1};
        plusAlignInfo = new AlignInfo();
        minusAlignInfo = new AlignInfo();
        // SMS superitem does not save read string yet
        if (("").equals(readStr)){
            return;
        }
        
        // add "N" to ref at prefix and suffix
        
        String protectStr = createRepeatStr();                
        StringBuilder sb = new StringBuilder(protectStr);       
        sb.append(refStr);
        sb.append(protectStr);
        char[] refArray = sb.toString().toCharArray(); 
        
        // Align read from left to right                                
//        System.out.println(sb.toString());
        char[] readArray = readStr.toCharArray();                       
        
        List<List<Integer>> baseMatchBothStrandAtRef = new ArrayList<>(maxMisMatch + 2);        
        
        findBaseAtRef(refArray, readArray[0], baseMatchBothStrandAtRef);
        
        String plusPrefix = Character.toString(readArray[0]);
        String minusPrefix = Character.toString(baseMap[(int)readArray[0]]);
                
        expandMatchApprox(plusPrefix, minusPrefix, baseMatchBothStrandAtRef, baseMatchBothStrandAtRef, 1, 
                readArray, refArray, true);
                
        
        // Align read from right to left
        sb = new StringBuilder(readStr);
        sb.reverse();
        char[] readRevArray = sb.toString().toCharArray();        
        baseMatchBothStrandAtRef.clear();
        findBaseAtRef(refArray, readRevArray[0], baseMatchBothStrandAtRef);
        
        plusPrefix = Character.toString(readRevArray[0]);
        minusPrefix = Character.toString(baseMap[(int)readRevArray[0]]);
        
        expandMatchApprox(plusPrefix, minusPrefix, baseMatchBothStrandAtRef, baseMatchBothStrandAtRef, 1, 
                readRevArray, refArray, false);
        
        // Matched string length >= 30bp
        if (plusAlignInfo.getReadForMaxMatchLength() >= 30 || minusAlignInfo.getReadForMaxMatchLength() >= 30){
            if (plusAlignInfo.getReadForMaxMatchLength() > minusAlignInfo.getReadForMaxMatchLength()){
                maxMatchAtRef = plusAlignInfo.getMaxMatchIdxAtRef();
                maxMatchAtRead = plusAlignInfo.getMaxMatchIdxAtRead();            
                isPlusMatch = true;
            }
            else{
                maxMatchAtRef = minusAlignInfo.getMaxMatchIdxAtRef();
                maxMatchAtRead = minusAlignInfo.getMaxMatchIdxAtRead();            
                isMinusMatch = true;
            }
        }
                
        
    }    
         
    
    private void expandMatchApprox(String plusPrefix, String minusPrefix, List<List<Integer>> plusPosList, List<List<Integer>> minusPosList,
            int baseAtReadIdx, char[] readArray, char[] refArray, boolean readForwardAlign){
        
        if (baseAtReadIdx >= readArray.length){
            return;
        }
        
        if (plusPosList.get(0).isEmpty() && plusPosList.get(1).isEmpty() && minusPosList.get(0).isEmpty() && minusPosList.get(1).isEmpty()){
            return;
        }
                
        
        char baseAtRead = readArray[baseAtReadIdx];
        char comBaseAtRead = baseMap[(int)readArray[baseAtReadIdx]];
        
        List<List<Integer>> newPlusPosList = new ArrayList<>();
        List<List<Integer>> newMinusPosList = new ArrayList<>();
        
        for (int i = 0; i < plusPosList.size(); i++){
            newPlusPosList.add(new ArrayList<>(plusPosList.get(i).size()));
            newMinusPosList.add(new ArrayList<>(minusPosList.get(i).size()));
        }
        
        // plus strand alignment
        for (int i = 0; i < plusPosList.size(); i++){                       
            if (readForwardAlign) {
                for (Integer idx : plusPosList.get(i)){
                    if (refArray[idx + 1] == baseAtRead){                        
                        newPlusPosList.get(i).add(idx + 1);
                    }else{
                        switch (i){
                            case 0:
                                newPlusPosList.get(1).add(idx + 1);
                                break;
                            case 1:
                                newPlusPosList.get(2).add(idx + 1); 
                                break;
                        }                        
                    }
                }
            }else{
                for (Integer idx : plusPosList.get(i)){
                    if (refArray[idx - 1] == baseAtRead){
                        newPlusPosList.get(i).add(idx - 1);
                    }else{
                        switch (i){
                            case 0:
                                newPlusPosList.get(1).add(idx - 1);
                                break;
                            case 1:
                                newPlusPosList.get(2).add(idx - 1);
                                break;
                        }                        
                    }
                }
            }
        }
        
        // minus strand alignment
        for (int i = 0; i < minusPosList.size(); i ++){
            if (readForwardAlign){
                for (Integer idx : minusPosList.get(i)){
                    if (baseMap[(int)refArray[idx - 1]] == comBaseAtRead){
                        newMinusPosList.get(i).add(idx - 1);
                    }else{
                        switch (i){
                            case 0:
                                newMinusPosList.get(1).add(idx - 1);
                                break;
                            case 1:
                                newMinusPosList.get(2).add(idx - 1);
                                break;
                        }
                    }
                }
            }else{
                for (Integer idx : minusPosList.get(i)){
                    if (baseMap[(int)refArray[idx + 1]] == comBaseAtRead){
                        newMinusPosList.get(i).add(idx + 1);
                    }else{
                        switch (i){
                            case 0:
                                newMinusPosList.get(1).add(idx + 1);
                                break;
                            case 1:
                                newMinusPosList.get(2).add(idx + 1);
                                break;
                        }                        
                    }
                }
            }
        }
        
               
        int levelZeroPlusMatch = newPlusPosList.get(0).size();
        int levelOnePlusMatch = newPlusPosList.get(1).size();
//        int levelTwoPlusMatch = newPlusPosList.get(2).size();
        
        boolean plusAbleToExtend = levelZeroPlusMatch + levelOnePlusMatch > 0 ? true : false;
        
        int levelZeroMinusMatch = newMinusPosList.get(0).size();
        int levelOneMinusMatch = newMinusPosList.get(1).size();
//        int levelTwoMinusMatch = newMinusPosList.get(2).size();
        
        boolean minusAbleToExtend = levelZeroMinusMatch + levelOneMinusMatch > 0 ? true : false;
        
        
        if (plusAbleToExtend && minusAbleToExtend){
            plusAlignInfo.updateMatchInfoAtRead(plusPosList, newPlusPosList, baseAtReadIdx, readForwardAlign);
            minusAlignInfo.updateMatchInfoAtRead(minusPosList, newMinusPosList, baseAtReadIdx, readForwardAlign);
            
            plusAlignInfo.updateMaxMatchedStrInfo(plusPrefix + baseAtRead, newPlusPosList, baseAtReadIdx, readForwardAlign);
            minusAlignInfo.updateMaxMatchedStrInfo(minusPrefix + comBaseAtRead, newMinusPosList, baseAtReadIdx, readForwardAlign);
            expandMatchApprox(plusPrefix + baseAtRead, minusPrefix + comBaseAtRead, newPlusPosList, newMinusPosList, 
                    baseAtReadIdx + 1, readArray, refArray, readForwardAlign);
            
        }
        
        else if (plusAbleToExtend && !minusAbleToExtend){
            
            plusAlignInfo.updateMatchInfoAtRead(plusPosList, newPlusPosList, baseAtReadIdx, readForwardAlign);
            
            plusAlignInfo.updateMaxMatchedStrInfo(plusPrefix + baseAtRead, newPlusPosList, baseAtReadIdx, readForwardAlign);
            minusAlignInfo.updateMaxMatchedStrInfo(minusPrefix, minusPosList, baseAtReadIdx, readForwardAlign);
            
            expandMatchApprox(plusPrefix + baseAtRead, minusPrefix, newPlusPosList, minusPosList, baseAtReadIdx + 1, 
                    readArray, refArray, readForwardAlign);
        }
        else if (!plusAbleToExtend && minusAbleToExtend){
            
            minusAlignInfo.updateMatchInfoAtRead(minusPosList, newMinusPosList, baseAtReadIdx, readForwardAlign);
            
            plusAlignInfo.updateMaxMatchedStrInfo(plusPrefix, plusPosList, baseAtRead, readForwardAlign);
            minusAlignInfo.updateMaxMatchedStrInfo(minusPrefix + comBaseAtRead, newMinusPosList, baseAtReadIdx, readForwardAlign);
            expandMatchApprox(plusPrefix, minusPrefix + comBaseAtRead, plusPosList, newMinusPosList, baseAtReadIdx + 1, 
                    readArray, refArray, readForwardAlign);
        }
        else{
            expandMatchApprox(plusPrefix, minusPrefix, newPlusPosList, newMinusPosList, baseAtReadIdx + 1, 
                    readArray, refArray, readForwardAlign);
        }
    }
    
    
    private int[] bpAdjust(int forMin, int forMax, int revMin, int revMax, int readLen){
        int[] bpPos = new int[]{forMax, revMax};
        // read reverse align is longer than forward align
        boolean isFound = false;
        for (int i = forMin;i<=forMax;i++){            
            for (int j = revMax;j>=revMin;j--){
                if (i + j == readLen){
                    bpPos[0] = i;
                    bpPos[1] = j;
                    isFound = true;
                    break;
                }
            }
            if (isFound){
                break;
            }           
        }
        return bpPos;
    }
    
//    private int bpLeftShift(int leftPosAtRef, int rightPosAtRef, String refStr){
//        int shiftedLeft = leftPosAtRef;
//        return shiftedLeft;
//    }
    
    private void findBaseAtRef(char[] refArray, char base, List<List<Integer>> baseMatchPos){       
        int len = refArray.length;
        for (int i = 0; i < maxMisMatch + 2; i++){
            baseMatchPos.add(new ArrayList<>());
        }
        for (int i = 0; i < len; i++){  
            if (refArray[i] == base){
                baseMatchPos.get(0).add(i);
            }
        }        
    }
                   
    public void strCrossMatch(List<String> mStrings, List<String> sForwardStrings, List<String> sReverseStrings){
        ReadStringDB mDB = new ReadStringDB(mStrings, "M");
        ReadStringDB sForwardDB = new ReadStringDB(sForwardStrings, "SF");
        ReadStringDB sReverseDB = new ReadStringDB(sReverseStrings, "SR");
        
        stringsDB = new HashMap<>();
        
        stringsDB.put("M", mDB.sequences);
        stringsDB.put("SF", sForwardDB.sequences);
        stringsDB.put("SR", sReverseDB.sequences);
        
        List<StringIdentity> initialStringIds = new ArrayList<>();
        for (Entry<String, List<char[]>> entry : stringsDB.entrySet()){
            List<char[]> seqList = entry.getValue();
            for (int i = 0; i < seqList.size(); i++){
                initialStringIds.add(new StringIdentity(i, entry.getKey()));
            }
        }
        
        List<PseudoRead> initialProjectDB = mDB.pseduoDB;
        initialProjectDB.addAll(sForwardDB.pseduoDB);
        initialProjectDB.addAll(sReverseDB.pseduoDB);
        
        char initialPrefix[] = new char[]{'A','T','C','G'};
        for (char prefix : initialPrefix){            
            
            String base = Character.toString(prefix);
            List<PseudoRead> projectedDB = buildProjectedDB(base, initialProjectDB);
            
                        
            depthFirstSearch(base, projectedDB, true, true, initialStringIds);
        }
    } 
    
    
    private void depthFirstSearch(String prefix, List<PseudoRead> database, boolean isForward, boolean isBoth, List<StringIdentity> strList){
        BaseMatchInfo commonCharInfo = mostCommonCharInProjectDB(database);
        char sharedBase = commonCharInfo.base;
        isForward &= commonCharInfo.isForward;
        isBoth &= commonCharInfo.isBoth;
        
        if (sharedBase != 'N'){
            String newPrefix = prefix + sharedBase;
            List<PseudoRead> projectDB = buildProjectedDB(newPrefix, database);
            strList = commonCharInfo.baseStringIdentitys;
            depthFirstSearch(newPrefix, projectDB, commonCharInfo.isForward, commonCharInfo.isBoth, strList);
            
        }else{       
            saveString(isForward, isBoth, prefix, strList);                        
        }               
    }
    
    private List<PseudoRead> buildProjectedDB(String growthStr, List<PseudoRead> pseudoDB){
        List<PseudoRead> newProjectDB = new ArrayList<>();
        if (growthStr.length() == 1){            
            for (PseudoRead psRead : pseudoDB){
                char[] oriSequence = stringsDB.get(psRead.dbId).get(psRead.seqId);
                for (int i = 0; i < oriSequence.length; i++){
                    int idxOfItemInSeq = psRead.indexOf(growthStr, oriSequence, i);
                    if (idxOfItemInSeq != -1 && idxOfItemInSeq != psRead.seqSize - 1){
                        newProjectDB.add(new PseudoRead(psRead, idxOfItemInSeq + 1));
                    }
                }
            }

        }else{
            char lastChar = growthStr.charAt(growthStr.length() - 1);
            String base = Character.toString(lastChar);
            for (PseudoRead psRead : pseudoDB){
                char[] oriSequence = stringsDB.get(psRead.dbId).get(psRead.seqId);
                int idxOfItemInSeq = psRead.indexOf(base, oriSequence, psRead.firstBaseId);
                if (idxOfItemInSeq != -1 && idxOfItemInSeq != oriSequence.length - 1){
                    newProjectDB.add(new PseudoRead(psRead, 1));
                }
            }
        }
        
        return newProjectDB;
    }
    
    
    private BaseMatchInfo mostCommonCharInProjectDB(List<PseudoRead> projectedDB){
        int[] aCount = new int[3];
        int[] tCount = new int[3];
        int[] cCount = new int[3];
        int[] gCount = new int[3];
        
        List<StringIdentity> aAppearList = new ArrayList<>();
        List<StringIdentity> tAppearList = new ArrayList<>();
        List<StringIdentity> cAppearList = new ArrayList<>();
        List<StringIdentity> gAppearList = new ArrayList<>();
        
        for (PseudoRead pseudo : projectedDB){
            char ch = pseudo.getBase(pseudo, stringsDB);
            String dbId = pseudo.dbId;
            if (dbId.equals("M")){
                switch(ch){
                    case 'A':
                        aCount[0] += 1;
                        aAppearList.add(new StringIdentity(pseudo.seqId, "M"));
                        break;
                    case 'T':
                        tCount[0] += 1;
                        tAppearList.add(new StringIdentity(pseudo.seqId, "M"));
                        break;
                    case 'C':
                        cCount[0] += 1;
                        cAppearList.add(new StringIdentity(pseudo.seqId, "M"));
                        break;
                    case 'G':
                        gCount[0] += 1;
                        gAppearList.add(new StringIdentity(pseudo.seqId, "M"));
                }
            }
            if (dbId.equals("SF")){
                switch(ch){
                    case 'A':
                        aCount[1] += 1;
                        aAppearList.add(new StringIdentity(pseudo.seqId, "SF"));
                        break;
                    case 'T':
                        tCount[1] += 1;
                        tAppearList.add(new StringIdentity(pseudo.seqId, "SF"));
                        break;
                    case 'C':
                        cCount[1] += 1;
                        cAppearList.add(new StringIdentity(pseudo.seqId, "SF"));
                        break;
                    case 'G':
                        gCount[1] += 1;
                        gAppearList.add(new StringIdentity(pseudo.seqId, "SF"));
                }
            }
            if (dbId.equals("SR")){
                switch(ch){
                    case 'A':
                        aCount[2] += 1;
                        aAppearList.add(new StringIdentity(pseudo.seqId, "SR"));
                        break;
                    case 'T':
                        tCount[2] += 1;
                        tAppearList.add(new StringIdentity(pseudo.seqId, "SR"));
                        break;
                    case 'C':
                        cCount[2] += 1;
                        cAppearList.add(new StringIdentity(pseudo.seqId, "SR"));
                        break;
                    case 'G':
                        gCount[2] += 1;
                        gAppearList.add(new StringIdentity(pseudo.seqId, "SR"));
                }
            }
            
        }
        boolean isForward = true;
        boolean isBoth = true;
        if (!(aCount[0] != 0 && (aCount[1] != 0|| aCount[2] != 0))){
            aCount = new int[3];
            
        }
        if (!(tCount[0] != 0 && (tCount[1] != 0|| tCount[2] != 0))){
            tCount = new int[3];
            
        }
        if (!(cCount[0] != 0 && (cCount[1] != 0|| cCount[2] != 0))){
            cCount = new int[3];
            
        }
        if (!(gCount[0] != 0 && (gCount[1] != 0|| gCount[2] != 0))){
            gCount = new int[3];
            
        }
        
        int aSum = aCount[0] + aCount[1] + aCount[2];
        int tSum = tCount[0] + tCount[1] + tCount[2];
        int cSum = cCount[0] + cCount[1] + cCount[2];
        int gSum = gCount[0] + gCount[1] + gCount[2];
        
        int maxSum = Math.max(Math.max(aSum, tSum), Math.max(cSum, gSum));
        
        char commonChar = 'N';
        BaseMatchInfo matchInfo = new BaseMatchInfo(commonChar, isForward, isBoth);
        if (maxSum == 0){
            return matchInfo;
        }
        if (maxSum == aSum){            
            commonChar = 'A';
            matchInfo.setBase(commonChar);
            matchInfo.setBaseStringIDs(aAppearList);
            if (aCount[1] == 0){
                isForward = false;                  
            }
            if (aCount[1] == 0 || aCount[2] == 0){
                isBoth = false;
            }
            matchInfo.setIsForward(isForward); 
            matchInfo.setIsBoth(isBoth);
        }
        if (maxSum == tSum){
            commonChar = 'T';
            matchInfo.setBase(commonChar);
            matchInfo.setBaseStringIDs(tAppearList);
            if (tCount[1] == 0){
                isForward = false;
            }
            if (tCount[1] == 0 || tCount[2] == 0){
                isBoth = false;
            }
            matchInfo.setIsForward(isForward); 
            matchInfo.setIsBoth(isBoth);
        }
        if (maxSum == cSum){
            commonChar = 'C';
            matchInfo.setBase(commonChar);
            matchInfo.setBaseStringIDs(cAppearList);
            if (cCount[1] == 0){
                isForward = false;
            }
            if (cCount[1] == 0 || cCount[2] == 0){
                isBoth = false;
            }
            matchInfo.setIsForward(isForward); 
            matchInfo.setIsBoth(isBoth);
        }
        if (maxSum == gSum){
            commonChar = 'G';
            matchInfo.setBase(commonChar);
            matchInfo.setBaseStringIDs(gAppearList);
            if (gCount[1] == 0){
                isForward = false;
            }
            if (gCount[1] == 0 || gCount[2] == 0){
                isBoth = false;
            }
            matchInfo.setIsForward(isForward); 
            matchInfo.setIsBoth(isBoth);
        }
        return matchInfo;
    }
    
    
    
    
    /**
     * estimate the maximum length of a common string shared by a set of strings.
     * @param numOfStrs
     * @param avgOfStrLength
     * @return 
     */
    public int estimateStrLength(int numOfStrs, int avgOfStrLength){
        int length = 0;
        
        
        
        return length;
    }
    
    public int[] isCrossLinked(){                       
//        List<Integer> mCount = new ArrayList<>();
//        List<Integer> sCount = new ArrayList<>();
        
        int[] linkInfo = new int[]{0,0,0};
        
        if (reverseSharedStrings != null && forwardSharedStrings != null){
            return linkInfo;
        }
        if (forwardSharedStrings != null){
            List<StringIdentity> identitys = forwardSharedStrings.strIdentitys;
            for (StringIdentity id : identitys){
                if (id.seqType.equals("M")){
//                    mCount.add(id.seqId);
                    linkInfo[0] += 1;
                }
                if (id.seqType.equals("SF")){
//                    sCount.add(id.seqId);
                    linkInfo[1] += 1;
                }
            }
//            System.out.println("M size: " + mCount.size() + " S size: " + sCount.size());
//            if (mCount.size() == sCount.size() && mCount.size() == 1){
//                if (mCount.get(0) != sCount.get(0)){                
//                    return true;
//                }
//            }
        }
        if (reverseSharedStrings != null){
            List<StringIdentity> identitys = reverseSharedStrings.strIdentitys;
            for (StringIdentity id : identitys){
                if (id.seqType.equals("M")){
//                    mCount.add(id.seqId);
                    linkInfo[0] += 1;
                }
                if (id.seqType.equals("SR")){
//                    sCount.add(id.seqId);
                    linkInfo[2] += 1;
                }
            }
//            System.out.println("M size: " + mCount.size() + " S size: " + sCount.size());
//            if (mCount.size() == sCount.size() && mCount.size() == 1){
//                if (mCount.get(0) != sCount.get(0)){                
//                    return true;
//                }
//            }
        }
                                           
        return linkInfo;
    }
    
    private void saveString(boolean isForward, boolean isBoth, String prefix, List<StringIdentity> strList){
        int preForwardStringLength = forwardSharedStrings == null? 0 : forwardSharedStrings.strLength;
        int preReverseStringLength = reverseSharedStrings == null? 0 : reverseSharedStrings.strLength;
        if (!isBoth){
            if (isForward){                        
                if (prefix.length() > preForwardStringLength){
                    forwardSharedStrings = new sharedString(strList, prefix);
                }
            }else{                
                if (prefix.length() > preReverseStringLength){
                    reverseSharedStrings = new sharedString(strList, prefix);
                }
            }
        } 
    }
    public List<StringIdentity> getForwardSharedStrIdentitys(){
        List<StringIdentity> ids = new ArrayList<>();
        if (forwardSharedStrings != null){
            ids = forwardSharedStrings.strIdentitys;
        }
        return ids;
    }
    public boolean isForwardExist(){
        boolean exist = false;
        if (forwardSharedStrings != null){
            exist = true;
        }
        return exist;
    }
    public List<StringIdentity> getReverseSharedStrIdentitys(){
        List<StringIdentity> ids = new ArrayList<>();
        if (reverseSharedStrings != null){
            ids = reverseSharedStrings.strIdentitys;
        }
        return ids;
    }
    public String getForwardSharedString(){        
        return forwardSharedStrings.commonString;
    }
    public int getForwardSharedStrLength(){
        String forwardStr = forwardSharedStrings.commonString;
        char[] charArray = forwardStr.toCharArray();
        int len = 1;
        char initialCh = charArray[0];
        for (char ch : charArray){
            if (ch != initialCh){
                len += 1;
                initialCh = ch;
            }
        }
        return len;
    }
    public String getReversedSharedString(){
        return reverseSharedStrings.commonString;
    }
    public int getReverseSharedStrLength(){
        String forwardStr = reverseSharedStrings.commonString;
        char[] charArray = forwardStr.toCharArray();
        int len = 1;
        char initialCh = charArray[0];
        for (char ch : charArray){
            if (ch != initialCh){
                len += 1;
                initialCh = ch;
            }
        }
        return len;
    }
    public int[] getEstimateBp(List<PseudoSuperItem> superItems, SequenceDatabase database){
        int[] pos = new int[2];
        if (forwardSharedStrings != null){
            List<StringIdentity> identitys = forwardSharedStrings.strIdentitys;
            Collections.sort(identitys, new Comparator<StringIdentity>(){
            @Override
                public int compare(StringIdentity o1, StringIdentity o2){
                    return o1.getSeqId() - o2.getSeqId();
                }        
            });
            
            PseudoSuperItem psLeftItem = superItems.get(identitys.get(0).seqId);
            Node nodeLeft = psLeftItem.getSuperItem(database);
            pos[0] = nodeLeft.getPos();
            
            PseudoSuperItem psRightItem = superItems.get(identitys.get(identitys.size() - 1).seqId);
            Node nodeRight = psRightItem.getSuperItem(database);
            pos[1] = nodeRight.getPos();
            
        }        
        if (reverseSharedStrings!= null){
            List<StringIdentity> identitys = reverseSharedStrings.strIdentitys;
            Collections.sort(identitys, new Comparator<StringIdentity>(){
            @Override
                public int compare(StringIdentity o1, StringIdentity o2){
                    return o1.getSeqId() - o2.getSeqId();
                }        
            });
            
            PseudoSuperItem psLeftItem = superItems.get(identitys.get(0).seqId);
            Node nodeLeft = psLeftItem.getSuperItem(database);
            pos[0] = nodeLeft.getPos();
            
            PseudoSuperItem psRightItem = superItems.get(identitys.get(identitys.size() - 1).seqId);
            Node nodeRight = psRightItem.getSuperItem(database);
            pos[1] = nodeRight.getPos();
        } 
        return pos;
    }
    public void printSharedString(List<PseudoSuperItem> superItems, SequenceDatabase database){
        if (forwardSharedStrings != null){
            StringBuilder sb = new StringBuilder();
            for (StringIdentity id : forwardSharedStrings.strIdentitys){
                PseudoSuperItem psItem = superItems.get(id.seqId);
                Node node = psItem.getSuperItem(database);
                sb.append(node.getType());
                sb.append(",");
                sb.append(node.getPos());
                sb.append("\t");
            }
            sb.append(forwardSharedStrings.commonString);
            sb.append("\t");            
//            sb.append(getForwardSharedStrLength());
            System.out.println(sb.toString());
        }        
        if (reverseSharedStrings!= null){
            StringBuilder sb = new StringBuilder();
            for (StringIdentity id : reverseSharedStrings.strIdentitys){
                PseudoSuperItem psItem = superItems.get(id.seqId);
                Node node = psItem.getSuperItem(database);
                sb.append(node.getType());
                sb.append(",");
                sb.append(node.getPos());
                sb.append("\t");
            }
            sb.append(reverseSharedStrings.commonString);
            sb.append("\t");
//            sb.append(getReverseSharedStrLength());
            System.out.println(sb.toString());
        }             
    }
    
    class sharedString{
        List<StringIdentity> strIdentitys = new ArrayList<>();
        String commonString;
        int strLength = 0;
        
        public sharedString(){
            
        }        
        public sharedString(List<StringIdentity> identity, String str) {
            strIdentitys = identity;
            commonString = str;
            strLength = str.length();
        }  
        public boolean inValid(){
            return commonString == null;
        }
        @Override
        public String toString(){
            StringBuilder sb = new StringBuilder();
            for (StringIdentity identity : strIdentitys){
                sb.append(identity.toString());
                sb.append(" ");
            }
            sb.append("\n");
            sb.append(commonString);
            return sb.toString();
        }
    }
}
