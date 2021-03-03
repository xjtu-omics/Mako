/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matcher;

import java.util.List;

/**
 *
 * @author jiadonglin
 */
public class AlignInfo {
    
    // keep the pos of maximum and minimum unique string match at read [forMin, forMax, revMin, revMax]
    private final int[] uMatchIdxAtRead = new int[4]; 
    
    // keep the alignment pos at ref [forMatchIdx, revMatchIdx]
    private final int[] matchIdxAtRef = new int[2];
    
    // keep the mis-match index at read [forIdx, revIdx]
    private final int[] misMatchIdxAtRead = new int[2];
    private String uMaxForMatchStr = "";
    private String uMaxRevMatchStr = "";   
    
        
    public AlignInfo(){
        
    }
        
    
    public void updateMaxMatchedStrInfo(String maxStr, List<List<Integer>> refPosList, int baseAtReadIdx, boolean readForwarAlign){
        
        if (refPosList.get(0).size() == 1 && refPosList.get(1).isEmpty()){
            if (readForwarAlign){
                matchIdxAtRef[0] = refPosList.get(0).get(0) - 10;
                uMatchIdxAtRead[1] = baseAtReadIdx;
                uMaxForMatchStr = maxStr.length() > uMaxForMatchStr.length() ? maxStr : uMaxForMatchStr;
            }else{
                matchIdxAtRef[1] = refPosList.get(0).get(0) - 10;
                uMatchIdxAtRead[3] = baseAtReadIdx;
                uMaxRevMatchStr = maxStr.length() > uMaxRevMatchStr.length() ? maxStr : uMaxRevMatchStr;
            }
        }
        if (refPosList.get(0).isEmpty() && refPosList.get(1).size() == 1){
            if (readForwarAlign){
                matchIdxAtRef[0] = refPosList.get(1).get(0) - 10;
                uMatchIdxAtRead[1] = baseAtReadIdx;
                uMaxForMatchStr = maxStr.length() > uMaxForMatchStr.length() ? maxStr : uMaxForMatchStr;
            }else{
                matchIdxAtRef[1] = refPosList.get(1).get(0) - 10;
                uMatchIdxAtRead[3] = baseAtReadIdx;
                uMaxRevMatchStr = maxStr.length() > uMaxRevMatchStr.length() ? maxStr : uMaxRevMatchStr;
            }            
        }        
    }
    
    public void updateMatchInfoAtRead(List<List<Integer>> refPosList, List<List<Integer>> newRefPosList, 
            int baseAtReadIdx, boolean readForwardAlign){
        
        if (refPosList.get(0).size() > 0 && newRefPosList.get(0).isEmpty() && newRefPosList.get(1).size() == 1){
            if (readForwardAlign){
                misMatchIdxAtRead[0] = baseAtReadIdx;
            }else{
                misMatchIdxAtRead[1] = baseAtReadIdx;             
            }
        }        
        
        if ((refPosList.get(0).size() > 1 && newRefPosList.get(0).size() == 1) 
                || (refPosList.get(1).size() > 1 && newRefPosList.get(1).size() == 1)){
            if (readForwardAlign){
                uMatchIdxAtRead[0] = baseAtReadIdx;
            }else{
                uMatchIdxAtRead[2] = baseAtReadIdx;
            }
        }
    }
    public int[] getMaxMatchIdxAtRef(){
        if (uMatchIdxAtRead[1] - misMatchIdxAtRead[0] < 5){
            matchIdxAtRef[0] -= (uMatchIdxAtRead[1] - misMatchIdxAtRead[0] + 1);
        }
        if (uMatchIdxAtRead[3] - misMatchIdxAtRead[0] < 5){
            matchIdxAtRef[1] += (uMatchIdxAtRead[3] - misMatchIdxAtRead[1] + 1);
        }
        return matchIdxAtRef;
    }
    public int[] getMaxMatchIdxAtRead(){
        if (uMatchIdxAtRead[1] - misMatchIdxAtRead[0] < 5){
            uMatchIdxAtRead[1] -= (uMatchIdxAtRead[1] - misMatchIdxAtRead[0] + 1);
        }
        if (uMatchIdxAtRead[3] - misMatchIdxAtRead[0] < 5){
            uMatchIdxAtRead[3] -= (uMatchIdxAtRead[3] - misMatchIdxAtRead[1] + 1);
        }
        return uMatchIdxAtRead;
    }
    public String getReadForMaxMatchStr(){
        return uMaxForMatchStr;
    }
    public int getReadForMaxMatchLength(){
        return uMaxForMatchStr.length();
    }
      
    
    public void printMatchStatus(String strand){
        System.out.println("\n"+strand);
        System.out.println("Max str match at ref: [" + matchIdxAtRef[0] + " : " + matchIdxAtRef[1] + "]");
        System.out.println("Left to right: " + uMaxForMatchStr);
        System.out.println("Right to left: " + uMaxRevMatchStr);
        System.out.println("Matched str at read: [" + uMatchIdxAtRead[0] + " " + uMatchIdxAtRead[1] + " " + uMatchIdxAtRead[2] + " " + uMatchIdxAtRead[3] +"]");
    }
}
