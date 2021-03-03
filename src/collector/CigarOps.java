/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package collector;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import java.util.List;

/**
 *
 * @author jiadonglin
 */
public class CigarOps {
    
    private String cigarString;
    private int qsPos;
    // indicate co-occurence of I and D on a single read
    private boolean coID = false;  
    private int[] clippedLength;
    
    public CigarOps(){
        
    }
        
    public CigarOps(SAMRecord record){
        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        CigarElement firstElement = cigarElements.get(0);
        CigarElement lastElement = cigarElements.get(cigarElements.size() - 1);
        clippedLength = new int[2];
        if (firstElement.getOperator().toString().equals("S") && !lastElement.getOperator().toString().equals("S")){
            clippedLength[0] = firstElement.getLength();
        }
        if (!firstElement.getOperator().toString().equals("S") && lastElement.getOperator().toString().equals("S")){
            clippedLength[1] = lastElement.getLength();
        }
        else if (firstElement.getOperator().toString().equals("S") && lastElement.getOperator().toString().equals("S")){
            clippedLength[0] = firstElement.getLength();
            clippedLength[1] = lastElement.getLength();
        }
    }
    public int[] getClippedStatus(){
        return clippedLength;
    }
    public String getCigarStr(){
        return this.cigarString;
    }
    public int getqsPos(){
        return this.qsPos;
    }

    public int getLongestD(List<CigarElement> cigarElements){
        int operators = cigarElements.size();
        int longestD = 0;
        for (int i = 0; i < operators; i++){
            CigarElement cigarElement = cigarElements.get(i);
            String operation = cigarElement.getOperator().toString();                                    
            
            int opLength = cigarElement.getLength();
            if (operation.equals("D") && opLength > longestD){
                longestD = opLength;
            }
        }
        return longestD;
    }
    
    public void calQueryPosFromCigar(List<CigarElement> cigarElements, int nIgnore, boolean isReverse, int readLen){

        int qsPos = 0;        
        String cigarStr = "";
        int IDnum = 0;  // number of 'I' 'D' operation in cigar
        int lastIDindex = -1;
        int firstIDindex = -1;
        int longestMatch = 0;
        int operators = cigarElements.size();
        for (int i = 0; i < operators; i++){
            CigarElement cigarElement = cigarElements.get(i);
            String operation = cigarElement.getOperator().toString();
                                    
            operation = operation.equals("H") ? "S":operation;
            int opLength = cigarElement.getLength();                         
            longestMatch = (operation.equals("M") && opLength > longestMatch) ? opLength : longestMatch; 
            if (i == 0){                
                cigarStr += operation;
                if (operation.equals("M")){
                    qsPos += opLength;
                }                
            }
            else if (operation.equals("I")||operation.equals("D")){
                
                lastIDindex = i;
                if (opLength <= nIgnore){
                    qsPos += opLength;
                }else{
                    cigarStr += operation;
                    IDnum += 1;                    
                }                    
                
                if (firstIDindex == -1) {
                    firstIDindex = i;
                }
            }else if (!cigarStr.isEmpty()){
                Character lastChar = cigarStr.charAt(cigarStr.length() - 1);
                Character curChar = operation.charAt(0);
                // = 0, two char same
                if (curChar.compareTo(lastChar) == 0){
                    continue;
                }else if (opLength > nIgnore){
                    cigarStr += operation;
                }
            }
        }
 
        // If alignment of a read has multiple 'I' and 'D', but still have long matched part.
        if (IDnum > 1 && longestMatch > 0.5 * readLen){
            coID = true;
        }
        

        this.cigarString = cigarStr; 
        this.qsPos = qsPos;
    }
    /**
     * 'I' 'D' operation occur together in Cigar, which is not good.
     * @return 
     */
    public boolean isCoIDread(){
        return coID;
    }    
    
}
