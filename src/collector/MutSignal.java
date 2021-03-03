    /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package collector;

import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

    /**
 *
 * @author jiadonglin
 */


public class MutSignal implements Comparable<MutSignal>{

    //  variable to keep the read information
    protected String contig1;
    protected String contig2;
    protected int pos1;
    protected int pos2;
    protected String ori1;
    protected String ori2;

    // attributes of each signal
    protected String queryName;
    protected int signalRefIndex;
    protected String mutSignalType;
    protected int mutSignalPos;
    protected String mutSignalOri;
    protected int insertSize;
    protected String cigarString;
    protected int mapQ;

    protected boolean isizeNormal = true;
    protected boolean isARP = false;
    protected boolean isSplitAlign = false;
    protected boolean isInterChrom = false;
    protected int splitAlignPos;
//    protected String splitAlignChr;
    protected String readSequence;
    protected int[] clippedStatus;
    protected int longestD;
    
    public MutSignal(){
        
    }

    @Override
    public int compareTo(MutSignal otherMutSignal){
        return mutSignalPos - otherMutSignal.getMutPos();
    }
    
    public MutSignal(SAMRecord record, String signalType, int pos, String ori){
        
        
        queryName = record.getReadName();        
        signalRefIndex = record.getReferenceIndex();
        contig1 = record.getReferenceName();
        contig2 = record.getMateReferenceName();

        isInterChrom = !contig1.equals(contig2);

        pos1 = record.getAlignmentStart();
        pos2 = record.getMateAlignmentStart();

        ori1 = record.getReadNegativeStrandFlag() ? "rev" : "fwd";
        ori2 = record.getMateNegativeStrandFlag() ? "rev" : "fwd";
        insertSize = Math.abs(record.getInferredInsertSize());

        mutSignalType = signalType;
        isARP = mutSignalType.contains("ARP");

        mutSignalPos = pos;
        mutSignalOri = ori;

        cigarString = record.getCigarString();
        readSequence = record.getReadString();
        mapQ = record.getMappingQuality();
    
        CigarOps cigarOps = new CigarOps(record);
        clippedStatus = cigarOps.getClippedStatus();
        longestD = cigarOps.getLongestD(record.getCigar().getCigarElements());


        if (record.isSecondaryOrSupplementary()){
            splitAlignPos = decodeSplitAlign(record.getAttribute("SA").toString(), mapQ);            
        }
    }

    @Override
    public boolean equals(Object obj){
        if (obj instanceof MutSignal){
            MutSignal mutSignal = (MutSignal) obj;
            return (mutSignal.getMutSignalType().equals(this.mutSignalType));
        }else{
            return false;
        }
    }
    @Override
    public int hashCode(){
        return mutSignalType.hashCode();
    }
    
    
    public boolean withinDistance(MutSignal otherMutSignal, int maxDist){
        int disDiff = Math.abs(otherMutSignal.getMutPos() - mutSignalPos);
//        int spanDiff = Math.abs(otherMutSignal.getSpanned() - getSpanned());
//        if (mutSignalType.contains("ARP")){
//            disDiff += spanDiff;
//        }

        int typePenalty = 0;
        if (!mutSignalType.equals(otherMutSignal.getMutSignalType())){
            typePenalty = 1;
        }
        return disDiff <= maxDist / (1 + typePenalty);
    }

    @Override
    public String toString(){     
        StringBuilder sb = new StringBuilder();
        sb.append(queryName);
        sb.append("\t");
        sb.append(mutSignalType);
        sb.append("\t");
        sb.append(mutSignalPos);
        sb.append("\t");
        sb.append(pos1);
        sb.append("\t");
        sb.append(contig1);
        sb.append("\t");
        sb.append(insertSize);
        return sb.toString();
    }
    public int getLongestD(){
        return longestD;
    }    
    public int[] getClippedStatus(){
        return clippedStatus;
    }
    public String getqName(){
        return new String(queryName);
    }
    public String getMutSignalType(){
        return mutSignalType;
    }
    public String getMutSignalOri(){
        return mutSignalOri;
    }
    public String getReadOri() {return ori1 + ori2;}

    public int getSignalChromIdx(){
        return signalRefIndex;
    }

    public String getSignalRef(){
        return contig1;
    }
    public int getRecordPos(){
        return pos1;
    }
    public int getMateRecordPos(){
        return pos2;
    }
    public int getInsertSize() {
        return insertSize;
    }
    public boolean isARPSignal(){
        return isARP;
    }
    public boolean isIsizeNormal(){
        return isizeNormal;
    }
    public boolean isSplitAlign(){
        return isSplitAlign;
    }
    public boolean isInterChrom(){
        return isInterChrom;
    }
    public int getMutPos(){
        return mutSignalPos;
    }
    public int getMapQ(){
        return mapQ;
    }

//    public String getSplitAlignInfo(){
//        return splitAlignInfo;
//    }
    public int getSplitAlignPos(){
        return splitAlignPos;
    }
//    public String getSplitAlignChr(){
//        return splitAlignChr;
//    }
    public String getReadSequence(){
        return readSequence;
    }
    public void setIsizeNormal(int isizeUpper, int isizeLower){
        if (this.insertSize > isizeUpper || this.insertSize < isizeLower){
            isizeNormal = false;
        }
    }

    /**
     * We only consider split alignment at the same chromosome
     * @param splitAlignString
     * @param mapQ
     * @return
     */
    private int decodeSplitAlign(String splitAlignString, int mapQ){
        
        String[] tokens = splitAlignString.split(",");        
        String splitAlignContig = tokens[0];
        int pos = Integer.parseInt(tokens[1]);
        
        if (splitAlignContig.equals(contig1) && mapQ >= 20){
            isSplitAlign = true;
            List<EasyCigar> cigarOps = getSplitAlignPosFromCigarString(tokens[3]);
            EasyCigar firstOp = cigarOps.get(0);
            if (firstOp.getOp().equals("M")){
                pos += firstOp.getOpLength();
            }      
        }else{
            pos = -1;
        }          
        return pos;
    }


    private List<EasyCigar> getSplitAlignPosFromCigarString(String cigar){
       List<EasyCigar> cigarOps = new ArrayList<>();
       Pattern cigarPattern = Pattern.compile("[0-9]+[MIDNSHP]");
       Matcher cigarOpStrings = cigarPattern.matcher(cigar);

       while (cigarOpStrings.find()){
           String opString = cigarOpStrings.group();
           int strLen = opString.length();
           String op = opString.substring(strLen - 1, strLen);
           int opLength = Integer.parseInt(opString.substring(0, strLen - 1));

           EasyCigar cigarop = new EasyCigar(op, opLength);
           cigarOps.add(cigarop);

       }

       return cigarOps;
    }

    class EasyCigar {
        String op;
        int opLength;

        EasyCigar(String op, int opLength) {
            this.op = op;
            this.opLength = opLength;
        }

        String getOp(){
            return op;
        }
        int getOpLength(){
            return opLength;
        }
    }
}

