package collector;

import htsjdk.samtools.SAMRecord;

/**
 * User: jiadonglin
 * Date: 2020/9/8
 */

public class DiscordantRp implements Comparable<DiscordantRp>{

    //  variable to keep the read information

    protected SAMRecord base;
    // attributes of each signal

    protected String signalType;
    protected int mutPos1;
    protected int mutPos2;
    protected String mutSignalOri;
    protected int insertSize;
    protected String cigarString;
    protected int mapQ;
    protected String bpType; // Breakpoint type: DEL, DUP, INV, INS
    
    protected boolean isARP = false;
    protected boolean isInterChrom = false;
    
    //    protected String splitAlignChr;

    @Override
    public int compareTo(DiscordantRp otherDiscRpSignal){
        return mutPos1 - otherDiscRpSignal.mutPos1;
    }

    public DiscordantRp(int pos1, int pos2, String signalType) {
        mutPos1 = pos1;
        mutPos2 = pos2;
        this.signalType = signalType;
        this.insertSize = pos2 - pos1;

    }

    public DiscordantRp(SAMRecord record, String signalType, String bpType, String ori, int pos1, int pos2){

        base = record;
        mutSignalOri = ori;

        isInterChrom = !record.getReferenceName().equals(record.getMateReferenceName());

        mutPos1 = pos1;
        mutPos2 = pos2;

        insertSize = Math.abs(record.getInferredInsertSize());

        this.signalType = signalType;
        this.bpType = bpType;
        isARP = signalType.contains("ARP");

        cigarString = record.getCigarString();
        mapQ = record.getMappingQuality();

    }

    @Override
    public boolean equals(Object obj){
        if (obj instanceof DiscordantRp){
            DiscordantRp discRpSignal = (DiscordantRp) obj;
            return (discRpSignal.getMutSignalType().equals(this.signalType));
        }else{
            return false;
        }
    }
    @Override
    public int hashCode(){
        return signalType.hashCode();
    }


    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(signalType);
        sb.append("\t");
        sb.append(mutPos1);
        sb.append("\t");
        sb.append(mutPos2);
        sb.append("\t");
        sb.append(insertSize);
        return sb.toString();
    }
    
    
    public String getqName(){
        return base.getReadName();
    }
    public String getMutSignalType(){
        return signalType;
    }

    public int getReadLen() {return base.getReadLength(); }

}
