/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;

/**
 * A class to keep identify the sequence.
 * @author jiadonglin
 */
public class NodeSeqIdentifier {
    
    private final int seqID;
    private final int subSeqStart;
    // Item index in the sequence in terms of itemset index.
    
    private final int itemSetIdx;
    private final int itemIdx;
    // Allow sequence identifier of strings.
//    private String sampleName;

    public NodeSeqIdentifier(int seqID, int subseqStart, int itemsetIdx, int itemIdx) {
        this.seqID = seqID;
        this.subSeqStart = subseqStart;
        this.itemSetIdx = itemsetIdx;
        this.itemIdx = itemIdx;
    }
    public int getSeqID(){
        return this.seqID;
    }
    public int getSubSeqID(){
        return this.subSeqStart;
    }
    public int getItemSetIdx(){
        return this.itemSetIdx;
    }
    public int getItemIdx(){
        return this.itemIdx;
    }
    public Node getSuperItem(SequenceDatabase database){
        return database.getSequenceByID(seqID).getItemsets().get(itemSetIdx).get(itemIdx);
    }
    
    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(seqID);
        sb.append(",");
        sb.append(itemSetIdx);
        sb.append(",");
        sb.append(itemIdx);
        return sb.toString();
    }
}
