/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matcher;

/**
 *
 * @author jiadonglin
 */
public class StringIdentity {
    int seqId;
    String seqType;        

    public StringIdentity(int seqId, String seqType) {
        this.seqId = seqId;
        this.seqType = seqType;
    }
    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(seqType);
        sb.append(":");
        sb.append(seqId);
        return sb.toString();
    } 
    public int getSeqId(){
        return seqId;
    }
}
