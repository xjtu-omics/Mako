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
public class BaseIdentity {
    String dbID;
    int seqID;
    int baseIdx;
    
    
    public BaseIdentity(String dbString, int seqIdx, int baseIdx){
        this.dbID = dbString;
        this.seqID = seqIdx;
        this.baseIdx = baseIdx;
    }
}
