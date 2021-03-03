/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import structures.Node;
import structures.SequenceDatabase;

/**
 *
 * @author jiadonglin
 */
public class PseudoSuperItem implements Comparable<PseudoSuperItem>{
    int sequenceId;
    protected int itemsetIdx;
    protected int itemIdx;   
    protected int superItemLeftPos;
    String superitemType;
    int[] interval;
    int[] mateInterval;

    public PseudoSuperItem(int seqId, int itemsetIndex, int itemIndex){
        sequenceId = seqId;
        itemsetIdx = itemsetIndex;
        itemIdx = itemIndex;
    }
    public Node getSuperItem(SequenceDatabase database){
        Node si = database.getSequenceByID(sequenceId).superItemAtPos(itemsetIdx, itemIdx);
        return si;
    }
    public void setPsSuperitemInfo(SequenceDatabase database){
        superItemLeftPos = getSuperItem(database).getPos();
        superitemType = getSuperItem(database).getType();
        interval = getSuperItem(database).getInterval();
        mateInterval = getSuperItem(database).getMateInterval();
    }
    @Override
    public int compareTo(PseudoSuperItem otherItem){
        return superItemLeftPos - otherItem.superItemLeftPos;
    }

    public int[] getInterval(){
        return interval;
    }
    public int[] getMateInterval(){
        return mateInterval;
    }

}
