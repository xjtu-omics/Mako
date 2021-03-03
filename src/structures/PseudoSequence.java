/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;

import java.util.*;
/**
 *
 * @author jiadonglin
 */
public class PseudoSequence {
    
    // keep the orignial sequence id of this pseudo-sequence
    private int seqID;
//    private int subSeqStart;
    // pseudo-sequence corresponding sequence size in terms of itemset.
    private int seqSize;
    /**
     * Keep the start of the pseudo-sequence while doing projection of length-1 prefix.
     * This will be used to check pattern span range while doing pattern growth.  
     */

    private int genomeCordStart;
    
    
    // properties of the pseudo-sequence
    private int firstItemset;
    private int firstItem;
    private int size;
    
    public PseudoSequence(){
        
    }
    /**
     * Create a new pseudo sequence from a pseudo sequence.
     * @param sequence 
     * @param indexItemset
     * @param indexItem 
     */
    public PseudoSequence(PseudoSequence sequence, int indexItemset, int indexItem){
        
        this.seqID = sequence.seqID;
        this.seqSize = sequence.seqSize;
//        this.subSeqStart = sequence.firstItemset;
        
        this.genomeCordStart = sequence.genomeCordStart;
        
        this.size = sequence.size - indexItemset;
        this.firstItemset = indexItemset + sequence.firstItemset;

        
        
        if(this.firstItemset == sequence.firstItemset){
            this.firstItem = indexItem + sequence.firstItem;
        }else{
            this.firstItem = indexItem; 
        }
        
    }
    /**
     * Create a new pseudo sequence from a sequence
     * @param sequence
     * @param indexItemset
     * @param indexItem 
     */
    public PseudoSequence(Sequence sequence, int indexItemset, int indexItem, int subseqStart){
        this.firstItem = indexItem;
        this.firstItemset = indexItemset;
        this.seqSize = sequence.size();
        this.seqID = sequence.getId();
//        this.subSeqStart = subseqStart;
        
        this.size = sequence.size() - firstItemset;        
        if (this.size == 1 && sequence.getItemsets().get(firstItemset).isEmpty()){
            this.size = 0;
        }
    }
    public int getId(){
        return this.seqID;
    }

    public int getSize(){
        return this.size;
    }
    public int getOriSeqSize(){
        return this.seqSize;
    }
    public void setGenomeStartPos(int posStart){
        genomeCordStart = posStart;
    }
    public int getFirstItemsetIdx(){
        return firstItemset;
    }
    
    public int getFirstItemIdx(){
        return firstItem;
    }
    public int getGenomeStartPos(){
        return genomeCordStart;
    }
    /**
     * Return true if an itemset is a postfix. Eg, <(_A)(C)(B)>. 
     * @param indexItemset
     * @return 
     */
    public boolean isPostfix(int indexItemset){
        return indexItemset == 0 && firstItem != 0;
    }
    public boolean isFirstItemset(int index){
        return index == 0;
    }
    public boolean isLastItemset(int index) {
        return (index + firstItemset) == this.seqSize -1;
    }
    
    public Node getItemAtItemsetAt(int indexItem, int indexItemset, Sequence sequence){
        if(isFirstItemset(indexItemset)){
            return sequence.get(indexItemset + firstItemset).get(indexItem + firstItem);
        }else{
            return sequence.get(indexItemset + firstItemset).get(indexItem);
        }
    }
    
    public int indexOf(int indexItemset, String item, Sequence sequence){
        List<Node> itemset = sequence.get(indexItemset );
        for(int i = 0; i < itemset.size() ; i++){
            Node superitemAtIdx = getItemAtItemsetAt(i, indexItemset, sequence);
            String superitemType = superitemAtIdx.type;
            if(superitemType.equals(item)){
                return i ;
            }
        }
        return -1;
    }
    public int getSizeOfItemsetAt(int idx, Sequence sequence){
        int size = sequence.getItemsets().get(idx + firstItemset).size();
        if (isFirstItemset(idx)){
            size -= firstItem;
        }
        return size;
    }
}
