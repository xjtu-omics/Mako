/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;

import java.util.*;
/**
 * This class represents a pair of value (1) item (2) a boolean indicating
 * if the item is contained in an itemset that was cut or not.
 * @author jiadonglin
 */

public class Pair {
    
    protected final String item;
    protected final boolean postfix;
    
    
    private List<NodeSeqIdentifier> itemAppears = new ArrayList<NodeSeqIdentifier>();
    
    public Pair(boolean postfix, String item) {
        this.item = item;
        this.postfix = postfix;
    }
        
    public boolean equals (Object object){
        Pair paire = (Pair) object;
        if((paire.postfix == this.postfix)&& (paire.item.equals(this.item))){
            return true;
        }
        return false;
    }
    /**
    * Method to calculate an hashcode (because pairs are stored in a map).
    */
    public int hashCode()
    {// Ex: 127333,P,X,1  127333,N,Z,2
           // transform it into a string
           StringBuilder r = new StringBuilder();
           r.append((postfix ? 'P' : 'N')); // the letters here have no meanings. they are just used for the hashcode
           r.append(item);
           // then use the hashcode method from the string class
           return r.toString().hashCode();
    }

   /**
    * Check if this is the case of the item appearing in a postfix
    * @return true if this is the case.
    */
    public boolean isPostfix() {
        return postfix;
    }

   /**
    * Get the item represented by this pair
    * @return the item.
    */
    public String getItem() {
           return item;
    }

   /**
    * Get the support of this item (the number of sequences 
    * containing it).
    * @return the support (an integer)
    */

    public int getCount(){
       return itemAppears.size();
    }
   /**
    * Get the list of sequence IDs associated with this item.
    * @return  the list of sequence IDs.
    */
    public List<NodeSeqIdentifier> getItemAppear() {
       return itemAppears;
    }

    public void addItemAppearIdx(int seqID, int subseqStart, int itemsetIdx, int itemIdx){
       itemAppears.add(new NodeSeqIdentifier(seqID, subseqStart, itemsetIdx, itemIdx));
    }
    
    public void printItemAppears(){
        for (NodeSeqIdentifier itemIdentity : itemAppears){
            System.out.println("type: " + item + " subseq: " + itemIdentity.getSubSeqID() 
            + " itemsetIdx: " + itemIdentity.getItemSetIdx());
        }
    }
}
