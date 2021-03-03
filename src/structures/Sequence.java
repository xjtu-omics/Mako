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
public class Sequence {
    
    /**
     * A sequence is like <(ARP_LARGE_INSERT)(MS)(SM)(ARP_LARGE_INSERT)>
     * There might be several superitems in a set
     */
    private List<List<Node>> itemsets = new ArrayList<List<Node>>();
    /** Each sequence has an identical id */
    private int id;
    
    private List<Node> smallIndelSuperitem = new ArrayList<Node>();
    
    /**
     * Sequence constructor
     * @param id 
     */
    public Sequence(int id){
        this.id = id;
    }
    /**
     * Add a new superitem-set to the sequence.
     * @param itemset 
     */
    public void addItemset(List<Node> itemset){
        itemsets.add(itemset);
    }
    /**
     * Return a string representation of this sequence
     */
    public String sequence2String(){
        StringBuilder r = new StringBuilder("");
        for (List<Node> itemset : itemsets){
            r.append('(');
            for (Node item : itemset){
                r.append(item.type);
                r.append(' ');
            }
            r.append(')');
        }
        return r.append("  ").toString();
    }
    public void printSequence(){
        System.out.println(sequence2String());
    }
    public int getId(){
        return id;
    }
    /** Return all superitem-sets in the sequence */
    public List<List<Node>> getItemsets(){
        return itemsets;
    }
    public Sequence getSubSequence(int itemsetIdx){
        Sequence newSequence = new Sequence(getId());
        List<List<Node>> newItemsets = itemsets.subList(itemsetIdx, itemsets.size());
        newSequence.itemsets = newItemsets;
        return newSequence;
    }
    
    
    /** Get an superitem-sets at a given position in the sequence*/
    public List<Node> get(int idx){
        return itemsets.get(idx);
    }
    /** Get the size of this sequence, number of superitem-sets*/
    public int size(){
        return itemsets.size();
    }
    
    public void setSmallIndelSuperitem(List<Node> alist){
        smallIndelSuperitem = alist;
    }
    public List<Node> getLastSet(){
        return itemsets.get(itemsets.size());
    }
    
    public Node superItemAtPos(int indexItemset, int indexItem){
        List<Node> superitemSet = get(indexItemset);
        return superitemSet.get(indexItem);
    }
    public Sequence removeInFreSuperitemFromSequence(Map<String, Integer> superitemTypeCount, double relativeMinSup){
        Sequence sequence = new Sequence(getId());
        for (List<Node> itemset : itemsets){
            List<Node> newItemset = removeInFreSuperitemFromItemset(itemset, superitemTypeCount, relativeMinSup);
            if(newItemset.size() != 0){
                sequence.addItemset(newItemset);
            }
        }
        return sequence;
    }
    public List<Node> removeInFreSuperitemFromItemset(List<Node> itemset, Map<String, Integer> superitemTypeCount, double relativeMinSup){
        List<Node> newItemset = new ArrayList<Node>();
        for (Node superitem : itemset){
            String superitemType = superitem.type;
            int support = superitemTypeCount.get(superitemType);
            if (support >= relativeMinSup){
                newItemset.add(superitem);
            }
        }
        return newItemset;
    }
    
}
