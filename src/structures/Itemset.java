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
public class Itemset {
    private final List<String> superitemTypes = new ArrayList<String>();
    
    public Itemset(String item){
        addItem(item);
    }
    
    public void addItem(String superitem){
        if (!superitemTypes.contains(superitem)){
            superitemTypes.add(superitem);
        }
    }
    
    public List<String> getItems(){
        return superitemTypes;
    }
	
    public String get(int index){
        return superitemTypes.get(index);
    }

    public void print(){
            System.out.print(toString());
    }

    public String toString(){
            StringBuilder r = new StringBuilder ();
            for(String superitem : superitemTypes){
                    r.append(superitem);
                    r.append(' ');
            }
            return r.toString();
    }


    public int size(){
            return superitemTypes.size();
    }
}
