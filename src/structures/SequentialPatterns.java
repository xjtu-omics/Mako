/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;

import java.util.*;
/**
 * Used for save discovered sequential patterns in all length.
 * @author jiadonglin
 */
public class SequentialPatterns {
    
    /** A list of list is used to store all patterns
     At position i, a list of sequential patterns contains all patterns of size i.
     */
    private final List<List<SequentialPattern>> levels = new ArrayList<List<SequentialPattern>>();
    private int patternCount = 0;
    
    public SequentialPatterns(){
        levels.add(new ArrayList<SequentialPattern>());
    }
    
    public void addPatterns(SequentialPattern pattern, int k){
        while(levels.size() <= k){
            levels.add(new ArrayList<SequentialPattern>());
        }
        levels.get(k).add(pattern);
        patternCount ++;
    }
    
}
