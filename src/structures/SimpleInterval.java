/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;

/**
 *
 * @author jiadonglin
 */
public class SimpleInterval implements Comparable<SimpleInterval>{
    int start;
    int end;
    public SimpleInterval(int s, int e){
        start = s;
        end = e;
    }
    @Override
    public int compareTo(SimpleInterval other){
        return this.start - other.start;
    }
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        sb.append(",");
        sb.append(end);
        return sb.toString();
    }
}
