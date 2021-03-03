package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

public class TestGraph {
    public static void main(String[] args){
        List<Integer> nodeNames = new ArrayList<>();

        nodeNames.add(5);
        nodeNames.add(6);
        nodeNames.add(7);
        nodeNames.add(8);
        nodeNames.add(9);

        Map<Integer, Integer> nameMaps = new HashMap<>();

        for (int i = 0; i < nodeNames.size(); i++) {
            nameMaps.put(nodeNames.get(i), i);
        }

        LinkGraph lg = new LinkGraph(nodeNames);
        lg.addEdge(nameMaps.get(5), nameMaps.get(6));
        lg.addEdge(nameMaps.get(5), nameMaps.get(5));
        lg.addEdge(nameMaps.get(7), nameMaps.get(8));
        lg.addEdge(nameMaps.get(8), nameMaps.get(9));

        lg.connectedComponents();
        List<List<Integer>> comps = lg.getComps();
        for (List<Integer> comp : comps) {
            for (Integer in : comp) {
                System.out.print(in +" ");
            }
            System.out.println();
        }
    }
}
