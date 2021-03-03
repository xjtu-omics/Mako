package utils;

import java.util.*;

public class LinkGraph {
    private List<Integer> nodes;
    private int V;
    LinkedList<Integer>[] adjListArray;
    private List<List<Integer>> comps;

    // constructor
    public LinkGraph(List<Integer> nodeNames) {
        nodes = nodeNames;
        V = nodes.size();
        // define the size of array as
        // number of vertices
        adjListArray = new LinkedList[V];
        comps = new ArrayList<>();

        // Create a new list for each vertex
        // such that adjacent nodes can be stored

        for(int i = 0; i < V ; i++){
            adjListArray[i] = new LinkedList<>();
        }
    }

    public void addEdge(int src, int dest) {
        // Add an edge from src to dest.

        adjListArray[src].add(dest);

        // Since graph is undirected, add an edge from dest
        // to src also
        adjListArray[dest].add(src);
    }

    private void depthFirstSearch(int v, boolean[] visited, List<Integer> node) {
        // Mark the current node as visited and print it
        visited[v] = true;
        node.add(nodes.get(v));
//        System.out.print(nodes[v]+" ");
        // Recur for all the vertices
        // adjacent to this vertex
        for (int x : adjListArray[v]) {
            if(!visited[x]) depthFirstSearch(x, visited, node);
        }
    }
    public void connectedComponents() {
        // Mark all the vertices as not visited
        boolean[] visited = new boolean[V];
        for(int v = 0; v < V; ++v) {
            List<Integer> currComp = new ArrayList<>();
            if(!visited[v]) {
                // print all reachable vertices
                // from v
                depthFirstSearch(v, visited, currComp);
//                System.out.println();
                comps.add(currComp);
            }
        }
    }

    public List<List<Integer>> getComps() {
        return comps;
    }

}
