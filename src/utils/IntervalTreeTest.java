package utils;
import java.util.List;
import java.util.ArrayList;

public class IntervalTreeTest {

    public static void main(String[] args){
        IntervalTree intervalSearchTree = new IntervalTree();
        intervalSearchTree.add(4, 10,3);
        intervalSearchTree.add(5, 8, 1);
        intervalSearchTree.add(7, 10,5);
        intervalSearchTree.add(15, 18,4);
        intervalSearchTree.add(15, 22,6);
        intervalSearchTree.add(17, 19, 0);

        intervalSearchTree.add(21, 24,2);

        List<Integer> overlappedNodes = new ArrayList<>();
        intervalSearchTree.overlap(6, 16, overlappedNodes);
        for (Integer idx : overlappedNodes) {
            System.out.println(idx);
        }
    }

}
