package utils;
import java.util.List;

public class IntervalTree {
    private IntervalNode root;

    private class IntervalNode {
        int nodeIdx;
        IntervalNode left;
        int start;
        int end;
        int maxEnd;
        IntervalNode right;

        public IntervalNode(int idx, IntervalNode left, int start, int end, int maxEnd, IntervalNode right) {
            this.nodeIdx = idx;
            this.left = left;
            this.start = start;
            this.end = end;
            this.maxEnd = maxEnd;
            this.right = right;
        }

        @Override
        public String toString(){
            return nodeIdx + "," + left + "-" + right;
        }
    }

    /**
     * Adds an interval to the the calendar
     *
     * @param start the start of interval
     * @param end   the end of the interval.
 *    @param idx    An index of the added interval
     */
    public void add (int start, int end, int idx) {
        if (start > end) throw new IllegalArgumentException("The end " + end + " should be greater than start " + start);

        IntervalNode inode = root;
        while (inode != null) {
            inode.maxEnd = (end > inode.maxEnd) ? end : inode.maxEnd;
            if (start < inode.start) {
                if (inode.left == null) {
                    inode.left = new IntervalNode(idx,null, start, end, end, null);
                    return;
                }
                inode = inode.left;
            } else {
                if (inode.right == null) {
                    inode.right = new IntervalNode(idx, null, start, end, end, null);
                    return;
                }
                inode = inode.right;
            }
        }
        root =  new IntervalNode(idx,null, start, end, end, null);
    }

    /**
     * Tests if the input interval overlaps with the existing intervals.
     *
     * Rules:
     * 1.  If interval intersects return true. obvious.
     * 2.  if (leftsubtree == null || leftsubtree.max <=  low) go right
     * 3.  else go left
     *
     * @param start     the start of the interval
     * @param end       the end of the interval
     * return           true if overlap, else false.
     */
    public void overlap(int start, int end, List<Integer> overlapped) {
        if (start > end) throw new IllegalArgumentException("The end " + end + " should be greater than start " + start);

        IntervalNode intervalNode = root;

        while (intervalNode != null) {
            boolean isOverlap = intersection(start, end, intervalNode.start, intervalNode.end);
            if (isOverlap) {
                overlapped.add(intervalNode.nodeIdx);
//                return true;
            }

            if (goLeft(start, end, intervalNode.left)) {
                intervalNode = intervalNode.left;
            } else {
                intervalNode = intervalNode.right;
            }
        }
//        return false;
    }

    /**
     * Returns if there is an intersection in the two intervals
     * Two intervals such that one of the points coincide:
     * eg: [10, 20] and [20, 40] are NOT considered as intersecting.
     */
    private boolean intersection (int start, int end, int intervalStart, int intervalEnd) {
        return start <= intervalEnd && end >= intervalStart;
    }

    private boolean goLeft(int start, int end, IntervalNode intervalLeftSubtree) {
        return intervalLeftSubtree != null && intervalLeftSubtree.maxEnd > start;
    }

}
