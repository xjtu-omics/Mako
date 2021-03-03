package utils;

public class Link {
    // Property for this link
//    private int linkedItemIdx;
//    private int linkedMateItemIdx;
    private int linkedSups;
    private int linkedLeftPos;
    private int linkedRightPos;
    private int patternIdx;
    private int matePatternIdx;
    private String linkType;
    private String svType;

    // Property for breakpoints derived from this link
//    private int[] bpMapQ;
//    private String[] bpType;
    private double[] bpRatio;
    private int[] bpWeight;

    // Property for split links
    private int splitStatus;

    public Link(int sup, int leftPos, int rightPos, int paIdx1, int paIdx2, String linkType, String svType){
//        linkedItemIdx = idx;
//        linkedMateItemIdx = mateIdx;
        this.svType = svType;
        this.linkType = linkType;
        if (leftPos < rightPos) {
            linkedLeftPos = leftPos;
            linkedRightPos = rightPos;
        }else{
            linkedLeftPos = rightPos;
            linkedRightPos = leftPos;
        }

        linkedSups = sup;
        patternIdx = paIdx1;
        matePatternIdx = paIdx2;
    }

    public void setBreakInfo(double[] ratio, int[] weight) {
        bpRatio = ratio;
        bpWeight = weight;
    }

    public void setSplitLinkStatus(int status){
        splitStatus = status;
    }

    public int getSups() {
        return linkedSups;
    }


    public boolean identicalLink(Link other) {
        int aSize = linkedRightPos - linkedLeftPos;
        int bSize = other.linkedRightPos - other.linkedLeftPos;
        int minDist = Math.min((linkedLeftPos - other.linkedLeftPos), (linkedRightPos-other.linkedRightPos));
        float sizeSim = (float) Math.min(aSize, bSize) / Math.max(aSize, bSize);

        return minDist <= 200 && sizeSim >= 0.7;
    }
    public int[] getLinkedItemPos() {
        if (linkedLeftPos < linkedRightPos) {
            return new int[] {linkedLeftPos, linkedRightPos};
        }
        return new int[] {linkedRightPos, linkedLeftPos};
    }
    public boolean isSelfLink() {
        return patternIdx == matePatternIdx;
    }
    public int getPatternIdx () {
        return patternIdx;
    }
    public int getMatePatternIdx() {
        return matePatternIdx;
    }

    public String getLinkType() { return linkType;}
    public double[] getBpRatio() {return bpRatio;}
    public int[] getBpWeight() {return bpWeight;}
    public String getBpType() { return svType; }

}
