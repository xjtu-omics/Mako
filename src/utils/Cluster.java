package utils;

import java.util.*;

public class Cluster implements Comparable<Cluster>{

    private double[] mean;
    private double[] sum;

    private List<double[]> vectors = new ArrayList<>();
    private List<Link> links;
    private Set<String> bpTypes;

    // cluster properties
    private String clusterType;
    private int clusterStd;
    private int clusterSa;
    private int clusterRp;

    public Cluster(Link link) {
        sum = new double[2];
        links = new ArrayList<>();
        bpTypes = new HashSet<>();
        clusterType = link.getLinkType();
    }

    @Override
    public int compareTo(Cluster other) {
        return (int)mean[0] - (int)other.mean[0];
    }

    public void setMean(double[] mean) {
        this.mean = mean;
    }

    public void addVector(Link link) {
        links.add(link);
        if (!link.getBpType().equals("")){
            bpTypes.add(link.getBpType());
        }

        double[] vector = new double[]{link.getLinkedItemPos()[0], link.getLinkedItemPos()[1]};
        vectors.add(vector);
        for(int i = 0; i < vector.length; i++){
            sum[i] += vector[i];
        }
    }

    public void recomputeClusterMean() {
        for (int i = 0; i < sum.length; i++){
            mean[i] = sum[i] / (double) vectors.size();
        }
    }

    public void computeClusterStats() {
        double sd = 0;
        int totalRpSup = 0;
        int totalSaSup = 0;
        double meanDiff = mean[1] - mean[0];

        for (Link link : links) {
            double size = Math.abs(link.getLinkedItemPos()[1] - link.getLinkedItemPos()[0]);
            if (link.getLinkType().equals("sa")) {
                totalSaSup += link.getSups();
            }
            if (link.getLinkType().equals("rp")){
                totalRpSup += link.getSups();
            }
            if (link.getLinkType().equals("bt")){
                totalRpSup += link.getSups();
            }
            sd += Math.pow(size - meanDiff, 2);
        }

        clusterStd = (int)Math.sqrt(sd / links.size());
        clusterRp = totalRpSup / links.size();
        clusterSa = totalSaSup / links.size();
    }


    public List<Link> getLinks() { return links; }

    public Set<String> getBpTypes() { return bpTypes; }

    public double[] getMean(){ return mean; }

    public String getSupType() { return clusterType; }

    public int[] getBp() { return new int[]{(int)mean[0], (int)mean[1]}; }

    public int getClusterStd() { return clusterStd; }

    public String toString() {
        List<String> list = new ArrayList<>(bpTypes);
        StringBuilder sb = new StringBuilder();
        if (list.size() == 0) {
            sb.append("BND,");
        }else{
            sb.append(list.get(0));
            sb.append(",");
        }
        sb.append((int)mean[0]);
        sb.append("-");
        sb.append((int)mean[1]);
        sb.append(",sa=");
        sb.append(clusterSa);
        sb.append(",rp=");
        sb.append(clusterRp);
        return sb.toString();
    }
}
