package utils;

import collector.MutSignal;

import java.util.List;

public class BreakEnd {

    private String sourceContig;
    private String destContig;
    private int sourePos;
    private int destPos;


    private double sourceStd;
    private double destStd;

    public BreakEnd(String contig1, String contig2, int pos1, int pos2,
                    double std1, double std2) {
        sourceContig = contig1;
        destContig = contig2;
        sourePos = pos1;
        destPos = pos2;

        sourceStd = std1;
        destStd = std2;

    }



    public String createSourceEntry() {
        StringBuilder sb = new StringBuilder();
        sb.append(sourceContig);
        sb.append("\t");
        sb.append(sourePos);
        sb.append("\t");
        sb.append(sourePos + 1);
        sb.append("\t");
        sb.append(sourceStd);
        sb.append("\t");
        String destStr = "bnd;>" + destContig + ":" + destPos + ":" + destStd;
        sb.append(destStr);

        return sb.toString();
    }

    public String createDestEntry(){
        StringBuilder sb = new StringBuilder();
        sb.append(destContig);
        sb.append("\t");
        sb.append(destPos);
        sb.append("\t");
        sb.append(destPos + 1);
        sb.append("\t");

        sb.append("\t");
        String sourceStr = "bnd;<" + sourceContig + ":" + sourePos;
        sb.append(sourceStr);

        return sb.toString();
    }
}
