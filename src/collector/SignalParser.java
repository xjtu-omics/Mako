/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package collector;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

import structures.Node;
import utils.BreakEnd;

/**
 *
 * @author jiadonglin
 */
public class SignalParser {
    private List<Node> superitems = new ArrayList<>();
    // number of normal reads is not set for some break superitems in previous window, will be set in the following window.
    private List<Node> unSetSuperitems = new ArrayList<>();
    // save abnormal mapped reads - mutation signals 
    private List<MutSignal> mutSignals = new ArrayList<>();
    // Read-pair signals are separated into forward and reverse sub-channels.
    private List<MutSignal> forwardMutSignals = new ArrayList<>();
    private List<MutSignal> reverseMutSignals = new ArrayList<>();
    // some break reads may indicate small indels are saved to a separate channel while dealing with break reads
    private List<MutSignal> smallIndelMutSignals = new ArrayList<>();

    // Trans signals
    Map<String, List<MutSignal>> traSignals = new HashMap<>();
    

    private final int maximumDeletionSize = 1000000;
    // number of superitems generated in this channel
    private int SuperItemCount = 0;


    // Abnormal signal writer
    private BufferedWriter writer = null;
    private boolean isARPchannel = false;
    private final int BufferSize = 2000;
    private String channelName;
    private int rpClusterMaxDist;
    
    public SignalParser(int maxDist, String name, String signalOutPath){
        rpClusterMaxDist = maxDist;
        this.channelName = name;   
        if (signalOutPath != null){
            String fileName = signalOutPath + channelName + ".abnormal.signals.txt";
            try {
                writer = new BufferedWriter(new FileWriter(fileName));
                writer.write("qname\tmut_type\tmut_pos\tse_pos\tref\tinsert_size\tcigar\n");
            } catch (IOException e) {
                System.out.println(e);
            }            
        }
        
    }
    /**
     * Only add mutational signals at current chrom
     * @param signal
     * @param curRefName 
     */
     public void addSignals(MutSignal signal, String curRefName){
        String mutSignalType = signal.getMutSignalType();    
        String signalRefName = signal.getSignalRef();
        if (writer != null){  
            try {
                writer.write(signal.toString());
                writer.newLine(); 
            } catch (IOException e) {
                e.printStackTrace();
            }
                                  
        }
        // Read-pairs in different chromosomes
        if(signal.isInterChrom) {
            String oriKey = signal.getReadOri();
            List<MutSignal> signals = traSignals.get(oriKey);
            if (signals == null) {
                signals = new ArrayList<>();
                signals.add(signal);
                traSignals.put(oriKey, signals);
            }else {
                signals.add(signal);
                traSignals.put(oriKey, signals);
            }
        }
        else{
            // Abnormal read-pairs
            if (signal.isARPSignal()){
                isARPchannel = true;

                if (signal.insertSize <= maximumDeletionSize){
                    int forwardDistToLast = 0;
                    int reverseDistToLast = 0;

                    if (forwardMutSignals.isEmpty() && signal.getMutSignalOri().equals("+") && signalRefName.equals(curRefName)){
                        forwardMutSignals.add(signal);
                    }
                    if (reverseMutSignals.isEmpty() && signal.getMutSignalOri().equals("-") && signalRefName.equals(curRefName)){
                        reverseMutSignals.add(signal);
                    }
                    else if (signalRefName.equals(curRefName)){
                        if (signal.getMutSignalOri().equals("+")){
                            MutSignal forwardLastSignal = forwardMutSignals.get(forwardMutSignals.size() - 1);

                            forwardMutSignals.add(signal);
                            forwardDistToLast = signal.getMutPos() - forwardLastSignal.getMutPos();
                        }else{
                            MutSignal reverseLastSignal = reverseMutSignals.get(reverseMutSignals.size() - 1);

                            reverseMutSignals.add(signal);
                            reverseDistToLast = signal.getMutPos() - reverseLastSignal.getMutPos();
                        }

                    }
                    if ((forwardMutSignals.size() > BufferSize && forwardDistToLast > rpClusterMaxDist) ){
                        signalLinearClustering(forwardMutSignals, rpClusterMaxDist);

//                    forwardMutSignals.subList(0, forwardMutSignals.size() - 1).clear();
                        MutSignal lastSignal = forwardMutSignals.get(forwardMutSignals.size() - 1);
                        mutSignalListClear(forwardMutSignals, lastSignal);
                    }

                    if ((reverseMutSignals.size() > BufferSize && reverseDistToLast > rpClusterMaxDist) ){
                        signalLinearClustering(reverseMutSignals, rpClusterMaxDist);

//                    reverseMutSignals.subList(0, reverseMutSignals.size() - 1).clear();
                        MutSignal lastSignal = reverseMutSignals.get(reverseMutSignals.size() - 1);
                        mutSignalListClear(reverseMutSignals, lastSignal);
                    }
                }
            }

            // Normal read pairs, signals from one read
            else {

                if ((mutSignalType.contains("I") || mutSignalType.contains("D")) && signal.isIsizeNormal()){
                    MutSignal lastSignal = signal;
                    if (!smallIndelMutSignals.isEmpty()){
                        lastSignal = smallIndelMutSignals.get(smallIndelMutSignals.size() - 1);
                    }
                    smallIndelMutSignals.add(signal);
                    int distToLastSignal = signal.getMutPos() - lastSignal.getMutPos();

                    if (smallIndelMutSignals.size() > BufferSize && distToLastSignal > 0){
                        signalLinearClustering(smallIndelMutSignals, 0);

                        // clear the list only keep the last element.
//                    smallIndelMutSignals.subList(0, smallIndelMutSignals.size() - 1).clear();
                        MutSignal restSignal = smallIndelMutSignals.get(smallIndelMutSignals.size() - 1);
                        mutSignalListClear(smallIndelMutSignals, restSignal);

                    }
                }else if (!mutSignalType.contains("I") && !mutSignalType.contains("D")){

                    MutSignal lastSignal = signal;
                    if (! mutSignals.isEmpty()){
                        lastSignal = mutSignals.get(mutSignals.size() - 1);
                    }

                    int distToLastSignal = signal.getMutPos() - lastSignal.getMutPos();
                    mutSignals.add(signal);
                    if (mutSignals.size() > BufferSize && distToLastSignal > 0 ){
                        signalLinearClustering(mutSignals, 0);

//                    mutSignals.subList(0, mutSignals.size() - 1).clear();
                        MutSignal restSignal = mutSignals.get(mutSignals.size() - 1);
                        mutSignalListClear(mutSignals, restSignal);
                    }
                }
            }
        }

        
    }

    private void mutSignalListClear(List<MutSignal> listToClear, MutSignal signalToAdd){
        listToClear.clear();
        listToClear.add(signalToAdd);
    }
    
    public int getSuperitemCount(){
        return SuperItemCount;
    }
    
    private void addSuperItem(Node superitem){
        SuperItemCount += 1;                                    
        superitems.add(superitem);                        
    }

    /**
     * Create BND candidates within each sliding windows
     */

    public void createBreakEndCandidates(BufferedWriter bndWriter) throws IOException{
        // We first form partition of signals with same orientation
        for(Entry<String, List<MutSignal>> entry : traSignals.entrySet()) {
            List<MutSignal> signals = entry.getValue();

            // First sort signals by contig and pos
            Collections.sort(signals, new Comparator<MutSignal>() {
                @Override
                public int compare(MutSignal o1, MutSignal o2) {
                    int t;
                    t = o1.contig1.compareTo(o2.contig1);
                    // equal contig
                    if (t == 0) {
                        t = o1.pos1 - o2.pos1;
                    }
                    return t;
                }
            });

            List<List<MutSignal>> thisOriSourcePartitions = breakEndSourcePartition(signals);

            int sourcePartitionSize = thisOriSourcePartitions.size();

            for (int i = 0; i < sourcePartitionSize; i++) {
                List<MutSignal> sourcePartition = thisOriSourcePartitions.get(i);
                String contig1 = sourcePartition.get(0).contig1;
                int pos1 = sourcePartition.get(0).pos1;

                List<List<MutSignal>> thisOriDestPartitions = breakEndDestPartition(sourcePartition);

                double sourceStd = std(sourcePartition);

                for (List<MutSignal> destPartition : thisOriDestPartitions) {

                    String contig2 = destPartition.get(0).contig2;
                    int pos2 = destPartition.get(0).pos2;

                    double destStd = std(destPartition);

                    if (destStd <= 50 && sourceStd <= 50) {
                        BreakEnd bnd = new BreakEnd(contig1, contig2, pos1, pos2, sourceStd, destStd);
                        bndWriter.write(bnd.createSourceEntry());
                        bndWriter.newLine();

                    }

//                    bndWriter.write(bnd.createDestEntry());
                }

            }
        }
        traSignals = new HashMap<>();

    }

    private List<List<MutSignal>> breakEndSourcePartition(List<MutSignal> signals) {
        List<List<MutSignal>> partitions = new ArrayList<>();
        List<MutSignal> currentPartition = new ArrayList<>();

        for (MutSignal signal : signals) {
            int currentPartSize = currentPartition.size();
            if (currentPartSize > 0 && sourceDistance(currentPartition.get(currentPartSize - 1), signal) > 200){
                if (currentPartSize >= 10) {
                    partitions.add(currentPartition);
                    currentPartition = new ArrayList<>();
                }

            }
            currentPartition.add(signal);
        }

        if (currentPartition.size() >= 10) {
            partitions.add(currentPartition);
        }

        return partitions;
    }

    private double std(List<MutSignal> signals) {
        int sum = 0;
        double sd = 0;
        for (MutSignal signal : signals) {
            int pos = signal.getMutPos();
            sum += pos;
        }
        double mean = (double) sum / signals.size();

        for (MutSignal signal : signals) {
            sd += Math.pow(signal.getMutPos() - mean, 2);
        }

        return (int)Math.sqrt(sd / signals.size());
    }

    /**
     * For a given parition of same orientaion, we need further find their correponding destination
     * @param signals
     */
    private List<List<MutSignal>> breakEndDestPartition(List<MutSignal> signals) {
        Collections.sort(signals, new Comparator<MutSignal>() {
            @Override
            public int compare(MutSignal o1, MutSignal o2) {
                int t;
                t = o1.contig2.compareTo(o2.contig2);
                // equal contig
                if (t == 0) {
                    t = o1.pos2 - o2.pos2;
                }
                return t;
            }
        });

        List<List<MutSignal>> partitions = new ArrayList<>();
        List<MutSignal> currentPartition = new ArrayList<>();

        for (MutSignal signal : signals) {
            if (currentPartition.size() < 1) {
                currentPartition.add(signal);
                continue;
            }
            int diff = destDistance(currentPartition.get(0), signal);
            if(diff > 1000) {
                if (currentPartition.size() >= 10) {
                    List<MutSignal> thisPartition = new ArrayList<>();
                    thisPartition.addAll(currentPartition);

                    partitions.add(thisPartition);
                }
                while (currentPartition.size() > 0 && destDistance(currentPartition.get(0), signal) > 1000) {
                    currentPartition.remove(0);
                }

            }
            currentPartition.add(signal);
        }

        if (currentPartition.size() >= 10) {
            partitions.add(currentPartition);
        }

        return partitions;
    }

    private int destDistance(MutSignal sig1, MutSignal sig2) {
        return sig1.contig2.equals(sig2.contig2) ? Math.abs(sig1.pos2 - sig2.pos2) : Integer.MAX_VALUE;
    }

    private int sourceDistance(MutSignal sig1, MutSignal sig2) {
        if (sig1.mutSignalType.equals(sig2.mutSignalType) && sig1.contig1.equals(sig2.contig1)) {
            return Math.abs(sig1.pos1 - sig2.pos1);
        }else{
            return Integer.MAX_VALUE;
        }
    }
    
    /**
     * get number of normal read per base.
     * @param readDepthContainer
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer 
     */
    public void setSuperitemWeightRatio(int[] readDepthContainer, int windowStart, int windowSize, int[] preReadDepthBuffer){
        List<Node> setNode = new ArrayList<>();
        for (Node si : superitems){
            int pos = si.getPos();
            int indexOfArray = pos - windowStart;
            if (indexOfArray < windowSize && indexOfArray >= 0){
                int readDepthAtPos = readDepthContainer[indexOfArray];
                si.setSuperitemReadDepth(readDepthAtPos);                   
                setNode.add(si);

            }else if (indexOfArray < 0){
                int indexInPreBuffer = preReadDepthBuffer.length + indexOfArray;
                int readDepthAtPos = preReadDepthBuffer[indexInPreBuffer];
                si.setSuperitemReadDepth(readDepthAtPos);                   
                setNode.add(si);
            }
            else{
                unSetSuperitems.add(si);
            }                                                
        }
        superitems.clear();
        superitems = setNode;
    }
    
    public void writeSuperItemsInChannel(BufferedWriter superitemWriter) throws IOException{
        
        if (isARPchannel){          
            for (Node si : superitems){
                String strOut = si.toString();
                superitemWriter.write(strOut);
                superitemWriter.newLine();   
                
            }       
            superitems.clear();
            SuperItemCount = 0;
        }else{
           for (Node si : superitems){
                String strOut = si.toString();
                superitemWriter.write(strOut);
                superitemWriter.newLine();            
            }       
            superitems.clear();
            for (Node si : unSetSuperitems){
                superitems.add(si);
            }
            unSetSuperitems.clear();
            SuperItemCount = 0;
        }
        
    }        
    
    public void processFinalSignals(){

        if (isARPchannel){
            if (!forwardMutSignals.isEmpty()){
                signalLinearClustering(forwardMutSignals, rpClusterMaxDist);

                forwardMutSignals.clear();
            }
            if (!reverseMutSignals.isEmpty()){
                signalLinearClustering(reverseMutSignals, rpClusterMaxDist);
 
                reverseMutSignals.clear();
            }
        }else{
            if (!mutSignals.isEmpty()){
                signalLinearClustering(mutSignals, 0);
               
                mutSignals.clear();
            }
            if (!smallIndelMutSignals.isEmpty()){
                signalLinearClustering(smallIndelMutSignals, 0);
                
                smallIndelMutSignals.clear();
            }
        }
        
    }

    /**
     * Clustering existing signals to node in the buffer
     * @param signals
     * @param maxDist
     */
    private void signalLinearClustering(List<MutSignal> signals, int maxDist){
                           
        if (signals.size() < 3){
            return;
        }
        // sort mutation signals in ascending order
        Collections.sort(signals);
        Iterator<MutSignal> iter = signals.iterator();
//        List<List<MutSignal>> clusters = new ArrayList<>();
        List<MutSignal> cluster = new ArrayList<>();
        cluster.add(iter.next());
//        clusters.add(cluster);
        
        while (iter.hasNext()){
            MutSignal mutSignal = iter.next();
            if (mutSignal.withinDistance(cluster.get(cluster.size() - 1), maxDist)){
                cluster.add(mutSignal);
            }else{
                if (cluster.size() < 3){
                    cluster.clear();
                }else{
                    Map<String, List<MutSignal>> combinedParts = signalCombineInCluster(cluster);
                    for(Entry<String, List<MutSignal>> entry : combinedParts.entrySet()){
                        List<MutSignal> thisCluster = entry.getValue();
                        Node node = new Node(thisCluster);
                        addSuperItem(node);
                        cluster.clear();
                    }

                }

                cluster.add(mutSignal);
            }            
            
        }
        
        if (cluster.size() >= 3) {
            Node superitem = new Node(cluster);
            addSuperItem(superitem);
        }
                
    }

    private Map<String, List<MutSignal>> signalCombineInCluster(List<MutSignal> cluster){

        Map<String, List<MutSignal>> combinedClusters = new HashMap<>();
        for (MutSignal signal : cluster) {
            String signalType = signal.getMutSignalType();
            List<MutSignal> thisPart = combinedClusters.get(signalType);
            if (thisPart == null) {
                thisPart = new ArrayList<>();
                thisPart.add(signal);
                combinedClusters.put(signalType, thisPart);
            }else{
                thisPart.add(signal);
                combinedClusters.put(signalType, thisPart);
            }
        }

        return combinedClusters;
    }

    public void setARPSuperItemRatio(int[] readDepthContainer, int windowStart, int windowSize, int[] preReadDepthBuffer){
        for (Node node : superitems){
            int superItemIntervalStart = node.getSuperitemRegion().start - windowStart - 1;
            int superItemIntervalEnd = node.getSuperitemRegion().end - windowStart;
            int nreadSum = 0;
            double superItemCov = 0;
            if (superItemIntervalEnd < windowSize && superItemIntervalStart >= 0){
                 for (int i = superItemIntervalStart; i< superItemIntervalEnd; i++){
                    nreadSum += readDepthContainer[i];
                }
                superItemCov = (double) nreadSum / (superItemIntervalEnd - superItemIntervalStart);
            }else if (superItemIntervalEnd < 0 && superItemIntervalEnd > -1000000){                
                int endIndexPreBuffer = preReadDepthBuffer.length + superItemIntervalEnd;
                int startIndexPreBuffer = preReadDepthBuffer.length + superItemIntervalStart;
                for (int i = startIndexPreBuffer; i < endIndexPreBuffer ; i ++){
                    try {
                        nreadSum += preReadDepthBuffer[i];
                    } catch (Exception e) {
                        System.err.println("error");
                    }
                    nreadSum += preReadDepthBuffer[i];
                }
                superItemCov = (double) nreadSum / (superItemIntervalEnd - superItemIntervalStart);
            }else{
                superItemCov = 0;
            }                                        
            node.setARPsuperitemRatio(superItemCov);
        }        
    }
}