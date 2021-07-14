/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mako;

import com.beust.jcommander.*;
import com.beust.jcommander.internal.Lists;
import fspm.CFSPM;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import structures.SequenceDatabase;
import utils.FileReader;
import collector.Filter;

import collector.Collect;


/**
 *
 * @author jiadonglin
 */
public class Mako {

      /**
     * @param args the command line arguments
     */

    static String VERSION="V1.2";
    public static void main(String[] args) throws IOException{
    // TODO code application logic here
        Options options = new Options();
        JCommander jc = new JCommander(options);

        if (args.length == 0) {
            printUsage(jc);
            System.exit(1);
        }

        try {
            jc.parse(args);
            options.processArguments();
            System.out.println("\n============ Running parameters ================");
            options.printParams();

            System.out.println("\n\n============ Step1: Collection ================ \n");
            Filter filter = new Filter(options.insertMu, options.insertStd, options.cutStd, options.readLen,
                    options.maxClusterRp, options.minMapQ);

            Set<String> PassedRps = filter.runFilter(options.bamFile, options.faiFile, options.chrom, options.regionStart, options.regionEnd, options.signalSummary);

            Collect collector = new Collect(options.insertMu, options.insertStd, options.cutStd, options.readLen,
                    options.maxClusterRp, options.minMapQ, options.faiFile);

            collector.runCollect(options.bamFile, PassedRps, options.chrom, options.regionStart, options.regionEnd, options.sigOut,
                    options.alignOut, options.bndOut);

            System.out.println("\n\n================ Step2: Detection ================ \n");

            SequenceDatabase sequenceDatabase = new SequenceDatabase();

            BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(options.svOut));

            sequenceDatabase.loadSequencesFromFile(options.sigOut, options.minAf, options.minWeight);

            CFSPM algoContiguousFSPM = new CFSPM(options.searchRange, collector.chromNameMap, null);

            algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, null, svRegionWriter, options.fastaFile);
            algoContiguousFSPM.printAlgoStatistics();

        }

        catch (ParameterException e) {
            System.out.println("Failed to parse parameteres, detailed message below: ");
            System.out.println(e.getMessage());
            System.out.println();
            System.out.println("See usage: -h");
            System.exit(1);
        }


        // Debug code

//        String[] chroms = new String[]{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
//        "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"};
//        String fastaIndexFile = "/Users/apple/Data/genome/GRCh38_autosomes.fa.fai";
//        String fastaFile = "/Users/apple/Data/genome/GRCh38_autosomes.fa";
//        String bamFile = "/Users/apple/IdeaProjects/Mako/demo/NA19240.30X.chr20.1000k-2000k.bam";
//        String superitemOut = "";
//        String bamSummaryFile = "/Users/apple/IdeaProjects/Mako/demo/NA19240.mako.cfg";
//        String svOut = "/Users/apple/IdeaProjects/Mako/demo/demo_out.txt";
//
//        FileReader fileReader = new FileReader();
//        Map<String, Float> bamSummary = fileReader.loadBamSummary(bamSummaryFile);
//        Filter filter = new Filter(559, 151, 3, 126, 120, 20);
//        Set<String> failedRps = filter.runFilter(bamFile, fastaIndexFile, "chr20", 1000000, 2000000, bamSummary);
//
//        Collect collector = new Collect(514, 151, 3, 101, 120, 20, fastaIndexFile);
//        collector.runCollect(bamFile, failedRps, "chr20", 1000000 , 2000000, superitemOut, null, null);
//
//        structures.SequenceDatabase sequenceDatabase = new structures.SequenceDatabase();
//
//        BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(svOut));
//
//        sequenceDatabase.loadSequencesFromFile(superitemOut, 0.2, 10);
//        CFSPM algoContiguousFSPM = new CFSPM(2000, collector.chromNameMap, null);
//
//        algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, null, svRegionWriter, fastaFile);
//        algoContiguousFSPM.printAlgoStatistics();



    }


    private static void printUsage(JCommander commander) {

        if (commander.getDescriptions() == null) {
            commander.createDescriptions();
        }

        System.out.println("============  Mako info  ==================");
        System.out.println("\nVersion: " + VERSION);
        System.out.println("Author: Jiadong Lin (Xi'an JiaoTong University)\nContact: jiadong66@stu.xjtu.edu.cn");
        System.out.println("Usage: java -jar Mako.jar <Options>");
        System.out.println("Example: java -jar Mako.jar -R /path/to/ref.fa -F /path/to/sample.mako.cfg");
        System.out.println("================================================\n");
        StringBuilder required = new StringBuilder();
        required.append("Required\n");
        StringBuilder optional = new StringBuilder();
        optional.append("Optional\n");

        List<ParameterDescription> params = Lists.newArrayList();
        params.addAll(commander.getFields().values());
        if (params.size() > 0) {

            for (ParameterDescription pd : params) {
                if (pd.getParameter().required()) {
                    required.append("\t");
                    required.append(pd.getNames());
                    required.append("\t");
                    required.append(pd.getDescription());
                    if (pd.getDefault() != null) {
                        required.append(" (");
                        required.append("Default=");
                        required.append(pd.getDefault());
                        required.append(")");
                    }
                    required.append("\n");
                }else{
                    optional.append("\t");
                    optional.append(pd.getNames());
                    optional.append("\t");
                    optional.append(pd.getDescription());
                    if (pd.getDefault() != null) {
                        optional.append("(");
                        optional.append("Default=");
                        optional.append(pd.getDefault());
                        optional.append(")");
                    }
                    optional.append("\n");
                }

            }
        }


        System.out.println(required.toString() + "\n" + optional.toString());
    }

    @Parameters(commandDescription="Detect Structural variants from short-read data")

    static class Options {

        @Parameter(names = {"--reference", "-R"}, description = "Path to the reference genome", required = true)
        String fastaFile;

        @Parameter(names = {"--config", "-F"}, description = "Path to the Mako configuration file", required = true)
        String configFile;

        @Parameter(names = {"--minAf", "-A"}, description = "Minimum allele fraction for nodes")
        double minAf = 0.2;

        @Parameter(names = {"--minWeight", "-W"}, description = "Minimum weigth for nodes produced from abnormal reads")
        int minWeight = 10;

        @Parameter(names = {"--insertCut", "-C"}, description = "Cutoff of determining abnormal insert size of read pairs")
        int cutStd = 3;

        @Parameter(names = {"--maxClusterRp", "-D"}, description = "Maximum distance allowed for clustering discordant read pairs")
        int maxClusterRp;

        @Parameter(names = {"--minMapQ", "-Q"}, description = "Minimum mapping quality of an alignment")
        int minMapQ = 20;

        @Parameter(names = {"--searchRange", "-S"}, description = "Maximum distance allowed to search along an edge")
        int searchRange;

        @Parameter(names = {"--genomeRegion", "-G"}, description = "Discover SVs at a given region or a chromosome (e.g chr1:1000-3000 or chr1)")
        String genomeRegion;

        String bamFile;
        String faiFile;

        int insertMu;
        int insertStd;
        int readLen;

        String workDir;
        String sampleName;

        String chrom;
        int regionStart = -1;
        int regionEnd = -1;
        Map<String, Float> signalSummary;

        // output files
        String sigOut;
        String svOut;
        String bndOut;

        // output file for debug usage
        String alignOut;


        private void printParams() {

            StringBuilder sb = new StringBuilder();
            sb.append("Files");
            sb.append("\n\tWorking directory: ");
            sb.append(workDir);
            sb.append("\n\tBAM path: ");
            sb.append(bamFile);
            sb.append("\n\tReference: ");
            sb.append(fastaFile);
            sb.append("\n\tReference index: ");
            sb.append(faiFile);
            sb.append("\n");

            sb.append("Parameters\n");
            sb.append("\tInsertMu: ");
            sb.append(insertMu);
            sb.append("\n\tInsertStd: ");
            sb.append(insertStd);
            sb.append("\n\tNode search range: ");
            sb.append(searchRange);
            sb.append("\n\tRead-pair clustering distance: ");
            sb.append(maxClusterRp);
            sb.append("\n\tMinimum node weight: ");
            sb.append(minWeight);
            sb.append("\n\tMinimum node allele fraction: ");
            sb.append(minAf);
            sb.append("\n\tMinimum mapping quality: ");
            sb.append(minMapQ);
            if (chrom != null) {
                sb.append("\n\tProcessing region: ");
                sb.append(genomeRegion);
            }

            System.out.println(sb.toString());
        }

        private void processArguments() throws IOException {

            faiFile = fastaFile + ".fai";
            signalSummary = new HashMap<>();
            if (genomeRegion != null) {
                if (!genomeRegion.contains(":")) {
                    chrom = genomeRegion;
                } else {
                    chrom = genomeRegion.split(":")[0];
                    regionStart = Integer.parseInt(genomeRegion.split(":")[1].split("-")[0]);
                    regionEnd = Integer.parseInt(genomeRegion.split(":")[1].split("-")[1]);
                }
            }


            if (configFile != null) {
                FileInputStream fin = new FileInputStream(new File(configFile));
                BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
                String thisLine;

                while ((thisLine = myInput.readLine()) != null) {
                    String[] tokens = thisLine.split(":");

                    if (tokens[0].equals("bam")) {
                        bamFile = tokens[1];
                    } else if (tokens[0].equals("mean")) {
                        insertMu = Integer.parseInt(tokens[1]);
                        searchRange = 4 * insertMu;
                    } else if (tokens[0].equals("stdev")) {
                        insertStd = Integer.parseInt(tokens[1]);
                    } else if (tokens[0].equals("readlen")) {
                        readLen = Integer.parseInt(tokens[1]);
                        maxClusterRp = readLen;
                    } else if (tokens[0].equals("workDir")) {
                        workDir = tokens[1];
                    } else if (tokens[0].equals("name")) {
                        sampleName = tokens[1];
                    } else if (tokens[0].equals("ARP_LARGE")) {
                        signalSummary.put(tokens[0], Float.parseFloat(tokens[1]));
                    } else if (tokens[0].equals("ARP_SMALL")) {
                        signalSummary.put(tokens[0], Float.parseFloat(tokens[1]));
                    } else if (tokens[0].equals("ARP_RF")) {
                        signalSummary.put(tokens[0], Float.parseFloat(tokens[1]));
                    } else if (tokens[0].equals("ARP_FF")) {
                        signalSummary.put(tokens[0], Float.parseFloat(tokens[1]));
                    } else if (tokens[0].equals("ARP_RR")) {
                        signalSummary.put(tokens[0], Float.parseFloat(tokens[1]));
                    }

                    sigOut = workDir + sampleName + "_signatures.txt";
                    bndOut = workDir + sampleName + "_BNDs.txt";
                    svOut = workDir + sampleName + "_mako_calls.txt";

                    if (chrom != null) {
                        sigOut = workDir + sampleName + "." + chrom + "_signatures.txt";
                        bndOut = workDir + sampleName + "." + chrom + "_BNDs.txt";
                        svOut = workDir + sampleName + "." + chrom + "_mako_calls.txt";
                    }
                }
            }
        }
    }
}
