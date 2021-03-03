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
            System.out.println("==== Step1: running Mako collection ==== ");

            options.printParams();

            Filter filter = new Filter(options.insertMu, options.insertStd, options.cutStd, options.readLen,
                    options.maxClusterRp, options.minMapQ);

            Set<String> failedRps = filter.runFilter(options.bamFile, options.faiFile,options.chrom, options.signalSummary);

            Collect collector = new Collect(options.insertMu, options.insertStd, options.cutStd, options.readLen,
                    options.maxClusterRp, options.minMapQ, options.faiFile);

            collector.runCollect(options.bamFile, failedRps, options.chrom, options.regionStart, options.regionEnd, options.sigOut,
                    options.alignOut, options.bndOut);

            System.out.println(" ==== Step2: running Mako detection ==== ");

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
//        String fastaIndexFile = "/Users/jiadonglin/Data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai";
//        String fastaFile = "/Users/jiadonglin/Data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";
//        String bamFile = "/Users/jiadonglin/Data/HG00514/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram.bam";
//        String superitemOut = "/Users/jiadonglin/SV_Research/Mako/Revision2/HG00514/items_exclude_failedRps.txt";
//        String bamSummaryFile = "/Users/jiadonglin/SV_Research/Mako/Revision2/HG00514/HG00514.mako.cfg";

//        String fastaIndexFile = "/Users/jiadonglin/Data/ref_genome/hs37d5.fa.fai";
//        String fastaFile = "/Users/jiadonglin/Data/ref_genome/hs37d5.fa";
//        String bamFile = "/Users/jiadonglin/Data/SKBR3/SKBR3_550bp_pcrFREE_S1_L001_AND_L002_R1_001.101bp.bwamem.ill.mapped.sort.bam";
//        String superitemOut = "/Users/jiadonglin/SV_Research/Mako/Revision2/SKBR3/items_exclude_failedRps.txt";
//        String bamSummaryFile = "/Users/jiadonglin/SV_Research/Mako/Revision2/SKBR3/SKBR3.mako.cfg";



//        String fastaIndexFile = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/grch38.reference.fa.fai";
//        String fastaFile = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/grch38.reference.fa";
//        String bamFile = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/p100_c30/random_sim/sim.srt.bam";
//        String bamSummaryFile = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/p100_c30/random_sim/wgs_sim.mako.cfg";
//        String superitemOut = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/p100_c30/random_sim/svHasType/nodes_exclude_failedRps.txt";
//        String svOut = "/Users/jiadonglin/SV_Research/onGoingProjects/Mako/Revision2/simulation/p100_c30/random_sim/svHasType/svout_test.txt";

//        for (String chrom: chroms) {
//
//            String svOut = "/Users/jiadonglin/SV_Research/Mako/Revision2/HG00514/" + chrom + "_svout.txt";
//            FileReader fileReader = new FileReader();
//            Map<String, Float> bamSummary = fileReader.loadBamSummary(bamSummaryFile);
//            Filter filter = new Filter(514, 151, 3, 101, 120, 20);
//            Set<String> failedRps = filter.runFilter(bamFile, fastaIndexFile, chrom, bamSummary);
//
//            Collect collector = new Collect(514, 151, 3, 101, 120, 20, fastaIndexFile);
//            collector.runCollect(bamFile, failedRps, chrom, 0 , 0, superitemOut, null, null);
//
//            structures.SequenceDatabase sequenceDatabase = new structures.SequenceDatabase();
//
//            BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(svOut));
//
//            sequenceDatabase.loadSequencesFromFile(superitemOut, 0.2, 10);
//            CFSPM algoContiguousFSPM = new CFSPM(2000, collector.chromNameMap, null);
//
//            algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, null, svRegionWriter, fastaFile);
//            algoContiguousFSPM.printAlgoStatistics();
//        }



    }


    private static void printUsage(JCommander commander) {

        if (commander.getDescriptions() == null) {
            commander.createDescriptions();
        }

        System.out.println("============  Mako info  ==================");
        System.out.println("\nVersion: V1.0\nAuthor: Jiadong Lin (Xi'an JiaoTong University)\nContact: jiadong324@gmail.com");
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
        int regionStart;
        int regionEnd;
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
            sb.append("\n================================\n");

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
