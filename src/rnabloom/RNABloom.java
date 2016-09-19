/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.File;
import java.io.IOException;
import static java.lang.Math.pow;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqPair;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.ReadPair;
import rnabloom.io.FastqReader;
import static rnabloom.util.GraphUtils.assembleFragment;
import static rnabloom.util.GraphUtils.correctMismatches;
import static rnabloom.util.GraphUtils.findBackbonePath;
import static rnabloom.util.GraphUtils.naiveExtend;
import static rnabloom.util.SeqUtils.*;
import static rnabloom.util.GraphUtils.assembleTranscript;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    private final static long NUM_PARSED_INTERVAL = 100000;
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    
    private final int k;
    private final Pattern qualPattern;
    private BloomFilterDeBruijnGraph graph;
        
    public RNABloom(boolean strandSpecific, long dbgbfNumBits, long cbfNumBytes, long pkbfNumBits, int dbgbfNumHash, int cbfNumHash, int pkbfNumHash, int seed, int k, int q) {
        this.k = k;
        this.qualPattern = getPhred33Pattern(q, k);
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            seed,
                                            k,
                                            strandSpecific);
    }
    
    public void saveGraph(File f) {
        try {
            graph.save(f);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void restoreGraph(File f) {
        try {
            graph = new BloomFilterDeBruijnGraph(f);
            
            //BloomFilterDeBruijnGraph graph2 = new BloomFilterDeBruijnGraph(f);
            //System.out.println(graph2.getDbgbf().equivalent(graph.getDbgbf()));
            //System.out.println(graph2.getCbf().equivalent(graph.getCbf()));
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void createDBG(String[] forwardFastqs, String[] reverseFastqs) {
        int lineNum = 0;
        
        try {
            FastqReader fr;
            
            for (String fastq : forwardFastqs) {
                System.out.println("Parsing `" + fastq + "`...");
                
                fr = new FastqReader(fastq, false);
                while (fr.hasNext()) {
                    if (++lineNum % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(lineNum) + " reads...");
                    }

                    for (String seq : filterFastq(fr.next(), qualPattern)) {
                        graph.addKmersFromSeq(seq);
                    }
                }
                fr.close();
            }
            
            for (String fastq : reverseFastqs) {
                System.out.println("Parsing `" + fastq + "`...");
                
                fr = new FastqReader(fastq, false);
                while (fr.hasNext()) {
                    if (++lineNum % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(lineNum) + " reads...");
                    }

                    for (String seq : filterFastq(fr.next(), qualPattern)) {
                        graph.addKmersFromSeq(reverseComplement(seq));
                    }
                }
                fr.close();
            }
        }
        catch (NoSuchElementException e) {
            /**@TODO handle invalid format*/
        }
        catch (IOException e) {
            /**@TODO Handle it!!! */
        }
        finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(lineNum) + " reads...");
        }
    }

    public BloomFilterDeBruijnGraph getGraph() {
        return graph;
    }
    
    private boolean isPolyA(ReadPair p) {
        return p.right.endsWith("AAAA") || (!graph.isStranded() && p.left.startsWith("TTTT"));
    }
    
    private boolean okToAssemble(ReadPair p) {
        if (p.left.length() >= k && p.right.length() >= k) {
            String[] leftKmers = kmerize(p.left, k);
            String[] rightKmers = kmerize(p.right, k);
            int threshold = 2*k-1;

            int numLeftAssembled = 0;
            for (String s : leftKmers) {
                if (graph.lookupFragmentKmer(s)) {
                    ++numLeftAssembled;
                }
            }

            if (numLeftAssembled <= threshold && leftKmers.length - numLeftAssembled > 0) {
                return true;
            }

            int numRightAssembled = 0;
            for (String s : rightKmers) {
                if (graph.lookupFragmentKmer(s)) {
                    ++numRightAssembled;
                }
            }

            if (numRightAssembled <= threshold && rightKmers.length - numRightAssembled > 0) {
                return true;
            }
        }

        return false;        
    }
    
    public void assembleFragments(FastqPair[] fastqs, String outFasta, int mismatchesAllowed, int bound, int lookahead, int minOverlap, int maxTipLen, int sampleSize) {
        System.out.println("Assembling fragments ...");
        
        graph.initializePairKmersBloomFilter();
        
        long readPairsParsed = 0;
        int correctionWorks = 0;
        
        try {
            FastqReader lin, rin;
            FastaWriter out = new FastaWriter(outFasta);
            ArrayList<String> sampleFragments = new ArrayList<>(sampleSize);
            int[] fragmentLengths = new int[sampleSize];
            int fid = 0;
            
            for (FastqPair fqPair: fastqs) {
                lin = new FastqReader(fqPair.leftFastq, true);
                rin = new FastqReader(fqPair.rightFastq, true);

                FastqPairReader fqpr = new FastqPairReader(lin, rin, qualPattern, fqPair.leftRevComp, fqPair.rightRevComp);
                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");
                
                ReadPair p;
                while (fqpr.hasNext()) {
                    p = fqpr.next();
                    
                    if (++readPairsParsed % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
                    }                    
                    
                    if (okToAssemble(p)) {
                        // correct individual reads
                        p.left = correctMismatches(p.left, graph, lookahead, mismatchesAllowed);
                        p.right = correctMismatches(p.right, graph, lookahead, mismatchesAllowed);
                        
                        if (okToAssemble(p)) {
                            String fragment = assembleFragment(p.left, p.right, graph, bound, lookahead, minOverlap);
                            
                            // correct fragment
                            fragment = correctMismatches(fragment, graph, lookahead, mismatchesAllowed);

                            int fragLen = fragment.length();

                            if (fragLen > k) {                            
                                /** extend on both ends unambiguously*/
                                fragment = naiveExtend(fragment, graph, maxTipLen);

                                ++fid;
                                //out.write(Integer.toString(fid), fragment);
                                out.write(Integer.toString(fid) + " " + p.left + " " + p.right, fragment);

                                if (fid > sampleSize) {
                                    /** store paired kmers */
                                    graph.addPairedKmersFromSeq(fragment);
                                }
                                else if (fid == sampleSize) {
                                    fragmentLengths[0] = fragLen;
                                    sampleFragments.add(fragment);

                                    /** Calculate median fragment length */
                                    Arrays.sort(fragmentLengths);
                                    int half = sampleSize/2;
                                    int medianFragLen = (fragmentLengths[half] + fragmentLengths[half - 1])/2;

                                    System.out.println("Median fragment length: " + fragLen);

                                    /** set kmer pair distance */
                                    graph.setPairedKmerDistance(medianFragLen - k);

                                    /** readjust bound to be based on 1.5*IQR */
                                    int whisker = (fragmentLengths[sampleSize*3/4] - fragmentLengths[sampleSize/4]) * 3/2;
                                    bound = medianFragLen + whisker;

                                    System.out.println("longest fragment allowed: " + bound);

                                    /** clear sample fragment lengths */
                                    fragmentLengths = null;

                                    /** store paired kmers of all sample fragments */
                                    for (String frag : sampleFragments) {
                                        graph.addPairedKmersFromSeq(frag);
                                    }

                                    /** clear sample fragments */
                                    sampleFragments = null;
                                }
                                else {
                                    /** store fragment length*/
                                    fragmentLengths[fid] = fragLen;
                                    sampleFragments.add(fragment);

                                    //System.out.println(fragLen);
                                }
                            }

                            /* assemble 3' UTR only
                            if (p.right.endsWith("AAA")) {
                                String fragment = assembleFragment(p.left, p.right, graph, mismatchesAllowed, bound, lookahead, minOverlap);
                                System.out.println("LEFT:  " + p.left);
                                System.out.println("RIGHT: " + p.right);
                                System.out.println(fragment.length() + ": " + fragment);
                            }
                            */
                        }
                    }
                    else {
                        ++correctionWorks;
                    }
                }

                lin.close();
                rin.close();
            }
            
            out.close();
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
            System.out.println(correctionWorks + " reads not needed for assembly due to correction");
        }
    }
    
    public void savePairedKmersBloomFilter(File graphFile) {
        try {
            graph.savePkbf(graphFile);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void restorePairedKmersBloomFilter(File graphFile) {
        try {
            graph.restorePkbf(graphFile);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void assembleTranscripts(String inFasta, String outFasta, int lookAhead) {
        long numFragmentsParsed = 0;
        
        try {
            System.out.println("Parsing `" + inFasta + "`...");
            
            FastaReader fin = new FastaReader(inFasta);
            FastaWriter fout = new FastaWriter(outFasta);
            
            int cid = 0;
            while (fin.hasNext()) {
                if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
                    System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
                }
                
                String fragment = fin.next();
                String transcript = assembleTranscript(fragment, graph, lookAhead);
                
                fout.write(Integer.toString(++cid), transcript);
            }
            
            fin.close();
            fout.close();
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
        }
    }
    
    public void test2() {
        String left = "ATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCG";
        String right = "CCTCCACCTTCGACGCTGGGGCTGGCATTGCC";
        
        String fragment = assembleFragment(left, right, graph, 500, 5, 10);
        
        
    }
    
    public void test() {
        //String f = "TCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCC";
        String f = "CCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCT";
        String transcript = assembleTranscript(f, graph, 5);
        System.out.print(transcript);
        for (Kmer kmer : graph.getKmers(transcript)) {
            System.out.print(kmer.count);
            System.out.print(" ");
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
                
        ///*
        String fastq1 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz"; //right
        String fastq2 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz"; //left        
        String fragsFasta = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/fragments.fa";
        String transcriptsFasta = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/transcripts.fa";
        String graphFile = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/graph";
        //*/
        /*
        String fastq1 = "/home/gengar/test_data/GAPDH/GAPDH_1.fq.gz";
        String fastq2 = "/home/gengar/test_data/GAPDH/GAPDH_2.fq.gz";
        String fragsFasta = "/home/gengar/test_assemblies/GAPDH/fragments.fa";
        String transcriptsFasta = "/home/gengar/test_assemblies/GAPDH/transcripts.fa";
        String graphFile = "/home/gengar/test_assemblies/GAPDH/graph";
        */        
        
        boolean revComp1 = true;
        boolean revComp2 = false;
        
        int mismatchesAllowed = 5;
        int bound = 500;
        int lookahead = 5;
        int minOverlap = 10;
        int maxTipLen = 10;
        
        String[] forwardFastqs = new String[]{fastq2};
        String[] backwardFastqs = new String[]{fastq1};
        
        boolean strandSpecific = true;
        long dbgbfSize = NUM_BITS_1GB;
        long cbfSize = NUM_BYTES_1GB;
        long pkbfSize = NUM_BITS_1GB;
        int dbgbfNumHash = 3;
        int cbfNumHash = 4;
        int pkbfNumHash = 1;
        int seed = 689;
        int k = 25;
        int q = 3;
        int sampleSize = 1000;
        
        ///*
        RNABloom assembler = new RNABloom(strandSpecific, dbgbfSize, cbfSize, pkbfSize, dbgbfNumHash, cbfNumHash, pkbfNumHash, seed, k, q);
        
        /*
        assembler.createDBG(forwardFastqs, backwardFastqs);
        
        
        System.out.println("Saving graph to file `" + graphFile + "`...");
        assembler.saveGraph(new File(graphFile));
        */
        
        System.out.println("Loading graph from file `" + graphFile + "`...");
        assembler.restoreGraph(new File(graphFile));
        
        ///*
        FastqPair fqPair = new FastqPair(fastq2, fastq1, revComp2, revComp1);
        FastqPair[] fqPairs = new FastqPair[]{fqPair};
        
        ///*
        assembler.assembleFragments(fqPairs, fragsFasta, mismatchesAllowed, bound, lookahead, minOverlap, maxTipLen, sampleSize);
        
        System.out.println("Saving paired kmers Bloom filter to file...");
        assembler.savePairedKmersBloomFilter(new File(graphFile));
        
        System.out.println("Restoring paired kmers Bloom filter from file...");
        assembler.restorePairedKmersBloomFilter(new File(graphFile));
        
        assembler.assembleTranscripts(fragsFasta, transcriptsFasta, lookahead);
        //*/
    }
    
}
