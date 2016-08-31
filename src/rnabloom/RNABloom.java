/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import static java.lang.Math.pow;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import rnabloom.io.FastqPair;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.ReadPair;
import rnabloom.io.FastqReader;
import static rnabloom.util.GraphUtils.assemble;
import static rnabloom.util.GraphUtils.assembleFragment;
import static rnabloom.util.GraphUtils.findBackbonePath;
import static rnabloom.util.GraphUtils.naiveExtend;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    private final static String GZIP_EXTENSION = ".gz";
    
    private final int k;
    private final Pattern qualPattern;
    private final BloomFilterDeBruijnGraph graph;
    
    public RNABloom(boolean strandSpecific, long dbgbfNumBits, long cbfNumBytes, long pkbfNumBits, int dbgbfNumHash, int cbfNumHash, int pkbfSingleKeyNumHash, int pkbfPairedKeysNumHash, int seed, int k, int q) {
        this.k = k;
        this.qualPattern = getPhred33Pattern(q, k);
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfSingleKeyNumHash,
                                            pkbfPairedKeysNumHash,
                                            seed,
                                            k,
                                            strandSpecific);
    }
    
    public void createDBG(String[] forwardFastqs, String[] reverseFastqs) {
        int lineNum = 0;
        
        try {
            BufferedReader br;
            FastqReader fr;
            
            for (String fastq : forwardFastqs) {
                System.out.println("Parsing `" + fastq + "`...");
                
                if (fastq.endsWith(GZIP_EXTENSION)) {
                    br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
                }
                else {
                    br = new BufferedReader(new InputStreamReader(new FileInputStream(fastq)));
                }
                fr = new FastqReader(br, false);
                while (fr.hasNext()) {
                    if (++lineNum % 100000 == 0) {
                        System.out.println("Parsed " + lineNum + " reads...");
                    }

                    for (String seq : filterFastq(fr.next(), qualPattern)) {
                        graph.addKmersFromSeq(seq);
                    }
                }
                br.close();
            }
            
            for (String fastq : reverseFastqs) {
                System.out.println("Parsing `" + fastq + "`...");
                
                if (fastq.endsWith(GZIP_EXTENSION)) {
                    br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
                }
                else {
                    br = new BufferedReader(new InputStreamReader(new FileInputStream(fastq)));
                }
                fr = new FastqReader(br, false);
                while (fr.hasNext()) {
                    if (++lineNum % 100000 == 0) {
                        System.out.println("Parsed " + lineNum + " reads...");
                    }

                    for (String seq : filterFastq(fr.next(), qualPattern)) {
                        graph.addKmersFromSeq(reverseComplement(seq));
                    }
                }
                br.close();
            }
        }
        catch (NoSuchElementException e) {
            /**@TODO handle invalid format*/
        }
        catch (IOException e) {
            /**@TODO Handle it!!! */
        }
        finally {
            System.out.println("Parsed " + lineNum + " reads...");
        }
    }

    public BloomFilterDeBruijnGraph getGraph() {
        return graph;
    }
        
    public void assembleFragments(FastqPair[] fastqs, String outFasta, int mismatchesAllowed, int bound, int lookahead, int minOverlap, int maxTipLen, int sampleSize) {
        try {
            BufferedReader lin, rin;
            FastqReader leftReader, rightReader;
            BufferedWriter out = new BufferedWriter(new FileWriter(outFasta));
            ArrayList<String> sampleFragments = new ArrayList<>(sampleSize);
            int[] fragmentLengths = new int[sampleSize];
            int fid = 0;
            
            for (FastqPair fqPair: fastqs) {
                if (fqPair.leftFastq.endsWith(GZIP_EXTENSION)) {
                    lin = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fqPair.leftFastq))));
                }
                else {
                    lin = new BufferedReader(new InputStreamReader(new FileInputStream(fqPair.leftFastq)));
                }
                leftReader = new FastqReader(lin, true);

                if (fqPair.rightFastq.endsWith(GZIP_EXTENSION)) {
                    rin = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fqPair.rightFastq))));
                }
                else {
                    rin = new BufferedReader(new InputStreamReader(new FileInputStream(fqPair.rightFastq)));
                }
                rightReader = new FastqReader(rin, true);

                FastqPairReader fqpr = new FastqPairReader(leftReader, rightReader, qualPattern, fqPair.leftRevComp, fqPair.rightRevComp);
                ReadPair p;
                long readPairsParsed = 0;
                
                while (fqpr.hasNext()) {
                    p = fqpr.next();
                                        
                    ++readPairsParsed;
                    if (++readPairsParsed % 100000 == 0) {
                        System.out.println("Parsed " + readPairsParsed + " read pairs...");
                    }                    
                    
                    if (p.left.length() >= k && p.right.length() >= k) {
                        /**@TODO check whether kmers in reads had been assembled */
                        
                        String fragment = assembleFragment(p.left, p.right, graph, mismatchesAllowed, bound, lookahead, minOverlap);
                        int fragLen = fragment.length();
                        
                        if (fragLen > 0) {                            
                            /** extend on both ends unambiguously*/
                            fragment = naiveExtend(fragment, graph, maxTipLen);
                            
                            ++fid;
                            
                            if (fid > sampleSize) {
                                /** store paired kmers */
                                graph.addPairedKmersFromSeq(fragment);
                                
                                /** write fragment */
                                out.write(">" + fid);
                                out.newLine();
                                out.write(fragment);
                                out.newLine();
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
                                fid = 0;
                                for (String frag : sampleFragments) {
                                    graph.addPairedKmersFromSeq(frag);

                                    /** write fragment */
                                    out.write(">" + ++fid);
                                    out.newLine();
                                    out.write(frag);
                                    out.newLine();
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
                        
                        /*
                        if (p.right.endsWith("AAA")) {
                            String fragment = assembleFragment(p.left, p.right, graph, mismatchesAllowed, bound, lookahead, minOverlap);
                            System.out.println("LEFT:  " + p.left);
                            System.out.println("RIGHT: " + p.right);
                            System.out.println(fragment.length() + ": " + fragment);
                        }
                        */

                        /**@TODO assemble with paired k-mers bloom filter*/
                    }
                }

                lin.close();
                rin.close();
            }
            
            out.close();
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
        
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        /*
        String kmer = "AAAAATTTTTCCCCCGGGGG";
        String kmer2 = "AACAATTTTTCCCCCGGGGG";
        
        BloomFilter f1 = new BloomFilter(NUM_BITS_1GB, 3, 689, kmer.length());
        f1.add(kmer);
        System.out.println(f1.lookup(kmer));
        System.out.println(f1.lookup(kmer2));
        
        CountingBloomFilter f2 = new CountingBloomFilter((int) NUM_BYTES_1GB, 3, 689, kmer.length());
        
        for (int i=0; i < 20; ++i) {
            f2.increment(kmer);
        }
        
        System.out.println(f2.getCount(kmer));
        */
        
        /*
        String seq  = "NATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTTTCAAAAACCTCCTGCCGATTCTTGGCAGGCAAACATAC";
        String qual = "#0<BFFBFFBFFFIFIIFFFFIIIFIIBFFFBFFFBFFF<BBBBB7F######################################################";
        FastqRecord fq = new FastqRecord(seq, qual);
        int k = 25;
        int q = 3;
        Pattern p = SequenceOperations.getPhred33Pattern(q, k);
        System.out.println( SequenceOperations.filterFastq(fq, p).toString() );
        */
        
        
        String fastq1 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz"; //right
        String fastq2 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz"; //left        
        String fragsFasta = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/fragments.fa";
        
        //String fastq1 = "/home/gengar/test_data/GAPDH/GAPDH_1.fq.gz";
        //String fastq2 = "/home/gengar/test_data/GAPDH/GAPDH_2.fq.gz";
        
        boolean revComp1 = true;
        boolean revComp2 = false;
        
        int mismatchesAllowed = 5;
        int bound = 500;
        int lookahead = 5;
        int minOverlap = 10;
        int maxTipLen = 5;
        
        String[] forwardFastqs = new String[]{fastq2};
        String[] backwardFastqs = new String[]{fastq1};
        
        boolean strandSpecific = true;
        long dbgbfSize = NUM_BITS_1GB;
        long cbfSize = NUM_BYTES_1GB;
        long pkbfSize = NUM_BITS_1GB;
        int dbgbfNumHash = 3;
        int cbfNumHash = 4;
        int pkbfSingleKeyNumHash = 1;
        int pkbfPairedKeysNumHash = 1;
        int seed = 689;
        int k = 25;
        int q = 3;
        int sampleSize = 1000;
                
                
        RNABloom assembler = new RNABloom(strandSpecific, dbgbfSize, cbfSize, pkbfSize, dbgbfNumHash, cbfNumHash, pkbfSingleKeyNumHash, pkbfPairedKeysNumHash, seed, k, q);
        assembler.createDBG(forwardFastqs, backwardFastqs);
        
        /*
        String left = "AAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCA";
        String right = "GGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACCCCCCCTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTT";
        String fragment = assembleFragment(left, right, assembler.getGraph(), mismatchesAllowed, bound, lookahead, minOverlap);
        System.out.println(fragment);
        */
        
        FastqPair fqPair = new FastqPair(fastq2, fastq1, revComp2, revComp1);
        FastqPair[] fqPairs = new FastqPair[]{fqPair};
        
        assembler.assembleFragments(fqPairs, fragsFasta, mismatchesAllowed, bound, lookahead, minOverlap, maxTipLen, sampleSize);
    }
    
}
