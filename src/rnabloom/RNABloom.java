/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.File;
import java.io.IOException;
import static java.lang.Math.pow;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
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
import static rnabloom.util.GraphUtils.extendWithPairedKmers;
import static rnabloom.util.GraphUtils.findBackbonePath;
import static rnabloom.util.GraphUtils.naiveExtend;
import static rnabloom.util.SeqUtils.*;

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
    private final BloomFilterDeBruijnGraph graph;
    
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
                String transcript = extendWithPairedKmers(fragment, graph, lookAhead);
                
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
    
    public static void test2() {
        String bfPath = "/projects/btl2/kmnip/rna-bloom/tests/test.bf";
        String bfDescPath = bfPath + ".desc";
        
        String kmer =  "AAAAATTTTTCCCCCGGGGG";
        String kmer2 = "AAAAAAAAAACCCCCGGGGG";
        
        CountingBloomFilter bf = new CountingBloomFilter(NUM_BYTES_1GB*4, 3, new HashFunction(3, 689, kmer.length()));
        bf.increment(kmer);
        System.out.println("true:" + bf.getCount(kmer));
        System.out.println("false:" + bf.getCount(kmer2));
        
        try {
            bf.write(new File(bfDescPath), new File(bfPath));
            bf.destroy();
            bf = null;
            
            bf = CountingBloomFilter.read(new File(bfDescPath), new File(bfPath));
            
            System.out.println("true:" + bf.getCount(kmer));
            System.out.println("false:" + bf.getCount(kmer2));
            
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void test() {
        //String f = "TCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCC";
        String f = "CCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCT";
        String transcript = extendWithPairedKmers(f, graph, 5);
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
        
        ///*
        String fastq1 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz"; //right
        String fastq2 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz"; //left        
        String fragsFasta = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/fragments.fa";
        String transcriptsFasta = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/transcripts.fa";
        //*/
        /*
        String fastq1 = "/home/gengar/test_data/GAPDH/GAPDH_1.fq.gz";
        String fastq2 = "/home/gengar/test_data/GAPDH/GAPDH_2.fq.gz";
        String fragsFasta = "/home/gengar/test_assemblies/GAPDH/fragments.fa";
        String transcriptsFasta = "/home/gengar/test_assemblies/GAPDH/transcripts.fa";
        */        
        
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
        int pkbfNumHash = 1;
        int seed = 689;
        int k = 25;
        int q = 3;
        int sampleSize = 1000;

        test2();
        
        /*        
        RNABloom assembler = new RNABloom(strandSpecific, dbgbfSize, cbfSize, pkbfSize, dbgbfNumHash, cbfNumHash, pkbfNumHash, seed, k, q);
        assembler.createDBG(forwardFastqs, backwardFastqs);
                
        FastqPair fqPair = new FastqPair(fastq2, fastq1, revComp2, revComp1);
        FastqPair[] fqPairs = new FastqPair[]{fqPair};
        
        assembler.assembleFragments(fqPairs, fragsFasta, mismatchesAllowed, bound, lookahead, minOverlap, maxTipLen, sampleSize);
        
        assembler.test();
        
        //assembler.assembleTranscripts(fragsFasta, transcriptsFasta, lookahead);
        */
    }
    
}
