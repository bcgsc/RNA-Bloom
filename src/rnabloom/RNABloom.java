/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import static java.lang.Math.pow;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.ReadPair;
import rnabloom.io.FastqReader;
import static rnabloom.util.GraphUtils.assembleFragment;
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
    
    public RNABloom(boolean strandSpecific, long dbgbfSize, long cbfSize, int dbgbfNumHash, int cbfNumHash, int seed, int k, int q) {
        this.k = k;
        this.qualPattern = getPhred33Pattern(q, k);
        
        graph = new BloomFilterDeBruijnGraph(dbgbfSize,
                                            cbfSize,
                                            dbgbfNumHash,
                                            cbfNumHash,
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
        
    public void assembleFragments(String leftFastq, String rightFastq, boolean revCompLeft, boolean revCompRight, int mismatchesAllowed, int bound, int lookahead, int minOverlap) {
        BufferedReader lbr, rbr;
        FastqReader leftReader, rightReader;

        try {
            if (leftFastq.endsWith(GZIP_EXTENSION)) {
                lbr = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(leftFastq))));
            }
            else {
                lbr = new BufferedReader(new InputStreamReader(new FileInputStream(leftFastq)));
            }
            leftReader = new FastqReader(lbr, true);

            if (rightFastq.endsWith(GZIP_EXTENSION)) {
                rbr = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(rightFastq))));
            }
            else {
                rbr = new BufferedReader(new InputStreamReader(new FileInputStream(rightFastq)));
            }
            rightReader = new FastqReader(rbr, true);
            
            FastqPairReader fqpr = new FastqPairReader(leftReader, rightReader, qualPattern, revCompLeft, revCompRight);
            ReadPair p;
            while (fqpr.hasNext()) {
                p = fqpr.next();
                
                System.out.println("LEFT:  " + p.left);
                System.out.println("RIGHT: " + p.right);
                
                String fragment = assembleFragment(p.left, p.right, graph, mismatchesAllowed, bound, lookahead, minOverlap);
                
                System.out.println(fragment);
                /**@TODO*/
            }
            
            lbr.close();
            rbr.close();
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
        
        boolean revCompLeft = false;
        boolean revCompRight = true;
        int mismatchesAllowed = 5;
        int bound = 500;
        int lookahead = 5;
        int minOverlap = 10;
        
        String[] forwardFastqs = new String[]{fastq2};
        String[] backwardFastqs = new String[]{fastq1};
        
        boolean strandSpecific = true;
        long dbgbfSize = NUM_BITS_1GB;
        long cbfSize = NUM_BYTES_1GB;
        int dbgbfNumHash = 3;
        int cbfNumHash = 4;
        int seed = 689;
        int k = 25;
        int q = 3;
                
                
        RNABloom assembler = new RNABloom(strandSpecific, dbgbfSize, cbfSize, dbgbfNumHash, cbfNumHash, seed, k, q);
        assembler.createDBG(forwardFastqs, backwardFastqs);
        
        //String left = "AAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCA";
        //String right = "GGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACCCCCCCTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTT";
        
        
        
        //assembler.assembleFragments(fastq2, fastq1, revCompLeft, revCompRight, mismatchesAllowed, bound, lookahead, minOverlap);
    }
    
}
