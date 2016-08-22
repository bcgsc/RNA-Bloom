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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import rnabloom.io.FastqReader;
import rnabloom.io.FastqRecord;
import static rnabloom.util.GraphUtils.assemble;
import static rnabloom.util.GraphUtils.correctMismatches;
import static rnabloom.util.GraphUtils.getMaxCoveragePath;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    
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
        
        
        long dbgbfSize = NUM_BITS_1GB;
        int cbfSize = (int) NUM_BYTES_1GB;
        int dbgbfNumHash = 3;
        int cbfNumHash = 4;
        int seed = 689;
        int k = 25;
        int q = 3;
        
        BloomFilterDeBruijnGraph graph = new BloomFilterDeBruijnGraph(dbgbfSize,
                                                                    cbfSize,
                                                                    dbgbfNumHash,
                                                                    cbfNumHash,
                                                                    seed,
                                                                    k,
                                                                    true);
        
        Pattern qualPattern = getPhred33Pattern(q, k);
        
        String fastq1 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz";
        String fastq2 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz";

        int lineNum = 0;
        
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq1))));
            FastqReader fr = new FastqReader(br);
            Stream<FastqRecord> s = Stream.generate(fr);
            Iterator<FastqRecord> itr = s.iterator();
            try {
                while (true) {
                    if (++lineNum % 100000 == 0) {
                        System.out.println("Parsed " + lineNum + " reads...");
                    }

                    for (String seq : filterFastq(itr.next(), qualPattern)) {
                        for (String kmer : kmerize(reverseComplement(seq), k)) {
                            graph.add(kmer);
                        }
                    }
                }
            }
            catch (NoSuchElementException e) {
                //do nothing
            }
            
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq2))));
            fr = new FastqReader(br);
            s = Stream.generate(fr);
            itr = s.iterator();
            try {
                while (true) {
                    if (++lineNum % 100000 == 0) {
                        System.out.println("Parsed " + lineNum + " reads...");
                    }

                    for (String seq : filterFastq(itr.next(), qualPattern)) {
                        for (String kmer : kmerize(seq, k)) {
                            graph.add(kmer);
                        }
                    }
                }
            }
            catch (NoSuchElementException e) {
                //do nothing
            }
        }
        catch (IOException e) {
            
        }
        finally {
            System.out.println("Parsed " + lineNum + " reads...");
        }
        

        String read1 = "TCGCTGTTGAAGTCAGAGGAGACCACCTGGTGCTCAGTGTAGCCCAGGATGCCCTTGAGGGGGCCCTCCGACGCCTGCTTCACCACCTTCTTGATGTCAT";
        String read2 = "ACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTAAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACAT";
        read1 = reverseComplement(read1);
        
        
        //System.out.println(read1);
        //String read1Corrected = assemble(correctMismatches(read1, graph, 5, 5));
        //System.out.println(read1Corrected);
        
        //System.out.println(read2);
        //String read2Corrected = assemble(correctMismatches(read2, graph, 5, 5));
        //System.out.println(read2Corrected);
        
        ArrayList<Kmer> rightKmers = correctMismatches(read1, graph, 5, 5);
        ArrayList<Kmer> leftKmers = correctMismatches(read2, graph, 5, 5);
        
        String path = assemble(getMaxCoveragePath(graph, leftKmers.get(leftKmers.size()-1), rightKmers.get(0), 500, 5));
        System.out.println(path);
    }
    
}
