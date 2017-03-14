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
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Option.Builder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.CanonicalNTHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.bloom.hash.ReverseComplementNTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqPair;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.FastqReadPair;
import rnabloom.io.FastqReader;
import rnabloom.io.FastqRecord;
import static rnabloom.util.GraphUtils.*;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    private final static long NUM_PARSED_INTERVAL = 100000;
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    
    private int k;
    private int kMinus1;
//    private boolean strandSpecific;
    private Pattern qualPatternDBG;
    private Pattern qualPatternFrag;
    private Pattern homoPolymerKmerPattern;
    private BloomFilterDeBruijnGraph graph = null;
    private BloomFilter screeningBf = null;

    private int maxTipLength;
    private int lookahead;
    private float maxCovGradient;
    private int maxIndelSize;
    private float percentIdentity;
    private float percentError;
    private int minNumKmerPairs = 5;
    
    private float dbgFPR = -1;
    private float covFPR = -1;
    private final static String[] COVERAGE_ORDER = {"e0", "e1", "e2", "e3", "e4", "e5"};
    
    public RNABloom(int k, int qDBG, int qFrag) {
        this.k = k;
        this.kMinus1 = k-1;
        this.qualPatternDBG = getPhred33Pattern(qDBG, k);
        this.qualPatternFrag = getPhred33Pattern(qFrag, k);
        this.homoPolymerKmerPattern = getHomoPolymerPattern(k);
    }
    
    public void setParams(int maxTipLength, int lookahead, float maxCovGradient, int maxIndelSize, float percentIdentity) {
        this.maxTipLength = maxTipLength;
        this.lookahead = lookahead;
        this.maxCovGradient = maxCovGradient;
        this.maxIndelSize = maxIndelSize;
        this.percentIdentity = percentIdentity;
        this.percentError = 1.0f - percentIdentity;
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
            if (graph != null) {
                graph.destroy();
            }
            
            graph = new BloomFilterDeBruijnGraph(f);

            dbgFPR = graph.getDbgbfFPR();
            covFPR = graph.getCbfFPR();
            
            //BloomFilterDeBruijnGraph graph2 = new BloomFilterDeBruijnGraph(f);
            //System.out.println(graph2.getDbgbf().equivalent(graph.getDbgbf()));
            //System.out.println(graph2.getCbf().equivalent(graph.getCbf()));
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public int addKmersFromFastq(String fastq, boolean stranded, boolean reverseComplement, int numHash) throws IOException {
        int numReads = 0;
        
        NTHashIterator itr;
        
        if (stranded) {
            if (reverseComplement) {
                itr = new ReverseComplementNTHashIterator(k, numHash);
            }
            else {
                itr = new NTHashIterator(k, numHash);
            }
        }
        else {
            itr = new CanonicalNTHashIterator(k, numHash);
        }
        
        FastqReader fr = new FastqReader(fastq, false);
        FastqRecord record = fr.record;
        Matcher m = qualPatternDBG.matcher("");
        long[] hashVals = itr.hVals;

        while (fr.hasNext()) {
            ++numReads;

            fr.nextWithoutNameFunction();
            m.reset(record.qual);

            while (m.find()) {
                itr.start(record.seq.substring(m.start(), m.end()));
                while (itr.hasNext()) {
                    itr.next();
                    if (screeningBf.lookupThenAdd(hashVals)) {
                        graph.add(hashVals);
                    }
                }
            }
        }
        fr.close();
        
        return numReads;
    }

    public class FastqParser implements Runnable {
        private final int id;
        
        private final String fastq;
        
        private final NTHashIterator itr;
        private int numReads = 0;
        private boolean successful = false;
        
        public FastqParser(int id, String fastq, boolean stranded, boolean reverseComplement, int numHash) {
            this.id = id;
            this.fastq = fastq;
            
            if (stranded) {
                if (reverseComplement) {
                    itr = new ReverseComplementNTHashIterator(k, numHash);
                }
                else {
                    itr = new NTHashIterator(k, numHash);
                }
            }
            else {
                itr = new CanonicalNTHashIterator(k, numHash);
            }
        }
        
        @Override
        public void run() {
            System.out.println("[" + id + "] Parsing `" + fastq + "`...");
            
            try {
                FastqReader fr = new FastqReader(fastq, false);
                FastqRecord record = fr.record;
                Matcher m = qualPatternDBG.matcher("");
                long[] hashVals = itr.hVals;
                
                while (fr.hasNext()) {
                    ++numReads;
                    
                    fr.nextWithoutNameFunction();
                    m.reset(record.qual);
                    
                    while (m.find()) {
                        itr.start(record.seq.substring(m.start(), m.end()));
                        while (itr.hasNext()) {
                            itr.next();
//                            if (screeningBf.lookupThenAdd(hashVals)) {
                                graph.add(hashVals);
//                            }
                        }
                    }
                }
                fr.close();
                
                successful = true;
                System.out.println("[" + id + "] Parsed " + NumberFormat.getInstance().format(numReads) + " reads.");
            } catch (IOException e) {
                /**@TODO */
                e.printStackTrace();
            }
            catch (NoSuchElementException e) {
                /**@TODO handle invalid format*/
                e.printStackTrace();
            }
        }
        
        public boolean isSuccessful() {
            return successful;
        }
        
        public long getReadCount() {
            return numReads;
        }
    }
    
    public void createGraph(String[] forwardFastqs,
                            String[] reverseFastqs,
                            boolean strandSpecific,
//                            long sbfNumBits,
                            long dbgbfNumBits,
                            long cbfNumBytes,
                            long pkbfNumBits,
//                            int sbfNumHash,
                            int dbgbfNumHash,
                            int cbfNumHash,
                            int pkbfNumHash,
                            int numThreads) {        
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            k,
                                            strandSpecific);
        
//        screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
        
        /** parse the reads */
        
        int numReads = 0;
        int numHash = graph.getMaxNumHash();
        
        ExecutorService service = Executors.newFixedThreadPool(numThreads);
        
        ArrayList<FastqParser> threadPool = new ArrayList<>();
        int threadId = 0;
           
        for (String fastq : forwardFastqs) {
            FastqParser t = new FastqParser(++threadId, fastq, strandSpecific, false, numHash);
            service.submit(t);
            threadPool.add(t);
        }

        for (String fastq : reverseFastqs) {
            FastqParser t = new FastqParser(++threadId, fastq, strandSpecific, true, numHash);
            service.submit(t);
            threadPool.add(t);
        }

        service.shutdown();
        
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            
            for (FastqParser t : threadPool) {
                numReads += t.getReadCount();
            }

            System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");
//            System.out.println("Screening Bloom filter FPR:  " + screeningBf.getFPR() * 100 + " %");
            
//            screeningBf.destroy();
            
        } catch (InterruptedException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        dbgFPR = graph.getDbgbfFPR();
        covFPR = graph.getCbfFPR();
    }
    
    public void insertIntoDeBruijnGraph(String[] fastas) throws IOException { 
        
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;
        
        for (String fasta : fastas) {
            FastaReader fin = new FastaReader(fasta);
            
            while (fin.hasNext()) {
                itr.start(fin.next());
                while (itr.hasNext()) {
                    itr.next();
                    graph.addDbgOnly(hashVals);
                }
            }
            
            fin.close();
        }
    }
    
    private boolean isPolyA(String left, String right) {
        return right.endsWith("AAAA") || (!graph.isStranded() && left.startsWith("TTTT"));
    }
    
//    private boolean okToConnectPair(String left, String right) {
//        NTHashIterator itr = graph.getHashIterator();
//        itr.start(left);
//        long[] hVals = itr.hVals;
//        float c;
//        
//        while (itr.hasNext()) {
//            if (!graph.lookupLeftKmer(hVals)) {
//                return true;
//            }
//        }
//        
//        itr.start(right);
//        
//        while (itr.hasNext()) {
//            if (!graph.lookupRightKmer(hVals)) {
//                return true;
//            }
//        }
//        
//        return false;
//        
////        int numKmersNotSeenLeft = 0;
////        float minCovLeft = Float.MAX_VALUE;
////        
////        // scan the entire read for c=0 kmers
////        while (itr.hasNext()) {
////            itr.next();
////
////            c = graph.getCount(hVals);
////            
////            if (c == 0) {
////                return false;
////            }
////            
////            if (!graph.lookupFragmentKmer(hVals)) {
////                ++numKmersNotSeenLeft;
////                
////                if (c < minCovLeft) {
////                    minCovLeft = c;
////                }
////            }
////        }
////        
////        itr.start(right);
////        
////        int numKmersNotSeenRight = 0;
////        float minCovRight = Float.MAX_VALUE;
////        
////        // scan the entire read for c=0 kmers
////        while (itr.hasNext()) {
////            itr.next();
////
////            c = graph.getCount(hVals);
////            
////            if (c == 0) {
////                return false;
////            }
////            
////            if (!graph.lookupFragmentKmer(hVals)) {
////                ++numKmersNotSeenRight;
////                
////                if (c < minCovRight) {
////                    minCovRight = c;
////                }
////            }
////        }
////        return numKmersNotSeenLeft >= k || numKmersNotSeenLeft >= minCovLeft || numKmersNotSeenLeft == getNumKmers(left, k ) ||
////                numKmersNotSeenRight >= k || numKmersNotSeenRight >= minCovRight || numKmersNotSeenRight == getNumKmers(right, k);
//    }
    
    private boolean okToConnectPair(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers) {
        if (leftKmers.isEmpty() || rightKmers.isEmpty()) {
            return false;
        }
        
        for (Kmer kmer : leftKmers) {
            if (!graph.lookupLeftKmer(kmer.hashVals)) {
                return true;
            }
        }

        for (Kmer kmer : rightKmers) {
            if (!graph.lookupRightKmer(kmer.hashVals)) {
                return true;
            }
        }
        
        return false;
        
//        int numKmersNotSeenLeft = 0;
//        for (Kmer kmer : leftKmers) {
//            if (!graph.lookupLeftKmer(kmer.hashVals)) {
//                if (++numKmersNotSeenLeft >= kMinus1) {
//                    return true;
//                }
//            }
//        }
//        
//        int numKmersNotSeenRight = 0;
//        for (Kmer kmer : rightKmers) {
//            if (!graph.lookupRightKmer(kmer.hashVals)) {
//                if (++numKmersNotSeenRight >= kMinus1) {
//                    return true;
//                }
//            }
//        }
//        
//        return numKmersNotSeenLeft + numKmersNotSeenRight >= kMinus1;
    }
    
//    
//    private int findBackboneIdNonStranded(String fragment) {
//        KmerIterator itr = graph.new KmerIterator(fragment);
//        Kmer seed = itr.next();
//        Kmer kmer;
//        String seq;
//        
//        while (itr.hasNext()) {
//            kmer = itr.next();
//            
//            seq = smallestStrand(kmer.seq);
//            
//            if (kmerToBackboneID.containsKey(seq)) {
//                return kmerToBackboneID.get(seq);
//            }
//            
//            if (kmer.count > seed.count) {
//                seed = kmer;
//            }
//        }
//        
//        ArrayList<Kmer> path = null;
//        boolean randomSeed = false;
//        for (int i=0; i<bbMaxIteration; ++i) {
//            if (i>0) {
//                if (randomSeed) {
//                    seed = path.get(random.nextInt(path.size()));
//                    randomSeed = false;
//                }
//                else {
//                    seed = findMaxCoverageWindowKmer(path, graph, bbWindowSize);
//                    randomSeed = true;
//                }
//            }
//
//            /* greedy extend on both sides */
//            HashSet<String> pathKmerStr = new HashSet<>(1000);
//            pathKmerStr.add(smallestStrand(seed.seq));
//            
//            /* extend on right side */
//            ArrayList<Kmer> rightPath = new ArrayList<>(1000);
//            Kmer best = seed;
//            while (true) {
//                best = greedyExtendRightOnce(graph, best, bbLookahead);
//                if (best != null) {
//                    seq = smallestStrand(best.seq);
//                    
//                    if (kmerToBackboneID.containsKey(seq)) {
//                        return kmerToBackboneID.get(seq);
//                    }
//                    
//                    if (pathKmerStr.contains(seq)) {
//                        break;
//                    }
//                    
//                    pathKmerStr.add(seq);
//                    rightPath.add(best);
//                }
//                else {
//                    break;
//                }
//            }
//
//            /* extend on left side */
//            ArrayList<Kmer> leftPath = new ArrayList<>(1000);
//            best = seed;
//            while (true) {
//                best = greedyExtendLeftOnce(graph, best, bbLookahead);
//                if (best != null) {
//                    seq = smallestStrand(best.seq);
//                    
//                    if (kmerToBackboneID.containsKey(seq)) {
//                        return kmerToBackboneID.get(seq);
//                    }
//                    
//                    if (pathKmerStr.contains(seq)) {
//                        break;
//                    }
//                    
//                    pathKmerStr.add(seq);
//                    leftPath.add(best);
//                }
//                else {
//                    break;
//                }
//            }
//
//            Collections.reverse(leftPath);
//            leftPath.add(seed);
//            leftPath.addAll(rightPath);
//
//            path = leftPath;
//        }
//        
//        /* new backbone */
//        //backbones.add(assemble(path));
//        int id = ++currentBackboneId;
//        
//        /* store kmers in path */
//        int numKmers = path.size();
//        for (int i=0; i<numKmers; ++i) {
//            if (i % backboneHashKmerDistance == 0) {
//                kmerToBackboneID.put(smallestStrand(path.get(i).seq), i);
//            }
//        }
//        
//        return id;
//    }
//    
//    private int findBackboneIdStranded(String fragment) {
//        KmerIterator itr = graph.new KmerIterator(fragment);
//        Kmer seed = itr.next();
//        Kmer kmer;
//        while (itr.hasNext()) {
//            kmer = itr.next();
//            
//            if (kmerToBackboneID.containsKey(kmer.seq)) {
//                return kmerToBackboneID.get(kmer.seq);
//            }
//            
//            if (kmer.count > seed.count) {
//                seed = kmer;
//            }
//        }
//        
//        ArrayList<Kmer> path = null;
//        boolean randomSeed = false;
//        for (int i=0; i<bbMaxIteration; ++i) {
//            if (i>0) {
//                if (randomSeed) {
//                    seed = path.get(random.nextInt(path.size()));
//                    randomSeed = false;
//                }
//                else {
//                    seed = findMaxCoverageWindowKmer(path, graph, bbWindowSize);
//                    randomSeed = true;
//                }
//            }
//
//            /* greedy extend on both sides */
//            HashSet<String> pathKmerStr = new HashSet<>(1000);
//            pathKmerStr.add(seed.seq);
//            
//            /* extend on right side */
//            ArrayList<Kmer> rightPath = new ArrayList<>(1000);
//            Kmer best = seed;
//            while (true) {
//                best = greedyExtendRightOnce(graph, best, bbLookahead);
//                if (best != null) {
//                    String seq = best.seq;
//                    
//                    if (kmerToBackboneID.containsKey(seq)) {
//                        return kmerToBackboneID.get(seq);
//                    }
//                    
//                    if (pathKmerStr.contains(seq)) {
//                        break;
//                    }
//                    
//                    pathKmerStr.add(seq);
//                    rightPath.add(best);
//                }
//                else {
//                    break;
//                }
//            }
//
//            /* extend on left side */
//            ArrayList<Kmer> leftPath = new ArrayList<>(1000);
//            best = seed;
//            while (true) {
//                best = greedyExtendLeftOnce(graph, best, bbLookahead);
//                if (best != null) {
//                    String seq = best.seq;
//                    
//                    if (kmerToBackboneID.containsKey(seq)) {
//                        return kmerToBackboneID.get(seq);
//                    }
//                    
//                    if (pathKmerStr.contains(seq)) {
//                        break;
//                    }
//                    
//                    pathKmerStr.add(seq);
//                    leftPath.add(best);
//                }
//                else {
//                    break;
//                }
//            }
//
//            Collections.reverse(leftPath);
//            leftPath.add(seed);
//            leftPath.addAll(rightPath);
//
//            path = leftPath;
//            
//            //System.out.println(">" + i + "\n" + assemble(path));
//        }
//        
//        /* new backbone */
//        //backbones.add(assemble(path));
//        int id = ++currentBackboneId;;
//        
//        /* store kmers in path */
//        
//        //System.out.println(">bb\n" + assemble(path));
//        
//        int numKmers = path.size();
//        for (int i=0; i<numKmers; ++i) {
//            if (i % backboneHashKmerDistance == 0) {
//                kmerToBackboneID.put(path.get(i).seq, id);
//            }
//        }
//        
//        return id;
//    }
    
    public static class ReadPair {
        ArrayList<Kmer> leftKmers;
        ArrayList<Kmer> rightKmers;
        boolean corrected = false;
        
        public ReadPair(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers, boolean corrected) {
            this.leftKmers = leftKmers;
            this.rightKmers = rightKmers;
            this.corrected = corrected;
        }
    }
    
    private class Fragment {
//        String left;
//        String right;
        String seq;
        int length;
        float minCov;
        boolean isUnconnectedRead;
        
        public Fragment(String seq, int length, float minCov, boolean isUnconnectedRead) {
//            this.left = left;
//            this.right = right;
            this.seq = seq;
            this.length = length;
            this.minCov = minCov;
            this.isUnconnectedRead = isUnconnectedRead;
        }
    }
    
    private class Transcript {
        String fragment;
        String transcript;
        
        public Transcript(String fragment, String transcript) {
            this.fragment = fragment;
            this.transcript = transcript;
        }
    }
    
    private class AssembledTranscriptsQueue {
        
        private ArrayBlockingQueue<Transcript> queue;
        
        public AssembledTranscriptsQueue(int size) {
            queue = new ArrayBlockingQueue<> (size);
        }
        
        public synchronized void add(String fragment, ArrayList<Kmer> kmers) {
            if (!represented(kmers,
                                graph,
                                screeningBf,
                                lookahead,
                                maxIndelSize,
                                percentIdentity)) {

                for (Kmer kmer : kmers) {
                    screeningBf.add(kmer.hashVals);
                }

                String transcript = assemble(kmers, k);

                try {
                    queue.put(new Transcript(fragment, transcript));
                } catch (InterruptedException ex) {
                    Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        public int remainingCapacity() {
            return queue.remainingCapacity();
        }
        
        public boolean isEmpty() {
            return queue.isEmpty();
        }
        
        public Transcript poll() {
            return queue.poll();
        }
    }
    
    private class TranscriptAssembler implements Runnable {
        
        private String fragment;
        private boolean beGreedy;
        private AssembledTranscriptsQueue outList;
        
        public TranscriptAssembler(String fragment,
                                    AssembledTranscriptsQueue outList,
                                    boolean beGreedy) {
            this.fragment = fragment;
            this.outList = outList;
            this.beGreedy = beGreedy;
        }

        @Override
        public void run() {
//            ArrayList<Kmer> fragKmers = graph.getKmers(correctMismatches(fragment, graph, lookahead, (int) Math.ceil(fragment.length()*percentError)));
            ArrayList<Kmer> fragKmers = graph.getKmers(fragment);
            
            if (!represented(fragKmers,
                                graph,
                                screeningBf,
                                lookahead,
                                maxIndelSize,
                                percentIdentity)) {
                
                extendWithPairedKmers2(fragKmers, graph, lookahead, maxTipLength, beGreedy, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs);
                
                outList.add(fragment, fragKmers);
            }
        }
    }
    
    private class FragmentAssembler implements Runnable {
        private FastqReadPair p;
        private ArrayBlockingQueue<Fragment> outList;
        private int bound;
        private int minOverlap;
        private boolean storeKmerPairs;
        private int errorCorrectionIterations;
        private int leftReadLengthThreshold;
        private int rightReadLengthThreshold;
        
        public FragmentAssembler(FastqReadPair p,
                                ArrayBlockingQueue<Fragment> outList,
                                int bound, 
                                int minOverlap, 
                                boolean storeKmerPairs, 
                                int errorCorrectionIterations,
                                int leftReadLengthThreshold,
                                int rightReadLengthThreshold) {
            
            this.p = p;
            this.outList = outList;
            this.bound = bound;
            this.minOverlap = minOverlap;
            this.storeKmerPairs = storeKmerPairs;
            this.errorCorrectionIterations = errorCorrectionIterations;
            this.leftReadLengthThreshold = leftReadLengthThreshold;
            this.rightReadLengthThreshold = rightReadLengthThreshold;
        }
        
        @Override
        public void run() {
            // connect segments of each read
            String left = connect(p.left, graph, lookahead);
            String right = connect(p.right, graph, lookahead);

            if (left.length() >= this.leftReadLengthThreshold 
                    && right.length() >= this.rightReadLengthThreshold) { 
                                
                ArrayList<Kmer> leftKmers = graph.getKmers(left);
                ArrayList<Kmer> rightKmers = graph.getKmers(right);
                
                if (okToConnectPair(leftKmers, rightKmers)) {
                    boolean corrected = false;
                    
                    if (this.errorCorrectionIterations > 0) {
                        
                        ReadPair correctedReadPair = correctErrors2(leftKmers,
                                                            rightKmers,
                                                            graph, 
                                                            lookahead, 
                                                            maxIndelSize, 
                                                            maxCovGradient, 
                                                            covFPR,
                                                            this.errorCorrectionIterations,
                                                            2,
                                                            percentIdentity);
                        
                        if (correctedReadPair.corrected) {
                            corrected = true;
                            leftKmers = correctedReadPair.leftKmers;
                            rightKmers = correctedReadPair.rightKmers;
                        }
                    }
                    
                    if (!corrected || okToConnectPair(leftKmers, rightKmers)) {
                        ArrayList<Kmer> fragmentKmers = null;
                        
                        if (!isLowComplexity2(leftKmers.get(leftKmers.size()-1).seq) &&  
                                !isLowComplexity2(rightKmers.get(0).seq)) {
                            fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound, lookahead, minOverlap);
                        }
                        
                        if (fragmentKmers != null) {
                            int fragLength = fragmentKmers.size() + k - 1;
                            
                            if (fragLength >= k + lookahead) {
                                if (this.storeKmerPairs) {
                                    graph.addPairedKmers(fragmentKmers);
                                }

                                float minCov = Float.MAX_VALUE;
                                for (Kmer kmer : fragmentKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }
                                }

                                try {
                                    outList.put(new Fragment(assemble(fragmentKmers, k), fragLength, minCov, false));
                                } catch (InterruptedException ex) {
                                    ex.printStackTrace();
                                }
                            }
                        }
                        else {
                            // write unconnected reads to file
                            try {
                                float minCov = Float.MAX_VALUE;
                                
                                if (leftKmers.size() >= lookahead) {
                                    for (Kmer kmer : leftKmers) {
                                        if (kmer.count < minCov) {
                                            minCov = kmer.count;
                                        }
                                    }
                                    outList.put(new Fragment(assemble(leftKmers, k), leftKmers.size()+k-1, minCov, true));
                                }
                                
                                if (rightKmers.size() >= lookahead) {
                                    minCov = Float.MAX_VALUE;
                                    for (Kmer kmer : rightKmers) {
                                        if (kmer.count < minCov) {
                                            minCov = kmer.count;
                                        }
                                    }
                                    outList.put(new Fragment(assemble(rightKmers, k), rightKmers.size()+k-1, minCov, true));
                                }
                            } catch (InterruptedException ex) {
                                ex.printStackTrace();
                            }
                        }
                    }
                }
            }
        }
    }
    
    public class MyExecutorService {
        private final BlockingQueue<Runnable> queue;
        private final ExecutorService service;
        
        public MyExecutorService(int numThreads, int queueSize) {
            queue = new ArrayBlockingQueue<>(queueSize);    
            service = new ThreadPoolExecutor(numThreads, numThreads,
                        0L, TimeUnit.MILLISECONDS,
                        queue);
        }
        
        public void submit(Runnable r) {
            while (true) {
                try {
                    service.submit(r);
                    break;
                }
                catch(RejectedExecutionException e) {
                    // do nothing
                }
            }
        }
        
        public void terminate() throws InterruptedException {
            service.shutdown();
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        }
    }
        
    public int[] getMinQ1MedianQ3Max(ArrayList<Integer> a) {
        Collections.sort(a);
        
        int arrLen = a.size();
        int halfLen = arrLen/2;
        int q1Index = arrLen/4;
        int q3Index = halfLen+q1Index;
        
        int q1, median, q3;
        
        if (arrLen % 2 == 0) {
            median = (a.get(halfLen-1) + a.get(halfLen))/2;
        }
        else {
            median = a.get(halfLen);
        }
        
        if (arrLen % 4 == 0) {
            q1 = (a.get(q1Index-1) + a.get(q1Index))/2;
            q3 = (a.get(q3Index-1) + a.get(q3Index))/2;
        }
        else {
            q1 = a.get(q1Index);
            q3 = a.get(q3Index);
        }
        
        return new int[]{a.get(0), q1, median, q3, a.get(arrLen-1)};
    }
    
    public static int getMinCoverageOrderOfMagnitude(float c) {
        if (c >= 1e5) {
            return 5;
        }
        else if (c >= 1e4) {
            return 4;
        }
        else if (c >= 1e3) {
            return 3;
        }
        else if (c >= 1e2) {
            return 2;
        }
        else if (c >= 1e1) {
            return 1;
        }
        else if (c >= 1e0) {
            return 0;
        }
        else {
            return -1;
        }
    }
    
    public void assembleFragmentsMultiThreaded(FastqPair[] fastqs, 
                                                String[] longFragmentsFastaPaths,
                                                String[] shortFragmentsFastaPaths,
                                                int bound,
                                                int minOverlap,
                                                int sampleSize, 
                                                int numThreads, 
                                                int maxErrCorrIterations) {
        
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }
        
        System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
        System.out.println("Counting Bloom filter FPR: " + covFPR * 100 + " %");
        
        
        System.out.println("Assembling fragments...");
        
        graph.initializePairKmersBloomFilter();
        
        long fragmentId = 0;
        long readPairsParsed = 0;
        
        int maxTasksQueueSize = 2;
        
        int newBound = bound;
        boolean pairedKmerDistanceIsSet = false;
        int longFragmentLengthThreshold = -1;
        int shortestFragmentLengthAllowed = k + lookahead;
        
        // set up thread pool
        MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);
                
        try {
            FastqReader lin, rin;
            FastqPairReader fqpr;
            
            FastaWriter[] longFragmentsOut = new FastaWriter[]{new FastaWriter(longFragmentsFastaPaths[0], true),
                                                                new FastaWriter(longFragmentsFastaPaths[1], true),
                                                                new FastaWriter(longFragmentsFastaPaths[2], true),
                                                                new FastaWriter(longFragmentsFastaPaths[3], true),
                                                                new FastaWriter(longFragmentsFastaPaths[4], true),
                                                                new FastaWriter(longFragmentsFastaPaths[5], true)};
            
            FastaWriter[] shortFragmentsOut = new FastaWriter[]{new FastaWriter(shortFragmentsFastaPaths[0], true),
                                                                new FastaWriter(shortFragmentsFastaPaths[1], true),
                                                                new FastaWriter(shortFragmentsFastaPaths[2], true),
                                                                new FastaWriter(shortFragmentsFastaPaths[3], true),
                                                                new FastaWriter(shortFragmentsFastaPaths[4], true),
                                                                new FastaWriter(shortFragmentsFastaPaths[5], true)};
            
            FastqReadPair p;
            ArrayBlockingQueue<Fragment> fragments = new ArrayBlockingQueue<>(sampleSize);
            
//            boolean readLengthThresholdIsSet = false;
//            ArrayList<Integer> leftReadLengths = new ArrayList<>(sampleSize);
//            ArrayList<Integer> rightReadLengths = new ArrayList<>(sampleSize);
            
//            int minNumKmersNotAssembled = 1;
            
            for (FastqPair fqPair: fastqs) {
                lin = new FastqReader(fqPair.leftFastq, true);
                rin = new FastqReader(fqPair.rightFastq, true);

                fqpr = new FastqPairReader(lin, rin, qualPatternFrag, fqPair.leftRevComp, fqPair.rightRevComp);
                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");

                int leftReadLengthThreshold = k;
                int rightReadLengthThreshold = k;
                
                if (!pairedKmerDistanceIsSet) {
                                        
                    // Assembled initial sample of fragments
                    while (fqpr.hasNext()) {
                        p = fqpr.next();
                        ++readPairsParsed;

                            
//                            if (!readLengthThresholdIsSet) {
//                                int best = 0;
//                                for (String s : p.left) {
//                                    if (s.length() > best) {
//                                        best = s.length();
//                                    }
//                                }
//                                leftReadLengths.add(best);
//
//                                best = 0;
//                                for (String s : p.right) {
//                                    if (s.length() > best) {
//                                        best = s.length();
//                                    }
//                                }
//                                rightReadLengths.add(best);
//                                    
//                                if (leftReadLengths.size() >= sampleSize) {
//                                    int[] leftReadLengthsStats = getMinQ1MedianQ3Max(leftReadLengths);
//                                    int[] rightReadLengthsStats = getMinQ1MedianQ3Max(rightReadLengths);
//                                   
//                                    System.out.println("Read Lengths Distribution (n=" + sampleSize + ")");
//                                    System.out.println("  \tmin\tQ1\tM\tQ3\tmax");
//                                    System.out.println("L:\t" + leftReadLengthsStats[0] + "\t" + leftReadLengthsStats[1] + "\t" + leftReadLengthsStats[2] + "\t" + leftReadLengthsStats[3] + "\t" + leftReadLengthsStats[4]);
//                                    System.out.println("R:\t" + rightReadLengthsStats[0] + "\t" + rightReadLengthsStats[1] + "\t" + rightReadLengthsStats[2] + "\t" + rightReadLengthsStats[3] + "\t" + rightReadLengthsStats[4]);
//                                    
//                                    
//                                    leftReadLengthThreshold = Math.max(k, leftReadLengthsStats[1] - (leftReadLengthsStats[3] - leftReadLengthsStats[1]) * 3/2);
//                                    rightReadLengthThreshold = Math.max(k, rightReadLengthsStats[1] - (rightReadLengthsStats[3] - rightReadLengthsStats[1]) * 3/2);
//                                    
//                                    System.out.println("Left read length threshold:  " + leftReadLengthThreshold);
//                                    System.out.println("Right read length threshold: " + rightReadLengthThreshold);
//                                    
//                                    readLengthThresholdIsSet = true;
//                                    
//                                    leftReadLengths = null;
//                                    rightReadLengths = null;
//                                }
//                            }
                            
                            service.submit(new FragmentAssembler(p,
                                                                fragments,
                                                                bound,
                                                                minOverlap,
                                                                false,
                                                                maxErrCorrIterations, 
                                                                leftReadLengthThreshold,
                                                                rightReadLengthThreshold));
                            
                            if (fragments.remainingCapacity() == 0) {
                                break;
                            }
//                            else {
//                                System.out.println(readPairsParsed + " read pairs parsed");
//                                System.out.println(fragments.size() + " fragments assembled");
//                            }
                    }
                    
                    // Calculate length stats
                    ArrayList<Integer> fragLengths = new ArrayList<>(sampleSize);
                    for (Fragment frag : fragments) {
                        if (!frag.isUnconnectedRead) {
                            fragLengths.add(frag.length);
                        }
                    }
                   
                    int[] fragLengthsStats = getMinQ1MedianQ3Max(fragLengths);
                    System.out.println("Fragment Lengths Distribution (n=" + fragLengths.size() + ")");
                    System.out.println("\tmin\tQ1\tM\tQ3\tmax");
                    System.out.println("\t" + fragLengthsStats[0] + "\t" + fragLengthsStats[1] + "\t" + fragLengthsStats[2] + "\t" + fragLengthsStats[3] + "\t" + fragLengthsStats[4]);

                    longFragmentLengthThreshold = fragLengthsStats[1];
                    graph.setPairedKmerDistance(longFragmentLengthThreshold - k - minNumKmerPairs);
                    
                    // Set new bound for graph search
                    int iqr15 = (fragLengthsStats[3] - fragLengthsStats[1]) * 3 / 2;
                    newBound = fragLengthsStats[3] + iqr15; // 1.5*IQR
//                    shortestFragmentLengthAllowed = Math.max(k, fragLengthsStats[1] - iqr15);
                    
                    System.out.println("Paired kmers distance:       " + (longFragmentLengthThreshold - k));
                    System.out.println("Max graph traversal depth:   " + newBound);
//                    System.out.println("Shortest fragment allowed:   " + shortestFragmentLengthAllowed);

                    // Store kmer pairs and write fragments to file
                    int m;
                    Fragment frag;
                    while (!fragments.isEmpty()) {
                        frag = fragments.poll();
                        if (frag.length >= shortestFragmentLengthAllowed) {
                            m = getMinCoverageOrderOfMagnitude(frag.minCov);

                            if (m >= 0) {
                                if (frag.length >= longFragmentLengthThreshold) {
                                    graph.addPairedKmersFromSeq(frag.seq);
                                    longFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                                }
                                else {
                                    shortFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                                }
                            }
                        }
                    }                    
                    
                    pairedKmerDistanceIsSet = true;
                }
                
                // assemble the remaining fragments in multi-threaded mode
                while (fqpr.hasNext()) {
                    p = fqpr.next();
                    ++readPairsParsed;
                    
                    // ignore read pairs when more than half of raw read length were trimmed for each read
                    if (p.originalLeftLength > 2*p.numLeftBasesTrimmed &&
                            p.originalRightLength > 2*p.numRightBasesTrimmed) {
                        
                        service.submit(new FragmentAssembler(p,
                                                            fragments,
                                                            newBound, 
                                                            minOverlap, 
                                                            true,
                                                            maxErrCorrIterations, 
                                                            leftReadLengthThreshold, 
                                                            rightReadLengthThreshold));

                        if (fragments.remainingCapacity() <= numThreads) {

                            // write fragments to file
                            int m;
                            Fragment frag;
                            while (!fragments.isEmpty()) {
                                frag = fragments.poll();
                                if (frag.length >= shortestFragmentLengthAllowed) {
                                    m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                    if (m >= 0) {
                                        if (frag.length >= longFragmentLengthThreshold) {
                                            longFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                                        }
                                        else {
                                            shortFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                service.terminate();

                // write fragments to file
                int m;
                Fragment frag;
                while (!fragments.isEmpty()) {
                    frag = fragments.poll();
                    if (frag.length >= shortestFragmentLengthAllowed) {
                        m = getMinCoverageOrderOfMagnitude(frag.minCov);

                        if (m >= 0) { 
                            if (frag.length >= longFragmentLengthThreshold) {
                                longFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                            }
                            else {
                                shortFragmentsOut[m].write(Long.toString(++fragmentId), frag.seq);
                            }
                        }
                    }
                }
                                
                lin.close();
                rin.close();
            }
            
//            service.terminate();
            
//            // store kmer pairs and write fragments to file
//            int m;
//            for (Fragment frag : fragments) {
//                m = getMinCoverageOrderOfMagnitude(frag.minCov);
//                
//                if (m >= 0) {
//                    outs[m].write(Long.toString(++fragmentId), frag.seq);
//                }
//            }
            
            for (FastaWriter out : longFragmentsOut) {
                out.close();
            }
            
            for (FastaWriter out : shortFragmentsOut) {
                out.close();
            }
            
        } catch (IOException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
            
            System.out.println("Paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100 + " %");
        }        
    }
    
//    public void assembleFragments(FastqPair[] fastqs, String outDir, int mismatchesAllowed, int bound, int lookahead, int minOverlap, int maxTipLen, int sampleSize) {
//        System.out.println("Assembling fragments...");
//        
//        graph.initializePairKmersBloomFilter();
//        
//        long readPairsParsed = 0;
//        
//        ArrayList<Integer> clustersMaxContigId = new ArrayList<>(10000);
//        
//        try {
//            FastqReader lin, rin;
//            FastaWriter out;
//            ArrayList<String> sampleFragments = new ArrayList<>(sampleSize);
//            //ArrayList<Float> coverageGradients = new ArrayList<>(2*sampleSize*lookahead);
//            int[] fragmentLengths = new int[sampleSize];
//            int fid = 0;
//            
//            for (FastqPair fqPair: fastqs) {
//                lin = new FastqReader(fqPair.leftFastq, true);
//                rin = new FastqReader(fqPair.rightFastq, true);
//
//                FastqPairReader fqpr = new FastqPairReader(lin, rin, qualPatternFrag, fqPair.leftRevComp, fqPair.rightRevComp);
//                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");
//                
//                ReadPair p;
//                String rawLeft, rawRight;
//                while (fqpr.hasNext()) {
//                    p = fqpr.next();
//                    
//                    ++readPairsParsed;
//                    
//                    /*
//                    if (++readPairsParsed % NUM_PARSED_INTERVAL == 0) {
//                        System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
//                    }
//                    */
//                    
//                    rawLeft = connect(p.left, graph, k+p.numLeftBasesTrimmed+1, lookahead);
//                    rawRight = connect(p.right, graph, k+p.numRightBasesTrimmed+1, lookahead);
//                    
//                    if (okToConnectPair(rawLeft, rawRight)) {
//                        //System.out.println(">left\n" + rawLeft);
//                        //System.out.println(">right\n" + rawRight);
//                        
//                        /*
//                        if (fid < sampleSize) {
//                            if (p.left.length() > k+lookahead*2*2) {
//                                for (Float g : coverageGradients(p.left, graph, lookahead)) {
//                                    coverageGradients.add(g);
//                                }
//                            }
//                            
//                            if (p.right.length() > k+lookahead*2*2) {
//                                for (Float g : coverageGradients(p.right, graph, lookahead)) {
//                                    coverageGradients.add(g);
//                                }
//                            }
//                        }
//                        */
//                        
//                        // correct individual reads
//                        String left = correctMismatches(rawLeft, graph, lookahead, mismatchesAllowed);
//                        String right = correctMismatches(rawRight, graph, lookahead, mismatchesAllowed);
//                        
//                        if (okToConnectPair(left, right)) {
//                            String fragment = overlapThenConnect(left, right, graph, bound, lookahead, minOverlap);
//                                                        
//                            // correct fragment
//                            fragment = correctMismatches(fragment, graph, lookahead, mismatchesAllowed);
//                            
//                            //System.out.println(">fragment\n" + fragment);
//                            
//                            int fragLen = fragment.length();
//
//                            if (fragLen > k) {
//                                int backboneId = findBackboneId.apply(fragment);
//                                
//                                int cid = 0;
//                                boolean append = true;
//                                if (backboneId >= clustersMaxContigId.size()) {
//                                    clustersMaxContigId.add(0);
//                                    append = false;
//                                }
//                                else {
//                                    cid = clustersMaxContigId.get(backboneId)+1;
//                                    clustersMaxContigId.set(backboneId, cid);
//                                }
//                                
//                                //float minCov = graph.getMinKmerCoverage(fragment);
//                                
//                                /** extend on both ends unambiguously*/
//                                fragment = naiveExtend(fragment, graph, maxTipLen);
//
//                                out = new FastaWriter(outDir + File.separator + backboneId + ".fa", append);
//                                out.write(Integer.toString(cid) + " " + rawLeft + " " + rawRight, fragment);
//                                out.close();
//                                
//                                ++fid;
//                                if (fid > sampleSize) {
//                                    /** store paired kmers */
//                                    graph.addPairedKmersFromSeq(fragment);
//                                }
//                                else if (fid == sampleSize) {
//                                    fragmentLengths[0] = fragLen;
//                                    sampleFragments.add(fragment);
//
//                                    /** Calculate median fragment length */
//                                    Arrays.sort(fragmentLengths);
//                                    int half = sampleSize/2;
//                                    int medianFragLen = (fragmentLengths[half] + fragmentLengths[half - 1])/2;
//
//                                    System.out.println("Median fragment length: " + medianFragLen);
//
//                                    /** set kmer pair distance */
//                                    graph.setPairedKmerDistance(medianFragLen - k);
//
//                                    /** readjust bound to be based on 1.5*IQR */
//                                    int whisker = (fragmentLengths[sampleSize*3/4] - fragmentLengths[sampleSize/4]) * 3/2;
//                                    bound = medianFragLen + whisker;
//
//                                    System.out.println("Max graph traversal depth: " + bound);
//
//                                    /** clear sample fragment lengths */
//                                    fragmentLengths = null;
//
//                                    /** store paired kmers of all sample fragments */
//                                    for (String frag : sampleFragments) {
//                                        graph.addPairedKmersFromSeq(frag);
//                                    }
//
//                                    /** clear sample fragments */
//                                    sampleFragments = null;
//                                    
//                                    /*
//                                    Collections.sort(coverageGradients);
//                                    int cgSize = coverageGradients.size();
//                                    float cgIqr15 = (coverageGradients.get(cgSize*3/4) - coverageGradients.get(cgSize/4)) * 3/2;
//                                    float cgMedian = coverageGradients.get(cgSize/2);
//                                    System.out.println("median cov gradient: " + cgMedian + "+/-" + cgIqr15);
//                                    coverageGradients = null;
//                                    */
//                                }
//                                else {
//                                    /** store fragment length*/
//                                    fragmentLengths[fid] = fragLen;
//                                    sampleFragments.add(fragment);
//
//                                    //System.out.println(fragLen);
//                                }
//                            }
//
//                            /* assemble 3' UTR only
//                            if (p.right.endsWith("AAA")) {
//                                String fragment = assembleFragment(p.left, p.right, graph, mismatchesAllowed, bound, lookahead, minOverlap);
//                                System.out.println("LEFT:  " + p.left);
//                                System.out.println("RIGHT: " + p.right);
//                                System.out.println(fragment.length() + ": " + fragment);
//                            }
//                            */
//                        }
//                    }
//                }
//
//                lin.close();
//                rin.close();
//            }
//            
//            kmerToBackboneID = null;
//            
//        } catch (IOException ex) {
//            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
//        } finally {
//            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
//            System.out.println("Assembled fragments in " + currentBackboneId + " clusters.");
//        }
//    }
    

    public void savePairedKmersBloomFilter(File graphFile) {
        try {
            graph.savePkbf(graphFile);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void restorePairedKmersBloomFilter(File graphFile) {
        try {
            graph.destroyPkbf();
            graph.restorePkbf(graphFile);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private boolean areUnassembledKmers(ArrayList<Kmer> kmers) {
        int numKmers = kmers.size();
        int length = numKmers + k - 1;
        
        int maxMismatchesAllowed = (int) (length * percentError);
        int numMismatchBases = 0;
        
        Kmer lastGoodKmer = null;
        
        int numKmersNotSeen = 0;
        
        Kmer kmer;
        for (int i=0; i<numKmers; ++i) {
            kmer = kmers.get(i);
            
            if (screeningBf.lookup(kmer.hashVals)) {
                if (numKmersNotSeen > 0) {
                    if (lastGoodKmer == null) {
                        // left edge is not assembled
                        if (numKmersNotSeen < k) {
                            int depth = k-numKmersNotSeen;
                            Kmer edgeKmer = kmers.get(0);
                            if (hasDepthLeft(edgeKmer, graph, depth)) {
                                ArrayDeque<Kmer> extension = greedyExtendLeft(graph, edgeKmer, lookahead, depth+maxIndelSize);

                                if (!extension.isEmpty() && hasValidPath(graph,
                                                                        extension.getFirst(),
                                                                        kmer,
                                                                        screeningBf,
                                                                        k-maxIndelSize,
                                                                        k+maxIndelSize)) {
                                    numMismatchBases += 1;
                                }
                                else {
                                    numMismatchBases += numKmersNotSeen;
                                }
                            }
                        }
                        else {
                            numMismatchBases += numKmersNotSeen;
                        }
                    }
                    else {
                        // unassembled gap
                        if (numKmersNotSeen <= 2*k && hasValidPath(graph, lastGoodKmer, kmer, screeningBf, numKmersNotSeen-maxIndelSize, numKmersNotSeen+maxIndelSize)) {
                            numMismatchBases += 1;
                        }
                        else {
                            return true;
                        }
                    }
                
                    if (numMismatchBases > maxMismatchesAllowed) {
                        return true;
                    }
                }
                
                numKmersNotSeen = 0;
                lastGoodKmer = kmer;
            }
            else {
                ++numKmersNotSeen;
            }
        }
        
        if (lastGoodKmer == null) {
            return true;
        }
        
        if (numKmersNotSeen > 0 && numKmersNotSeen < k) {
            // right edge is not assembled
            int depth = k-numKmersNotSeen;
            Kmer edgeKmer = kmers.get(numKmers-1);
            if (hasDepthRight(edgeKmer, graph, depth)) {
                ArrayDeque<Kmer> extension = greedyExtendRight(graph, edgeKmer, lookahead, depth+maxIndelSize);

                if (!extension.isEmpty() && hasValidPath(graph,
                                                        lastGoodKmer,
                                                        extension.getLast(),
                                                        screeningBf,
                                                        k-maxIndelSize,
                                                        k+maxIndelSize)) {
                    numMismatchBases += 1;
                }
                else {
                    numMismatchBases += numKmersNotSeen;
                }
            }
        }
        else {
            numMismatchBases += numKmersNotSeen;
        }
        
        return numMismatchBases > maxMismatchesAllowed;
    }
    
    public void assembleTranscriptsMultiThreaded(String[] longFragmentsFastas, 
                                                String[] shortFragmentsFastas,
                                                String outFasta,
                                                String outFastaShort,
                                                long sbfNumBits, 
                                                int sbfNumHash,
                                                int numThreads,
                                                int sampleSize,
                                                int minTranscriptLength) {
        
        long numFragmentsParsed = 0;

        try {
            System.out.println("Creating graph from fragment kmers...");
            graph.getDbgbf().empty();
            insertIntoDeBruijnGraph(longFragmentsFastas);
            insertIntoDeBruijnGraph(shortFragmentsFastas);
        
            dbgFPR = graph.getDbgbf().getFPR();
            System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
            
            if (covFPR <= 0) {
                covFPR = graph.getCbf().getFPR();
            }

            System.out.println("Assembling transcripts...");

            boolean append = false;

            // set up thread pool
            int maxTasksQueueSize = 2;
        

            FastaWriter fout = new FastaWriter(outFasta, append);
            FastaWriter foutShort = new FastaWriter(outFastaShort, append);
            
            screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());

            long cid = 0;

            FastaReader fin;
            
            String tag = ".L.";
            
            for (String[] fragmentsFastas : new String[][]{longFragmentsFastas, shortFragmentsFastas}) {
                for (int mag=longFragmentsFastas.length-1; mag>=0; --mag) {
                    boolean beGreedy = false;

                    String prefix = "E" + mag + tag;

                    String fragmentsFasta = fragmentsFastas[mag];

                    System.out.println("Parsing fragments in `" + fragmentsFasta + "`...");

                    fin = new FastaReader(fragmentsFasta);

                    MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);
                    AssembledTranscriptsQueue transcripts = new AssembledTranscriptsQueue(sampleSize);
                    
                    String fragment;
                    Transcript t;
//                    float minCoverageThreshold = (float) Math.pow(10, mag);
                    while (fin.hasNext()) {
                        ++numFragmentsParsed;
//                        if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
//                            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
//                        }

                        fragment = fin.next();

                        int numFragKmers = getNumKmers(fragment, k);
                        
                        if (numFragKmers > 0 && !isHomoPolymer(fragment)) {
                            
                            service.submit(new TranscriptAssembler(fragment,
                                    transcripts,
                                    beGreedy
                            ));
                            
                            if (transcripts.remainingCapacity() <= numThreads) {

                                // write fragments to file
                                int len;
                                while (!transcripts.isEmpty()) {
                                    t = transcripts.poll();
                                    len = t.transcript.length();
                                    
                                    if (len >= minTranscriptLength) {
                                        fout.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + t.fragment + "]", t.transcript);
                                    }
                                    else {
                                        foutShort.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + t.fragment + "]", t.transcript);
                                    }
                                }
                            }
                        }
                    }

                    fin.close();
                    
                    service.terminate();

                    // write fragments to file
                    int len;
                    while (!transcripts.isEmpty()) {
                        t = transcripts.poll();
                        len = t.transcript.length();

                        if (len >= minTranscriptLength) {
                            fout.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + t.fragment + "]", t.transcript);
                        }
                        else {
                            foutShort.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + t.fragment + "]", t.transcript);
                        }
                    }
                }
                
                tag = ".S.";
            }
            
            fout.close();
            foutShort.close();
            
            screeningBf.destroy();
        } catch (IOException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
        }
    }
    
//    public void assembleTranscripts(String[] longFragmentsFastas, String[] shortFragmentsFastas, String outFasta, int lookAhead, int maxTipLength, float maxCovGradient, long sbfNumBits, int sbfNumHash) {
//        if (dbgFPR <= 0) {
//            dbgFPR = graph.getDbgbf().getFPR();
//        }
//        
//        if (covFPR <= 0) {
//            covFPR = graph.getCbf().getFPR();
//        }
//        
//        int maxIndelSize = 1;
//        float percentIdentity = 0.95f;
//        float percentError = 1-percentIdentity;
//        
//        System.out.println("Assembling transcripts...");
//        long numFragmentsParsed = 0;
//        boolean append = false;
//        
//        try {
//            FastaWriter fout = new FastaWriter(outFasta, append);
//            
//            screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
//
//            long cid = 0;
//
//            FastaReader fin;
//            
//            for (int mag=longFragmentsFastas.length-1; mag>=0; --mag) {
//                boolean beGreedy = mag > 0;
//                
//                String prefix = "E" + mag + ".L.";
//                
//                for (String fragmentsFasta : new String[]{longFragmentsFastas[mag], shortFragmentsFastas[mag]}) {
//                    
//                    System.out.println("Parsing fragments in `" + fragmentsFasta + "`...");
//
//                    fin = new FastaReader(fragmentsFasta);
//
//                    String fragment;
////                    float minCoverageThreshold = (float) Math.pow(10, mag);
//
//                    while (fin.hasNext()) {
//                        ++numFragmentsParsed;
////                        if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
////                            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
////                        }
//
//                        fragment = fin.next();
//
//                        int numFragKmers = getNumKmers(fragment, k);
//                        
//                        if (numFragKmers > 0 && !isHomoPolymer(fragment)) {
//                            
//                            ArrayList<Kmer> fragKmers = graph.getKmers(correctMismatches(fragment, graph, lookAhead, (int) Math.ceil(fragment.length()*percentError)));
//                            
//                            if (hasNotYetAssembled(fragKmers, maxIndelSize, percentIdentity)) {
//
//                                /** check whether sequence-wide coverage differences are too large */
//                                int numFalsePositivesAllowed = (int) Math.round(numFragKmers * covFPR);
//                                float[] covs = new float[numFragKmers];
//                                for (int i=0; i<numFragKmers; ++i) {
//                                    covs[i] = fragKmers.get(i).count;
//                                }
//
//                                boolean covDiffTooLarge = false;
//                                Arrays.sort(covs);
//                                float covLow = covs[0];
//                                float covHigh;
//                                for (int i=1; i<numFragKmers-numFalsePositivesAllowed; ++i) {
//                                    covHigh = covs[i];
//                                    if (covHigh * maxCovGradient > covLow) {
//                                        covDiffTooLarge = true;
//                                        break;
//                                    }
//                                    covLow = covHigh;
//                                }
//
//                                extendWithPairedKmers(fragKmers, graph, lookAhead, maxTipLength, !covDiffTooLarge && beGreedy, screeningBf);
//
//                                if (hasNotYetAssembled(fragKmers, maxIndelSize, percentIdentity)) {
//
//                                    for (Kmer kmer : fragKmers) {
//                                        screeningBf.add(kmer.hashVals);
//                                    }
//
//                                    String transcript = assemble(fragKmers, k);
//                                    
//                                    fout.write(prefix +  Long.toString(++cid) + " l=" + transcript.length() + " F=[" + fragment + "]", transcript);
//                                }
//                            }
//                        }
//                    }
//
//                    fin.close();
//                    
//                    prefix = "E" + mag + ".S.";
//                }
//            }
//            
//            fout.close();
//            
//            screeningBf.destroy();
//        } catch (IOException ex) {
//            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
//        } finally {
//            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
//        }
//    }
        
    public void printFloatArray(float[] arr){
        StringBuilder sb = new StringBuilder();
        sb.append("[ ");
        for (float f : arr) {
            sb.append(f);
            sb.append(" ");
        }
        sb.append("]");
        System.out.println(sb.toString());
    }
    
    public void testErrorCorrection() {
        int lookahead = 5;
        float maxCovGradient = 0.5f;
        int bound = 500; 
        int maxIndelSize = 1; 
        float covFPR = graph.getCbf().getFPR();
        int minOverlap = 10;
        int errorCorrectionIterations = 2;

        String leftRead =  "GGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCTTCGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTTGTCAAGCTCATTTCCTGG";
        String rightRead = "AGCAAGAGCACAAGAGGAAGAGGGAGCCCCTCCCTGCTGGGGAGTCCCTGCCACACTAAGTCCCCCACCACACTGAATCTCCC";
        
//        String leftRead =  "CAGTCAGCCGCATCTTCTTTTGCGTCGCCAGCCGAGCCACATCGCTCAGACACCATGGGGAAGGTGAAGGTCGGAGTCAACGGGTGAGTTCGCGGGTGGC";
//        String rightRead = "CTGGGGGGCCCTGGGCTGCGACCGCCCCCGAACCGCGTCTACGAGCCTTGCGGGCTCCGGGTCTTTGCAGTCGTATGGGGGCAGGGTAGCTGTTCCCCGC";
        
        System.out.println(leftRead);
        System.out.println(rightRead);
        
    }
    
    public void testInsertionCorrection() {
        String seq =  "TGAAGCAGGCGTCGGAGGGCCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGAC";
        String note = "                        ^                                                           ";
        
        int lookahead = 5;
        int errorsallowed = 5;
        
        String corrected = correctErrors(seq, graph, lookahead, errorsallowed);
        
        System.out.println(seq);
        System.out.println(note);
        System.out.println(corrected);
    }

    
    public void testMismatchCorrection() {
        String seq =  "AGAGGGGTGATGTGGGGAGTACGCTGCAGGGCCTCACTCCTTTTGCAGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCTTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAAC";
        
        int lookahead = 7;
        float pid = 0.95f;
        
        String corrected = correctMismatches(seq, graph, lookahead, (int) Math.ceil((1.0f - pid) * seq.length()));

        
        System.out.println(seq);
        System.out.println(corrected);
    }
    
    public static void touch(File f) throws IOException {
        f.getParentFile().mkdirs();
        if (!f.createNewFile()){
            f.setLastModified(System.currentTimeMillis());
        }
    }
    
    public static void clearDirectory(File dir) {
        for(File file: dir.listFiles()) 
            if (!file.isDirectory()) 
                file.delete();
    }    

    public static void printHelp(Options options, boolean error) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(null);
        formatter.printHelp( "java -jar RNA-Bloom.jar", options, true);
        
        if (error) {
            System.exit(1);
        }
        else {
            System.exit(0);
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        final String STARTED = "STARTED";
        final String DBG_DONE = "DBG.DONE";
        final String FRAGMENTS_DONE = "FRAGMENTS.DONE";
        final String TRANSCRIPTS_DONE = "TRANSCRIPTS.DONE";
        
        long globalStartTime = System.nanoTime();
        
        System.out.println("args: " + Arrays.toString(args));
        
        // -left /home/gengar/test_data/GAPDH/GAPDH_2.fq.gz -right /home/gengar/test_data/GAPDH/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /home/gengar/test_assemblies/GAPDH
        // -dm 1 -cm 2.5 -pm 0.5 -left /home/gengar/test_data/SRR1360926/SRR1360926_2.fastq.gz -right /home/gengar/test_data/SRR1360926/SRR1360926_1.fastq.gz -revcomp-right -stranded -name SRR1360926 -outdir /home/gengar/test_assemblies/SRR1360926

        
        // -left /projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/gapdh
        // -dm 1 -cm 2.5 -pm 0.5 -left /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_2.fq -right /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_1.fq -revcomp-right -stranded -name SRR1360926 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/SRR1360926
        // -dm 5 -cm 20 -pm 3 -sm 2 -left /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_2.fq.gz -right /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_1.fq.gz -revcomp-right -stranded -name mcf7 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/mcf7        
        // Based on: http://commons.apache.org/proper/commons-cli/usage.html
        CommandLineParser parser = new DefaultParser();

        Options options = new Options();

        Builder builder;
                
        builder = Option.builder("l");
        builder.longOpt("left");
        builder.desc("left reads file");
        builder.hasArg(true);
        builder.argName("FILE");
//        builder.required(true);
        Option optLeftReads = builder.build();
        options.addOption(optLeftReads);
        
        builder = Option.builder("r");
        builder.longOpt("right");
        builder.desc("right reads file");
        builder.hasArg(true);
        builder.argName("FILE");
//        builder.required(true);
        Option optRightReads = builder.build();
        options.addOption(optRightReads);
        
        builder = Option.builder("rcl");
        builder.longOpt("revcomp-left");
        builder.desc("reverse-complement left reads");
        builder.hasArg(false);
        Option optRevCompLeft = builder.build();
        options.addOption(optRevCompLeft);

        builder = Option.builder("rcr");
        builder.longOpt("revcomp-right");
        builder.desc("reverse-complement right reads");
        builder.hasArg(false);
        Option optRevCompRight = builder.build();
        options.addOption(optRevCompRight);
        
        builder = Option.builder("ss");
        builder.longOpt("stranded");
        builder.desc("strand specific");
        builder.hasArg(false);
        Option optStranded = builder.build();
        options.addOption(optStranded);
        
        builder = Option.builder("n");
        builder.longOpt("name");
        builder.desc("assembly name");
        builder.hasArg(true);
        builder.argName("STR");
        //builder.required(true);
        Option optName = builder.build();
        options.addOption(optName);
        
        builder = Option.builder("t");
        builder.longOpt("threads");
        builder.desc("run in INT threads");
        builder.hasArg(true);
        builder.argName("INT");
        Option optThreads = builder.build();
        options.addOption(optThreads);
        
        builder = Option.builder("o");
        builder.longOpt("outdir");
        builder.desc("output directory");
        builder.hasArg(true);
        builder.argName("PATH");
        //builder.required(true);
        Option optOutdir = builder.build();
        options.addOption(optOutdir);
        
        builder = Option.builder("f");
        builder.longOpt("force");
        builder.desc("force overwrite existing files");
        builder.hasArg(false);
        Option optForce = builder.build();
        options.addOption(optForce);
        
        builder = Option.builder("k");
        builder.longOpt("kmer");
        builder.desc("kmer size");
        builder.hasArg(true);
        builder.argName("INT");
        Option optKmerSize = builder.build();
        options.addOption(optKmerSize);
        
        builder = Option.builder("q");
        builder.longOpt("qual-dbg");
        builder.desc("min base quality for constructing DBG");
        builder.hasArg(true);
        builder.argName("INT");
        Option optBaseQualDbg = builder.build();
        options.addOption(optBaseQualDbg);

        builder = Option.builder("Q");
        builder.longOpt("qual-frag");
        builder.desc("min base quality for fragment assembly");
        builder.hasArg(true);
        builder.argName("INT");
        Option optBaseQualFrag = builder.build();
        options.addOption(optBaseQualFrag);        
                
        builder = Option.builder("sh");
        builder.longOpt("sbf-hash");
        builder.desc("number of hash functions for screening Bloom filter");
        builder.hasArg(true);
        builder.argName("INT");
        Option optSbfHash = builder.build();
        options.addOption(optSbfHash); 
        
        builder = Option.builder("dh");
        builder.longOpt("dbgbf-hash");
        builder.desc("number of hash functions for de Bruijn graph Bloom filter");
        builder.hasArg(true);
        builder.argName("INT");
        Option optDbgbfHash = builder.build();
        options.addOption(optDbgbfHash);

        builder = Option.builder("ch");
        builder.longOpt("cbf-hash");
        builder.desc("number of hash functions for kmer counting Bloom filter");
        builder.hasArg(true);
        builder.argName("INT");
        Option optCbfHash = builder.build();
        options.addOption(optCbfHash);
        
        builder = Option.builder("ph");
        builder.longOpt("pkbf-hash");
        builder.desc("number of hash functions for paired kmers Bloom filter");
        builder.hasArg(true);
        builder.argName("INT");
        Option optPkbfHash = builder.build();
        options.addOption(optPkbfHash);        

        builder = Option.builder("sm");
        builder.longOpt("sbf-mem");
        builder.desc("allocate DECIMAL-gigabyte for screening Bloom filter");
        builder.hasArg(true);
        builder.argName("DECIMAL");
        Option optSbfMem = builder.build();
        options.addOption(optSbfMem);
        
        builder = Option.builder("dm");
        builder.longOpt("dbgbf-mem");
        builder.desc("allocate DECIMAL-gigabyte for de Bruijn graph Bloom filter");
        builder.hasArg(true);
        builder.argName("DECIMAL");
        Option optDbgbfMem = builder.build();
        options.addOption(optDbgbfMem);

        builder = Option.builder("cm");
        builder.longOpt("cbf-mem");
        builder.desc("allocate DECIMAL-gigabyte for kmer counting Bloom filter");
        builder.hasArg(true);
        builder.argName("DECIMAL");
        Option optCbfMem = builder.build();
        options.addOption(optCbfMem);
        
        builder = Option.builder("pm");
        builder.longOpt("pkbf-mem");
        builder.desc("allocate DECIMAL-gigabyte for paired kmers Bloom filter");
        builder.hasArg(true);
        builder.argName("DECIMAL");
        Option optPkbfMem = builder.build();
        options.addOption(optPkbfMem);
        
//        builder = Option.builder("fpr");
//        builder.longOpt("max-fpr");
//        builder.desc("maximum FPR (%) allowed for Bloom filters");
//        builder.hasArg(true);
//        builder.argName("DECIMAL");
//        Option optMaxFPR = builder.build();
//        options.addOption(optMaxFPR);
                
//        builder = Option.builder("m");
//        builder.longOpt("mismatch");
//        builder.desc("max number of mismatch bases allowed per read");
//        builder.hasArg(true);
//        builder.argName("INT");
//        Option optMismatch = builder.build();
//        options.addOption(optMismatch);

        builder = Option.builder("tl");
        builder.longOpt("tiplength");
        builder.desc("max tip length allowed");
        builder.hasArg(true);
        builder.argName("INT");
        Option optTipLength = builder.build();
        options.addOption(optTipLength);  
        
        builder = Option.builder("la");
        builder.longOpt("lookahead");
        builder.desc("number of kmers to look ahead during graph traversal");
        builder.hasArg(true);
        builder.argName("INT");
        Option optLookahead = builder.build();
        options.addOption(optLookahead);        
        
        builder = Option.builder("ol");
        builder.longOpt("overlap");
        builder.desc("min number of overlapping bases between mates");
        builder.hasArg(true);
        builder.argName("INT");
        Option optOverlap = builder.build();
        options.addOption(optOverlap);
        
        builder = Option.builder("b");
        builder.longOpt("bound");
        builder.desc("max distance between mates");
        builder.hasArg(true);
        builder.argName("INT");
        Option optBound = builder.build();
        options.addOption(optBound);

        builder = Option.builder("S");
        builder.longOpt("sample");
        builder.desc("sample size for estimating median fragment length");
        builder.hasArg(true);
        builder.argName("INT");
        Option optSample = builder.build();
        options.addOption(optSample);
        
        builder = Option.builder("mcg");
        builder.longOpt("maxcovgrad");
        builder.desc("max coverage gradient for error correction");
        builder.hasArg(true);
        builder.argName("DECIMAL");
        Option optMaxCovGrad = builder.build();
        options.addOption(optMaxCovGrad);
        
        builder = Option.builder("i");
        builder.longOpt("indel");
        builder.desc("maximum indel size allowed");
        builder.hasArg(true);
        builder.argName("INT");
        Option optIndelSize = builder.build();
        options.addOption(optIndelSize);  

        builder = Option.builder("p");
        builder.longOpt("percent");
        builder.desc("percent identity allowed");
        builder.hasArg(true);
        builder.argName("FLOAT");
        Option optPercentIdentity = builder.build();
        options.addOption(optPercentIdentity); 
        
        builder = Option.builder("e");
        builder.longOpt("errcorritr");
        builder.desc("max number of iterations of read error correction");
        builder.hasArg(true);
        builder.argName("INT");
        Option optErrCorrItr = builder.build();
        options.addOption(optErrCorrItr);        

        builder = Option.builder("len");
        builder.longOpt("length");
        builder.desc("min transcript length in final assembly");
        builder.hasArg(true);
        builder.argName("INT");
        Option optMinLength = builder.build();
        options.addOption(optMinLength);  
        
        builder = Option.builder("h");
        builder.longOpt("help");
        builder.desc("print this message and exits");
        Option optHelp = builder.build();
        options.addOption(optHelp);
        

        try {
            CommandLine line = parser.parse(options, args);

            if (line.getOptions().length == 0 || line.hasOption(optHelp.getOpt())) {
                printHelp(options, false);
            }
            
            int numThreads = Integer.parseInt(line.getOptionValue(optThreads.getOpt(), "2"));
            boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            String name = line.getOptionValue(optName.getOpt(), "rna-bloom");
            String outdir = line.getOptionValue(optOutdir.getOpt(), System.getProperty("user.dir") + File.separator + name + "_assembly");
            /**@TODO evaluate whether out dir is a valid dir */
            
            String longFragmentsFastaPrefix = outdir + File.separator + name + ".fragments.long.";
            String shortFragmentsFastaPrefix = outdir + File.separator + name + ".fragments.short.";
            String transcriptsFasta = outdir + File.separator + name + ".transcripts.fa";
            String shortTranscriptsFasta = outdir + File.separator + name + ".transcripts.short.fa";
//            String tmpFasta = outdir + File.separator + name + ".tmp.fa";
            String graphFile = outdir + File.separator + name + ".graph";
            
            File startedStamp = new File(outdir + File.separator + STARTED);
            File dbgDoneStamp = new File(outdir + File.separator + DBG_DONE);
            File fragsDoneStamp = new File(outdir + File.separator + FRAGMENTS_DONE);
            File txptsDoneStamp = new File(outdir + File.separator + TRANSCRIPTS_DONE);
            
            if (forceOverwrite) {
                if (startedStamp.exists()) {
                    startedStamp.delete();
                }
                
                if (dbgDoneStamp.exists()) {
                    dbgDoneStamp.delete();
                }
                
                if (fragsDoneStamp.exists()) {
                    fragsDoneStamp.delete();
                }
                
                if (txptsDoneStamp.exists()) {
                    txptsDoneStamp.delete();
                }
            }
            
            String fastqLeft = line.getOptionValue(optLeftReads.getOpt());
            String fastqRight = line.getOptionValue(optRightReads.getOpt());
            
            boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            boolean strandSpecific = line.hasOption(optStranded.getOpt());
            
            int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), "25"));
            int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), "3"));
            int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), "10"));
            
            long sbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), "1")));
            long dbgbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), "1")));
            long cbfSize = (long) (NUM_BYTES_1GB * Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), "1")));
            long pkbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), "1")));
            
            int sbfNumHash = Integer.parseInt(line.getOptionValue(optSbfHash.getOpt(), "2"));
            int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), "2"));
            int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), "2"));
            int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), "2"));
            
            /**@TODO ensure that sbfNumHash and pkbfNumHash <= max(dbgbfNumHash, cbfNumHash) */
                        
//            int mismatchesAllowed = Integer.parseInt(line.getOptionValue(optMismatch.getOpt(), "5"));
            int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), "10"));
            int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), "1000"));
            int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), "500"));
            int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), "5"));
            int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), "3"));
            float maxCovGradient = Float.parseFloat(line.getOptionValue(optMaxCovGrad.getOpt(), "0.5"));
            float percentIdentity = Float.parseFloat(line.getOptionValue(optPercentIdentity.getOpt(), "0.9"));
            int maxIndelSize = Integer.parseInt(line.getOptionValue(optIndelSize.getOpt(), "1"));
            int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), "2"));
            int minTranscriptLength = Integer.parseInt(line.getOptionValue(optMinLength.getOpt(), "200"));
            
            boolean saveGraph = true;
            boolean saveKmerPairs = true;

            System.out.println("name: " + name);
            System.out.println("outdir: " + outdir);
            
            File f = new File(outdir);
            if (!f.exists()) {
                f.mkdirs();
            }

            RNABloom assembler = new RNABloom(k, qDBG, qFrag);
            assembler.setParams(maxTipLen, lookahead, maxCovGradient, maxIndelSize, percentIdentity);

            try {
                touch(startedStamp);
            } catch (IOException ex) {
                Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            if (!forceOverwrite && dbgDoneStamp.exists()) {
                System.out.println("Loading graph from file `" + graphFile + "`...");
                assembler.restoreGraph(new File(graphFile));
            }
            else {
                String[] forwardFastqs = new String[]{fastqLeft};
                String[] backwardFastqs = new String[]{fastqRight};
                
                long startTime = System.nanoTime();
                
                assembler.createGraph(forwardFastqs, backwardFastqs, 
                        strandSpecific, 
                        dbgbfSize, cbfSize, pkbfSize, 
                        dbgbfNumHash, cbfNumHash, pkbfNumHash,
                        numThreads);

                System.out.println("Time elapsed: " + (System.nanoTime() - startTime) / Math.pow(10, 9) + " seconds");
                
                if (saveGraph) {
                    System.out.println("Saving graph to file `" + graphFile + "`...");
                    assembler.saveGraph(new File(graphFile));
                }
                
                try {
                    touch(dbgDoneStamp);
                } catch (IOException ex) {
                    Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
//            /**@TODO */
//            assembler.testMismatchCorrection();
//            System.exit(0);
            
            FastqPair fqPair = new FastqPair(fastqLeft, fastqRight, revCompLeft, revCompRight);
            FastqPair[] fqPairs = new FastqPair[]{fqPair};        

            String[] longFragmentsFastaPaths = {longFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
            
            String[] shortFragmentsFastaPaths = {shortFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
            
            if (!forceOverwrite && fragsDoneStamp.exists()) {
                System.out.println("Restoring paired kmers Bloom filter from file...");
                assembler.restorePairedKmersBloomFilter(new File(graphFile));            
            }
            else {
                for (String fragmentsFasta : longFragmentsFastaPaths) {
                    File fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }
                
                for (String fragmentsFasta : shortFragmentsFastaPaths) {
                    File fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }
                
                long startTime = System.nanoTime();
                
                assembler.assembleFragmentsMultiThreaded(fqPairs, 
                        longFragmentsFastaPaths, 
                        shortFragmentsFastaPaths, 
                        bound, 
                        minOverlap,
                        sampleSize,
                        numThreads,
                        maxErrCorrItr);

                System.out.println("Time elapsed: " + (System.nanoTime() - startTime) / Math.pow(10, 9) + " seconds");
                
                if (saveKmerPairs) {
                    System.out.println("Saving paired kmers Bloom filter to file...");
                    assembler.savePairedKmersBloomFilter(new File(graphFile));            
                }
                
                try {
                    touch(fragsDoneStamp);
                } catch (IOException ex) {
                    Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            if (forceOverwrite || !txptsDoneStamp.exists()) {
                
                File transcriptsFile = new File(transcriptsFasta);
                if (transcriptsFile.exists()) {
                    transcriptsFile.delete();
                }
                
                File shortTranscriptsFile = new File(shortTranscriptsFasta);
                if (shortTranscriptsFile.exists()) {
                    shortTranscriptsFile.delete();
                }
                
                long startTime = System.nanoTime();
                
                assembler.assembleTranscriptsMultiThreaded(longFragmentsFastaPaths, 
                                                            shortFragmentsFastaPaths, 
                                                            transcriptsFasta, 
                                                            shortTranscriptsFasta,
                                                            sbfSize,
                                                            sbfNumHash,
                                                            numThreads,
                                                            sampleSize,
                                                            minTranscriptLength);

                System.out.println("Transcripts assembled in `" + transcriptsFasta + "`");
                System.out.println("Time elapsed: " + (System.nanoTime() - startTime) / Math.pow(10, 9) + " seconds");
                
                try {
                    touch(txptsDoneStamp);
                } catch (IOException ex) {
                    Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
            
        }
        catch (ParseException exp) {
            System.out.println("ERROR:" + exp.getMessage() );
        }
        
        System.out.println("Total Runtime: " + (System.nanoTime() - globalStartTime) / Math.pow(10, 9) + " seconds");
    }
}
