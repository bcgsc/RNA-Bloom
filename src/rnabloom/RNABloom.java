/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.pow;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.CanonicalNTHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.bloom.hash.PairedNTHashIterator;
import rnabloom.bloom.hash.ReverseComplementNTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer2;
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
    
//    private final static long NUM_PARSED_INTERVAL = 100000;
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    
    private int k;
//    private boolean strandSpecific;
    private Pattern seqPattern;
    private Pattern qualPatternDBG;
    private Pattern qualPatternFrag;
//    private Pattern homoPolymerKmerPattern;
    private BloomFilterDeBruijnGraph graph = null;
    private BloomFilter screeningBf = null;

    private int maxTipLength;
    private int lookahead;
    private float maxCovGradient;
    private int maxIndelSize;
    private float percentIdentity;
    private float percentError;
    private int minNumKmerPairs;
    private int longFragmentLengthThreshold = -1;
    
    private int qDBG = -1;
    private int qFrag = -1;
    
    private float dbgFPR = -1;
    private float covFPR = -1;
    private final static String[] COVERAGE_ORDER = {"e0", "e1", "e2", "e3", "e4", "e5"};
    
    public RNABloom(int k, int qDBG, int qFrag) {
        this.qDBG = qDBG;
        this.qFrag = qFrag;
        
        this.setK(k);
    }
    
    public final void setK(int k) {
        this.k = k;
        
        this.seqPattern = getNucleotideCharsPattern(k);
        this.qualPatternDBG = getPhred33Pattern(qDBG, k);
        this.qualPatternFrag = getPhred33Pattern(qFrag, k);
//        this.homoPolymerKmerPattern = getHomoPolymerPattern(k);

        if (graph != null) {
            graph.setK(k);
        }
    }
    
    public void setParams(int maxTipLength, 
            int lookahead, 
            float maxCovGradient, 
            int maxIndelSize, 
            float percentIdentity, 
            int minNumKmerPairs) {
        
        this.maxTipLength = maxTipLength;
        this.lookahead = lookahead;
        this.maxCovGradient = maxCovGradient;
        this.maxIndelSize = maxIndelSize;
        this.percentIdentity = percentIdentity;
        this.percentError = 1.0f - percentIdentity;
        this.minNumKmerPairs = minNumKmerPairs;
    }
    
    private void handleException(Exception ex) {
        System.out.println("ERROR:" + ex.getMessage() );
        ex.printStackTrace();
        System.exit(1);
    }
    
    public void saveGraph(File f) {
        try {
            graph.save(f);
        } catch (Exception ex) {
            handleException(ex);
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
        } catch (Exception ex) {
            handleException(ex);
        }
    }
    
    public void clearGraph() {
        graph.clearDbgbf();
        graph.clearPkbf();
        
        dbgFPR = 0;
        covFPR = 0;
    }

    public void restoreDbg(File f) {
        try {
            dbgFPR = graph.getDbgbfFPR();
        } catch (Exception ex) {
            handleException(ex);
        }
    }
    
//    public int addKmersFromFastq(String fastq, boolean stranded, boolean reverseComplement, int numHash) throws IOException {
//        int numReads = 0;
//        
//        NTHashIterator itr;
//        
//        if (stranded) {
//            if (reverseComplement) {
//                itr = new ReverseComplementNTHashIterator(k, numHash);
//            }
//            else {
//                itr = new NTHashIterator(k, numHash);
//            }
//        }
//        else {
//            itr = new CanonicalNTHashIterator(k, numHash);
//        }
//        
//        FastqReader fr = new FastqReader(fastq, false);
//        FastqRecord record = fr.record;
//        Matcher m = qualPatternDBG.matcher("");
//        long[] hashVals = itr.hVals;
//
//        while (fr.hasNext()) {
//            ++numReads;
//
//            fr.nextWithoutNameFunction();
//            m.reset(record.qual);
//
//            while (m.find()) {
//                itr.start(record.seq.substring(m.start(), m.end()));
//                while (itr.hasNext()) {
//                    itr.next();
//                    if (screeningBf.lookupThenAdd(hashVals)) {
//                        graph.add(hashVals);
//                    }
//                }
//            }
//        }
//        fr.close();
//        
//        return numReads;
//    }

    public class SeqToGraphWorker implements Runnable {
//        public final static int FASTA_FILE = 1;
//        public final static int FASTQ_FILE = 2;
        
        private final int id;
        private final String path;
        private final NTHashIterator itr;
        private int numReads = 0;
        private boolean successful = false;
        private boolean insertIfAbsent = false;
        
        public SeqToGraphWorker(int id, String path, boolean stranded, boolean reverseComplement, int numHash) {
            /*@TODO: Support FASTA files*/
            
            this.id = id;
            this.path = path;
            this.insertIfAbsent = insertIfAbsent;
            
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
            System.out.println("[" + id + "] Parsing `" + path + "`...");
            
            try {
                FastqReader fr = new FastqReader(path, false);
                FastqRecord record = fr.record;
                Matcher mQual = qualPatternDBG.matcher("");
                Matcher mSeq = seqPattern.matcher("");
                long[] hashVals = itr.hVals;
                
                if (insertIfAbsent) {
                    while (fr.hasNext()) {
                        ++numReads;

                        fr.nextWithoutNameFunction();
                        mQual.reset(record.qual);

                        while (mQual.find()) {
                            mSeq.reset(record.seq.substring(mQual.start(), mQual.end()));
                            while (mSeq.find()) {
                                itr.start(mSeq.group());
                                while (itr.hasNext()) {
                                    itr.next();
                                    graph.addIfAbsent(hashVals); // Only insert if kmer is absent in the graph
                                }
                            }
                        }
                    }
                }
                else {
                    while (fr.hasNext()) {
                        ++numReads;

                        fr.nextWithoutNameFunction();
                        mQual.reset(record.qual);

                        while (mQual.find()) {
                            mSeq.reset(record.seq.substring(mQual.start(), mQual.end()));
                            while (mSeq.find()) {
                                itr.start(mSeq.group());
                                while (itr.hasNext()) {
                                    itr.next();
                                    graph.add(hashVals);
                                }
                            }
                        }
                    }
                }
                
                fr.close();
                successful = true;
                System.out.println("[" + id + "] Parsed " + NumberFormat.getInstance().format(numReads) + " sequences.");
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        
        public boolean isSuccessful() {
            return successful;
        }
        
        public long getReadCount() {
            return numReads;
        }
    }
    
    public void initializeGraph(boolean strandSpecific,
                            long dbgbfNumBits,
                            long cbfNumBytes,
                            long pkbfNumBits,
                            int dbgbfNumHash,
                            int cbfNumHash,
                            int pkbfNumHash) {
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            k,
                                            strandSpecific);
    }
    
    public void populateGraph(Collection<String> forwardFastqs,
                            Collection<String> reverseFastqs,
                            boolean strandSpecific,
                            int numThreads) {        
        
//        screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
        
        /** parse the reads */
        
        int numReads = 0;
        int numHash = graph.getMaxNumHash();
        
        ExecutorService service = Executors.newFixedThreadPool(numThreads);
        
        ArrayList<SeqToGraphWorker> threadPool = new ArrayList<>();
        int threadId = 0;
           
        for (String fastq : forwardFastqs) {
            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, fastq, strandSpecific, false, numHash);
            service.submit(t);
            threadPool.add(t);
        }

        for (String fastq : reverseFastqs) {
            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, fastq, strandSpecific, true, numHash);
            service.submit(t);
            threadPool.add(t);
        }

        service.shutdown();
        
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            
            for (SeqToGraphWorker t : threadPool) {
                numReads += t.getReadCount();
            }

            System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");
//            System.out.println("Screening Bloom filter FPR:  " + screeningBf.getFPR() * 100 + " %");
            
//            screeningBf.destroy();
            
        } catch (Exception ex) {
            handleException(ex);
        }
        
        dbgFPR = graph.getDbgbfFPR();
        covFPR = graph.getCbfFPR();
    }

    public void repopulateGraph(Collection<String> fastas, boolean strandSpecific) {
        /** insert into graph if absent */
        
        /** parse the reads */
                          
        NTHashIterator itr = graph.getHashIterator(graph.getMaxNumHash());
        long[] hashVals = itr.hVals;

        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;

        try {
            for (String path : fastas) {
                FastaReader fin = new FastaReader(path);

                String seq;
                while (fin.hasNext()) {
                    seq = fin.next();

                    itr.start(seq);
                    while (itr.hasNext()) {
                        itr.next();
                        graph.addIfAbsent(hashVals);
                    }

                    pItr.start(seq);
                    while (pItr.hasNext()) {
                        pItr.next();
                        graph.addPairedKmers(hashVals1, hashVals2, hashVals3);
                    }
                }

                fin.close();
            }
        } catch (Exception ex) {
            handleException(ex);
        }
                
        dbgFPR = graph.getDbgbfFPR();
        covFPR = graph.getCbfFPR();
    }
    
    public void insertIntoDeBruijnGraph(String fasta) throws IOException { 
        
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;
        
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

    public void insertIntoDeBruijnGraphAndPairedKmers(String fasta) throws IOException {
                
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;

        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;
        
        FastaReader fin = new FastaReader(fasta);
        
        String seq;
        while (fin.hasNext()) {
            seq = fin.next();
            
            itr.start(seq);
            while (itr.hasNext()) {
                itr.next();
                graph.addDbgOnly(hashVals);
            }
            
            pItr.start(seq);
            while (pItr.hasNext()) {
                pItr.next();
                graph.addPairedKmers(hashVals1, hashVals2, hashVals3);
            }
        }

        fin.close();
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
    
    private void insertIntoDeBruijnGraphMultiThreaded(ArrayDeque<String> fragmentsFasta, int sampleSize, int numThreads) throws InterruptedException, IOException {
        ArrayBlockingQueue<String> fragmentsQueue = new ArrayBlockingQueue<>(sampleSize, true);

        FragmentDbgWorker[] workers = new FragmentDbgWorker[numThreads];
        Thread[] threads = new Thread[numThreads];
        for (int i=0; i<numThreads; ++i) {
            workers[i] = new FragmentDbgWorker(fragmentsQueue);
            threads[i] = new Thread(workers[i]);
            threads[i].start();
        }

        for (String fa : fragmentsFasta) {
            FastaReader fin = new FastaReader(fa);

            while (fin.hasNext()) {
                fragmentsQueue.put(fin.next());
            }

            fin.close();
        }

        for (FragmentDbgWorker w : workers) {
            w.stopWhenEmpty();
        }

        for (Thread t : threads) {
            t.join();
        }
    }
    
    private class FragmentDbgWorker implements Runnable {
        private final ArrayBlockingQueue<String> fragments;
        private boolean keepGoing = true;
        
        public FragmentDbgWorker(ArrayBlockingQueue<String> fragments) {
            this.fragments = fragments;
        }

        public void stopWhenEmpty () {
            keepGoing = false;
        }
        
        @Override
        public void run() {
            NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
            long[] hashVals = itr.hVals;
            String fragment = null;
            
            try {
                while (true) {
                    fragment = fragments.poll(50, TimeUnit.MILLISECONDS);
                                        
                    if (fragment == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
                        itr.start(fragment);
                        while (itr.hasNext()) {
                            itr.next();
                            graph.addDbgOnly(hashVals);
                        }
                    }
                }
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
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
    
//    private boolean okToConnectPair(ArrayList<Kmer2> leftKmers, ArrayList<Kmer2> rightKmers) {
//        if (leftKmers.isEmpty() || rightKmers.isEmpty()) {
//            return false;
//        }
//        
//        for (Kmer2 kmer : leftKmers) {
//            if (!graph.lookupLeftKmer(kmer.hashVals)) {
//                return true;
//            }
//        }
//
//        for (Kmer2 kmer : rightKmers) {
//            if (!graph.lookupRightKmer(kmer.hashVals)) {
//                return true;
//            }
//        }
//        
//        return false;
//        
////        int numKmersNotSeenLeft = 0;
////        for (Kmer kmer : leftKmers) {
////            if (!graph.lookupLeftKmer(kmer.hashVals)) {
////                if (++numKmersNotSeenLeft >= kMinus1) {
////                    return true;
////                }
////            }
////        }
////        
////        int numKmersNotSeenRight = 0;
////        for (Kmer kmer : rightKmers) {
////            if (!graph.lookupRightKmer(kmer.hashVals)) {
////                if (++numKmersNotSeenRight >= kMinus1) {
////                    return true;
////                }
////            }
////        }
////        
////        return numKmersNotSeenLeft + numKmersNotSeenRight >= kMinus1;
//    }
    
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
        ArrayList<Kmer2> leftKmers;
        ArrayList<Kmer2> rightKmers;
        boolean corrected = false;
        
        public ReadPair(ArrayList<Kmer2> leftKmers, ArrayList<Kmer2> rightKmers, boolean corrected) {
            this.leftKmers = leftKmers;
            this.rightKmers = rightKmers;
            this.corrected = corrected;
        }
    }
    
    private class Fragment {
        String left;
        String right;
        String seq;
        int length;
        float minCov;
        boolean isUnconnectedRead;
        
        public Fragment(String left, String right, String seq, int length, float minCov, boolean isUnconnectedRead) {
            this.left = left;
            this.right = right;
            this.seq = seq;
            this.length = length;
            this.minCov = minCov;
            this.isUnconnectedRead = isUnconnectedRead;
        }
    }
    
    private class Transcript {
        String fragment;
        ArrayList<Kmer2> transcriptKmers;
        
        public Transcript(String fragment, ArrayList<Kmer2> transcriptKmers) {
            this.fragment = fragment;
            this.transcriptKmers = transcriptKmers;
        }
    }
    
//    private class AssembledTranscriptsQueue {
//        
//        private ArrayBlockingQueue<Transcript> queue;
//        
//        public AssembledTranscriptsQueue(int size) {
//            queue = new ArrayBlockingQueue<> (size);
//        }
//        
//        public synchronized void add(String fragment, ArrayList<Kmer> kmers) {
//            if (!represented(kmers,
//                                graph,
//                                screeningBf,
//                                lookahead,
//                                maxIndelSize,
//                                percentIdentity)) {
//
//                for (Kmer kmer : kmers) {
//                    screeningBf.add(kmer.hashVals);
//                }
//
//                String transcript = assemble(kmers, k);
//
//                try {
//                    queue.put(new Transcript(fragment, transcript));
//                } catch (Exception ex) {
//                    ex.printStackTrace();
//                }
//            }
//        }
//        
//        public int remainingCapacity() {
//            return queue.remainingCapacity();
//        }
//        
//        public boolean isEmpty() {
//            return queue.isEmpty();
//        }
//        
//        public Transcript poll() {
//            return queue.poll();
//        }
//    }
    
    private class TranscriptWriter {
        private final FastaWriter fout;
        private final FastaWriter foutShort;
        private final int minTranscriptLength;
        private String prefix = "";
        private long cid = 0;
        
        public TranscriptWriter(FastaWriter fout, 
                                FastaWriter foutShort,
                                int minTranscriptLength) {
            this.fout = fout;
            this.foutShort = foutShort;
            this.minTranscriptLength = minTranscriptLength;
        }
        
        public void setOutputPrefix(String prefix) {
            this.prefix = prefix;
        }
        
        public void write(ArrayList<Kmer2> transcriptKmers) throws IOException {
            if (!represented(transcriptKmers,
                                graph,
                                screeningBf,
                                lookahead,
                                maxIndelSize,
                                maxTipLength,
                                percentIdentity)) {

                for (Kmer2 kmer : transcriptKmers) {
                    screeningBf.add(kmer.getHash());
                }

                String transcript = graph.assemble(transcriptKmers);
                int len = transcript.length();
                
                if (len >= minTranscriptLength) {
                    fout.write(prefix +  Long.toString(++cid) + " l=" + len, transcript);
                }
                else {
                    foutShort.write(prefix +  Long.toString(++cid) + " l=" + len, transcript);
                }
            }
        }
        
        public void write(String fragment, ArrayList<Kmer2> transcriptKmers) throws IOException {
            if (!represented(transcriptKmers,
                                graph,
                                screeningBf,
                                lookahead,
                                maxIndelSize,
                                maxTipLength,
                                percentIdentity)) {

                for (Kmer2 kmer : transcriptKmers) {
                    screeningBf.add(kmer.getHash());
                }

                String transcript = graph.assemble(transcriptKmers);
                int len = transcript.length();
                
                if (len >= minTranscriptLength) {
                    fout.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + fragment + "]", transcript);
                }
                else {
                    foutShort.write(prefix +  Long.toString(++cid) + " l=" + len + " F=[" + fragment + "]", transcript);
                }
            }
        }
    }
    
    private class TranscriptAssemblyWorker implements Runnable {
        
        private final ArrayBlockingQueue<String> fragments;
        private final ArrayBlockingQueue<Transcript> transcripts;
        private boolean keepGoing = true;
        private boolean includeNaiveExtensions = false;
        private boolean extendBranchFreeFragmentsOnly = false;
        
        public TranscriptAssemblyWorker(ArrayBlockingQueue<String> fragments,
                                        ArrayBlockingQueue<Transcript> transcripts,
                                        boolean includeNaiveExtensions,
                                        boolean extendBranchFreeFragmentsOnly) {
            this.fragments = fragments;
            this.transcripts = transcripts;
            this.includeNaiveExtensions = includeNaiveExtensions;
            this.extendBranchFreeFragmentsOnly = extendBranchFreeFragmentsOnly;
        }

        public void stopWhenEmpty () {
            keepGoing = false;
        }
        
        @Override
        public void run() {
            try {
                while (true) {
                    String fragment = fragments.poll(10, TimeUnit.MICROSECONDS);
                    
                    if (fragment == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
//                        if (fragment.equals("CGTTAGGAAAGCCTGCCGGTGACTAACCCTGCGCTCCTGCCTCGATGGGTGGAGTCGCGTGTGGCGGGGAAGTCAGGTGGAGCGAGGCTAGCTGGCCCGATTTCTCCTCCGGGTGATGCTTTTCCTAGATTATTCTCTGGTAAATCAAAGAAGTGGGTTTATGGAGGTCCTCTTGTGTCCCCTCCCCGCAGAGGTGTGGTGGCTGTGGCATGGTGCCAAG")) {
//                            System.out.println(fragment);
//                        }
                        
                        ArrayList<Kmer2> fragKmers = graph.getKmers(fragment);

                        if (!extendBranchFreeFragmentsOnly || isBranchFree(fragKmers, graph, maxTipLength)) {
                            if (!represented(fragKmers,
                                                graph,
                                                screeningBf,
                                                lookahead,
                                                maxIndelSize,
                                                maxTipLength,
                                                percentIdentity)) {

                                if (includeNaiveExtensions) {
                                    extendWithPairedKmers(fragKmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs, 0.1f);
                                }
                                else {
                                    extendWithPairedKmersDFS(fragKmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs, 0.1f);
                                }

                                transcripts.put(new Transcript(fragment, fragKmers));
                            }
                        }
                    }
                }
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
    }
    
    private class TranscriptWriterWorker implements Runnable {
        
        private final ArrayBlockingQueue<Transcript> transcripts;
        private final TranscriptWriter writer;
        private boolean keepGoing = true;
        
        public TranscriptWriterWorker(ArrayBlockingQueue<Transcript> transcripts,
                                        TranscriptWriter writer) {
            this.transcripts = transcripts;
            this.writer = writer;
        }

        public void stopWhenEmpty () {
            keepGoing = false;
        }
        
        @Override
        public void run() {
            try {
                while (true) {
                    Transcript t = transcripts.poll(10, TimeUnit.MICROSECONDS);
                                        
                    if (t == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
                        writer.write(t.fragment, t.transcriptKmers);
                    }
                }
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
    }
//    
//    private class TranscriptAssembler implements Runnable {
//        
//        private String fragment;
//        private AssembledTranscriptsQueue outList;
//        
//        public TranscriptAssembler(String fragment,
//                                    AssembledTranscriptsQueue outList) {
//            this.fragment = fragment;
//            this.outList = outList;
//        }
//
//        @Override
//        public void run() {
//            ArrayList<Kmer> fragKmers = graph.getKmers(fragment);
//            
////            if (fragment.equals("TGCAGATCGAGAGCCTGAATGAAGAGCTAGCCTACATGAAGAAGAACCATGAAGAGGAGATGAAGGAATTTAGCAACCAGGTGGTCGGCCAGGTCAACGTGGAGATGGATGCCACCCCAGGCATTGACCTGACCCGCGTGCTGGCAGAGATGAGGGAGCAGTACGAGGCCATGGCAGAGAGGAA")) {
////                System.out.println("here");
////            }
//            
//            ArrayList<Kmer> correctedFragKmers = correctErrorsSE(fragKmers,
//                                                                graph, 
//                                                                lookahead,
//                                                                maxIndelSize,
//                                                                maxCovGradient, 
//                                                                covFPR,
//                                                                percentIdentity);
//            if (correctedFragKmers != null) {
//                fragKmers = correctedFragKmers;
//            }
//            
//            if (!represented(fragKmers,
//                                graph,
//                                screeningBf,
//                                lookahead,
//                                maxIndelSize,
//                                percentIdentity)) {
//                                
////                extendWithPairedKmers(fragKmers, graph, lookahead, maxTipLength, beGreedy, screeningBf, maxIndelSize, percentIdentity);
//                extendWithPairedKmers2(fragKmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs);
//
//                
//                outList.add(fragment, fragKmers);
//            }
//        }
//    }
    
    private class ReadConnector implements Runnable {
        private String left;
        private String right;
        private ArrayBlockingQueue<Fragment> outList;
        private int bound;
        private int minOverlap;
        private boolean storeKmerPairs;
        private boolean extendFragments;
        private int errorCorrectionIterations = 0;
        
        public ReadConnector(String left,
                                String right,
                                ArrayBlockingQueue<Fragment> outList,
                                int bound, 
                                int minOverlap,
                                int errorCorrectionIterations,
                                boolean storeKmerPairs,
                                boolean extendFragments) {
            
            this.left = left;
            this.right = right;
            this.outList = outList;
            this.bound = bound;
            this.minOverlap = minOverlap;
            this.storeKmerPairs = storeKmerPairs;
            this.extendFragments = extendFragments;
            this.errorCorrectionIterations = errorCorrectionIterations;
        }
        
        @Override
        public void run() {
            try {
                
                ArrayList<Kmer2> leftKmers = graph.getKmers(left);
                ArrayList<Kmer2> rightKmers = graph.getKmers(right);

                if (this.errorCorrectionIterations > 0) {

                    ReadPair correctedReadPair = correctErrorsPE(leftKmers,
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
                        leftKmers = correctedReadPair.leftKmers;
                        rightKmers = correctedReadPair.rightKmers;
                    }
                }
                
                if (!leftKmers.isEmpty() && !rightKmers.isEmpty()) {

                    ArrayList<Kmer2> fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound, lookahead, minOverlap, maxCovGradient);

                    if (fragmentKmers != null) {
                        int fragLength = fragmentKmers.size() + k - 1;

                        if (fragLength >= k + lookahead) {
                            boolean hasComplexKmer = false;

                            float minCov = Float.MAX_VALUE;
                            for (Kmer2 kmer : fragmentKmers) {
                                if (kmer.count < minCov) {
                                    minCov = kmer.count;
                                }

                                if (!hasComplexKmer) {
                                    if (!graph.isLowComplexity(kmer)) {
                                        hasComplexKmer = true;
                                    }
                                }
                            }

                            if (hasComplexKmer) {
                                if (extendFragments) {
                                    fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength);
                                }

                                if (this.storeKmerPairs) {
                                    graph.addPairedKmers(fragmentKmers);
                                }

                                outList.put(new Fragment(left, right, graph.assemble(fragmentKmers), fragLength, minCov, false));
                            }
                        }
                    }
                    else {
                        // this is an unconnected read pair
                        float minCov = Float.MAX_VALUE;

                        boolean hasComplexLeftKmer = false;

                        if (leftKmers.size() >= lookahead) {
                            for (Kmer2 kmer : leftKmers) {
                                if (kmer.count < minCov) {
                                    minCov = kmer.count;
                                }

                                if (!hasComplexLeftKmer && !graph.isLowComplexity(kmer)) {
                                    hasComplexLeftKmer = true;
                                }
                            }
                        }

                        boolean hasComplexRightKmer = false;

                        if (rightKmers.size() >= lookahead) {
                            for (Kmer2 kmer : rightKmers) {
                                if (kmer.count < minCov) {
                                    minCov = kmer.count;
                                }

                                if (!hasComplexRightKmer && !graph.isLowComplexity(kmer)) {
                                    hasComplexRightKmer = true;
                                }
                            }
                        }

                        if (hasComplexLeftKmer && hasComplexRightKmer) {
                            outList.put(new Fragment(graph.assemble(leftKmers), graph.assemble(rightKmers), null, 0, minCov, true));
                        }
                    }
                }
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
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
        private int polyXMinLen;
        private int polyXMaxMismatches;
        private boolean extendFragments;
        
        public FragmentAssembler(FastqReadPair p,
                                ArrayBlockingQueue<Fragment> outList,
                                int bound, 
                                int minOverlap, 
                                boolean storeKmerPairs, 
                                int errorCorrectionIterations,
                                int leftReadLengthThreshold,
                                int rightReadLengthThreshold,
                                int polyXMinLen,
                                int polyXMaxMismatches,
                                boolean extendFragments) {
            
            this.p = p;
            this.outList = outList;
            this.bound = bound;
            this.minOverlap = minOverlap;
            this.storeKmerPairs = storeKmerPairs;
            this.errorCorrectionIterations = errorCorrectionIterations;
            this.leftReadLengthThreshold = leftReadLengthThreshold;
            this.rightReadLengthThreshold = rightReadLengthThreshold;
            this.polyXMinLen = polyXMinLen;
            this.polyXMaxMismatches = polyXMaxMismatches;
            this.extendFragments = extendFragments;
        }
        
        @Override
        public void run() {
            try {
                // connect segments of each read
                String left = connect(p.left, graph, lookahead);
                
                if (left.length() < this.leftReadLengthThreshold) {
                    return;
                } 
                
                String right = connect(p.right, graph, lookahead);
                
                if (right.length() < this.rightReadLengthThreshold) {
                    return;
                }
                
//                right = chompRightPolyX(right, polyXMinLen, polyXMaxMismatches);
                
                if (right.length() < this.rightReadLengthThreshold) {
                    return;
                }
                
//                left = "GCAACAGGGTGGTGGACCTCATGGCCCACATGGCCTCCAAGGAGTAAGACCCCTGGACCACCAGCCCCAGCAAGAGCACAAGAGGAAGAGAGAGACCCTC";
//                right = "GCTGGGGAGTCCCTGCCACACTAAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTTCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGC";
                
                ArrayList<Kmer2> leftKmers = graph.getKmers(left);
                ArrayList<Kmer2> rightKmers = graph.getKmers(right);

                if (!leftKmers.isEmpty() && !rightKmers.isEmpty()) {
//                    if (okToConnectPair(leftKmers, rightKmers)) {
//                        boolean corrected = false;

                    if (this.errorCorrectionIterations > 0) {

//                            if (left.equals("CTCACGTATTCCCCCAGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGGAAGCTTGTCATCAATGG") &&
//                                    right.equals("TGGAAGAAATGTGCTTTGGGGAGGCAACTAGGATGGTGTGGCTCCCTTGGGTATATGGTAACCTTGTGTCCCTCAATATGGTCCTGTCCCCATCTCCCCC")) {
//                                System.out.println("here");
//                            }
//
//                            System.out.println(left + " " + right);

                        ReadPair correctedReadPair = correctErrorsPE(leftKmers,
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
//                                corrected = true;
                            leftKmers = correctedReadPair.leftKmers;
                            rightKmers = correctedReadPair.rightKmers;
                        }
                    }

//                        if (!corrected || okToConnectPair(leftKmers, rightKmers)) {
                        ArrayList<Kmer2> fragmentKmers = null;

//                        if (!graph.isLowComplexity(leftKmers.get(leftKmers.size()-1)) &&  
//                                !graph.isLowComplexity(rightKmers.get(0))) {
                            fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound, lookahead, minOverlap, maxCovGradient);
//                        }

                        if (fragmentKmers != null) {
                            int fragLength = fragmentKmers.size() + k - 1;

                            if (fragLength >= k + lookahead) {
                                boolean hasComplexKmer = false;

                                float minCov = Float.MAX_VALUE;
                                for (Kmer2 kmer : fragmentKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexKmer) {
                                        if (!graph.isLowComplexity(kmer)) {
                                            hasComplexKmer = true;
                                        }
                                    }
                                }

                                if (hasComplexKmer) {
                                    if (extendFragments) {
                                        fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength);
                                    }
                                    
                                    if (this.storeKmerPairs) {
                                        graph.addPairedKmers(fragmentKmers);
                                    }

                                    outList.put(new Fragment(left, right, graph.assemble(fragmentKmers), fragLength, minCov, false));
                                }
                            }
                        }
                        else {
                            // this is an unconnected read pair
                            float minCov = Float.MAX_VALUE;

                            boolean hasComplexLeftKmer = false;

                            if (leftKmers.size() >= lookahead) {
                                for (Kmer2 kmer : leftKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexLeftKmer && !graph.isLowComplexity(kmer)) {
                                        hasComplexLeftKmer = true;
                                    }
                                }
                            }

                            boolean hasComplexRightKmer = false;

                            if (rightKmers.size() >= lookahead) {
                                for (Kmer2 kmer : rightKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexRightKmer && !graph.isLowComplexity(kmer)) {
                                        hasComplexRightKmer = true;
                                    }
                                }
                            }

                            if (hasComplexLeftKmer && hasComplexRightKmer) {
                                outList.put(new Fragment(graph.assemble(leftKmers), graph.assemble(rightKmers), null, 0, minCov, true));
                            }
                        }
//                        }
                }
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
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
                if (queue.remainingCapacity() > 0) {
                    service.submit(r);
                    break;                    
                }
                
//                try {
//                    service.submit(r);
//                    break;
//                }
//                catch(RejectedExecutionException e) {
//                    // do nothing
//                }
            }
        }
        
        public void terminate() throws InterruptedException {
            service.shutdown();
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        }
        
        public int getQueueRemainingCapacity() {
            return queue.remainingCapacity();
        }
    }
        
    public int[] getMinQ1MedianQ3Max(ArrayList<Integer> a) {
        if (a.size() == 1) {
            int[] result = new int[5];
            Arrays.fill(result, a.get(0));
            return result;
        }
        
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
    
    public void setupKmerScreeningBloomFilter(long sbfNumBits, int sbfNumHash) {
        if (screeningBf == null) {
            screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
        }
        else {
            screeningBf.empty();
        }
    }
    
    public void rescueUnconnectedMultiThreaded(String[] fastas, 
                                                String[] longFragmentsFastaPaths,
                                                String[] shortFragmentsFastaPaths,
                                                String[] unconnectedReadsFastaPaths,
                                                String longSingletonsFasta,
                                                String shortSingletonsFasta,
                                                String unconnectedSingletonsFasta,
                                                int bound,
                                                int minOverlap,
                                                int sampleSize, 
                                                int numThreads, 
                                                int maxErrCorrIterations,
                                                boolean extendFragments) {
        
        /* make sure paired kmer distance is set */
        
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }
        
        System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
        System.out.println("Counting Bloom filter FPR: " + covFPR * 100 + " %");
        
        
        System.out.println("Assembling fragments...");
                
        long fragmentId = 0;
        long unconnectedReadId = 0;
        long readPairsParsed = 0;
        
        int maxTasksQueueSize = numThreads;
        
        int newBound = bound;
        int shortestFragmentLengthAllowed = k + lookahead;
        
        // set up thread pool
        MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);
                
        try {
            FastaReader in;
            
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

            FastaWriter[] unconnectedReadsOut = new FastaWriter[]{new FastaWriter(unconnectedReadsFastaPaths[0], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[1], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[2], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[3], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[4], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[5], true)};
            
            FastaWriter longSingletonsOut = new FastaWriter(longSingletonsFasta, true);
            FastaWriter shortSingletonsOut = new FastaWriter(shortSingletonsFasta, true);
            FastaWriter unconnectedSingletonsOut = new FastaWriter(unconnectedSingletonsFasta, true);
            
            ArrayBlockingQueue<Fragment> fragments = new ArrayBlockingQueue<>(sampleSize);
                        
            for (String fasta : fastas) {
                in = new FastaReader(fasta);
                
                System.out.println("Parsing `" + fasta + "`...");
                                
                // assemble the remaining fragments in multi-threaded mode
                while (in.hasNext()) {
                    String left = in.next();
                    String right = in.next();
                    ++readPairsParsed;
                    
                    if (left.length() >= k && right.length() >= k) {
                        service.submit(new ReadConnector(left,
                                                        right,
                                                        fragments,
                                                        newBound, 
                                                        minOverlap,
                                                        maxErrCorrIterations,
                                                        false, // don't store the kmer pairs
                                                        extendFragments
                        ));
                    }

                    if (fragments.remainingCapacity() <= numThreads) {

                        // write fragments to file
                        int m;
                        Fragment frag;
                        for (int i=0; i<sampleSize; ++i) {
                            frag = fragments.poll();

                            if (frag == null) {
                                break;
                            }

                            if (frag.isUnconnectedRead) {
                                if (frag.minCov == 1) {
                                    unconnectedSingletonsOut.write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                                    unconnectedSingletonsOut.write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                                }
                                else if (frag.minCov > 1) {
                                    m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                    if (m >= 0) {
                                        unconnectedReadsOut[m].write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                                        unconnectedReadsOut[m].write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                                    }
                                }
                            }
                            else {
                                if (frag.length >= shortestFragmentLengthAllowed) {
                                    ArrayList<Kmer2> fragKmers = graph.getKmers(frag.seq);

                                    if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                        if (frag.minCov == 1) {
                                            graph.addPairedKmers(fragKmers);

                                            if (frag.length >= longFragmentLengthThreshold) {
                                                for (Kmer2 kmer : fragKmers) {
                                                    screeningBf.add(kmer.getHash());
                                                } 

                                                longSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                            }
                                            else {
                                                shortSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                            }
                                        }
                                        else if (frag.minCov > 1) {
                                            m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                            if (m >= 0) {
                                                graph.addPairedKmers(fragKmers);

                                                if (frag.length >= longFragmentLengthThreshold) {
                                                    for (Kmer2 kmer : fragKmers) {
                                                        screeningBf.add(kmer.getHash());
                                                    }

                                                    longFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                }
                                                else {
                                                    shortFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                                                
                in.close();
            }

            service.terminate();

            // write fragments to file
            int m;
            Fragment frag;
            while (!fragments.isEmpty()) {
                frag = fragments.poll();

                if (frag.isUnconnectedRead) {
                    if (frag.minCov == 1) {
                        unconnectedSingletonsOut.write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                        unconnectedSingletonsOut.write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                    }
                    else if (frag.minCov > 1) {
                        m = getMinCoverageOrderOfMagnitude(frag.minCov);

                        if (m >= 0) {
                            unconnectedReadsOut[m].write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                            unconnectedReadsOut[m].write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                        }
                    }
                }
                else {
                    if (frag.length >= shortestFragmentLengthAllowed) {
                        ArrayList<Kmer2> fragKmers = graph.getKmers(frag.seq);

                        if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                            if (frag.minCov == 1) {
                                graph.addPairedKmers(fragKmers);

                                if (frag.length >= longFragmentLengthThreshold) {
                                    for (Kmer2 kmer : fragKmers) {
                                        screeningBf.add(kmer.getHash());
                                    }

                                    longSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                }
                                else {
                                    shortSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                }
                            }
                            else if (frag.minCov > 1)  {
                                m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                if (m >= 0) {
                                    graph.addPairedKmers(fragKmers);

                                    if (frag.length >= longFragmentLengthThreshold) {
                                        for (Kmer2 kmer : fragKmers) {
                                            screeningBf.add(kmer.getHash());
                                        }

                                        longFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                    }
                                    else {
                                        shortFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            unconnectedSingletonsOut.close();
            longSingletonsOut.close();
            shortSingletonsOut.close();

            for (FastaWriter out : longFragmentsOut) {
                out.close();
            }
            
            for (FastaWriter out : shortFragmentsOut) {
                out.close();
            }
            
            for (FastaWriter out : unconnectedReadsOut) {
                out.close();
            }
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
            
            System.out.println("Paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
            System.out.println("Screening Bloom filter FPR:    " + screeningBf.getFPR() * 100 + " %");
        }        
    }
    
    private static final String LABEL_SEPARATOR = ":";
    private final static String LABEL_MIN = "min";
    private final static String LABEL_Q1 = "Q1";
    private final static String LABEL_MEDIAN = "M";
    private final static String LABEL_Q3 = "Q3";
    private final static String LABEL_MAX = "max";
    
    public void writeFragStatsToFile(int[] fragStats, String path) {
        try {
            FileWriter writer = new FileWriter(path, false);

            writer.write(LABEL_MIN + LABEL_SEPARATOR + fragStats[0] + "\n" +
                        LABEL_Q1 + LABEL_SEPARATOR + fragStats[1] + "\n" +
                        LABEL_MEDIAN + LABEL_SEPARATOR + fragStats[2] + "\n" +
                        LABEL_Q3 + LABEL_SEPARATOR + fragStats[3] + "\n" +
                        LABEL_MAX + LABEL_SEPARATOR + fragStats[4] + "\n"
                    );
            writer.close();
        } catch (Exception ex) {
            handleException(ex);
        }
    }
    
    public int[] restoreFragStatsFromFile(String path) {
        int[] fragStats = new int[5];
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while ((line = br.readLine()) != null) {
                String[] entry = line.split(LABEL_SEPARATOR);
                String key = entry[0];
                String val = entry[1];
                switch(key) {
                    case LABEL_MIN:
                        fragStats[0] = Integer.parseInt(val);
                        break;
                    case LABEL_Q1:
                        fragStats[1] = Integer.parseInt(val);
                        break;
                    case LABEL_MEDIAN:
                        fragStats[2] = Integer.parseInt(val);
                        break;
                    case LABEL_Q3:
                        fragStats[3] = Integer.parseInt(val);
                        break;
                    case LABEL_MAX:
                        fragStats[4] = Integer.parseInt(val);
                        break;
                }
            }
            br.close();
        
        } catch (Exception ex) {
            handleException(ex);
        }
        
        return fragStats;
    }
    
    public void setPairedKmerDistance(int fragStatQ1) {
        graph.setPairedKmerDistance(fragStatQ1 - k - minNumKmerPairs);
    }
    
    public int getPairedReadsMaxDistance(int[] fragStats) {
        return fragStats[3] + ((fragStats[3] - fragStats[1]) * 3 / 2); // 1.5*IQR
    }
    
    public int[] assembleFragmentsMultiThreaded(FastqPair[] fastqPairs, 
                                                String[] longFragmentsFastaPaths,
                                                String[] shortFragmentsFastaPaths,
                                                String[] unconnectedReadsFastaPaths,
                                                String longSingletonsFasta,
                                                String shortSingletonsFasta,
                                                String unconnectedSingletonsFasta,
                                                int bound,
                                                int minOverlap,
                                                int sampleSize, 
                                                int numThreads, 
                                                int maxErrCorrIterations,
                                                boolean extendFragments) {
        
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
        long unconnectedReadId = 0;
        long readPairsParsed = 0;
        
        int maxTasksQueueSize = numThreads;
        
        int newBound = bound;
        int[] fragLengthsStats = null;
        boolean pairedKmerDistanceIsSet = false;
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

            FastaWriter[] unconnectedReadsOut = new FastaWriter[]{new FastaWriter(unconnectedReadsFastaPaths[0], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[1], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[2], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[3], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[4], true),
                                                                new FastaWriter(unconnectedReadsFastaPaths[5], true)};
            
            FastaWriter longSingletonsOut = new FastaWriter(longSingletonsFasta, true);
            FastaWriter shortSingletonsOut = new FastaWriter(shortSingletonsFasta, true);
            FastaWriter unconnectedSingletonsOut = new FastaWriter(unconnectedSingletonsFasta, true);
            
            FastqReadPair p;
            ArrayBlockingQueue<Fragment> fragments = new ArrayBlockingQueue<>(sampleSize);
            
//            boolean readLengthThresholdIsSet = false;
//            ArrayList<Integer> leftReadLengths = new ArrayList<>(sampleSize);
//            ArrayList<Integer> rightReadLengths = new ArrayList<>(sampleSize);
            
//            int minNumKmersNotAssembled = 1;
            
            for (FastqPair fqPair: fastqPairs) {
                lin = new FastqReader(fqPair.leftFastq, true);
                rin = new FastqReader(fqPair.rightFastq, true);

                fqpr = new FastqPairReader(lin, rin, qualPatternFrag, seqPattern, fqPair.leftRevComp, fqPair.rightRevComp);
                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");

                int leftReadLengthThreshold = k;
                int rightReadLengthThreshold = k;
                int polyXMinLen = k;
                int polyXMaxMismatches = Math.round(k*percentError);
                
                if (!pairedKmerDistanceIsSet) {
                                        
                    // Assembled initial sample of fragments
                    while (fqpr.hasNext()) {
                        p = fqpr.next();
                        ++readPairsParsed;
                            
                            service.submit(new FragmentAssembler(p,
                                                                fragments,
                                                                bound,
                                                                minOverlap,
                                                                false, // do not store paired kmers
                                                                maxErrCorrIterations, 
                                                                leftReadLengthThreshold,
                                                                rightReadLengthThreshold,
                                                                polyXMinLen,
                                                                polyXMaxMismatches,
                                                                extendFragments
                            ));
                            
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
                   
                    fragLengthsStats = getMinQ1MedianQ3Max(fragLengths);
                    System.out.println("Fragment Lengths Distribution (n=" + fragLengths.size() + ")");
                    System.out.println("\tmin\tQ1\tM\tQ3\tmax");
                    System.out.println("\t" + fragLengthsStats[0] + "\t" + fragLengthsStats[1] + "\t" + fragLengthsStats[2] + "\t" + fragLengthsStats[3] + "\t" + fragLengthsStats[4]);

                    longFragmentLengthThreshold = fragLengthsStats[1];
                    setPairedKmerDistance(longFragmentLengthThreshold);
                    
                    // Set new bound for graph search
                    newBound = getPairedReadsMaxDistance(fragLengthsStats);
                    
                    System.out.println("Paired kmers distance:       " + (longFragmentLengthThreshold - k - minNumKmerPairs));
                    System.out.println("Max graph traversal depth:   " + newBound);
//                    System.out.println("Shortest fragment allowed:   " + shortestFragmentLengthAllowed);

                    // Store kmer pairs and write fragments to file
                    int m;
                    Fragment frag;
                    while (!fragments.isEmpty()) {
                        frag = fragments.poll();
                        
                        if (frag.isUnconnectedRead) {
                            if (frag.minCov == 1) {
                                unconnectedSingletonsOut.write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                                unconnectedSingletonsOut.write(Long.toString(unconnectedReadId) + "R", frag.right);
                            }
                            else if (frag.minCov > 1)  {
                                m = getMinCoverageOrderOfMagnitude(frag.minCov);
                                
                                if (m >= 0) {
                                    unconnectedReadsOut[m].write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                                    unconnectedReadsOut[m].write(Long.toString(unconnectedReadId) + "R", frag.right);
                                }
                            }
                        }
                        else {
                            if (frag.length >= shortestFragmentLengthAllowed) {
                                ArrayList<Kmer2> fragKmers = graph.getKmers(frag.seq);

                                if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                    if (frag.minCov == 1) {
                                        graph.addPairedKmers(fragKmers);

                                        if (frag.length >= longFragmentLengthThreshold) {
                                            for (Kmer2 kmer : fragKmers) {
                                                screeningBf.add(kmer.getHash());
                                            }

                                            longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                        }
                                        else {
                                            shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                        }
                                    }
                                    else if (frag.minCov > 1)  {
                                        m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                        if (m >= 0) {
                                            graph.addPairedKmers(fragKmers);

                                            if (frag.length >= longFragmentLengthThreshold) {
                                                for (Kmer2 kmer : fragKmers) {
                                                    screeningBf.add(kmer.getHash());
                                                }

                                                longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                            }
                                            else {
                                                shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                            }
                                        }
                                    }
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
//                    if (p.originalLeftLength > 2*p.numLeftBasesTrimmed &&
//                            p.originalRightLength > 2*p.numRightBasesTrimmed) {
                    
                    if (!p.left.isEmpty() && !p.right.isEmpty()) {
                        
                        service.submit(new FragmentAssembler(p,
                                                            fragments,
                                                            newBound, 
                                                            minOverlap, 
                                                            false, // do not store paired k-mers
                                                            maxErrCorrIterations, 
                                                            leftReadLengthThreshold, 
                                                            rightReadLengthThreshold,
                                                            polyXMinLen,
                                                            polyXMaxMismatches,
                                                            extendFragments
                        ));

                        if (fragments.remainingCapacity() <= numThreads) {

                            // write fragments to file
                            int m;
                            Fragment frag;
                            for (int i=0; i<sampleSize; ++i) {
                                frag = fragments.poll();
                                
                                if (frag == null) {
                                    break;
                                }
                                
                                if (frag.isUnconnectedRead) {
                                    if (frag.minCov == 1) {
                                        unconnectedSingletonsOut.write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                                        unconnectedSingletonsOut.write(Long.toString(unconnectedReadId) + "R", frag.right);
                                    }
                                    else if (frag.minCov > 1) {
                                        m = getMinCoverageOrderOfMagnitude(frag.minCov);
                                        
                                        if (m >= 0) {
                                            unconnectedReadsOut[m].write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                                            unconnectedReadsOut[m].write(Long.toString(unconnectedReadId) + "R", frag.right);
                                        }
                                    }
                                }
                                else {
                                    if (frag.length >= shortestFragmentLengthAllowed) {
                                        ArrayList<Kmer2> fragKmers = graph.getKmers(frag.seq);

                                        if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                            if (frag.minCov == 1) {
                                                graph.addPairedKmers(fragKmers);

                                                if (frag.length >= longFragmentLengthThreshold) {
                                                    for (Kmer2 kmer : fragKmers) {
                                                        screeningBf.add(kmer.getHash());
                                                    } 

                                                    longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                }
                                                else {
                                                    shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                }
                                            }
                                            else if (frag.minCov > 1) {
                                                m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                                if (m >= 0) {
                                                    graph.addPairedKmers(fragKmers);

                                                    if (frag.length >= longFragmentLengthThreshold) {
                                                        for (Kmer2 kmer : fragKmers) {
                                                            screeningBf.add(kmer.getHash());
                                                        }

                                                        longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                    }
                                                    else {
                                                        shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                                    }
                                                }
                                            }
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
                    
                    if (frag.isUnconnectedRead) {
                        if (frag.minCov == 1) {
                            unconnectedSingletonsOut.write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                            unconnectedSingletonsOut.write(Long.toString(unconnectedReadId) + "R", frag.right);
                        }
                        else if (frag.minCov > 1) {
                            m = getMinCoverageOrderOfMagnitude(frag.minCov);
                            
                            if (m >= 0) {
                                unconnectedReadsOut[m].write(Long.toString(++unconnectedReadId) + "L ", frag.left);
                                unconnectedReadsOut[m].write(Long.toString(unconnectedReadId) + "R", frag.right);
                            }
                        }
                    }
                    else {
                        if (frag.length >= shortestFragmentLengthAllowed) {
                            ArrayList<Kmer2> fragKmers = graph.getKmers(frag.seq);

                            if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                if (frag.minCov == 1) {
                                    graph.addPairedKmers(fragKmers);

                                    if (frag.length >= longFragmentLengthThreshold) {
                                        for (Kmer2 kmer : fragKmers) {
                                            screeningBf.add(kmer.getHash());
                                        }

                                        longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                    }
                                    else {
                                        shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                    }
                                }
                                else if (frag.minCov > 1)  {
                                    m = getMinCoverageOrderOfMagnitude(frag.minCov);

                                    if (m >= 0) {
                                        graph.addPairedKmers(fragKmers);

                                        if (frag.length >= longFragmentLengthThreshold) {
                                            for (Kmer2 kmer : fragKmers) {
                                                screeningBf.add(kmer.getHash());
                                            }

                                            longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                        }
                                        else {
                                            shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", frag.seq);
                                        }
                                    }
                                }
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
            
            unconnectedSingletonsOut.close();
            longSingletonsOut.close();
            shortSingletonsOut.close();

            for (FastaWriter out : longFragmentsOut) {
                out.close();
            }
            
            for (FastaWriter out : shortFragmentsOut) {
                out.close();
            }
            
            for (FastaWriter out : unconnectedReadsOut) {
                out.close();
            }
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
            
            System.out.println("Paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
            System.out.println("Screening Bloom filter FPR:    " + screeningBf.getFPR() * 100 + " %");
        }

        return fragLengthsStats;
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
        } catch (Exception ex) {
            handleException(ex);
        }
    }
    
    public void restorePairedKmersBloomFilter(File graphFile) {
        try {
            graph.destroyPkbf();
            graph.restorePkbf(graphFile);
        } catch (Exception ex) {
            handleException(ex);
        }
    }
    
    private long extendFragmentsMultiThreadedHelper(String fragmentsFasta, 
                                                    TranscriptWriter writer, 
                                                    int sampleSize, 
                                                    int numThreads, 
                                                    boolean includeNaiveExtensions,
                                                    boolean extendBranchFreeFragmentsOnly) throws InterruptedException, IOException {
        
        long numFragmentsParsed = 0;
        FastaReader fin = new FastaReader(fragmentsFasta);

        ArrayBlockingQueue<String> fragmentsQueue = new ArrayBlockingQueue<>(sampleSize, true);
        ArrayBlockingQueue<Transcript> transcriptsQueue = new ArrayBlockingQueue<>(numThreads*2, true);
        
        TranscriptAssemblyWorker[] workers = new TranscriptAssemblyWorker[numThreads];
        Thread[] threads = new Thread[numThreads];
        for (int i=0; i<numThreads; ++i) {
            workers[i] = new TranscriptAssemblyWorker(fragmentsQueue, transcriptsQueue, includeNaiveExtensions, extendBranchFreeFragmentsOnly);
            threads[i] = new Thread(workers[i]);
            threads[i].start();
        }
        
        TranscriptWriterWorker writerWorker = new TranscriptWriterWorker(transcriptsQueue, writer);
        Thread writerThread = new Thread(writerWorker);
        writerThread.start();

        while (fin.hasNext()) {
            ++numFragmentsParsed;
            fragmentsQueue.put(fin.next());
        }

        fin.close();

        for (TranscriptAssemblyWorker w : workers) {
            w.stopWhenEmpty();
        }

        for (Thread t : threads) {
            t.join();
        }

        writerWorker.stopWhenEmpty();
        writerThread.join();
        
        return numFragmentsParsed;
    }
    
    public void assembleTransfragsMultiThreaded(String[] longFragmentsFastas, 
                                                String[] shortFragmentsFastas,
                                                String[] unconnectedReadsFastas,
                                                String[] outFastasLong,
                                                String[] outFastasShort,
                                                int numThreads,
                                                int sampleSize) {
        
        long numFragmentsParsed = 0;
        
        int minTransfragLength = graph.getPairedKmerDistance() + k + minNumKmerPairs;

        try {
            System.out.println("Extending fragments...");
            
            for (int mag=longFragmentsFastas.length-1; mag>=0; --mag) {
                graph.clearDbgbf();
                
                String longFragsFasta = longFragmentsFastas[mag];
                String shortFragsFasta = shortFragmentsFastas[mag];
                String unconnectedReadsFasta = unconnectedReadsFastas[mag];
                
                graph.clearPkbf();
                insertIntoDeBruijnGraphAndPairedKmers(longFragsFasta);
                insertIntoDeBruijnGraph(shortFragsFasta);
                insertIntoDeBruijnGraph(unconnectedReadsFasta);
                
                FastaWriter fout = new FastaWriter(outFastasLong[mag], false);
                FastaWriter foutShort = new FastaWriter(outFastasShort[mag], false);
                TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTransfragLength);
                
                System.out.println("Parsing fragments in `" + longFragsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(longFragsFasta, writer, sampleSize, numThreads, false, false);
                
                System.out.println("Parsing fragments in `" + shortFragsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(shortFragsFasta, writer, sampleSize, numThreads, false, false);
                
                fout.close();
                foutShort.close();
            }
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
            System.out.println("Screening Bloom filter FPR:      " + screeningBf.getFPR() * 100 + " %");
        }
    }
    
    public void assembleTranscriptsMultiThreaded(String[] longFragmentsFastas, 
                                                String[] shortFragmentsFastas,
                                                String[] unconnectedReadsFastas,
                                                String longSingletonsFasta,
                                                String shortSingletonsFasta,
                                                String unconnectedSingletonsFasta,
                                                String outFasta,
                                                String outFastaShort,
                                                String graphFile,
                                                int numThreads,
                                                int sampleSize,
                                                int minTranscriptLength,
                                                boolean createFragmentDBG,
                                                boolean useSingletonFragments,
                                                String txptNamePrefix) {
        
        long numFragmentsParsed = 0;

        try {

//            if (createFragmentDBG) {
//                System.out.println("Creating graph from fragment kmers...");
//                
//                graph.getDbgbf().empty();
//
//                if (numThreads > 1) {
//                    ArrayDeque<String> fragmentsFastas = new ArrayDeque<>(longFragmentsFastas.length + shortFragmentsFastas.length + 2);
//                    fragmentsFastas.addAll(Arrays.asList(longFragmentsFastas));
//                    fragmentsFastas.addAll(Arrays.asList(shortFragmentsFastas));
//                    fragmentsFastas.addAll(Arrays.asList(unconnectedReadsFastas));
//                    fragmentsFastas.add(longSingletonsFasta);
//                    fragmentsFastas.add(shortSingletonsFasta);
//                    fragmentsFastas.add(unconnectedSingletonsFasta);
//                    insertIntoDeBruijnGraphMultiThreaded(fragmentsFastas, sampleSize, numThreads);
//                }
//                else {
//                    insertIntoDeBruijnGraph(longFragmentsFastas);
//                    insertIntoDeBruijnGraph(shortFragmentsFastas);
//                    insertIntoDeBruijnGraph(unconnectedReadsFastas);
//                    insertIntoDeBruijnGraph(longSingletonsFasta);
//                    insertIntoDeBruijnGraph(shortSingletonsFasta);
//                    insertIntoDeBruijnGraph(unconnectedSingletonsFasta);
//                }
//            }
//            else {
//                System.out.println("Restoring graph from file...");
//                
//                restoreDbg(new File(graphFile));
//            }
//
//            graph.restorePkbf(new File(graphFile));
            
            dbgFPR = graph.getDbgbf().getFPR();
            System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
            
            if (covFPR <= 0) {
                covFPR = graph.getCbf().getFPR();
            }

            System.out.println("Assembling transcripts...");
        

            FastaWriter fout = new FastaWriter(outFasta, false);
            FastaWriter foutShort = new FastaWriter(outFastaShort, false);
            TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength);

            String fragmentsFasta;
            
            // extend LONG fragments
            
            for (int mag=longFragmentsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
                fragmentsFasta = longFragmentsFastas[mag];
                System.out.println("Extending fragments in `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads, true, false);
            }          
            
            // extend SHORT fragments
            
            for (int mag=shortFragmentsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
                fragmentsFasta = shortFragmentsFastas[mag];
                System.out.println("Extending fragments in `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads, true, false);
            }

            // extend UNCONNECTED reads
            
            for (int mag=unconnectedReadsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
                fragmentsFasta = unconnectedReadsFastas[mag];
                System.out.println("Extending fragments in `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads, true, false);
            }
            
            if (useSingletonFragments) {

                // extend LONG singleton fragments
            
                writer.setOutputPrefix(txptNamePrefix + "01.L.");
                System.out.println("Extending fragments in `" + longSingletonsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(longSingletonsFasta, writer, sampleSize, numThreads, true, true);

                // extend SHORT singleton fragments
                
                writer.setOutputPrefix(txptNamePrefix + "01.S.");
                System.out.println("Extending fragments in `" + shortSingletonsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(shortSingletonsFasta, writer, sampleSize, numThreads, true, true);
                
                // extend UNCONNECTED reads
                
                writer.setOutputPrefix(txptNamePrefix + "01.U.");
                System.out.println("Extending fragments in `" + unconnectedSingletonsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(unconnectedSingletonsFasta, writer, sampleSize, numThreads, true, true);
            }
            
            fout.close();
            foutShort.close();
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
            System.out.println("Screening Bloom filter FPR:      " + screeningBf.getFPR() * 100 + " %");
        }
    }
    
    private static class MyTimer {
        private final long globalStartTime;
        private long startTime;
        
        public MyTimer() {
            globalStartTime = System.currentTimeMillis();
            startTime = globalStartTime;
        }
        
        public void start() {
            startTime = System.currentTimeMillis();
        }
        
        public long elapsedMillis() {
            return System.currentTimeMillis() - startTime;
        }
                
        public long totalElapsedMillis() {
            return System.currentTimeMillis() - globalStartTime;
        }
        
        public static String hmsFormat(long millis) {
            long seconds = millis / 1000;
            
            long hours = seconds / 3600;
            
            seconds = seconds % 3600;
            
            long minutes = seconds / 60;
            
            seconds = seconds % 60;
            
            StringBuilder sb = new StringBuilder();
            
            if (hours > 0) {
                sb.append(hours);
                sb.append("h ");
            }
            
            if (minutes > 0 || hours > 0) {
                sb.append(minutes);
                sb.append("m ");
            }
                        
            if (hours == 0 && minutes == 0) {
                sb.append(millis/1000f);
                sb.append("s");
            }
            else {
                sb.append(seconds);
                sb.append("s");
            }
            
            return sb.toString();
        }
    }
    
    public static void touch(File f) throws IOException {
        f.getParentFile().mkdirs();
        if (!f.createNewFile()){
            f.setLastModified(System.currentTimeMillis());
        }
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
        final String FRAGMENTS_2K_DONE = "FRAGMENTS.2K.DONE";
//        final String TRANSFRAGS_DONE = "TRANSFRAGS.DONE";
        final String TRANSCRIPTS_DONE = "TRANSCRIPTS.DONE";
        
        MyTimer timer = new MyTimer();
        
        System.out.println("args: " + Arrays.toString(args));
        
        // -mem 0.5 -left /home/gengar/test_data/GAPDH/GAPDH_2.fq.gz -right /home/gengar/test_data/GAPDH/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /home/gengar/test_assemblies/GAPDH
        // -mem 4 -left /home/gengar/test_data/SRR1360926/SRR1360926_2.fastq.gz -right /home/gengar/test_data/SRR1360926/SRR1360926_1.fastq.gz -revcomp-right -stranded -name SRR1360926 -outdir /home/gengar/test_assemblies/SRR1360926
        // -mem 0.5  -left /home/gengar/test_data/SRR1360926/SRR1360926.RNF213.2.fq.gz -right /home/gengar/test_data/SRR1360926/SRR1360926.RNF213.1.fq.gz -revcomp-right -stranded -name RNF213 -outdir /home/gengar/test_assemblies/RNF213
        
        // -t 1 -pair 10 -mem 2 -left /projects/btl2/kmnip/rna-bloom/example/SRX983106/SRR1957705_2.fastq.gz.trim.fq.gz -right /projects/btl2/kmnip/rna-bloom/example/SRX983106/SRR1957705_1.fastq.gz.trim.fq.gz -revcomp-right -stranded -name SRR953119 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/c.elegans/SRR953119/apr07.debug
        // -mem 0.5 -left /projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/gapdh
        // -mem 4 -left /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_1.fq.gz -revcomp-right -stranded -name SRR1360926 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/SRR1360926
        // -mem 30 -left /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_2.fq.gz -right /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_1.fq.gz -revcomp-right -stranded -name mcf7 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/mcf7        
        // Based on: http://commons.apache.org/proper/commons-cli/usage.html
        CommandLineParser parser = new DefaultParser();

        Options options = new Options();

                
        Option optLeftReads = Option.builder("l")
                                    .longOpt("left")
                                    .desc("left reads file")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optLeftReads);
        
        Option optRightReads = Option.builder("r")
                                    .longOpt("right")
                                    .desc("right reads file")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optRightReads);
        
        Option optRevCompLeft = Option.builder("rcl")
                                    .longOpt("revcomp-left")
                                    .desc("reverse-complement left reads")
                                    .hasArg(false)
                                    .build();
        options.addOption(optRevCompLeft);

        Option optRevCompRight = Option.builder("rcr")
                                    .longOpt("revcomp-right")
                                    .desc("reverse-complement right reads")
                                    .hasArg(false)
                                    .build();
        options.addOption(optRevCompRight);
        
        Option optStranded = Option.builder("ss")
                                    .longOpt("stranded")
                                    .desc("strand specific")
                                    .hasArg(false)
                                    .build();
        options.addOption(optStranded);
        
        Option optName = Option.builder("n")
                                    .longOpt("name")
                                    .desc("assembly name")
                                    .hasArg(true)
                                    .argName("STR")
                                    .build();
        options.addOption(optName);
        
        Option optPrefix = Option.builder("prefix")
                                    .longOpt("prefix")
                                    .desc("assembled transcript name prefix")
                                    .hasArg(true)
                                    .argName("STR")
                                    .build();
        options.addOption(optPrefix);
        
        Option optThreads = Option.builder("t")
                                    .longOpt("threads")
                                    .desc("run in INT threads")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optThreads);
        
        Option optOutdir = Option.builder("o")
                                    .longOpt("outdir")
                                    .desc("output directory")
                                    .hasArg(true)
                                    .argName("PATH")
                                    .build();
        options.addOption(optOutdir);
        
        Option optForce = Option.builder("f")
                                    .longOpt("force")
                                    .desc("force overwrite existing files")
                                    .hasArg(false)
                                    .build();
        options.addOption(optForce);
        
        Option optKmerSize = Option.builder("k")
                                    .longOpt("kmer")
                                    .desc("kmer size")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optKmerSize);
        
        Option optBaseQualDbg = Option.builder("q")
                                    .longOpt("qual-dbg")
                                    .desc("min base quality for constructing DBG")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBaseQualDbg);

        Option optBaseQualFrag = Option.builder("Q")
                                    .longOpt("qual-frag")
                                    .desc("min base quality for fragment assembly")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBaseQualFrag);        
        
        Option optAllHash = Option.builder("hash")
                                    .longOpt("all-hash")
                                    .desc("number of hash functions for all Bloom filters")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optAllHash); 
        
        Option optSbfHash = Option.builder("sh")
                                    .longOpt("sbf-hash")
                                    .desc("number of hash functions for screening Bloom filter")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optSbfHash); 
        
        Option optDbgbfHash = Option.builder("dh")
                                    .longOpt("dbgbf-hash")
                                    .desc("number of hash functions for de Bruijn graph Bloom filter")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optDbgbfHash);

        Option optCbfHash = Option.builder("ch")
                                    .longOpt("cbf-hash")
                                    .desc("number of hash functions for kmer counting Bloom filter")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optCbfHash);
        
        Option optPkbfHash = Option.builder("ph")
                                    .longOpt("pkbf-hash")
                                    .desc("number of hash functions for paired kmers Bloom filter")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optPkbfHash);        

        Option optAllMem = Option.builder("mem")
                                    .longOpt("all-mem")
                                    .desc("total amount of memory (GB) for all Bloom filters")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optAllMem);
        
        Option optSbfMem = Option.builder("sm")
                                    .longOpt("sbf-mem")
                                    .desc("amount of memory (GB) for screening Bloom filter")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optSbfMem);
        
        Option optDbgbfMem = Option.builder("dm")
                                    .longOpt("dbgbf-mem")
                                    .desc("amount of memory (GB) for de Bruijn graph Bloom filter")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optDbgbfMem);

        Option optCbfMem = Option.builder("cm")
                                    .longOpt("cbf-mem")
                                    .desc("amount of memory (GB) for kmer counting Bloom filter")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optCbfMem);
        
        Option optPkbfMem = Option.builder("pm")
                                    .longOpt("pkbf-mem")
                                    .desc("amount of memory (GB) for paired kmers Bloom filter")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optPkbfMem);
 
        Option optTipLength = Option.builder("tiplength")
                                    .longOpt("tiplength")
                                    .desc("max tip length allowed")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optTipLength);  
        
        Option optLookahead = Option.builder("lookahead")
                                    .longOpt("lookahead")
                                    .desc("number of kmers to look ahead during graph traversal")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optLookahead);        
        
        Option optOverlap = Option.builder("overlap")
                                    .longOpt("overlap")
                                    .desc("min number of overlapping bases between mates")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optOverlap);
        
        Option optBound = Option.builder("bound")
                                    .longOpt("bound")
                                    .desc("max distance between mates")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBound);

        Option optSample = Option.builder("sample")
                                    .longOpt("sample")
                                    .desc("sample size for estimating median fragment length")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optSample);
        
        Option optMaxCovGrad = Option.builder("grad")
                                    .longOpt("maxcovgrad")
                                    .desc("max coverage gradient for error correction")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optMaxCovGrad);
        
        Option optIndelSize = Option.builder("indel")
                                    .longOpt("indel")
                                    .desc("maximum indel size allowed")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optIndelSize);  

        Option optPercentIdentity = Option.builder("p")
                                    .longOpt("percent")
                                    .desc("minimum percent identity allowed")
                                    .hasArg(true)
                                    .argName("FLOAT")
                                    .build();
        options.addOption(optPercentIdentity); 
        
        Option optErrCorrItr = Option.builder("e")
                                    .longOpt("errcorritr")
                                    .desc("max number of iterations of read error correction")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optErrCorrItr);        

        Option optExtend = Option.builder("extend")
                                    .longOpt("extend")
                                    .desc("extend assembled fragments during fragment assembly")
                                    .hasArg(false)
                                    .build();
        options.addOption(optExtend);
        
//        Option optFdbg = Option.builder("fdbg")
//                                    .longOpt("fdbg")
//                                    .desc("used only fragment kmers during transcript assembly")
//                                    .hasArg(false)
//                                    .build();
//        options.addOption(optFdbg);

        Option optSingleton = Option.builder("1")
                                    .longOpt("singleton")
                                    .desc("assemble transcripts from singleton fragments")
                                    .hasArg(false)
                                    .build();
        options.addOption(optSingleton);
        
        Option optMinKmerPairs = Option.builder("pair")
                                    .longOpt("pair")
                                    .desc("minimum number of consecutive kmer pairs for assembling transcripts")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinKmerPairs);  
        
        Option optMinLength = Option.builder("length")
                                    .longOpt("length")
                                    .desc("min transcript length in final assembly")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinLength);  
        
        Option optHelp = Option.builder("h")
                                    .longOpt("help")
                                    .desc("print this message and exits")
                                    .build();
        options.addOption(optHelp);
        

        try {
            CommandLine line = parser.parse(options, args);

            if (line.getOptions().length == 0 || line.hasOption(optHelp.getOpt())) {
                printHelp(options, false);
            }
            
            int numThreads = Integer.parseInt(line.getOptionValue(optThreads.getOpt(), "2"));
            boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            String name = line.getOptionValue(optName.getOpt(), "rnabloom");
            String outdir = line.getOptionValue(optOutdir.getOpt(), System.getProperty("user.dir") + File.separator + name + "_assembly");
            /**@TODO evaluate whether out dir is a valid dir */
            
            String longFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.long.";
            String shortFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.short.";
            String unconnectedReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.";
//            String long2kFragmentsFastaPrefix =    outdir + File.separator + name + ".fragments.2k.long.";
//            String short2kFragmentsFastaPrefix =   outdir + File.separator + name + ".fragments.2k.short.";
            String unconnected2kReadsFastaPrefix = outdir + File.separator + name + ".unconnected.2k.";
            String transcriptsFasta =              outdir + File.separator + name + ".transcripts.fa";
            String shortTranscriptsFasta =         outdir + File.separator + name + ".transcripts.short.fa";
//            String tmpFasta = outdir + File.separator + name + ".tmp.fa";
            String graphFile = outdir + File.separator + name + ".graph";
            String graph2kFile = outdir + File.separator + name + ".2k.graph";
            String fragStatsFile = outdir + File.separator + name + ".fragstats";
            
            File startedStamp = new File(outdir + File.separator + STARTED);
            File dbgDoneStamp = new File(outdir + File.separator + DBG_DONE);
            File fragsDoneStamp = new File(outdir + File.separator + FRAGMENTS_DONE);
            File frags2kDoneStamp = new File(outdir + File.separator + FRAGMENTS_2K_DONE);
//            File txfgsDoneStamp = new File(outdir + File.separator + TRANSFRAGS_DONE);
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
                
                if (frags2kDoneStamp.exists()) {
                    frags2kDoneStamp.delete();
                }
                
//                if (txfgsDoneStamp.exists()) {
//                    txfgsDoneStamp.delete();
//                }
                
                if (txptsDoneStamp.exists()) {
                    txptsDoneStamp.delete();
                }
            }
            
            String[] fastqsLeft = line.getOptionValues(optLeftReads.getOpt());
            String[] fastqsRight = line.getOptionValues(optRightReads.getOpt());
            
            if (fastqsLeft.length == 0) {
                System.out.println("ERROR: Please specify left read files!");
                System.exit(1);
            }

            if (fastqsRight.length == 0) {
                System.out.println("ERROR: Please specify right read files!");
                System.exit(1);
            }
            
            if (fastqsLeft.length != fastqsRight.length) {
                System.out.println("ERROR: Read files are not paired properly!");
                System.exit(1);
            }
            
            boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            boolean strandSpecific = line.hasOption(optStranded.getOpt());
            
            int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), "25"));
            int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), "3"));
            int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), "3"));
            
            double leftReadFilesTotalBytes = 0;
            for (String fq : fastqsLeft) {
                leftReadFilesTotalBytes += new File(fq).length();
            }
            
            float maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (leftReadFilesTotalBytes * 0.8d / NUM_BYTES_1GB))));
            float sbfGB = Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 7f)));
            float dbgGB = Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), Float.toString(maxBfMem * 1f / 7f)));
            float cbfGB = Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), Float.toString(maxBfMem * 5f / 7f)));
            float pkbfGB = Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 7f)));
            
            long sbfSize = (long) (NUM_BITS_1GB * sbfGB);
            long dbgbfSize = (long) (NUM_BITS_1GB * dbgGB);
            long cbfSize = (long) (NUM_BYTES_1GB * cbfGB);
            long pkbfSize = (long) (NUM_BITS_1GB * pkbfGB);
            
            int allNumHash = Integer.parseInt(line.getOptionValue(optAllHash.getOpt(), "2"));
            String allNumHashStr = Integer.toString(allNumHash);
            int sbfNumHash = Integer.parseInt(line.getOptionValue(optSbfHash.getOpt(), allNumHashStr));
            int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), allNumHashStr));
            int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), allNumHashStr));
            int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), allNumHashStr));
            
            /**@TODO ensure that sbfNumHash and pkbfNumHash <= max(dbgbfNumHash, cbfNumHash) */
                        
//            int mismatchesAllowed = Integer.parseInt(line.getOptionValue(optMismatch.getOpt(), "5"));
            int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), "10"));
            int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), "1000"));
            int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), "500"));
            int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), "3"));
            int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), "5"));
            float maxCovGradient = Float.parseFloat(line.getOptionValue(optMaxCovGrad.getOpt(), "0.5"));
            float percentIdentity = Float.parseFloat(line.getOptionValue(optPercentIdentity.getOpt(), "0.9"));
            int maxIndelSize = Integer.parseInt(line.getOptionValue(optIndelSize.getOpt(), "1"));
            int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), "1"));
            int minTranscriptLength = Integer.parseInt(line.getOptionValue(optMinLength.getOpt(), "200"));
            boolean createFragmentDBG = false;
//            boolean createFragmentDBG = line.hasOption(optFdbg.getOpt());
            boolean useSingletonFragments = line.hasOption(optSingleton.getOpt());
            boolean extendFragments = line.hasOption(optExtend.getOpt());
            int minNumKmerPairs = Integer.parseInt(line.getOptionValue(optMinKmerPairs.getOpt(), "10"));
            String txptNamePrefix = line.getOptionValue(optOverlap.getOpt(), "");
            
//            boolean saveGraph = true;
//            boolean saveKmerPairs = true;

            System.out.println("Bloom filters     Memory (GB)");
            System.out.println("=============================");
            System.out.println("de Bruijn graph:  " + dbgGB);
            System.out.println("kmer counting:    " + cbfGB);
            System.out.println("paired kmers:     " + pkbfGB);
            System.out.println("screening:        " + sbfGB);
            System.out.println("=============================");
            System.out.println("Total:            " + (dbgGB+cbfGB+pkbfGB+sbfGB));
            
            System.out.println("name:    " + name);
            System.out.println("outdir:  " + outdir);
            
            File f = new File(outdir);
            if (!f.exists()) {
                f.mkdirs();
            }

            RNABloom assembler = new RNABloom(k, qDBG, qFrag);
            assembler.setParams(maxTipLen, lookahead, maxCovGradient, maxIndelSize, percentIdentity, minNumKmerPairs);

            try {
                touch(startedStamp);
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(1);
            }
            
            if (!forceOverwrite && dbgDoneStamp.exists()) {
                System.out.println("Loading graph from file `" + graphFile + "`...");
                assembler.restoreGraph(new File(graphFile));
            }
            else {                
                ArrayList<String> forwardFilesList = new ArrayList<>();
                ArrayList<String> backwardFilesList = new ArrayList<>();
                
                if (revCompLeft) {
                    backwardFilesList.addAll(Arrays.asList(fastqsLeft));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(fastqsLeft));
                }
                
                if (revCompRight) {
                    backwardFilesList.addAll(Arrays.asList(fastqsRight));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(fastqsRight));
                }
                                
                timer.start();
                
                assembler.initializeGraph(strandSpecific, 
                        dbgbfSize, cbfSize, pkbfSize, 
                        dbgbfNumHash, cbfNumHash, pkbfNumHash);
                
                assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads);

                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                System.out.println("Saving graph to file `" + graphFile + "`...");
                assembler.saveGraph(new File(graphFile));
                
                try {
                    touch(dbgDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
                        

            FastqPair[] fqPairs = new FastqPair[fastqsLeft.length];
            for (int i=0; i<fastqsLeft.length; ++i) {
                fqPairs[i] = new FastqPair(fastqsLeft[i], fastqsRight[i], revCompLeft, revCompRight);
            }

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
            
            String[] unconnectedReadsFastaPaths = {unconnectedReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
            
            String longSingletonsFastaPath = longFragmentsFastaPrefix + "01.fa";
            String shortSingletonsFastaPath = shortFragmentsFastaPrefix + "01.fa";
            String unconnectedSingletonsFastaPath = unconnectedReadsFastaPrefix + "01.fa";
                        
            if (!forceOverwrite && fragsDoneStamp.exists()) {
                int[] fragStats = assembler.restoreFragStatsFromFile(fragStatsFile);
                assembler.setPairedKmerDistance(fragStats[1]);
            }
            else {
                File fragmentsFile;
                
                for (String fragmentsFasta : longFragmentsFastaPaths) {
                    fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }
                
                for (String fragmentsFasta : shortFragmentsFastaPaths) {
                    fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }
                
                for (String fragmentsFasta : unconnectedReadsFastaPaths) {
                    fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }
                
                fragmentsFile = new File(longSingletonsFastaPath);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
                
                fragmentsFile = new File(shortSingletonsFastaPath);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
                
                fragmentsFile = new File(unconnectedSingletonsFastaPath);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
                
                timer.start();
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                
                int[] fragStats = assembler.assembleFragmentsMultiThreaded(fqPairs, 
                                                longFragmentsFastaPaths, 
                                                shortFragmentsFastaPaths,
                                                unconnectedReadsFastaPaths,
                                                longSingletonsFastaPath,
                                                shortSingletonsFastaPath,
                                                unconnectedSingletonsFastaPath,
                                                bound, 
                                                minOverlap,
                                                sampleSize,
                                                numThreads,
                                                maxErrCorrItr,
                                                extendFragments);
                
                assembler.writeFragStatsToFile(fragStats, fragStatsFile);
                
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                                
                try {
                    touch(fragsDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }

            /* Connect unconnected reads by doubling the kmer size */
            
            
            // double the kmer size
            assembler.setK(2*k);
            
            String[] unconnected2kReadsFastaPaths = {unconnected2kReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            unconnected2kReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            unconnected2kReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            unconnected2kReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            unconnected2kReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            unconnected2kReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            String unconnected2kSingletonsFastaPath = unconnected2kReadsFastaPrefix + "01.fa";

            if (!forceOverwrite && frags2kDoneStamp.exists()) {
                System.out.println("Loading graph from file `" + graph2kFile + "`...");
                assembler.restoreGraph(new File(graph2kFile));
                
                System.out.println("Restoring paired kmers Bloom filter from file...");
                assembler.restorePairedKmersBloomFilter(new File(graph2kFile));
            }
            else {
                // clear DBG-Bf, c-Bf, pk-Bf, a-Bf
                assembler.clearGraph();
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);

                // populate graph with kmers from reads
                ArrayList<String> forwardFilesList = new ArrayList<>();
                ArrayList<String> backwardFilesList = new ArrayList<>();

                if (revCompLeft) {
                    backwardFilesList.addAll(Arrays.asList(fastqsLeft));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(fastqsLeft));
                }

                if (revCompRight) {
                    backwardFilesList.addAll(Arrays.asList(fastqsRight));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(fastqsRight));
                }

                assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads);

                // repopulate with NEW kmers from fragments

                ArrayList<String> fragmentPaths = new ArrayList<>(longFragmentsFastaPaths.length + shortFragmentsFastaPaths.length);
                fragmentPaths.addAll(Arrays.asList(longFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(shortFragmentsFastaPaths));
                fragmentPaths.add(longSingletonsFastaPath);
                fragmentPaths.add(shortSingletonsFastaPath);

                assembler.repopulateGraph(fragmentPaths, strandSpecific);

                File fragmentsFile = new File(unconnected2kSingletonsFastaPath);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }

                for (String fragmentsFasta : unconnected2kReadsFastaPaths) {
                    fragmentsFile = new File(fragmentsFasta);
                    if (fragmentsFile.exists()) {
                        fragmentsFile.delete();
                    }
                }

                String[] allUnconnectedReads = new String[unconnectedReadsFastaPaths.length + 1];
                for (int i=0; i<unconnectedReadsFastaPaths.length; ++i) {
                    allUnconnectedReads[i] = unconnectedReadsFastaPaths[unconnectedReadsFastaPaths.length - i - 1];
                }
                allUnconnectedReads[unconnectedReadsFastaPaths.length] = unconnectedSingletonsFastaPath;

                assembler.rescueUnconnectedMultiThreaded(allUnconnectedReads, 
                                                    longFragmentsFastaPaths,
                                                    shortFragmentsFastaPaths,
                                                    unconnected2kReadsFastaPaths,
                                                    longSingletonsFastaPath,
                                                    shortSingletonsFastaPath,
                                                    unconnected2kSingletonsFastaPath,
                                                    bound,
                                                    minOverlap,
                                                    sampleSize, 
                                                    numThreads, 
                                                    maxErrCorrItr,
                                                    extendFragments);

                /* Save DBG-Bf and pk-Bf to disk */
                assembler.saveGraph(new File(graph2kFile));
                assembler.savePairedKmersBloomFilter(new File(graph2kFile));  
                
                /* Touch stamp */
                try {
                    touch(frags2kDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
            
//            String[] longTransfragsFastaPaths = {longTransfragsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
//                                            longTransfragsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
//                                            longTransfragsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
//                                            longTransfragsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
//                                            longTransfragsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
//                                            longTransfragsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
//            
//            String[] shortTransfragsFastaPaths = {shortTransfragsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
//                                            shortTransfragsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
//                                            shortTransfragsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
//                                            shortTransfragsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
//                                            shortTransfragsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
//                                            shortTransfragsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
//            
//            if (forceOverwrite || !txfgsDoneStamp.exists()) {
//                
//                for (String transfragsFasta : longTransfragsFastaPaths) {
//                    File transfragsFile = new File(transfragsFasta);
//                    if (transfragsFile.exists()) {
//                        transfragsFile.delete();
//                    }
//                }
//                
//                for (String transfragsFasta : shortTransfragsFastaPaths) {
//                    File transfragsFile = new File(transfragsFasta);
//                    if (transfragsFile.exists()) {
//                        transfragsFile.delete();
//                    }
//                }
//                                
//                timer.start();
//                
//                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
//                
//                assembler.assembleTransfragsMultiThreaded(longFragmentsFastaPaths, 
//                                                        shortFragmentsFastaPaths,
//                                                        unconnectedReadsFastaPaths,
//                                                        longTransfragsFastaPaths,
//                                                        shortTransfragsFastaPaths,
//                                                        numThreads,
//                                                        sampleSize);
//
//                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
//                
//                try {
//                    touch(txfgsDoneStamp);
//                } catch (Exception ex) {
//                    ex.printStackTrace();
//                }
//            }
            
            if (forceOverwrite || !txptsDoneStamp.exists()) {
                
                File transcriptsFile = new File(transcriptsFasta);
                if (transcriptsFile.exists()) {
                    transcriptsFile.delete();
                }
                
                File shortTranscriptsFile = new File(shortTranscriptsFasta);
                if (shortTranscriptsFile.exists()) {
                    shortTranscriptsFile.delete();
                }
                
                timer.start();
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                
                assembler.assembleTranscriptsMultiThreaded(longFragmentsFastaPaths, 
                                                            shortFragmentsFastaPaths,
                                                            unconnected2kReadsFastaPaths,
                                                            longSingletonsFastaPath,
                                                            shortSingletonsFastaPath,
                                                            unconnected2kSingletonsFastaPath,
                                                            transcriptsFasta, 
                                                            shortTranscriptsFasta,
                                                            graphFile,
                                                            numThreads,
                                                            sampleSize,
                                                            minTranscriptLength,
                                                            createFragmentDBG,
                                                            useSingletonFragments,
                                                            txptNamePrefix);

                System.out.println("Transcripts assembled in `" + transcriptsFasta + "`");
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                try {
                    touch(txptsDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
            
            
        }
        catch (ParseException exp) {
            System.out.println("ERROR:" + exp.getMessage() );
            System.exit(1);
        }
        
        System.out.println("Total Runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
    }
}
