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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
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
import rnabloom.graph.BloomFilterDeBruijnGraph.KmerIterator;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqPair;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.ReadPair;
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
    private boolean strandSpecific;
    private Pattern qualPatternDBG;
    private Pattern qualPatternFrag;
    private BloomFilterDeBruijnGraph graph = null;
    private BloomFilter screeningBf = null;
    
    private Random random;
    
    private final int bbLookahead = 5;
    private final int bbWindowSize = 10;
    private final int bbMaxIteration = 2;
    //private Function<String, Integer> findBackboneId;
    
    //private ArrayList<String> backbones = new ArrayList<>(10000);
    private int currentBackboneId = 0;
    private ConcurrentHashMap<String, Integer> kmerToBackboneID = new ConcurrentHashMap<>(10);
    private final int backboneHashKmerDistance = 100;

    private float dbgFPR = -1;
    private float covFPR = -1;
    
    public RNABloom(int k, int qDBG, int qFrag) {
        this.k = k;
        this.qualPatternDBG = getPhred33Pattern(qDBG, k);
        this.qualPatternFrag = getPhred33Pattern(qFrag, k);
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
            
            random = new Random(graph.getSeed());

            this.strandSpecific = graph.isStranded();

            dbgFPR = graph.getDbgbf().getFPR();
            covFPR = graph.getCbf().getFPR();
            
            /*
            if (strandSpecific) {
                findBackboneId = this::findBackboneIdStranded;
            }
            else {
                findBackboneId = this::findBackboneIdNonStranded;
            }
            */
            
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
                            if (screeningBf.lookupThenAdd(hashVals)) {
                                graph.add(hashVals);
                            }
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
                            long sbfNumBits,
                            long dbgbfNumBits,
                            long cbfNumBytes,
                            long pkbfNumBits,
                            int sbfNumHash,
                            int dbgbfNumHash,
                            int cbfNumHash,
                            int pkbfNumHash,
                            int seed,
                            int numThreads) {        
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            seed,
                                            k,
                                            strandSpecific);
        
        screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
        
        random = new Random(seed);
        
        this.strandSpecific = strandSpecific;
        
        /*
        if (strandSpecific) {
            findBackboneId = this::findBackboneIdStranded;
        }
        else {
            findBackboneId = this::findBackboneIdNonStranded;
        }
        */
        
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
            
            screeningBf.destroy();
            
        } catch (InterruptedException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        dbgFPR = graph.getDbgbf().getFPR();
        covFPR = graph.getCbf().getFPR();
    }

    public BloomFilterDeBruijnGraph getGraph() {
        return graph;
    }
    
    private boolean isPolyA(String left, String right) {
        return right.endsWith("AAAA") || (!graph.isStranded() && left.startsWith("TTTT"));
    }
        
    private boolean okToConnectPair(String left, String right, int minNumKmersPerReadNotAssembled) {
        NTHashIterator itr = graph.getHashIterator();
        itr.start(left);
        long[] hVals = itr.hVals;
        int numKmersNotSeenLeft = 0;

        while (itr.hasNext()) {
            itr.next();

            if (graph.getCount(hVals) == 0) {
                return false;
            }

            if (!graph.lookupFragmentKmer(hVals)) {
                ++numKmersNotSeenLeft;
            }
        }

        itr.start(right);
        int numKmersNotSeenRight = 0;

        while (itr.hasNext()) {
            itr.next();

            if (graph.getCount(hVals) == 0) {
                return false;
            }

            if (!graph.lookupFragmentKmer(hVals)) {
                ++numKmersNotSeenRight;
            }
        }

        if (numKmersNotSeenLeft >= minNumKmersPerReadNotAssembled || numKmersNotSeenRight >= minNumKmersPerReadNotAssembled) {
            return true;
        }

        return false;        
    }
    
    private int findBackboneIdNonStranded(String fragment) {
        KmerIterator itr = graph.new KmerIterator(fragment);
        Kmer seed = itr.next();
        Kmer kmer;
        String seq;
        
        while (itr.hasNext()) {
            kmer = itr.next();
            
            seq = smallestStrand(kmer.seq);
            
            if (kmerToBackboneID.containsKey(seq)) {
                return kmerToBackboneID.get(seq);
            }
            
            if (kmer.count > seed.count) {
                seed = kmer;
            }
        }
        
        ArrayList<Kmer> path = null;
        boolean randomSeed = false;
        for (int i=0; i<bbMaxIteration; ++i) {
            if (i>0) {
                if (randomSeed) {
                    seed = path.get(random.nextInt(path.size()));
                    randomSeed = false;
                }
                else {
                    seed = findMaxCoverageWindowKmer(path, graph, bbWindowSize);
                    randomSeed = true;
                }
            }

            /* greedy extend on both sides */
            HashSet<String> pathKmerStr = new HashSet<>(1000);
            pathKmerStr.add(smallestStrand(seed.seq));
            
            /* extend on right side */
            ArrayList<Kmer> rightPath = new ArrayList<>(1000);
            Kmer best = seed;
            while (true) {
                best = greedyExtendRightOnce(graph, best, bbLookahead);
                if (best != null) {
                    seq = smallestStrand(best.seq);
                    
                    if (kmerToBackboneID.containsKey(seq)) {
                        return kmerToBackboneID.get(seq);
                    }
                    
                    if (pathKmerStr.contains(seq)) {
                        break;
                    }
                    
                    pathKmerStr.add(seq);
                    rightPath.add(best);
                }
                else {
                    break;
                }
            }

            /* extend on left side */
            ArrayList<Kmer> leftPath = new ArrayList<>(1000);
            best = seed;
            while (true) {
                best = greedyExtendLeftOnce(graph, best, bbLookahead);
                if (best != null) {
                    seq = smallestStrand(best.seq);
                    
                    if (kmerToBackboneID.containsKey(seq)) {
                        return kmerToBackboneID.get(seq);
                    }
                    
                    if (pathKmerStr.contains(seq)) {
                        break;
                    }
                    
                    pathKmerStr.add(seq);
                    leftPath.add(best);
                }
                else {
                    break;
                }
            }

            Collections.reverse(leftPath);
            leftPath.add(seed);
            leftPath.addAll(rightPath);

            path = leftPath;
        }
        
        /* new backbone */
        //backbones.add(assemble(path));
        int id = ++currentBackboneId;
        
        /* store kmers in path */
        int numKmers = path.size();
        for (int i=0; i<numKmers; ++i) {
            if (i % backboneHashKmerDistance == 0) {
                kmerToBackboneID.put(smallestStrand(path.get(i).seq), i);
            }
        }
        
        return id;
    }
    
    private int findBackboneIdStranded(String fragment) {
        KmerIterator itr = graph.new KmerIterator(fragment);
        Kmer seed = itr.next();
        Kmer kmer;
        while (itr.hasNext()) {
            kmer = itr.next();
            
            if (kmerToBackboneID.containsKey(kmer.seq)) {
                return kmerToBackboneID.get(kmer.seq);
            }
            
            if (kmer.count > seed.count) {
                seed = kmer;
            }
        }
        
        ArrayList<Kmer> path = null;
        boolean randomSeed = false;
        for (int i=0; i<bbMaxIteration; ++i) {
            if (i>0) {
                if (randomSeed) {
                    seed = path.get(random.nextInt(path.size()));
                    randomSeed = false;
                }
                else {
                    seed = findMaxCoverageWindowKmer(path, graph, bbWindowSize);
                    randomSeed = true;
                }
            }

            /* greedy extend on both sides */
            HashSet<String> pathKmerStr = new HashSet<>(1000);
            pathKmerStr.add(seed.seq);
            
            /* extend on right side */
            ArrayList<Kmer> rightPath = new ArrayList<>(1000);
            Kmer best = seed;
            while (true) {
                best = greedyExtendRightOnce(graph, best, bbLookahead);
                if (best != null) {
                    String seq = best.seq;
                    
                    if (kmerToBackboneID.containsKey(seq)) {
                        return kmerToBackboneID.get(seq);
                    }
                    
                    if (pathKmerStr.contains(seq)) {
                        break;
                    }
                    
                    pathKmerStr.add(seq);
                    rightPath.add(best);
                }
                else {
                    break;
                }
            }

            /* extend on left side */
            ArrayList<Kmer> leftPath = new ArrayList<>(1000);
            best = seed;
            while (true) {
                best = greedyExtendLeftOnce(graph, best, bbLookahead);
                if (best != null) {
                    String seq = best.seq;
                    
                    if (kmerToBackboneID.containsKey(seq)) {
                        return kmerToBackboneID.get(seq);
                    }
                    
                    if (pathKmerStr.contains(seq)) {
                        break;
                    }
                    
                    pathKmerStr.add(seq);
                    leftPath.add(best);
                }
                else {
                    break;
                }
            }

            Collections.reverse(leftPath);
            leftPath.add(seed);
            leftPath.addAll(rightPath);

            path = leftPath;
            
            //System.out.println(">" + i + "\n" + assemble(path));
        }
        
        /* new backbone */
        //backbones.add(assemble(path));
        int id = ++currentBackboneId;;
        
        /* store kmers in path */
        
        //System.out.println(">bb\n" + assemble(path));
        
        int numKmers = path.size();
        for (int i=0; i<numKmers; ++i) {
            if (i % backboneHashKmerDistance == 0) {
                kmerToBackboneID.put(path.get(i).seq, id);
            }
        }
        
        return id;
    }
    
    public class Fragment {
        String left;
        String right;
        String seq;
        
        public Fragment(String left, String right, String seq) {
            this.left = left;
            this.right = right;
            this.seq = seq;
        }
    }
    
    public class FragmentAssembler implements Runnable {
        private ReadPair p;
        private List<Fragment> outList;
        //private int maxCorrectionsAllowed;
        private int bound;
        private int lookahead;
        private int minOverlap;
        private int maxTipLen;
        private boolean storeKmerPairs;
        private int maxIndelSize = 1; /** @TODO turn this into an option */
        private int minNumKmersPerReadNotAssembled; /** @TODO turn this into an option */
        private float maxCovGradient;
        private float covFPR;
        private int errorCorrectionIterations;
        private int leftReadLengthThreshold;
        private int rightReadLengthThreshold;
        
        public FragmentAssembler(ReadPair p,
                                List<Fragment> outList,
                                int bound, 
                                int lookahead, 
                                int minOverlap, 
                                int maxTipLen, 
                                boolean storeKmerPairs, 
                                float maxCovGradient, 
                                float covFPR,
                                int errorCorrectionIterations,
                                int leftReadLengthThreshold,
                                int rightReadLengthThreshold,
                                int minNumKmersPerReadNotAssembled) {
            
            this.p = p;
            this.outList = outList;
            this.bound = bound;
            this.lookahead = lookahead;
            this.minOverlap = minOverlap;
            this.maxTipLen = maxTipLen;
            this.storeKmerPairs = storeKmerPairs;
            this.maxCovGradient = maxCovGradient;
            this.covFPR = covFPR;
            this.errorCorrectionIterations = errorCorrectionIterations;
            this.leftReadLengthThreshold = leftReadLengthThreshold;
            this.rightReadLengthThreshold = rightReadLengthThreshold;
            this.minNumKmersPerReadNotAssembled = minNumKmersPerReadNotAssembled;
        }
        
        @Override
        public void run() {
            // connect segments of each read
            String left = connect(p.left, graph, k+p.numLeftBasesTrimmed+this.maxIndelSize, this.lookahead);
            String right = connect(p.right, graph, k+p.numRightBasesTrimmed+this.maxIndelSize, this.lookahead);

            if (left.length() >= this.leftReadLengthThreshold 
                    && right.length() >= this.rightReadLengthThreshold 
                    && okToConnectPair(left, right, this.minNumKmersPerReadNotAssembled)) {
                                
                String[] readPair = correctErrors2(left,
                                                    right,
                                                    graph, 
                                                    this.lookahead, 
                                                    this.maxIndelSize, 
                                                    this.maxCovGradient, 
                                                    this.covFPR,
                                                    this.errorCorrectionIterations);
                
                String leftCorrected = readPair[0];
                String rightCorrected = readPair[1];
                
                if (okToConnectPair(leftCorrected, rightCorrected, this.minNumKmersPerReadNotAssembled)) {
                    String fragment = overlapThenConnect(leftCorrected, rightCorrected, graph, this.bound, this.lookahead, this.minOverlap);
                    int fraglen = fragment.length();
                    if (fraglen > k) {
                        fragment = naiveExtend(fragment, graph, this.maxTipLen);

                        // mark fragment kmers as assembled
                        graph.addFragmentKmersFromSeq(fragment);
                        if (this.storeKmerPairs) {
                            graph.addPairedKmersFromSeq(fragment);
                        }

                        outList.add(new Fragment(leftCorrected, rightCorrected, fragment));
                    }
                }
                
                /*
                // correct each read
                String left = correctErrors(connectedLeft, graph, lookahead, maxCorrectionsAllowed);
                String right = correctErrors(connectedRight, graph, lookahead, maxCorrectionsAllowed);
                
                if (okToConnectPair(left, right, minNumKmersPerReadNotAssembled)) {

                    // assemble fragment from read pair
                    String fragment = overlapThenConnect(left, right, graph, bound, lookahead, minOverlap);

                    int fraglen = fragment.length();
                    if (fraglen > k) {
                        fragment = naiveExtend(correctErrors(fragment, graph, lookahead, Math.max(Math.round(fraglen*0.05f), maxCorrectionsAllowed)), graph, maxTipLen);
                        //fragment = naiveExtend(fragment, graph, maxTipLen);
                        
                        // mark fragment kmers as assembled
                        graph.addFragmentKmersFromSeq(fragment);
                        if (storeKmerPairs) {
                            graph.addPairedKmersFromSeq(fragment);
                        }

                        outList.add(fragment);
                    }
                }
                */
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
        
        int min, q1, median, q3, max;
        
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
    
    public void assembleFragmentsMultiThreaded(FastqPair[] fastqs, 
                                                String outfile, 
                                                int bound,
                                                int lookahead, 
                                                float maxCovGradient, 
                                                int minOverlap, 
                                                int maxTipLen, 
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
        
        int maxTasksQueueSize = numThreads*10;
        
        int newBound = bound;
        boolean pairedKmerDistanceIsSet = false;
        
        // set up thread pool
        MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);
        
        try {
            FastqReader lin, rin;
            FastqPairReader fqpr;
            
            FastaWriter out = new FastaWriter(outfile, true);
            
            ReadPair p;
            List<Fragment> fragments = Collections.synchronizedList(new ArrayList<>(sampleSize));
            
            boolean readLengthThresholdIsSet = false;
            ArrayList<Integer> leftReadLengths = new ArrayList<>(sampleSize);
            ArrayList<Integer> rightReadLengths = new ArrayList<>(sampleSize);
            
            int minNumKmersNotAssembled = k;
            
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
                        
                        // ignore read pairs when more than half of raw read length were trimmed for each read
                        if (p.originalLeftLength > 2*p.numLeftBasesTrimmed &&
                                p.originalRightLength > 2*p.numRightBasesTrimmed) {
                            
                            if (!readLengthThresholdIsSet) {
                                int best = 0;
                                for (String s : p.left) {
                                    if (s.length() > best) {
                                        best = s.length();
                                    }
                                }
                                leftReadLengths.add(best);

                                best = 0;
                                for (String s : p.right) {
                                    if (s.length() > best) {
                                        best = s.length();
                                    }
                                }
                                rightReadLengths.add(best);
                                    
                                if (leftReadLengths.size() >= sampleSize) {
                                    int[] leftReadLengthsStats = getMinQ1MedianQ3Max(leftReadLengths);
                                    int[] rightReadLengthsStats = getMinQ1MedianQ3Max(rightReadLengths);
                                    
                                    leftReadLengthThreshold = Math.max(k, leftReadLengthsStats[2] - (leftReadLengthsStats[3] - leftReadLengthsStats[1]) * 3/2);
                                    rightReadLengthThreshold = Math.max(k, rightReadLengthsStats[2] - (rightReadLengthsStats[3] - rightReadLengthsStats[1]) * 3/2);
                                    
                                    System.out.println("Left read length threshold:  " + leftReadLengthThreshold);
                                    System.out.println("Right read length threshold: " + rightReadLengthThreshold);
                                    
                                    readLengthThresholdIsSet = true;
                                    
                                    leftReadLengths = null;
                                    rightReadLengths = null;
                                }
                            }
                            
                            service.submit(new FragmentAssembler(p,
                                                                fragments,
                                                                bound,
                                                                lookahead,
                                                                minOverlap, 
                                                                maxTipLen, 
                                                                false, 
                                                                maxCovGradient, 
                                                                covFPR, 
                                                                maxErrCorrIterations, 
                                                                leftReadLengthThreshold,
                                                                rightReadLengthThreshold,
                                                                minNumKmersNotAssembled));
                            
                            if (fragments.size() >= sampleSize) {
                                break;
                            }
//                            else {
//                                System.out.println(readPairsParsed + " read pairs parsed");
//                                System.out.println(fragments.size() + " fragments assembled");
//                            }
                        }
                    }

                    service.terminate();

                    
                    // Calculate length stats
                    int numFragments = fragments.size();
                    ArrayList<Integer> fragLengths = new ArrayList<>(numFragments);
                    
                    for (Fragment frag : fragments) {
                        fragLengths.add(frag.seq.length());
                    }
                    
                    int[] fragLengthsStats = getMinQ1MedianQ3Max(fragLengths);
                    
                    // Set new bound for graph search
                    int medianFragLen = fragLengthsStats[2];
                    graph.setPairedKmerDistance(medianFragLen - k);
                    newBound = medianFragLen + (fragLengthsStats[3] - fragLengthsStats[1]) * 3 / 2; // 1.5*IQR
                    
                    System.out.println("Median fragment length:      " + medianFragLen);
                    System.out.println("Max graph traversal depth:   " + newBound);

                    // Store kmer pairs and write fragments to file
                    for (Fragment frag : fragments) {
                        if (frag.left.length() >= leftReadLengthThreshold && frag.right.length() >= rightReadLengthThreshold) {
                            graph.addPairedKmersFromSeq(frag.seq);
                            out.write(Long.toString(++fragmentId), frag.seq);
                        }
                    }                    
                    
                    pairedKmerDistanceIsSet = true;
                }
                
                service = new MyExecutorService(numThreads, maxTasksQueueSize);
                
                // assemble the remaining fragments in multi-threaded mode
                for (fragments = Collections.synchronizedList(new ArrayList<>(sampleSize)); fqpr.hasNext();) {
                    p = fqpr.next();
                    ++readPairsParsed;
                    
                    // ignore read pairs when more than half of raw read length were trimmed for each read
                    if (p.originalLeftLength > 2*p.numLeftBasesTrimmed &&
                            p.originalRightLength > 2*p.numRightBasesTrimmed) {
                        
                        service.submit(new FragmentAssembler(p,
                                                            fragments,
                                                            newBound, 
                                                            lookahead, 
                                                            minOverlap, 
                                                            maxTipLen, 
                                                            true, 
                                                            maxCovGradient, 
                                                            covFPR, 
                                                            maxErrCorrIterations, 
                                                            leftReadLengthThreshold, 
                                                            rightReadLengthThreshold,
                                                            minNumKmersNotAssembled));

                        if (fragments.size() >= sampleSize) {
                            service.terminate();

                            // write fragments to file
                            for (Fragment frag : fragments) {
                                out.write(Long.toString(++fragmentId), frag.seq);
                            }

                            // reset thread pool
                            fragments = Collections.synchronizedList(new ArrayList<>(sampleSize));
                            service = new MyExecutorService(numThreads, maxTasksQueueSize);
                        }
                    }
                }
                
                service.terminate();

                // write fragments to file
                for (Fragment frag : fragments) {
                    out.write(Long.toString(++fragmentId), frag.seq);
                }
                                
                lin.close();
                rin.close();
            }
            
            service.terminate();
            
            // store kmer pairs and write fragments to file
            for (Fragment frag : fragments) {
                out.write(Long.toString(++fragmentId), frag.seq);
            }
            
            out.close();
            
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
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
        
    public void assembleTranscripts(String fragmentsFasta, String outFasta, String tmpFasta, int lookAhead, float maxCovGradient, long sbfNumBits, int sbfNumHash) {
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }
        
        System.out.println("Assembling transcripts...");
        long numFragmentsParsed = 0;
        boolean append = false;
        
        try {
            FastaWriter fout = new FastaWriter(outFasta, append);

            System.out.println("Parsing fragments in `" + fragmentsFasta + "`...");
            
            FastaWriter tmpFout = new FastaWriter(tmpFasta, append);
                
            FastaReader fin = new FastaReader(fragmentsFasta);

            //System.out.println("Parsing `" + fileName + "`...");
            
            screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());
            
            NTHashIterator fragHashItr = graph.getHashFunction().getHashIterator(graph.getMaxNumHash());
            long[] fragHvals = fragHashItr.hVals;
            
            NTHashIterator txptHashItr = graph.getHashFunction().getHashIterator(sbfNumHash);
            long[] txptHvals = txptHashItr.hVals;
            
            long cid = 0;
            long tmpCid = 0;
            int minNumKmersNotAssembled = k;

            while (fin.hasNext()) {
                if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
                    System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
                }

                String fragment = fin.next();
                
                /** count assembled kmers */
                int numFragKmers = getNumKmers(fragment, k);
                float[] covs = new float[numFragKmers];
                fragHashItr.start(fragment);
                
                int numKmersNotAssembled = 0;
                for (int i=0; i<numFragKmers; ++i) {
                    covs[i] = graph.getCount(fragHvals);
                    fragHashItr.next();
                    if (!screeningBf.lookup(fragHvals)) {
                        ++numKmersNotAssembled;
                    }
                }
                
                if (numKmersNotAssembled >= minNumKmersNotAssembled) {
                
                    /** check whether sequence-wide coverage differences are too large */
                    boolean covDiffTooLarge = false;
                    Arrays.sort(covs);
                    float covLow = covs[0];
                    float covHigh;
                    for (int i=1; i<numFragKmers; ++i) {
                        covHigh = covs[i];
                        if (covHigh * maxCovGradient > covLow) {
                            covDiffTooLarge = true;
                            break;
                        }
                        covLow = covHigh;
                    }

                    if (covDiffTooLarge) {
                        // postpone extension to next round
                        tmpFout.write(Long.toString(++tmpCid), fragment);
                        continue;
                    }

                    String transcript = extendWithPairedKmersGreedily(fragment, graph, lookAhead);

    //                System.out.println(">f\n" + fragment + "\n>t\n" + transcript);

                    /** store assembled kmers */
                    txptHashItr.start(transcript);
                    while (txptHashItr.hasNext()) {
                        txptHashItr.next();
                        screeningBf.add(txptHvals);
                    }

                    fout.write(Long.toString(++cid) + " " + transcript.length() + " " + fragment, transcript);
                }
            }

            fin.close();
            tmpFout.close();
            
            System.out.println("Parsing fragments in `" + tmpFasta + "`...");
            
            /** evaluate and extend postponed fragments */
            numFragmentsParsed = 0;
            fin = new FastaReader(tmpFasta);
            while (fin.hasNext()) {
                if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
                    System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
                }
                
                String fragment = fin.next();
                
                /** count assembled kmers */
                int numFragKmers = getNumKmers(fragment, k);
                fragHashItr.start(fragment);

                boolean notAssembledEnough = false;
                int numKmersNotAssembled = 0;
                for (int i=0; i<numFragKmers; ++i) {
                    fragHashItr.next();
                    if (!screeningBf.lookup(fragHvals)) {
                        if (++numKmersNotAssembled >= minNumKmersNotAssembled) {
                            notAssembledEnough = true;
                            break;
                        }
                    }
                }
                
                if (notAssembledEnough) {
                    String bestAltPath = correctMismatches(fragment, graph, lookAhead, (int) Math.ceil(0.05 * fragment.length()));
                    
                    boolean bestAltPathAssembled = true;
                    numKmersNotAssembled = 0;
                    fragHashItr.start(bestAltPath);
                    for (int i=0; i<getNumKmers(bestAltPath, k); ++i) {                    
                        fragHashItr.next();
                        if (!screeningBf.lookup(fragHvals)) {
                            if (++numKmersNotAssembled >= minNumKmersNotAssembled){
                                bestAltPathAssembled = false;
                                break;
                            }
                        }
                    }

                    if (!bestAltPathAssembled) {
                        String transcript = extendWithPairedKmersNonGreedily(fragment, graph, lookAhead);

                        /** store assembled kmers */
                        txptHashItr.start(transcript);
                        while (txptHashItr.hasNext()) {
                            txptHashItr.next();
                            screeningBf.add(txptHvals);
                        }

                        fout.write(Long.toString(++cid) + " " + transcript.length() + " " + fragment, transcript);
                    }
                }
            }
            
            fin.close();
            fout.close();
            
            screeningBf.destroy();
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
        }
    }
        
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
        
        String[] readPair = correctErrors2(leftRead,
                                            rightRead,
                                            graph, 
                                            lookahead, 
                                            maxIndelSize, 
                                            maxCovGradient, 
                                            covFPR, 
                                            errorCorrectionIterations);
                                            
        System.out.println(overlapThenConnect(readPair[0], readPair[1], graph, bound, lookahead, minOverlap));
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
        String seq =  "GCCTGATCACCAGGGCTGCTTGTAACTCTGGGAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTAACCTCAACTACATGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATCCCTCCAAAATCAAGTGGGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAGAAGGCTGGGGCTCAT";
        
        int lookahead = 10;
        int errorsallowed = 5;
        
        String corrected = correctErrors(seq, graph, lookahead, errorsallowed);
        
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
        
        // -left /projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/gapdh
        // -left /home/gengar/test_data/GAPDH/GAPDH_2.fq.gz -right /home/gengar/test_data/GAPDH/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /home/gengar/test_assemblies/GAPDH

        // -dm 1 -cm 2.5 -pm 0.5 -left /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_2.fq -right /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_1.fq -revcomp-right -stranded -name SRR1360926 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/SRR1360926
        // -dm 3 -cm 7.5 -pm 1.5 -left /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_2.fq.gz -right /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_1.fq.gz -revcomp-right -stranded -name mcf7 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/mcf7
        
        // Based on: http://commons.apache.org/proper/commons-cli/usage.html
        CommandLineParser parser = new DefaultParser();

        Options options = new Options();

        Builder builder;
                
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
        
        builder = Option.builder("l");
        builder.longOpt("left");
        builder.desc("left reads file");
        builder.hasArg(true);
        builder.argName("FILE");
        builder.required(true);
        Option optLeftReads = builder.build();
        options.addOption(optLeftReads);
        
        builder = Option.builder("r");
        builder.longOpt("right");
        builder.desc("right reads file");
        builder.hasArg(true);
        builder.argName("FILE");
        builder.required(true);
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
        
        builder = Option.builder("s");
        builder.longOpt("seed");
        builder.desc("seed for random number generator and hash function");
        builder.hasArg(true);
        builder.argName("INT");
        Option optSeed = builder.build();
        options.addOption(optSeed);
        
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
        
//        builder = Option.builder("m");
//        builder.longOpt("mismatch");
//        builder.desc("max number of mismatch bases allowed per read");
//        builder.hasArg(true);
//        builder.argName("INT");
//        Option optMismatch = builder.build();
//        options.addOption(optMismatch);

        builder = Option.builder("tl");
        builder.longOpt("tiplength");
        builder.desc("max length of tip allowed");
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
        
        builder = Option.builder("e");
        builder.longOpt("errcorritr");
        builder.desc("max number of iterations of read error correction");
        builder.hasArg(true);
        builder.argName("INT");
        Option optErrCorrItr = builder.build();
        options.addOption(optErrCorrItr);        
        

        try {
            CommandLine line = parser.parse(options, args);
            
            int numThreads = Integer.parseInt(line.getOptionValue(optThreads.getOpt(), "2"));
            boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            String name = line.getOptionValue(optName.getOpt(), "rna-bloom");
            String outdir = line.getOptionValue(optOutdir.getOpt(), System.getProperty("user.dir") + File.separator + name + "_assembly");
            /**@TODO evaluate whether out dir is a valid dir */
            
            String fragmentsFasta = outdir + File.separator + name + ".fragments.fa";
            String transcriptsFasta = outdir + File.separator + name + ".transcripts.fa";
            String tmpFasta = outdir + File.separator + name + ".tmp.fa";
            String graphFile = outdir + File.separator + name + ".graph";
            
            File startedStamp = new File(outdir + File.separator + STARTED);
            File dbgDoneStamp = new File(outdir + File.separator + DBG_DONE);
            File fragsDoneStamp = new File(outdir + File.separator + FRAGMENTS_DONE);
            File txptsDoneStamp = new File(outdir + File.separator + TRANSCRIPTS_DONE);
            
            String fastqLeft = line.getOptionValue(optLeftReads.getOpt());
            String fastqRight = line.getOptionValue(optRightReads.getOpt());
            
            boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            boolean strandSpecific = line.hasOption(optStranded.getOpt());
            
            int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), "25"));
            int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), "3"));
            int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), "3"));
            int seed = Integer.parseInt(line.getOptionValue(optSeed.getOpt(), "689"));
            
            long sbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), "1")));
            long dbgbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), "1")));
            long cbfSize = (long) (NUM_BYTES_1GB * Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), "1")));
            long pkbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), "1")));
            
            int sbfNumHash = Integer.parseInt(line.getOptionValue(optSbfHash.getOpt(), "1"));
            int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), "3"));
            int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), "4"));
            int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), "1"));
            
            /**@TODO ensure that sbfNumHash and pkbfNumHash <= max(dbgbfNumHash, cbfNumHash) */
                        
//            int mismatchesAllowed = Integer.parseInt(line.getOptionValue(optMismatch.getOpt(), "5"));
            int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), "10"));
            int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), "1000"));
            int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), "500"));
            int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), "7"));
            int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), "10"));
            float maxCovGradient = Float.parseFloat(line.getOptionValue(optMaxCovGrad.getOpt(), "0.5"));
            int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), "2"));
            
            boolean saveGraph = true;
            boolean saveKmerPairs = true;

            System.out.println("name: " + name);
            System.out.println("outdir: " + outdir);
            
            File f = new File(outdir);
            if (!f.exists()) {
                f.mkdirs();
            }

            RNABloom assembler = new RNABloom(k, qDBG, qFrag);

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
                        sbfSize, dbgbfSize, cbfSize, pkbfSize, 
                        sbfNumHash, dbgbfNumHash, cbfNumHash, pkbfNumHash,
                        seed,
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
            
            /**@TODO */
//            assembler.testErrorCorrection();
//            System.exit(0);
            
            FastqPair fqPair = new FastqPair(fastqLeft, fastqRight, revCompLeft, revCompRight);
            FastqPair[] fqPairs = new FastqPair[]{fqPair};        

            if (!forceOverwrite && fragsDoneStamp.exists()) {
                System.out.println("Restoring paired kmers Bloom filter from file...");
                assembler.restorePairedKmersBloomFilter(new File(graphFile));            
            }
            else {
                File fragmentsFile = new File(fragmentsFasta);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
                
                long startTime = System.nanoTime();
                
                assembler.assembleFragmentsMultiThreaded(fqPairs, fragmentsFasta, bound, lookahead, maxCovGradient, minOverlap, maxTipLen, sampleSize, numThreads, maxErrCorrItr);

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
                
                File tmpFile = new File(tmpFasta);
                if (tmpFile.exists()) {
                    tmpFile.delete();
                }
                
                long startTime = System.nanoTime();
                
                assembler.assembleTranscripts(fragmentsFasta, transcriptsFasta, tmpFasta, lookahead, maxCovGradient, sbfSize, sbfNumHash);

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
