/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import static java.lang.Math.pow;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
import rnabloom.bloom.hash.CanonicalNTHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.bloom.hash.ReverseComplementNTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
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
    
    private Random random;
    
    private final int bbLookahead = 5;
    private final int bbWindowSize = 10;
    private final int bbMaxIteration = 2;
    private Function<String, Integer> findBackboneId;
    
    //private ArrayList<String> backbones = new ArrayList<>(10000);
    private int currentBackboneId = 0;
    private HashMap<String, Integer> kmerToBackboneID = new HashMap<>(100000);
    private final int backboneHashKmerDistance = 100;
        
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

            if (strandSpecific) {
                findBackboneId = this::findBackboneIdStranded;
            }
            else {
                findBackboneId = this::findBackboneIdNonStranded;
            }
            
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
                itr = new NTHashIterator(k, numHash);
            }
            else {
                itr = new ReverseComplementNTHashIterator(k, graph.getMaxNumHash());
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
                    graph.add(hashVals);
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
                    itr = new NTHashIterator(k, numHash);
                }
                else {
                    itr = new ReverseComplementNTHashIterator(k, graph.getMaxNumHash());
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
                            graph.add(hashVals);
                        }
                    }
                }
                fr.close();
                
                successful = true;
                System.out.println("[" + id + "] Parsed " + NumberFormat.getInstance().format(numReads) + " reads...");
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
    
    public void createGraph(String[] forwardFastqs, String[] reverseFastqs, boolean strandSpecific, long dbgbfNumBits, long cbfNumBytes, long pkbfNumBits, int dbgbfNumHash, int cbfNumHash, int pkbfNumHash, int seed) {        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            seed,
                                            k,
                                            strandSpecific);
        
        random = new Random(seed);
        
        this.strandSpecific = strandSpecific;
        
        if (strandSpecific) {
            findBackboneId = this::findBackboneIdStranded;
        }
        else {
            findBackboneId = this::findBackboneIdNonStranded;
        }
        
        /** parse the reads */
        
        int numReads = 0;
        int numHash = graph.getMaxNumHash();
        
        ExecutorService service = Executors.newFixedThreadPool(2);
        
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
        } catch (InterruptedException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public BloomFilterDeBruijnGraph getGraph() {
        return graph;
    }
    
    private boolean isPolyA(String left, String right) {
        return right.endsWith("AAAA") || (!graph.isStranded() && left.startsWith("TTTT"));
    }
    
    private boolean okToConnectPair(String left, String right) {
        if (left.length() >= k && right.length() >= k) {
            int threshold = 2*k-1;
            int highCov = 16;

            int numCov1Kmers = 0;
            int numHighCovKmers = 0; // >
            float c;
            
            int numLeftAssembled = 0;
            
            KmerIterator leftItr = new KmerIterator(left, k);
            while (leftItr.hasNext()) {
                String s = leftItr.next();
                if (graph.lookupFragmentKmer(s)) {
                    ++numLeftAssembled;
                }
                
                c = graph.getCount(s);
                if (c == 1) {
                    ++numCov1Kmers;
                }
                else if (c > highCov) {
                    ++numHighCovKmers;
                }
                
                if (numCov1Kmers > 0 && numHighCovKmers >= k) {
                    return false;
                }
            }
                        
            int numRightAssembled = 0;
            KmerIterator rightItr = new KmerIterator(right, k);
            while (rightItr.hasNext()) {
                String s = rightItr.next();
                
                if (graph.lookupFragmentKmer(s)) {
                    ++numRightAssembled;
                }
                
                c = graph.getCount(s);
                if (c == 1) {
                    ++numCov1Kmers;
                }
                else if (c > highCov) {
                    ++numHighCovKmers;
                }
                
                if (numCov1Kmers > 0 && numHighCovKmers >= k) {
                    return false;
                }
            }

            if ((numLeftAssembled <= threshold && leftItr.numKmers - numLeftAssembled > 0) ||
                    (numRightAssembled <= threshold && rightItr.numKmers - numRightAssembled > 0)) {
                return true;
            }
        }

        return false;        
    }
    
    private int findBackboneIdNonStranded(String fragment) {
        ArrayList<Kmer> fragKmers = graph.getKmers(fragment);
        Kmer seed = fragKmers.get(0);
        for (Kmer kmer : fragKmers) {
            String seq = smallestStrand(kmer.seq);
            
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
                    String seq = smallestStrand(best.seq);
                    
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
                    String seq = smallestStrand(best.seq);
                    
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
        ArrayList<Kmer> fragKmers = graph.getKmers(fragment);
        Kmer seed = fragKmers.get(0);
        for (Kmer kmer : fragKmers) {
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
        }
        
        /* new backbone */
        //backbones.add(assemble(path));
        int id = ++currentBackboneId;;
        
        /* store kmers in path */
        int numKmers = path.size();
        for (int i=0; i<numKmers; ++i) {
            if (i % backboneHashKmerDistance == 0) {
                kmerToBackboneID.put(path.get(i).seq, id);
            }
        }
        
        return id;
    }
    
    public void assembleFragments(FastqPair[] fastqs, String outDir, int mismatchesAllowed, int bound, int lookahead, int minOverlap, int maxTipLen, int sampleSize) {
        System.out.println("Assembling fragments...");
        
        graph.initializePairKmersBloomFilter();
        
        long readPairsParsed = 0;
        
        ArrayList<Integer> clustersMaxContigId = new ArrayList<>(10000);
        
        try {
            FastqReader lin, rin;
            FastaWriter out;
            ArrayList<String> sampleFragments = new ArrayList<>(sampleSize);
            //ArrayList<Float> coverageGradients = new ArrayList<>(2*sampleSize*lookahead);
            int[] fragmentLengths = new int[sampleSize];
            int fid = 0;
            
            for (FastqPair fqPair: fastqs) {
                lin = new FastqReader(fqPair.leftFastq, true);
                rin = new FastqReader(fqPair.rightFastq, true);

                FastqPairReader fqpr = new FastqPairReader(lin, rin, qualPatternFrag, fqPair.leftRevComp, fqPair.rightRevComp);
                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");
                
                ReadPair p;
                String rawLeft, rawRight;
                while (fqpr.hasNext()) {
                    p = fqpr.next();
                    
                    if (++readPairsParsed % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
                    }                    

                    rawLeft = connect(p.left, graph, k+p.numLeftBasesTrimmed+1, lookahead);
                    rawRight = connect(p.right, graph, k+p.numRightBasesTrimmed+1, lookahead);
                    
                    if (okToConnectPair(rawLeft, rawRight)) {
                        /*
                        if (fid < sampleSize) {
                            if (p.left.length() > k+lookahead*2*2) {
                                for (Float g : coverageGradients(p.left, graph, lookahead)) {
                                    coverageGradients.add(g);
                                }
                            }
                            
                            if (p.right.length() > k+lookahead*2*2) {
                                for (Float g : coverageGradients(p.right, graph, lookahead)) {
                                    coverageGradients.add(g);
                                }
                            }
                        }
                        */
                        
                        // correct individual reads
                        String left = correctMismatches(rawLeft, graph, lookahead, mismatchesAllowed);
                        String right = correctMismatches(rawRight, graph, lookahead, mismatchesAllowed);
                        
                        if (okToConnectPair(left, right)) {
                            String fragment = overlapThenConnect(left, right, graph, bound, lookahead, minOverlap);
                                                        
                            // correct fragment
                            fragment = correctMismatches(fragment, graph, lookahead, mismatchesAllowed);
                            int fragLen = fragment.length();

                            if (fragLen > k) {
                                int backboneId = findBackboneId.apply(fragment);
                                
                                int cid = 0;
                                boolean append = true;
                                if (backboneId >= clustersMaxContigId.size()) {
                                    clustersMaxContigId.add(0);
                                    append = false;
                                }
                                else {
                                    cid = clustersMaxContigId.get(backboneId)+1;
                                    clustersMaxContigId.set(backboneId, cid);
                                }
                                
                                //float minCov = graph.getMinKmerCoverage(fragment);
                                
                                /** extend on both ends unambiguously*/
                                fragment = naiveExtend(fragment, graph, maxTipLen);

                                out = new FastaWriter(outDir + File.separator + backboneId + ".fa", append);
                                out.write(Integer.toString(cid) + " " + rawLeft + " " + rawRight, fragment);
                                out.close();
                                
                                ++fid;
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

                                    System.out.println("Max graph traversal depth: " + bound);

                                    /** clear sample fragment lengths */
                                    fragmentLengths = null;

                                    /** store paired kmers of all sample fragments */
                                    for (String frag : sampleFragments) {
                                        graph.addPairedKmersFromSeq(frag);
                                    }

                                    /** clear sample fragments */
                                    sampleFragments = null;
                                    
                                    /*
                                    Collections.sort(coverageGradients);
                                    int cgSize = coverageGradients.size();
                                    float cgIqr15 = (coverageGradients.get(cgSize*3/4) - coverageGradients.get(cgSize/4)) * 3/2;
                                    float cgMedian = coverageGradients.get(cgSize/2);
                                    System.out.println("median cov gradient: " + cgMedian + "+/-" + cgIqr15);
                                    coverageGradients = null;
                                    */
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
                }

                lin.close();
                rin.close();
                kmerToBackboneID = null;
            }
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
            System.out.println("Assembled fragments in " + currentBackboneId + " clusters...");
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
            graph.destroyPkbf();
            graph.restorePkbf(graphFile);
        } catch (IOException ex) {
            Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public class Fragment {
        public String seq;
        public float minCov;
        public float medCov;
        public float maxCov;
        
        public Fragment(String seq) {
            this.seq = seq;
            float[] c = graph.getMinMedianMaxKmerCoverage(seq);
            this.minCov = c[0];
            this.medCov = c[1];
            this.maxCov = c[2];
        }
    }

    public class FragmentComparator implements Comparator<Fragment> {
        @Override
        public int compare(Fragment f1, Fragment f2) {
            int compareMin = Float.compare(f1.minCov, f2.minCov);
            
            if (compareMin < 0) {
                return 1;
            }
            else if (compareMin == 0) {
                int compareMed = Float.compare(f1.medCov, f2.medCov);
                
                if (compareMed < 0) {
                    return 1;
                }
                else if (compareMed == 0) {
                    int compareMax = Float.compare(f1.maxCov, f2.maxCov);
                    
                    if (compareMax < 0) {
                        return 1;
                    }
                    else if (compareMax == 0) {
                        return 0;
                    }
                    else {
                        return -1;
                    }
                }
                else {
                    return -1;
                }
            }
            else {
                return -1;
            }
        }
    }    
    
    public class KmerSet {
        private final HashSet<String> set;
        private final Consumer<String> addFunction;
        private final Function<String, Boolean> containsFunction;
        
        public KmerSet(boolean stranded) {
            set = new HashSet<>();
            
            if (stranded) {
                this.addFunction = this::addStranded;
                this.containsFunction = this::containsStranded;
            }
            else {
                this.addFunction = this::addNonStranded;
                this.containsFunction = this::containsNonStranded;
            }
        }

        public void add(String s) {
            addFunction.accept(s);
        }
        
        private void addStranded(String s) {
            set.add(s);
        }

        private void addNonStranded(String s) {
            set.add(smallestStrand(s));
        }
        
        public boolean contains(String s) {
            return containsFunction.apply(s);
        }
        
        private boolean containsStranded(String s) {
            return set.contains(s);
        }
        
        private boolean containsNonStranded(String s) {
            return set.contains(smallestStrand(s));
        }
        
        public void clear() {
            set.clear();
        }
    }
        
    public void assembleTranscripts(String fragsDirPath, String outFasta, int lookAhead, float covGradient) {
        System.out.println("Assembling transcripts...");
        long numFragmentsParsed = 0;
        boolean append = false;
        FilenameFilter fragsFilenameFilter = (File dir, String name) -> name.endsWith(".fa");
        
        try {
            File fragsDir = new File(fragsDirPath);
            FastaWriter fout = new FastaWriter(outFasta, append);
            
            FragmentComparator fragComp = new FragmentComparator();
            
            File[] fragFiles = fragsDir.listFiles(fragsFilenameFilter);
            System.out.println("Found " + fragFiles.length + " fragment clusters...");
            System.out.println("Parsing fragments *.fa in `" + fragsDir.getName() + "`...");
            
            for (File inFasta : fragFiles) {
                ArrayList<Fragment> frags = new ArrayList<>();
                KmerSet assembledKmers = new KmerSet(strandSpecific);
                
                FastaReader fin = new FastaReader(inFasta);
                String fileName = inFasta.getName();
                
                //System.out.println("Parsing `" + fileName + "`...");
                
                String clusterId = fileName.substring(0, fileName.indexOf("."));
                int cid = 0;
                
                /** sort fragments by min kmer count */
                while (fin.hasNext()) {
                    if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
                    }

                    frags.add(new Fragment(fin.next()));
                }

                fin.close();
                
                Collections.sort(frags, fragComp);
                KmerIterator fragKmerItr = new KmerIterator(k);
                KmerIterator txptKmerItr = new KmerIterator(k);
                for (Fragment fragment : frags) {
                    String fragmentSeq = fragment.seq;
                    fragKmerItr.initialize(fragmentSeq);
                    
                    while (fragKmerItr.hasNext()) {
                        if (!assembledKmers.contains(fragKmerItr.next())) {
                            
                            String transcript = assembleTranscript(fragmentSeq, graph, lookAhead, covGradient, assembledKmers);
                            fout.write("c" + clusterId + "_t" + Integer.toString(++cid), transcript);
                            
                            txptKmerItr.initialize(transcript);
                            while (txptKmerItr.hasNext()) {
                                assembledKmers.add(txptKmerItr.next());
                            }

                            break;
                        }
                    }
                }
                
            }
            
            fout.close();
        } catch (IOException ex) {
            //Logger.getLogger(RNABloom.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
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
    
    public void test2() {
        String seq = "CTGCTGACGCCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCAT";
        
        String corrected = correctMismatches(seq, graph, 5, 5);
        
        System.out.println(seq);
        System.out.println(corrected);
    }
    
    public void test() {
        String f = "GCCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTTGTCAAGCTCATTTCCTGGTATGACAACGAATTTGGCTACAGCAACAGGGTGGTGGACCTCATGGCCCACATGGCCTCCAAGGAGTAAGACCCCTGGACCACCAGCCCCAGCAAGAGCACAAGAGGAAGAGAGAGACCCTCACTGCTGGGGAGTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAAAAAAAAAAAAAAAAAAAAAAAAA";
        String transcript = assembleTranscript(f, graph, 10, 0.1f, new KmerSet(strandSpecific));
        System.out.print(transcript);
        for (Kmer kmer : graph.getKmers(transcript)) {
            System.out.print(kmer.count);
            System.out.print(" ");
        }
    }
    
    public static void touch(File f) throws IOException {
        f.getParentFile().mkdirs();
        if (!f.createNewFile()){
            f.setLastModified(System.currentTimeMillis());
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
        
        builder = Option.builder("m");
        builder.longOpt("mismatch");
        builder.desc("max number of mismatch bases allowed per read");
        builder.hasArg(true);
        builder.argName("INT");
        Option optMismatch = builder.build();
        options.addOption(optMismatch);

        builder = Option.builder("t");
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
        
        

        try {
            CommandLine line = parser.parse(options, args);
            
            boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            String name = line.getOptionValue(optName.getOpt(), "rna-bloom");
            String outdir = line.getOptionValue(optOutdir.getOpt(), System.getProperty("user.dir") + File.separator + name + "_assembly");
            
            String fragsDirPath = outdir + File.separator + "fragments";
            String transcriptsFasta = outdir + File.separator + name + ".transcripts.fa";
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
            
            long dbgbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), "1")));
            long cbfSize = (long) (NUM_BYTES_1GB * Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), "1")));
            long pkbfSize = (long) (NUM_BITS_1GB * Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), "1")));
            
            int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), "3"));
            int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), "4"));
            int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), "1"));
                        
            int mismatchesAllowed = Integer.parseInt(line.getOptionValue(optMismatch.getOpt(), "5"));
            int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), "10"));
            int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), "1000"));
            int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), "500"));
            int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), "5"));
            int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), "10"));
            
            boolean saveGraph = false;
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
                                
                assembler.createGraph(forwardFastqs, backwardFastqs, strandSpecific, dbgbfSize, cbfSize, pkbfSize, dbgbfNumHash, cbfNumHash, pkbfNumHash, seed);
                
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
            System.exit(0);
            
            FastqPair fqPair = new FastqPair(fastqLeft, fastqRight, revCompLeft, revCompRight);
            FastqPair[] fqPairs = new FastqPair[]{fqPair};        

            File fragsDir = new File(fragsDirPath);
            if (!fragsDir.exists()) {
                fragsDir.mkdirs();
            }

            if (!forceOverwrite && fragsDoneStamp.exists()) {
                System.out.println("Restoring paired kmers Bloom filter from file...");
                assembler.restorePairedKmersBloomFilter(new File(graphFile));            
            }
            else {
                long startTime = System.nanoTime();
                
                assembler.assembleFragments(fqPairs, fragsDirPath, mismatchesAllowed, bound, lookahead, minOverlap, maxTipLen, sampleSize);

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
                assembler.assembleTranscripts(fragsDirPath, transcriptsFasta, 10, 0.1f);

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
    }
    
}
