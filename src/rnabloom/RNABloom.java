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
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqPair;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastqPairReader.ReadPair;
import rnabloom.io.FastqReader;
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
    private Pattern qualPattern;
    private BloomFilterDeBruijnGraph graph = null;
    
    private Random random;
    
    private final int bbLookahead = 5;
    private final int bbWindowSize = 10;
    private final int bbMaxIteration = 2;
    private Function<String, Integer> findBackboneId;
    
    private ArrayList<String> backbones = new ArrayList<>(10000);
    private HashMap<String, Integer> kmerToBackboneID = new HashMap<>(10000);
    private int backboneHashKmerDistance = 200;
        
    public RNABloom(int k, int q) {
        this.k = k;
        this.qualPattern = getPhred33Pattern(q, k);
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
            int highCov = 16;

            int numCov1Kmers = 0;
            int numHighCovKmers = 0; // >
            
            int numLeftAssembled = 0;
            float c;
            for (String s : leftKmers) {
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
            }

            if (numCov1Kmers > 0 && numHighCovKmers >= k) {
                return false;
            }
            
            int numRightAssembled = 0;
            for (String s : rightKmers) {
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
            }
            
            if (numCov1Kmers > 0 && numHighCovKmers >= k) {
                return false;
            }

            if ((numLeftAssembled <= threshold && leftKmers.length - numLeftAssembled > 0) || (numRightAssembled <= threshold && rightKmers.length - numRightAssembled > 0)) {
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
        backbones.add(assemble(path));
        int id = backbones.size() - 1;
        
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
        backbones.add(assemble(path));
        int id = backbones.size() - 1;
        
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
        System.out.println("Assembling fragments ...");
        
        graph.initializePairKmersBloomFilter();
        
        long readPairsParsed = 0;
        int correctionWorks = 0;
        
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

                FastqPairReader fqpr = new FastqPairReader(lin, rin, qualPattern, fqPair.leftRevComp, fqPair.rightRevComp);
                System.out.println("Parsing `" + fqPair.leftFastq + "` and `" + fqPair.rightFastq + "`...");
                
                ReadPair p;
                String rawLeft, rawRight;
                while (fqpr.hasNext()) {
                    p = fqpr.next();
                    
                    if (++readPairsParsed % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs...");
                    }                    
                    
                    if (okToAssemble(p)) {
                        rawLeft = p.left;
                        rawRight = p.right;
                        
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
                        p.left = correctMismatches(rawLeft, graph, lookahead, mismatchesAllowed);
                        p.right = correctMismatches(rawRight, graph, lookahead, mismatchesAllowed);
                        
                        if (okToAssemble(p)) {
                            String fragment = assembleFragment(p.left, p.right, graph, bound, lookahead, minOverlap);
                                                        
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
                                
                                out = new FastaWriter(outDir + File.separator + backboneId + ".fa", append);
                                
                                //float minCov = graph.getMinKmerCoverage(fragment);
                                
                                /** extend on both ends unambiguously*/
                                fragment = naiveExtend(fragment, graph, maxTipLen);

                                out.write(Integer.toString(cid) + " " + rawLeft + " " + rawRight, fragment);

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

                                    System.out.println("longest fragment allowed: " + bound);

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
                                
                                out.close();
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
    }
    
    public void assembleTranscripts(String fragsDirPath, String outFasta, int lookAhead, float covGradient) {
        long numFragmentsParsed = 0;
        boolean append = false;
        FilenameFilter fragsFilenameFilter = (File dir, String name) -> name.endsWith(".fa");
        
        try {
            File fragsDir = new File(fragsDirPath);
            System.out.println("Looking for fragments in `" + fragsDir.getName() + "`...");
            
            FastaWriter fout = new FastaWriter(outFasta, append);
            
            KmerSet assembledKmers = new KmerSet(strandSpecific);
            FragmentComparator fragComp = new FragmentComparator();
            
            for (File inFasta : fragsDir.listFiles(fragsFilenameFilter)) {
                FastaReader fin = new FastaReader(inFasta);
                String fileName = inFasta.getName();
                
                System.out.println("Parsing `" + fileName + "`...");
                
                String clusterId = fileName.substring(0, fileName.indexOf("."));
                int cid = 0;
                
                /** sort fragments by min kmer count */
                
                ArrayList<Fragment> frags = new ArrayList<>(100);
                
                while (fin.hasNext()) {
                    if (++numFragmentsParsed % NUM_PARSED_INTERVAL == 0) {
                        System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments...");
                    }

                    frags.add(new Fragment(fin.next()));
                }

                fin.close();
                
                Collections.sort(frags, fragComp);
                
                for (Fragment fragment : frags) {
                    String seq = fragment.seq;
                    
                    for (String kmer : kmerize(seq, k)) {
                        if (!assembledKmers.contains(kmer)) {
                            
                            String transcript = assembleTranscript(seq, graph, lookAhead, covGradient, assembledKmers);
                            fout.write("c" + clusterId + "_t" + Integer.toString(++cid) + " " + seq, transcript);

                            for (String kmer2 : kmerize(transcript, k)) {
                                assembledKmers.add(kmer2);
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
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
                
        ///*
        String fastq1 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz"; //right
        String fastq2 = "/projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz"; //left        
        String fragsDirPath = "/projects/btl2/kmnip/rna-bloom/tests/java_assemblies/fragments";
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
        RNABloom assembler = new RNABloom(k, q);
        
        /*
        assembler.createGraph(forwardFastqs, backwardFastqs, strandSpecific, dbgbfSize, cbfSize, pkbfSize, dbgbfNumHash, cbfNumHash, pkbfNumHash, seed);
        
        
        System.out.println("Saving graph to file `" + graphFile + "`...");
        assembler.saveGraph(new File(graphFile));
        */
        
        System.out.println("Loading graph from file `" + graphFile + "`...");
        assembler.restoreGraph(new File(graphFile));
        
        //assembler.test2();
        
        ///*
        FastqPair fqPair = new FastqPair(fastq2, fastq1, revComp2, revComp1);
        FastqPair[] fqPairs = new FastqPair[]{fqPair};
        
        //assembler.test3();
        
        File fragsDir = new File(fragsDirPath);
        if (!fragsDir.exists()) {
            fragsDir.mkdirs();
        }
        
        /*
        assembler.assembleFragments(fqPairs, fragsDirPath, mismatchesAllowed, bound, lookahead, minOverlap, maxTipLen, sampleSize);
        
        System.out.println("Saving paired kmers Bloom filter to file...");
        assembler.savePairedKmersBloomFilter(new File(graphFile));
        */
        
        System.out.println("Restoring paired kmers Bloom filter from file...");
        assembler.restorePairedKmersBloomFilter(new File(graphFile));
        
        //assembler.test2();
        assembler.assembleTranscripts(fragsDirPath, transcriptsFasta, 10, 0.1f);
        //*/
    }
    
}
