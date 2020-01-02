/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.pow;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.PairedKeysBloomFilter;
import rnabloom.bloom.hash.CanonicalNTHashIterator;
import rnabloom.bloom.hash.CanonicalPairedNTHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.bloom.hash.PairedNTHashIterator;
import rnabloom.bloom.hash.ReverseComplementNTHashIterator;
import rnabloom.bloom.hash.ReverseComplementPairedNTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer;
import rnabloom.io.FastaFilteredSequenceIterator;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqFilteredSequenceIterator;
import rnabloom.io.FastxFilePair;
import rnabloom.io.PairedReadSegments;
import rnabloom.io.FastqReader;
import rnabloom.io.FastqRecord;
import rnabloom.io.FastxPairSequenceIterator;
import rnabloom.io.FastxSequenceIterator;
import rnabloom.io.FileFormatException;
import static rnabloom.olc.OverlapLayoutConcensus.hasMinimap2;
import static rnabloom.olc.OverlapLayoutConcensus.hasRacon;
import static rnabloom.olc.OverlapLayoutConcensus.overlapLayout;
import static rnabloom.olc.OverlapLayoutConcensus.overlapLayoutConcensus;
import rnabloom.util.GraphUtils;
import static rnabloom.util.GraphUtils.*;
import rnabloom.util.NTCardHistogram;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author Ka Ming Nip
 */
public class RNABloom {
    public final static String VERSION = "1.2.1";
    
//    private final static long NUM_PARSED_INTERVAL = 100000;
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    public final static long NUM_BYTES_1MB = (long) pow(1024, 2);
    public final static long NUM_BYTES_1KB = (long) 1024;
    
    private final static String FASTA_EXT = ".fa";
    
    private boolean debug = false;
    private int k;
    private boolean strandSpecific;
    private Pattern seqPattern;
    private Pattern qualPatternDBG;
    private Pattern qualPatternFrag;
    private Pattern polyATailPattern;
    private Pattern polyATailOnlyPattern;
    private Pattern polyATailOnlyMatchingPattern;
    private Pattern polyTHeadOnlyPattern;
    private Pattern polyASignalPattern;
    private BloomFilterDeBruijnGraph graph = null;
    private BloomFilter screeningBf = null;

    private int maxTipLength;
    private int lookahead;
    private float maxCovGradient;
    private int maxIndelSize;
    private int minPolyATailLengthRequired;
    private float percentIdentity;
    private int minNumKmerPairs;
    private int longFragmentLengthThreshold = -1;
    
    private int qDBG = -1;
    private int qFrag = -1;
    
    private float dbgFPR = -1;
    private float covFPR = -1;
    
    private final static String STRATUM_01 = "01";
    private final static String STRATUM_E0 = "e0";
    private final static String STRATUM_E1 = "e1";
    private final static String STRATUM_E2 = "e2";
    private final static String STRATUM_E3 = "e3";
    private final static String STRATUM_E4 = "e4";
    private final static String STRATUM_E5 = "e5";
    private final static String[] STRATA = new String[]{STRATUM_01, STRATUM_E0, STRATUM_E1, STRATUM_E2, STRATUM_E3, STRATUM_E4, STRATUM_E5};
    private final static String[] COVERAGE_ORDER = {STRATUM_E0, STRATUM_E1, STRATUM_E2, STRATUM_E3, STRATUM_E4, STRATUM_E5};
    
    private static boolean isValidStratumName(String name) {
        switch (name) {
            case(STRATUM_01):
            case(STRATUM_E0):
            case(STRATUM_E1):
            case(STRATUM_E2):
            case(STRATUM_E3):
            case(STRATUM_E4):
            case(STRATUM_E5):
                return true;
            default:
                return false;
        }
    }
    
    private static int getStratumIndex(String name) {
        int index = -1;
        int numStrata = STRATA.length;
        
        for (int i=0; i<numStrata; ++i) {
            String s = STRATA[i];
            if (s.equals(name)) {
                return i;
            }
        }
        
        return index;
    }
    
    private static boolean isLowerStratum(String name1, String name2) {
        return getStratumIndex(name1) < getStratumIndex(name2);
    }
    
    private static boolean isHigherStratum(String name1, String name2) {
        return getStratumIndex(name1) > getStratumIndex(name2);
    }
    
    public RNABloom(int k, int qDBG, int qFrag, boolean debug) {
        this.qDBG = qDBG;
        this.qFrag = qFrag;
        this.debug = debug;
        this.setK(k);
    }
    
    public final void setK(int k) {
        this.k = k;
        
        this.seqPattern = getNucleotideCharsPattern(k);
        this.qualPatternDBG = getPhred33Pattern(qDBG, k);
        this.qualPatternFrag = getPhred33Pattern(qFrag, k);
        
        if (graph != null) {
            graph.setK(k);
        }
    }
    
    public void setParams(boolean strandSpecific,
            int maxTipLength, 
            int lookahead, 
            float maxCovGradient, 
            int maxIndelSize, 
            float percentIdentity, 
            int minNumKmerPairs,
            int minPolyATail) {
        
        this.strandSpecific = strandSpecific;
        this.maxTipLength = maxTipLength;
        this.lookahead = lookahead;
        this.maxCovGradient = maxCovGradient;
        this.maxIndelSize = maxIndelSize;
        this.percentIdentity = percentIdentity;
        this.minNumKmerPairs = minNumKmerPairs;
        this.minPolyATailLengthRequired = minPolyATail;
        
        if (minPolyATail > 0) {
            polyATailOnlyPattern = getPolyATailPattern(minPolyATail);
            polyTHeadOnlyPattern = getPolyTHeadPattern(minPolyATail);
            
            polyASignalPattern = getPolyASignalPattern();
            polyATailOnlyMatchingPattern = getPolyATailMatchingPattern(minPolyATail);
            
            if (strandSpecific) {
                polyATailPattern = polyATailOnlyPattern;
            }
            else {
                polyATailPattern = getPolyTHeadOrPolyATailPattern(minPolyATail);
            }
        }
    }
    
    private static void exitOnError(String msg) {
        System.out.println("ERROR: " + msg);
        System.exit(1);
    }
    
    private static void handleException(Exception ex) {
        System.out.println("ERROR: " + ex.getMessage() );
        ex.printStackTrace();
        System.exit(1);
    }
    
    public void saveGraph(File f) throws IOException {
        graph.save(f);
    }
    
    public void restoreGraph(File f, boolean loadDbgBits) throws IOException {
        if (graph != null) {
            graph.destroy();
        }

        graph = new BloomFilterDeBruijnGraph(f, loadDbgBits);

        if (loadDbgBits) {
            dbgFPR = graph.getDbgbfFPR();
            System.out.println("DBG Bloom filter FPR:                " + dbgFPR * 100 + " %");
        }

        covFPR = graph.getCbfFPR();
        System.out.println("Counting Bloom filter FPR:           " + covFPR * 100 + " %");
        System.out.println("Read paired k-mers Bloom filter FPR: " + graph.getRpkbf().getFPR() * 100 + " %");
    }
    
    public boolean isGraphInitialized() {
        return graph != null;
    }
    
    public void clearDbgBf() {
        graph.clearDbgbf();
        dbgFPR = 0;
    }

    public void clearCBf() {
        graph.clearCbf();
        covFPR = 0;
    }
    
    public void clearSBf() {
        if (screeningBf != null) {
            screeningBf.empty();
        }
    }
    
    public void clearPkBf() {
        graph.clearPkbf();
    }
    
    public void clearRpkBf() {
        graph.clearRpkbf();
    }
    
    public void clearAllBf() {
        graph.clearAllBf();
        
        dbgFPR = 0;
        covFPR = 0;
    }

    public void destroyAllBf() {
        graph.destroy();
        
        if (screeningBf != null) {
            screeningBf.destroy();
            screeningBf = null;
        }
        
        dbgFPR = 0;
        covFPR = 0;
    }
    
    public static void checkInputFileFormat(String[] paths) throws FileFormatException {
        for (String p : paths) {
            if (!FastaReader.isCorrectFormat(p) && !FastqReader.isCorrectFormat(p)) {
                throw new FileFormatException("Unsupported file format detected in input file `" + p + "`. Only FASTA and FASTQ formats are supported.");
            }
        }
    }
    
    public class PairedKmersToGraphWorker implements Runnable {
        private final int id;
        private final String path;
        private PairedNTHashIterator pitr = null;
        private int kmerPairDistance = 0;
        private long numReads = 0;
        private boolean successful = false;
        
        public PairedKmersToGraphWorker(int id, String path, boolean stranded, int numHash) {
            this.id = id;
            this.path = path;
            this.kmerPairDistance = graph.getReadPairedKmerDistance();
            if (stranded) {
                pitr = new PairedNTHashIterator(k, numHash, kmerPairDistance);
            } else {
                pitr = new CanonicalPairedNTHashIterator(k, numHash, kmerPairDistance);
            }
        }
        
        @Override
        public void run() {
            System.out.println("[" + id + "] Parsing `" + path + "`...");
            
            try {
                Matcher mSeq = seqPattern.matcher("");

                if (FastaReader.isCorrectFormat(path)) {
                    FastaReader fr = new FastaReader(path);

                    String seq;

                    long[] phashVals = pitr.hVals3;

                    while (fr.hasNext()) {
                        seq = fr.next();
                        mSeq.reset(seq);

                        while (mSeq.find()) {
                            int start = mSeq.start();
                            int end = mSeq.end();

                            if (end - start - k + 1 >= kmerPairDistance) {
                                pitr.start(seq, start, end);
                                while (pitr.hasNext()) {
                                    pitr.next();
                                    graph.addSingleReadPairedKmer(phashVals);
                                }
                            }
                        }

                        ++numReads;
                    }
                    
                    fr.close();
                }
                else {
                    throw new RuntimeException("Unsupported file format detected in input file `" + path + "`. Only FASTA and FASTQ formats are supported.");
                }
                
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
    
    public class SeqToGraphWorker implements Runnable {
        private final int id;
        private final String path;
        private NTHashIterator itr;
        private PairedNTHashIterator pitr = null;
        private int kmerPairDistance = 0;
        private long numReads = 0;
        private boolean successful = false;
        private final Consumer<long[]> addFunction;
        private boolean storeReadPairedKmers = false;
        
        public SeqToGraphWorker(int id, String path, boolean stranded, boolean reverseComplement, int numHash, boolean incrementIfPresent, boolean storeReadPairedKmers) {            
            this.id = id;
            this.path = path;
            this.storeReadPairedKmers = storeReadPairedKmers;
            
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
            
            if (storeReadPairedKmers) {
                kmerPairDistance = graph.getReadPairedKmerDistance();
                if (stranded) {
                    if (reverseComplement) {
                        pitr = new ReverseComplementPairedNTHashIterator(k, numHash, kmerPairDistance);
                    } else {
                        pitr = new PairedNTHashIterator(k, numHash, kmerPairDistance);
                    }
                } else {
                    pitr = new CanonicalPairedNTHashIterator(k, numHash, kmerPairDistance);
                }
            }
            
            if (incrementIfPresent) {
                addFunction = graph::addCountIfPresent;
            }
            else {
                addFunction = graph::add;
            }
        }
                
        @Override
        public void run() {
            System.out.println("[" + id + "] Parsing `" + path + "`...");
            
            try {
                Matcher mSeq = seqPattern.matcher("");
                long[] hashVals = itr.hVals;
                
                if (FastqReader.isCorrectFormat(path)) {
                    FastqReader fr = new FastqReader(path);
                    Matcher mQual = qualPatternDBG.matcher("");

                    FastqRecord record = new FastqRecord();

                    if (storeReadPairedKmers) {
                        long[] phashVals = pitr.hVals3;

                        while (fr.hasNext()) {
                            fr.nextWithoutName(record);
                            mQual.reset(record.qual);
                            mSeq.reset(record.seq);

                            while (mQual.find()) {
                                mSeq.region(mQual.start(), mQual.end());
                                while (mSeq.find()) {
                                    int start = mSeq.start();
                                    int end = mSeq.end();

                                    itr.start(record.seq, start, end);
                                    while (itr.hasNext()) {
                                        itr.next();
                                        addFunction.accept(hashVals);
                                    }

                                    if (end - start - k + 1 >= kmerPairDistance) {
                                        pitr.start(record.seq, start, end);
                                        while (pitr.hasNext()) {
                                            pitr.next();
                                            graph.addSingleReadPairedKmer(phashVals);
                                        }
                                    }
                                }
                            }

                            ++numReads;
                        }
                    }
                    else {
                        while (fr.hasNext()) {
                            fr.nextWithoutName(record);
                            mQual.reset(record.qual);
                            mSeq.reset(record.seq);

                            while (mQual.find()) {
                                mSeq.region(mQual.start(), mQual.end());
                                while (mSeq.find()) {
                                    itr.start(record.seq, mSeq.start(), mSeq.end());
                                    while (itr.hasNext()) {
                                        itr.next();
                                        addFunction.accept(hashVals);
                                    }
                                }
                            }

                            ++numReads;
                        }
                    }
                    
                    fr.close();
                }
                else if (FastaReader.isCorrectFormat(path)) {
                    FastaReader fr = new FastaReader(path);

                    String seq;

                    if (storeReadPairedKmers) {
                        long[] phashVals = pitr.hVals3;

                        while (fr.hasNext()) {
                            seq = fr.next();
                            mSeq.reset(seq);

                            while (mSeq.find()) {
                                int start = mSeq.start();
                                int end = mSeq.end();

                                itr.start(seq, start, end);
                                while (itr.hasNext()) {
                                    itr.next();
                                    addFunction.accept(hashVals);
                                }

                                if (end - start - k + 1 >= kmerPairDistance) {
                                    pitr.start(seq, start, end);
                                    while (pitr.hasNext()) {
                                        pitr.next();
                                        graph.addSingleReadPairedKmer(phashVals);
                                    }
                                }
                            }

                            ++numReads;
                        }
                    }
                    else {
                        while (fr.hasNext()) {
                            seq = fr.next();
                            mSeq.reset(seq);

                            while (mSeq.find()) {
                                itr.start(seq, mSeq.start(), mSeq.end());
                                while (itr.hasNext()) {
                                    itr.next();
                                    addFunction.accept(hashVals);
                                }
                            }

                            ++numReads;
                        }
                    }
                    
                    fr.close();
                }
                else {
                    throw new RuntimeException("Unsupported file format detected in input file `" + path + "`. Only FASTA and FASTQ formats are supported.");
                }
                
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
        
    public int getMaxReadLength(String path) throws FileFormatException, IOException{
        int max = -1;

        if (FastqReader.isCorrectFormat(path)) {
            FastqRecord r = new FastqRecord();
            FastqReader fr = new FastqReader(path);
            for (int i=0; i< 100 && fr.hasNext(); ++i) {
                fr.nextWithoutName(r);
                max = Math.max(max, r.seq.length());
            }
            fr.close();
        }
        else if (FastaReader.isCorrectFormat(path)) {
            FastaReader fr = new FastaReader(path);
            for (int i=0; i< 100 && fr.hasNext(); ++i) {
                max = Math.max(max, fr.next().length());
            }
            fr.close();                
        }
        else {
            throw new FileFormatException("Incompatible file format for `" + path + "`");
        }
        
        return max;
    }
    
    public void initializeGraph(boolean strandSpecific,
                            long dbgbfNumBits,
                            long cbfNumBytes,
                            long pkbfNumBits,
                            int dbgbfNumHash,
                            int cbfNumHash,
                            int pkbfNumHash,
                            boolean initPkbf,
                            boolean useReadPairedKmers) {
        
        graph = new BloomFilterDeBruijnGraph(dbgbfNumBits,
                                            cbfNumBytes,
                                            pkbfNumBits,
                                            dbgbfNumHash,
                                            cbfNumHash,
                                            pkbfNumHash,
                                            k,
                                            strandSpecific,
                                            useReadPairedKmers);
        
        if (initPkbf) {
            graph.initializePairKmersBloomFilter(pkbfNumBits, pkbfNumHash);
        }
    }
        
    public void setReadKmerDistance(Collection<String> forwardReadPaths,
                                    Collection<String> reverseReadPaths) throws IOException {
        
        int readLength = -1;
        
        for (String path : forwardReadPaths) {
            int len = this.getMaxReadLength(path);
            if (len > 0) {
                if (readLength < 0) {
                    readLength = len;
                }
                else {
                    readLength = Math.min(readLength, len);
                }
            }
            else {
                exitOnError("Cannot determine read length from read files.");
            }
        }
        
        for (String path : reverseReadPaths) {
            int len = this.getMaxReadLength(path);
            if (len > 0) {
                if (readLength < 0) {
                    readLength = len;
                }
                else {
                    readLength = Math.min(readLength, len);
                }
            }
            else {
                exitOnError("Cannot determine read length from read files.");
            }
        }
        
        if (readLength < 0) {
            exitOnError("Cannot determine read length from read files.");
        }

        if (readLength < k) {
            exitOnError("The read length (" + readLength + ") is too short for k-mer size (" + k + ").");
        }
        
        System.out.println("Read length: " + readLength);
        
        /*
            |<--d-->|
            ==------==     paired k-mers
             ==------==    
              ==------==   
               ==------==  
                ==------== 
            ============== single read
        */
        
        graph.setReadPairedKmerDistance(readLength - k - minNumKmerPairs);
    }
    
    public int getReadLength() {
        int d = graph.getReadPairedKmerDistance();
        return d <= 0 ? d : d + k + minNumKmerPairs;
    }
    
    public void populateGraph(Collection<String> forwardReadPaths,
                            Collection<String> reverseReadPaths,
                            Collection<String> longReadPaths,
                            Collection<String> refTranscriptsPaths,
                            boolean strandSpecific,
                            boolean reverseComplementLong,
                            int numThreads,
                            boolean addCountsOnly,
                            boolean storeReadKmerPairs) throws IOException, InterruptedException {        
        
//        screeningBf = new BloomFilter(sbfNumBits, sbfNumHash, graph.getHashFunction());

        if (storeReadKmerPairs) {
            setReadKmerDistance(forwardReadPaths, reverseReadPaths);
        }

        /** parse the reads */
        
        long numReads = 0;
        int numHash = graph.getMaxNumHash();
        
        ExecutorService service = Executors.newFixedThreadPool(numThreads);
        
        ArrayList<SeqToGraphWorker> threadPool = new ArrayList<>();
        int threadId = 0;
           
        for (String path : forwardReadPaths) {
            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, path, strandSpecific, false, numHash, addCountsOnly, storeReadKmerPairs);
            service.submit(t);
            threadPool.add(t);
        }

        for (String path : reverseReadPaths) {
            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, path, strandSpecific, true, numHash, addCountsOnly, storeReadKmerPairs);
            service.submit(t);
            threadPool.add(t);
        }
        
        for (String path : longReadPaths) {
            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, path, strandSpecific, reverseComplementLong, numHash, addCountsOnly, false);
            service.submit(t);
            threadPool.add(t);
        }

        for (String path : refTranscriptsPaths) {
            PairedKmersToGraphWorker t = new PairedKmersToGraphWorker(++threadId, path, strandSpecific, numHash);
            service.submit(t);
        }
        
        service.shutdown();
        service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        for (SeqToGraphWorker t : threadPool) {
            numReads += t.getReadCount();
        }

        System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");            

        
        dbgFPR = graph.getDbgbf().getFPR();
        covFPR = graph.getCbf().getFPR();
        
        System.out.println(    "DBG Bloom filter FPR:                 " + dbgFPR * 100 + " %");
        System.out.println(    "Counting Bloom filter FPR:            " + covFPR * 100 + " %");
        
        if (graph.getReadPairedKmerDistance() > 0) {
            System.out.println("Reads paired k-mers Bloom filter FPR: " + graph.getRpkbf().getFPR() * 100 + " %");
        }
//        System.out.println("Screening Bloom filter FPR:  " + screeningBf.getFPR() * 100 + " %");
    }
    
    public boolean withinMaxFPR(float fpr) {
        float maxFPR = fpr * 1.5f;
        return (graph.getDbgbf() == null || graph.getDbgbfFPR() <= maxFPR) && 
                (graph.getCbf() == null || graph.getCbfFPR() <= maxFPR) && 
                (graph.getRpkbf() == null || graph.getRpkbfFPR() <= maxFPR);
    }
    
    public long[] getOptimalBloomFilterSizes(float maxFPR) {
        long maxPopCount = Math.max(screeningBf.getPopCount(),
                                Math.max(graph.getDbgbf().getPopCount(),
                                    Math.max(graph.getCbf().getPopCount(), graph.getRpkbf().getPopCount())));
        
        long sbfSize = BloomFilter.getExpectedSize(maxPopCount, maxFPR, screeningBf.getNumHash());
        long dbgbfSize = BloomFilter.getExpectedSize(maxPopCount, maxFPR, graph.getDbgbfNumHash());
        long cbfSize = CountingBloomFilter.getExpectedSize(maxPopCount, maxFPR, graph.getCbfNumHash());
        long pkbfSize = PairedKeysBloomFilter.getExpectedSize(maxPopCount, maxFPR, graph.getPkbfNumHash());
                
        return new long[]{sbfSize, dbgbfSize, cbfSize, pkbfSize};
    }
    
    public void addPairedKmersFromSequences(String[] fastas) throws IOException {
        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;

        for (String path : fastas) {
            FastaReader fin = new FastaReader(path);

            System.out.println("Parsing `" + path + "`...");

            String seq;
            while (fin.hasNext()) {
                seq = fin.next();

                if (pItr.start(seq)) {
                    while (pItr.hasNext()) {
                        pItr.next();
                        graph.addPairedKmers(hashVals1, hashVals2, hashVals3);
                    }
                }
            }

            fin.close();
        }
    }
    
    public void populateGraphFromFragments(Collection<String> fastas, boolean strandSpecific, boolean loadPairedKmers) throws IOException {
        /** insert into graph if absent */
        
        /** parse the fragments */
                          
        NTHashIterator itr = graph.getHashIterator(graph.getMaxNumHash());
        long[] hashVals = itr.hVals;

        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;

        for (String path : fastas) {
            FastaReader fin = new FastaReader(path);

            System.out.println("Parsing `" + path + "`...");

            String seq;
            if (loadPairedKmers) {
                while (fin.hasNext()) {
                    seq = fin.next();

                    if (itr.start(seq)) {
                        while (itr.hasNext()) {
                            itr.next();
                            graph.addDbgOnly(hashVals);
                        }
                    }

                    if (pItr.start(seq)) {
                        while (pItr.hasNext()) {
                            pItr.next();
                            graph.addPairedKmers(hashVals1, hashVals2, hashVals3);
                        }
                    }
                }
            }
            else {
                while (fin.hasNext()) {
                    seq = fin.next();

                    if (itr.start(seq)) {
                        while (itr.hasNext()) {
                            itr.next();
                            graph.addDbgOnly(hashVals);
                        }
                    }
                }
            }

            fin.close();
        }
                
        dbgFPR = graph.getDbgbfFPR();
        System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
//        covFPR = graph.getCbfFPR();
    }
        
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
        String left;
        String right;
        ArrayList<Kmer> kmers;
        int length;
        float minCov;
        boolean isUnconnectedRead;
        
        public Fragment(String left, String right, ArrayList<Kmer> kmers, int length, float minCov, boolean isUnconnectedRead) {
            this.left = left;
            this.right = right;
            this.kmers = kmers;
            this.length = length;
            this.minCov = minCov;
            this.isUnconnectedRead = isUnconnectedRead;
        }
    }
    
    private class Transcript {
        String fragment;
        ArrayList<Kmer> transcriptKmers;
        
        public Transcript(String fragment, ArrayList<Kmer> transcriptKmers) {
            this.fragment = fragment;
            this.transcriptKmers = transcriptKmers;
        }
    }
    
    private class TranscriptWriter {
        private final FastaWriter fout;
        private final FastaWriter foutShort;
        private final int minTranscriptLength;
        private final int maxTipLength;
        private final boolean writeUracil;
        private String prefix = "";
        private long cid = 0;
        
        public TranscriptWriter(FastaWriter fout,
                                FastaWriter foutShort,
                                int minTranscriptLength,
                                int maxTipLength,
                                boolean writeUracil) {
            this.fout = fout;
            this.foutShort = foutShort;
            this.minTranscriptLength = minTranscriptLength;
            this.maxTipLength = maxTipLength;
            this.writeUracil = writeUracil;
        }
        
        public void setOutputPrefix(String prefix) {
            this.prefix = prefix;
        }
                
        public void write(String fragment, ArrayList<Kmer> transcriptKmers) throws IOException {            
            if (!represented(transcriptKmers,
                                graph,
                                screeningBf,
                                lookahead,
                                maxIndelSize,
                                maxTipLength,
                                percentIdentity)) {

                String transcript = graph.assemble(transcriptKmers);
                ArrayDeque<Integer> pasPositions = null;
                
                boolean txptReverseComplemented = false;
                if (minPolyATailLengthRequired > 0) {
                    boolean hasPolyATail = polyATailOnlyPattern.matcher(transcript).matches();
                    boolean hasPolyTHead = false;
                    
                    if (strandSpecific) {
                        pasPositions = getPolyASignalPositions(transcript, polyASignalPattern, polyATailOnlyMatchingPattern);
                    }
                    else {
                        hasPolyTHead = polyTHeadOnlyPattern.matcher(transcript).matches();
                        
                        if (hasPolyATail) {
                            pasPositions = getPolyASignalPositions(transcript, polyASignalPattern, polyATailOnlyMatchingPattern);
                        }
                        
                        if (hasPolyTHead && (pasPositions == null || pasPositions.isEmpty())) {
                            String transcriptRC = reverseComplement(transcript);
                            
                            pasPositions = getPolyASignalPositions(transcriptRC, polyASignalPattern, polyATailOnlyMatchingPattern);
                            
                            if (!pasPositions.isEmpty()) {
                                transcript = transcriptRC;
                                txptReverseComplemented = true;
                            }
                        }
                    }
                }
                                
                for (Kmer kmer : transcriptKmers) {
                    screeningBf.add(kmer.getHash());
                }
                
                int len = transcript.length();
                
                StringBuilder headerBuilder = new StringBuilder();
                headerBuilder.append(prefix);
                headerBuilder.append(++cid);
                headerBuilder.append(" l=");
                headerBuilder.append(len);
                if (!fragment.isEmpty()) {
                    headerBuilder.append(" F=[");
                    headerBuilder.append(fragment);
                    headerBuilder.append("]");
                }
                
                if (minPolyATailLengthRequired > 0) {                    
                    if (pasPositions != null && !pasPositions.isEmpty()) {
                        // add PAS and their positions to header
                        headerBuilder.append(" PAS=[");
                        
                        Iterator<Integer> itr = pasPositions.iterator();

                        int pasPos = itr.next();
                        int numTanscriptKmers = transcriptKmers.size();

                        headerBuilder.append(pasPos);
                        headerBuilder.append(':');
                        
                        if (txptReverseComplemented) {
                            // `transcriptKmers` has a poly-T-head
                            
                            headerBuilder.append(GraphUtils.getMinimumKmerCoverage(transcriptKmers, 0, Math.min(numTanscriptKmers-1, numTanscriptKmers-1-pasPos-6+k)));
                            headerBuilder.append(':');
                            headerBuilder.append(transcript.substring(pasPos, pasPos+6));

                            while (itr.hasNext()) {
                                headerBuilder.append(",");

                                pasPos = itr.next();

                                headerBuilder.append(pasPos);
                                headerBuilder.append(':');
                                headerBuilder.append(GraphUtils.getMinimumKmerCoverage(transcriptKmers, 0, Math.min(numTanscriptKmers-1, numTanscriptKmers-1-pasPos-6+k)));
                                headerBuilder.append(':');
                                headerBuilder.append(transcript.substring(pasPos, pasPos+6));
                            }                               
                        }
                        else {
                            // `transcriptKmers` has a poly-A-tail
                            
                            headerBuilder.append(GraphUtils.getMinimumKmerCoverage(transcriptKmers, Math.max(0, pasPos+6-k), numTanscriptKmers-1));
                            headerBuilder.append(':');
                            headerBuilder.append(transcript.substring(pasPos, pasPos+6));

                            while (itr.hasNext()) {
                                headerBuilder.append(",");

                                pasPos = itr.next();

                                headerBuilder.append(pasPos);
                                headerBuilder.append(':');
                                headerBuilder.append(GraphUtils.getMinimumKmerCoverage(transcriptKmers, Math.max(0, pasPos+6-k), numTanscriptKmers-1));
                                headerBuilder.append(':');
                                headerBuilder.append(transcript.substring(pasPos, pasPos+6));
                            }    
                        }
                        
                        // mask PAS in transcript
                        StringBuilder transcriptSB = new StringBuilder(transcript);

                        for (int p : pasPositions) {
                            for (int i=0; i<6; ++i) {
                                char c = transcript.charAt(p+i);
                                if (writeUracil && c=='T') {
                                    transcriptSB.setCharAt(p+i, 'u');
                                }
                                else {
                                    transcriptSB.setCharAt(p+i, Character.toLowerCase(c));
                                }
                            }
                        }
                        
                        transcript = transcriptSB.toString();
                        
                        headerBuilder.append("]");
                    }
                }
                
                if (writeUracil) {
                    transcript = transcript.replace('T', 'U');
                }
                
                if (len >= minTranscriptLength) {
                    fout.write(headerBuilder.toString(), transcript);
                }
                else {
                    foutShort.write(headerBuilder.toString(), transcript);
                }
            }
        }
    }
    
    private static int intervalOverlapSize(int[] x, int[] y) {
        int start = Math.max(x[0], y[0]);
        int end = Math.min(x[1], y[1]);
        
        return Math.max(0, end-start);
    }
    
    private class TranscriptAssemblyWorker implements Runnable {
        
        private final ArrayBlockingQueue<String> fragments;
        private final ArrayBlockingQueue<Transcript> transcripts;
        private boolean keepGoing = true;
        private boolean extendBranchFreeFragmentsOnly = false;
        private boolean keepArtifact = false;
        private boolean keepChimera = false;
        private boolean haveFragKmers = false;
        private final float minKmerCov;
        
        public TranscriptAssemblyWorker(ArrayBlockingQueue<String> fragments,
                                        ArrayBlockingQueue<Transcript> transcripts,
                                        boolean includeNaiveExtensions,
                                        boolean extendBranchFreeFragmentsOnly,
                                        boolean keepArtifact,
                                        boolean keepChimera,
                                        boolean haveFragKmers,
                                        float minKmerCov) {
            this.fragments = fragments;
            this.transcripts = transcripts;
            this.extendBranchFreeFragmentsOnly = extendBranchFreeFragmentsOnly;
            this.keepArtifact = keepArtifact;
            this.keepChimera = keepChimera;
            this.haveFragKmers = haveFragKmers;
            this.minKmerCov = minKmerCov;
        }

        public void stopWhenEmpty () {
            keepGoing = false;
        }
        
        @Override
        public void run() {
            try {
                int fragKmersDist = graph.getFragPairedKmerDistance();
                int maxEdgeClipLength = minPolyATailLengthRequired > 0 ? 0 : maxTipLength;
                boolean keepBluntEndArtifact = keepArtifact;
                
                while (true) {
//                    String seq = "";
//                    ArrayList<Kmer> kmers2 = graph.getKmers(seq);
//                    printPairedKmersPositions(kmers2, graph);
                    
                    String fragment = fragments.poll(10, TimeUnit.MICROSECONDS);
                    
                    if (fragment == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
                        ArrayList<Kmer> kmers = graph.getKmers(fragment);
                        
                        if (!kmers.isEmpty()) {
                            if ( (!extendBranchFreeFragmentsOnly || isBranchFree(kmers, graph, maxTipLength)) &&
                                 !represented(kmers,
                                                graph,
                                                screeningBf,
                                                lookahead,
                                                maxIndelSize,
                                                maxEdgeClipLength,
                                                percentIdentity) &&
                                 (keepChimera || !isChimera(kmers, graph, screeningBf, lookahead)) &&
                                 (keepBluntEndArtifact || !isBluntEndArtifact(kmers, graph, screeningBf, maxEdgeClipLength)) ) {
                                
                                int[] originalFragRange;
                                
                                if (haveFragKmers) {
                                    originalFragRange = extendPE(kmers, graph, maxTipLength, minKmerCov);
                                }
                                else {
                                    originalFragRange = extendSE(kmers, graph, maxTipLength, minKmerCov);
                                }

                                int[] currentRange = new int[]{0, kmers.size()};
                                
                                if (haveFragKmers) {
                                    if (kmers.size() >= fragKmersDist) {
                                        ArrayDeque<int[]> ranges = breakWithFragPairedKmers(kmers, graph, lookahead);
                                        int numFragSegs = ranges.size();

                                        if (numFragSegs == 1) {
                                            currentRange = ranges.peekFirst();
                                        }
                                        else if (numFragSegs > 1) {
                                            int bestOverlap = 0;
                                            int[] bestRange = null;
                                            for (int[] range : ranges) {
                                                int overlap = intervalOverlapSize(range, originalFragRange);
                                                if (overlap > bestOverlap) {
                                                    bestOverlap = overlap;
                                                    bestRange = range;
                                                }
                                            }

                                            currentRange = bestRange;
                                        }
                                        else {
                                            currentRange = null;
                                        }
                                    }
                                    else {
                                        currentRange = null;
                                    }
                                }
                                
                                if (currentRange != null) {
                                    ArrayDeque<int[]> ranges = breakWithReadPairedKmers(kmers, graph, lookahead, currentRange[0], currentRange[1]);
                                    int numFragSegs = ranges.size();

                                    if (numFragSegs > 0) {
                                        if (numFragSegs == 1) {
                                            currentRange = ranges.peekFirst();
                                        }
                                        else if (numFragSegs > 1) {
                                            int bestOverlap = 0;
                                            int[] bestRange = null;
                                            for (int[] range : ranges) {
                                                int overlap = intervalOverlapSize(range, originalFragRange);
                                                if (overlap > bestOverlap) {
                                                    bestOverlap = overlap;
                                                    bestRange = range;
                                                }
                                            }

                                            currentRange = bestRange;
                                        }
                                        
                                        if (currentRange != null) {
                                            ArrayList<Kmer> txptKmers = new ArrayList<>(kmers.subList(currentRange[0], currentRange[1]));

                                            String fragInfo = debug ? fragment : "";
                                            
                                            if (!keepArtifact) {
                                                txptKmers = trimReverseComplementArtifact(txptKmers, maxTipLength, maxIndelSize, percentIdentity, graph);
                                            }
                                            
//                                            if (!keepArtifact) {
//                                                if (!isTemplateSwitch2(txptKmers, graph, screeningBf, lookahead, percentIdentity)) {
//                                                    txptKmers = trimHairpinBySequenceMatching(txptKmers, k, percentIdentity, graph);
//                                                    transcripts.put(new Transcript(fragInfo, txptKmers));
//                                                }
//                                            }
//                                            else {
                                                transcripts.put(new Transcript(fragInfo, txptKmers));
//                                            }     
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
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
    
    private class ReadConnector implements Runnable {
        private String left;
        private String right;
        private ArrayBlockingQueue<Fragment> outList;
        private int bound;
        private int minOverlap;
        private boolean storeKmerPairs;
        private boolean extendFragments;
        private int errorCorrectionIterations = 0;
        private float minKmerCov;
        
        public ReadConnector(String left,
                                String right,
                                ArrayBlockingQueue<Fragment> outList,
                                int bound, 
                                int minOverlap,
                                int errorCorrectionIterations,
                                boolean storeKmerPairs,
                                boolean extendFragments,
                                float minKmerCov) {
            
            this.left = left;
            this.right = right;
            this.outList = outList;
            this.bound = bound;
            this.minOverlap = minOverlap;
            this.storeKmerPairs = storeKmerPairs;
            this.extendFragments = extendFragments;
            this.errorCorrectionIterations = errorCorrectionIterations;
            this.minKmerCov = minKmerCov;
        }
        
        @Override
        public void run() {
            try {
                
//                System.out.println("L: " + left);
//                System.out.println("R: " + right);
                
                ArrayList<Kmer> leftKmers = graph.getKmers(left);
                ArrayList<Kmer> rightKmers = graph.getKmers(right);

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
                                                        percentIdentity,
                                                        minKmerCov);

                    if (correctedReadPair.corrected) {
                        leftKmers = correctedReadPair.leftKmers;
                        rightKmers = correctedReadPair.rightKmers;
                    }
                }
                
                if (!leftKmers.isEmpty() && !rightKmers.isEmpty()) {

                    ArrayList<Kmer> fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound-k+1-leftKmers.size()-rightKmers.size(),
                            lookahead, minOverlap, maxCovGradient, maxTipLength, maxIndelSize, percentIdentity, minKmerCov);

                    ArrayDeque<int[]> ranges = breakWithReadPairedKmers(fragmentKmers, graph, lookahead);
                    
                    if (ranges.size() != 1) {
                        fragmentKmers = null;
                    }
                    
                    if (fragmentKmers != null) {
                        int fragLength = fragmentKmers.size() + k - 1;

                        if (fragLength >= k + lookahead) {
                            boolean hasComplexKmer = false;

                            float minCov = Float.MAX_VALUE;
                            for (Kmer kmer : fragmentKmers) {
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
                                    fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength, minKmerCov);
                                }

                                if (this.storeKmerPairs) {
                                    graph.addPairedKmers(fragmentKmers);
                                }

                                outList.put(new Fragment(left, right, fragmentKmers, fragLength, minCov, false));
                            }
                        }
                    }
                    else {
                        // this is an unconnected read pair
                        float minCov = Float.MAX_VALUE;

                        boolean hasComplexLeftKmer = false;

                        if (leftKmers.size() >= lookahead) {
                            for (Kmer kmer : leftKmers) {
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
                            for (Kmer kmer : rightKmers) {
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
        private ArrayBlockingQueue<PairedReadSegments> inList;
        private ArrayBlockingQueue<Fragment> outList;
        private int bound;
        private int minOverlap;
        private boolean storeKmerPairs;
        private int errorCorrectionIterations;
        private int leftReadLengthThreshold;
        private int rightReadLengthThreshold;
        private boolean extendFragments;
        private int minKmerCov;
        private boolean trimArtifact;
        private boolean keepGoing = true;
        
        public FragmentAssembler(ArrayBlockingQueue<PairedReadSegments> inList,
                                ArrayBlockingQueue<Fragment> outList,
                                int bound, 
                                int minOverlap, 
                                boolean storeKmerPairs, 
                                int errorCorrectionIterations,
                                int leftReadLengthThreshold,
                                int rightReadLengthThreshold,
                                boolean extendFragments,
                                int minKmerCov,
                                boolean keepArtifact) {
            
            this.inList = inList;
            this.outList = outList;
            this.bound = bound;
            this.minOverlap = minOverlap;
            this.storeKmerPairs = storeKmerPairs;
            this.errorCorrectionIterations = errorCorrectionIterations;
            this.leftReadLengthThreshold = leftReadLengthThreshold;
            this.rightReadLengthThreshold = rightReadLengthThreshold;
            this.extendFragments = extendFragments;
            this.minKmerCov = minKmerCov;
            this.trimArtifact = !keepArtifact;
        }
        
        @Override
        public void run() {
            try {
                while(true) {
                    PairedReadSegments p  = inList.poll(10, TimeUnit.MICROSECONDS);
                    
                    if (p == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {                        
                        ArrayList<Kmer> leftKmers = null;
                        ArrayList<Kmer> rightKmers = null;
                        
                        // connect segments of each read
                        //String left = getBestSegment(p.left, graph);
                        String left = connect(p.left, graph, lookahead);

                        if (left.length() >= this.leftReadLengthThreshold) {
                            if (!isRepeat(left)) {
                                if (minKmerCov > 1) {
                                    leftKmers = graph.getKmers(left, minKmerCov);
                                }
                                else {
                                    leftKmers = graph.getKmers(left);
                                }
                            }
                        }

                        //String right = getBestSegment(p.right, graph);
                        String right = connect(p.right, graph, lookahead);

                        if (right.length() >= this.rightReadLengthThreshold) {
                            if (!isRepeat(right)) {
                                if (minKmerCov > 1) {
                                    rightKmers = graph.getKmers(right, minKmerCov);
                                }
                                else {
                                    rightKmers = graph.getKmers(right);
                                }
                            }
                        }

                        boolean leftBad = leftKmers == null || leftKmers.isEmpty();
                        boolean rightBad = rightKmers == null || rightKmers.isEmpty();

                        if (leftBad && rightBad) {
                            continue;
                        }

                        ArrayList<Kmer> fragmentKmers = null;

                        if (!leftBad && !rightBad) {
                            if (errorCorrectionIterations > 0) {
                                ReadPair correctedReadPair = correctErrorsPE(leftKmers,
                                                                            rightKmers,
                                                                            graph, 
                                                                            lookahead, 
                                                                            maxIndelSize, 
                                                                            maxCovGradient, 
                                                                            covFPR,
                                                                            this.errorCorrectionIterations,
                                                                            2,
                                                                            percentIdentity,
                                                                            minKmerCov);

                                if (correctedReadPair.corrected) {
                                    leftKmers = correctedReadPair.leftKmers;
                                    rightKmers = correctedReadPair.rightKmers;
                                }
                            }

                            fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound,
                                    lookahead, minOverlap, maxCovGradient, maxTipLength, maxIndelSize, percentIdentity, minKmerCov);
                        }
                        else if (!leftBad) {
                            if (errorCorrectionIterations > 0) {
                                ArrayList<Kmer> corrected = correctErrorsSE(leftKmers,
                                                                            graph, 
                                                                            lookahead, 
                                                                            maxIndelSize, 
                                                                            maxCovGradient, 
                                                                            covFPR,
                                                                            percentIdentity,
                                                                            minKmerCov);
                                if (corrected != null && !corrected.isEmpty()) {
                                    leftKmers = corrected;
                                }
                            }
                        }
                        else if (!rightBad) {
                            if (errorCorrectionIterations > 0) {
                                ArrayList<Kmer> corrected = correctErrorsSE(rightKmers,
                                                                            graph, 
                                                                            lookahead, 
                                                                            maxIndelSize, 
                                                                            maxCovGradient, 
                                                                            covFPR,
                                                                            percentIdentity,
                                                                            minKmerCov);
                                if (corrected != null && !corrected.isEmpty()) {
                                    rightKmers = corrected;
                                }
                            }
                        }

                        // check for read consistency if fragment is long enough
                        int preExtensionFragLen = -1;
                        
                        if (fragmentKmers != null) {
                            if (graph.getReadPairedKmerDistance() < fragmentKmers.size()) {
                                ArrayDeque<int[]> ranges = breakWithReadPairedKmers(fragmentKmers, graph, lookahead);

                                if (ranges.size() != 1) {
                                    fragmentKmers = null;
                                }
                                else {
                                    int[] range = ranges.peek();
                                    if (range[0] >= leftKmers.size() || range[1] < fragmentKmers.size() - rightKmers.size()) {
                                        // fragment is a spurious connection between left and right reads
                                        fragmentKmers = null;
                                    }
                                    else if (range[0] > 0 || range[1] < fragmentKmers.size()) {
                                        // trim the fragment
                                        fragmentKmers = new ArrayList<>(fragmentKmers.subList(range[0], range[1]));
                                    }
                                }
                            }
                            
                            if (fragmentKmers != null && trimArtifact) {
                                fragmentKmers = trimReverseComplementArtifact(fragmentKmers, maxTipLength, maxIndelSize, percentIdentity, graph);
                                //fragmentKmers = trimHairpinBySequenceMatching(fragmentKmers, k, percentIdentity, graph);
                            }
                            
                            if (fragmentKmers != null) {
                                preExtensionFragLen = fragmentKmers.size() + k - 1;
                                
                                if (extendFragments) {
                                    ArrayList<Kmer> extendedFragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength, minKmerCov);
                                    
                                    if (extendedFragmentKmers.size() != preExtensionFragLen) {
                                        // fragment was extended; check consistency with reads
                                        ArrayDeque<int[]> ranges = breakWithReadPairedKmers(extendedFragmentKmers, graph, lookahead);
                                        
                                        if (ranges.size() == 1) {
                                            // there is one consistent section in the extended fragment
                                            int[] range = ranges.peek();
                                            if (range[0] > 0 || range[1] < extendedFragmentKmers.size()) {
                                                // trim extended fragment
                                                fragmentKmers = new ArrayList<>(extendedFragmentKmers.subList(range[0], range[1]));
                                            }
                                            else {
                                                // trimming not needed; replace original fragment with extended fragment
                                                fragmentKmers = extendedFragmentKmers;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (fragmentKmers != null) {
                            if (fragmentKmers.size() + k - 1 >= k + lookahead) {
                                boolean hasComplexKmer = false;

                                float minCov = Float.MAX_VALUE;
                                for (Kmer kmer : fragmentKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexKmer) {
                                        if (!graph.isRepeatKmer(kmer)) {
                                            hasComplexKmer = true;
                                        }
                                    }
                                }

                                if (hasComplexKmer) {
                                    if (this.storeKmerPairs) {
                                        graph.addPairedKmers(fragmentKmers);
                                    }

                                    if (!debug) {
                                       left = "";
                                       right = "";
                                    }

                                    outList.put(new Fragment(left, right, fragmentKmers, preExtensionFragLen, minCov, false));
                                }
                            }
                        }
                        else {
                            // this is an unconnected read pair
                            float minCov = Float.MAX_VALUE;

                            boolean hasComplexLeftKmer = false;

                            if (!leftBad && leftKmers.size() >= lookahead) {
                                for (Kmer kmer : leftKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexLeftKmer && !graph.isRepeatKmer(kmer)) {
                                        hasComplexLeftKmer = true;
                                    }
                                }
                            }

                            boolean hasComplexRightKmer = false;

                            if (!rightBad && rightKmers.size() >= lookahead) {
                                for (Kmer kmer : rightKmers) {
                                    if (kmer.count < minCov) {
                                        minCov = kmer.count;
                                    }

                                    if (!hasComplexRightKmer && !graph.isRepeatKmer(kmer)) {
                                        hasComplexRightKmer = true;
                                    }
                                }
                            }

                            if (hasComplexLeftKmer || hasComplexRightKmer) {
                                left = leftBad ? "" : left;
                                right = rightBad ? "" : right;

                                outList.put(new Fragment(left, right, null, 0, minCov, true));
                            }
                        }
                    }
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
                throw new RuntimeException(ex);
            }
        }
        
        public void updateBound(int bound) {
            this.bound = bound;
        }
        
        public void stopWhenEmpty () {
            keepGoing = false;
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
        
    private int[] getMinQ1MedianQ3Max(ArrayList<Integer> a) {
        if (a.isEmpty()) {
            int[] result = new int[5];
            Arrays.fill(result, 0);
            return result;
        }
        
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
    
    private static int getCoverageOrderOfMagnitude(float c) {
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
        else {
            return 0;
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
    
    public void setupFragmentPairedKmersBloomFilter(long numBits, int numHash) {
        graph.initializePairKmersBloomFilter(numBits, numHash);
    }
    
    public void updateFragmentKmerDistance(String graphFile) throws IOException {
        graph.updateFragmentKmerDistance(new File(graphFile));
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
                                                boolean extendFragments,
                                                float minKmerCov) throws IOException, InterruptedException {
        
        /* make sure paired kmer distance is set */
        
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }
        
        System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
        System.out.println("Counting Bloom filter FPR: " + covFPR * 100 + " %");
        
        
        System.out.println("Rescuing unconnected read pairs...");
                
        long fragmentId = 0;
        long unconnectedReadId = 0;
        long readPairsParsed = 0;
        long rescuedReadPairs = 0;
        
        int maxTasksQueueSize = numThreads;
        
        int newBound = bound;
        int shortestFragmentLengthAllowed = k + lookahead;
        
        // set up thread pool
        MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);

        int maxConcurrentSubmissions = numThreads + maxTasksQueueSize;
        
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
                                                    extendFragments,
                                                    minKmerCov
                    ));
                }

                if (fragments.remainingCapacity() <= maxConcurrentSubmissions) {

                    // write fragments to file
                    int m;
                    Fragment frag;
                    for (int i=0; i<sampleSize; ++i) {
                        frag = fragments.poll();

                        if (frag == null) {
                            break;
                        }

                        if (frag.minCov >= 2) {
                            /*
                            When reads were parsed at k2, kmers common to both fragments and reads have counts incremented by 1.
                            Kmers unique to fragments would have a count of 1.
                            So, min kmer counts >= 2 need to be decremented by 1 when assigning fragments.
                            */ 
                            --frag.minCov;
                        }

                        if (frag.isUnconnectedRead) {
                            if (frag.minCov == 1) {
                                unconnectedSingletonsOut.write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                                unconnectedSingletonsOut.write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                            }
                            else if (frag.minCov > 1) {
                                m = getCoverageOrderOfMagnitude(frag.minCov);

                                if (m >= 0) {
                                    unconnectedReadsOut[m].write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                                    unconnectedReadsOut[m].write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                                }
                            }
                        }
                        else {
                            ++rescuedReadPairs;

                            if (frag.length >= shortestFragmentLengthAllowed) {
                                ArrayList<Kmer> fragKmers = frag.kmers;

                                if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                    if (frag.minCov == 1) {
                                        graph.addPairedKmers(fragKmers);

                                        if (frag.length >= longFragmentLengthThreshold) {
                                            for (Kmer kmer : fragKmers) {
                                                screeningBf.add(kmer.getHash());
                                            } 

                                            longSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                        }
                                        else {
                                            shortSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                        }
                                    }
                                    else if (frag.minCov > 1) {
                                        m = getCoverageOrderOfMagnitude(frag.minCov);

                                        if (m >= 0) {
                                            graph.addPairedKmers(fragKmers);

                                            if (frag.length >= longFragmentLengthThreshold) {
                                                for (Kmer kmer : fragKmers) {
                                                    screeningBf.add(kmer.getHash());
                                                }

                                                longFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                            }
                                            else {
                                                shortFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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
                    m = getCoverageOrderOfMagnitude(frag.minCov);

                    if (m >= 0) {
                        unconnectedReadsOut[m].write("r" + Long.toString(++unconnectedReadId) + "L ", frag.left);
                        unconnectedReadsOut[m].write("r" + Long.toString(unconnectedReadId) + "R", frag.right);
                    }
                }
            }
            else {
                ++rescuedReadPairs;

                if (frag.length >= shortestFragmentLengthAllowed) {
                    ArrayList<Kmer> fragKmers = frag.kmers;

                    if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                        if (frag.minCov == 1) {
                            graph.addPairedKmers(fragKmers);

                            if (frag.length >= longFragmentLengthThreshold) {
                                for (Kmer kmer : fragKmers) {
                                    screeningBf.add(kmer.getHash());
                                }

                                longSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                            }
                            else {
                                shortSingletonsOut.write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                            }
                        }
                        else if (frag.minCov > 1)  {
                            m = getCoverageOrderOfMagnitude(frag.minCov);

                            if (m >= 0) {
                                graph.addPairedKmers(fragKmers);

                                if (frag.length >= longFragmentLengthThreshold) {
                                    for (Kmer kmer : fragKmers) {
                                        screeningBf.add(kmer.getHash());
                                    }

                                    longFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                }
                                else {
                                    shortFragmentsOut[m].write("r" + Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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

        System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
        System.out.println("Rescued " + NumberFormat.getInstance().format(rescuedReadPairs) + " read pairs.");
        System.out.println("Paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
        System.out.println("Screening Bloom filter FPR:    " + screeningBf.getFPR() * 100 + " %");    
    }
    
    private final static String LABEL_SEPARATOR = ":";
    private final static String LABEL_MIN = "min";
    private final static String LABEL_Q1 = "Q1";
    private final static String LABEL_MEDIAN = "M";
    private final static String LABEL_Q3 = "Q3";
    private final static String LABEL_MAX = "max";
    
    public void writeFragStatsToFile(int[] fragStats, String path) throws IOException {
        FileWriter writer = new FileWriter(path, false);

        writer.write(LABEL_MIN + LABEL_SEPARATOR + fragStats[0] + "\n" +
                    LABEL_Q1 + LABEL_SEPARATOR + fragStats[1] + "\n" +
                    LABEL_MEDIAN + LABEL_SEPARATOR + fragStats[2] + "\n" +
                    LABEL_Q3 + LABEL_SEPARATOR + fragStats[3] + "\n" +
                    LABEL_MAX + LABEL_SEPARATOR + fragStats[4] + "\n"
                );
        writer.close();
    }
    
    public int[] restoreFragStatsFromFile(String path) throws IOException {
        int[] fragStats = new int[5];
        
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

        longFragmentLengthThreshold = fragStats[1];
          
        return fragStats;
    }
    
    private void setPairedKmerDistance(int fragStatQ1) {
        graph.setFragPairedKmerDistance(fragStatQ1 - k - minNumKmerPairs);
    }
    
    private int getPairedReadsMaxDistance(int[] fragStats) {
        return fragStats[3] + ((fragStats[3] - fragStats[1]) * 3 / 2); // 1.5*IQR
    }
    
    public class ContainmentCalculator implements Runnable {
        private final long[] queryHashVals;
        private final ArrayBlockingQueue<Integer> sketchIDsQueue;
        private int bestOverlap = -1;
        private float bestOverlapProportion = -1;
        private int bestTargetIndex = -1;
        private final ArrayList<long[]> targetSketches;
        private boolean terminateWhenInputExhausts = false;
        private boolean successful = false;
        private Exception exception = null;
        private int minSketchOverlap = -1;
        private float minSketchOverlapProportion = -1;
        private ArrayDeque<Integer> overlappingSketchIndexes = new ArrayDeque<>();
        
        public ContainmentCalculator(long[] queryHashVals, ArrayList<long[]> targetSketches, 
                                    ArrayBlockingQueue<Integer> sketchIDsQueue,
                                    int minSketchOverlap,
                                    float minSketchOverlapProportion) {
            this.queryHashVals = queryHashVals;
            this.targetSketches = targetSketches;
            this.sketchIDsQueue = sketchIDsQueue;
            this.minSketchOverlap = minSketchOverlap;
            this.minSketchOverlapProportion = minSketchOverlapProportion;
        }
        
        @Override
        public void run() {
            try {
                while(!(terminateWhenInputExhausts && sketchIDsQueue.isEmpty())) {
                    Integer sketchID = sketchIDsQueue.poll(1, TimeUnit.MILLISECONDS);

                    if (sketchID == null) {
                        continue;
                    }
                    
                    long[] sketch = targetSketches.get(sketchID);
                    int overlap = getNumIntersection(queryHashVals, sketch);
                    float overlapProportion = Math.max((float)overlap / (float)sketch.length, (float)overlap / (float)queryHashVals.length);

                    if (overlapProportion >= minSketchOverlapProportion || overlap >= minSketchOverlap) {
                        if (bestTargetIndex < 0) {
                            bestTargetIndex = sketchID;
                            bestOverlap = overlap;
                            bestOverlapProportion = overlapProportion;
                        }
                        else if (overlapProportion >= bestOverlapProportion || overlap >= bestOverlap) {
                            overlappingSketchIndexes.add(bestTargetIndex);
                            bestTargetIndex = sketchID;
                            bestOverlap = overlap;
                            bestOverlapProportion = overlapProportion;
                        }
                    }
                }
                
                successful = true;
                
            } catch (InterruptedException ex) {
                exception = ex;
                successful = false;
            }
        }
        
        public int getBestTargetIndex() {
            return bestTargetIndex;
        }
        
        public int getMaxIntersection() {
            return bestOverlap;
        }
        
        public ArrayDeque<Integer> getOverlappingSketchIndexes() {
            return overlappingSketchIndexes;
        }
        
        public void terminateWhenInputExhausts() {
            terminateWhenInputExhausts = true;
        }
        
        public boolean isSucessful() {
            return successful;
        }
        
        public Exception getExceptionCaught() {
            return exception;
        }
    }
    
    private static void mergeClusterFastas(int target, int maxAltID, String clusteredLongReadsDirectory) throws IOException {
        String path = clusteredLongReadsDirectory + File.separator + target + FASTA_EXT;
        
        FastaWriter writer = new FastaWriter(path, true);
        for (int i=0; i<=maxAltID; ++i) {
            path = clusteredLongReadsDirectory + File.separator + target + "_" + i + FASTA_EXT;
            FastaReader reader = new FastaReader(path);
            while (reader.hasNext()) {
                String[] nameSeqPair = reader.nextWithName();
                writer.write(nameSeqPair[0], nameSeqPair[1]);
            }
            reader.close();
            
            Files.deleteIfExists(FileSystems.getDefault().getPath(path));
        }
        writer.close();
    }
    
    private static int aggregateClusterFastas(int target, HashSet<Integer> overlaps, ArrayList<Integer> altIDs, String clusteredLongReadsDirectory) throws IOException {
        int altID = altIDs.get(target);
        
        for (int i : overlaps) {
            Path source = Paths.get(clusteredLongReadsDirectory + File.separator + i + FASTA_EXT);
            Path dest = Paths.get(clusteredLongReadsDirectory + File.separator + target + "_" + ++altID + FASTA_EXT);
            Files.move(source, dest);
            
            int overlapAltID = altIDs.get(i);
            for (int j=0; j<=overlapAltID; ++j) {
                source = Paths.get(clusteredLongReadsDirectory + File.separator + i + "_" + j + FASTA_EXT);
                dest = Paths.get(clusteredLongReadsDirectory + File.separator + target + "_" + ++altID + FASTA_EXT);
                Files.move(source, dest);
            }
            
        }
        
        return altID;
    }
    
    public void clusterLongReads(String[][] correctedLongReadFileNames, String clusteredLongReadsDirectory,
            int sketchSize, int numThreads, float minCoverage, boolean useCompressedMinimizers) throws IOException, InterruptedException {

        final int minimizerWindowSize = k;
        
        int maxQueueSize = 100;
        ArrayBlockingQueue<Integer> sketchIndexQueue = new ArrayBlockingQueue<>(maxQueueSize);
        
        ArrayList<long[]> targetSketches = new ArrayList<>();
        ArrayList<Integer> targetSketchesAltIDs = new ArrayList<>();
        ArrayDeque<Integer> targetSketchesNullIndexes = new ArrayDeque<>();
        NTHashIterator itr = graph.getHashIterator();
        
        // skip l=0 because their lengths are too short to have a sketch
        for (int c=COVERAGE_ORDER.length-1; c>=0; --c) {
            for (int l=LENGTH_STRATUM_NAMES.length-1; l>0; --l) {
                
                String readFile = correctedLongReadFileNames[c][l];
                FastaReader fr = new FastaReader(readFile);
                System.out.println("Parsing file `" + readFile + "`...");
                
                while (fr.hasNext()) {
                    String[] nameSeqPair = fr.nextWithName();
                    String name = nameSeqPair[0];
                    String seq = nameSeqPair[1];
                    
                    int numKmers = getNumKmers(seq, k);
                    long[] sortedHashVals = useCompressedMinimizers ? 
                                            getAscendingHashValuesWithCompressedHomoPolymers(seq, itr, k) : 
                                            getAscendingHashValues(seq, itr, graph, numKmers, minCoverage);
                    
                    if (sortedHashVals.length < sketchSize) {
                        // not enough good kmers
                        continue;
                    }
                    
//                    if (sortedHashVals.length < minSketchOverlap) {
//                        // not enough good kmers
//                        continue;
//                    }                    

                    int bestTargetSketchID = -1;
                    
                    if (targetSketches.isEmpty()) {
                        long[] sketch = useCompressedMinimizers ?
                                        getMinimizersWithCompressedHomoPolymers(seq, k, itr, minimizerWindowSize) :
                                        getMinimizers(seq, numKmers, itr, minimizerWindowSize, graph, minCoverage);
                        targetSketches.add(sketch);
                        targetSketchesAltIDs.add(-1);
                        bestTargetSketchID = 0;
                    }
                    else {
                        float minSketchOverlapProportion = (float)sortedHashVals.length / (float)numKmers;
                        int minSketchOverlap = Math.max(1, (int) Math.floor(minSketchOverlapProportion * (sortedHashVals.length-minimizerWindowSize+1)/minimizerWindowSize));
                        
                        /** start thread pool*/
                        MyExecutorService service = new MyExecutorService(numThreads, numThreads);
                        
                        ContainmentCalculator[] workers = new ContainmentCalculator[numThreads];
                        for (int i=0; i<numThreads; ++i) {
                            ContainmentCalculator worker = new ContainmentCalculator(sortedHashVals, targetSketches, sketchIndexQueue, minSketchOverlap, minSketchOverlapProportion);
                            workers[i] = worker;
                            service.submit(worker);
                        }
                        
                        int numTargets = targetSketches.size();
                        for (int i=0; i<numTargets; ++i) {
                            if (targetSketches.get(i) != null) {
                                sketchIndexQueue.put(i);
                            }
                        }
                        
                        for (ContainmentCalculator worker : workers) {
                            worker.terminateWhenInputExhausts();
                        }
                        
                        service.terminate();
                        
                        float bestIntersectionSize = -1;
                        HashSet<Integer> overlapSketchIdSet = new HashSet<>();
                        for (ContainmentCalculator worker : workers) {
                            float mc = worker.getMaxIntersection();
                            if (mc > 0) {
                                ArrayDeque<Integer> overlaps = worker.getOverlappingSketchIndexes();
                                if (!overlaps.isEmpty()) {
                                    overlapSketchIdSet.addAll(overlaps);
                                }
                                
                                if (mc > bestIntersectionSize) {
                                    if (bestTargetSketchID >= 0) {
                                        overlapSketchIdSet.add(bestTargetSketchID);
                                    }
                                    bestIntersectionSize = mc;
                                    bestTargetSketchID = worker.getBestTargetIndex();
                                }
                            }
                        }
                        
                        if (bestIntersectionSize <= 0) {
                            // start a new cluster
                            long[] sketch = useCompressedMinimizers ?
                                            getMinimizersWithCompressedHomoPolymers(seq, k, itr, minimizerWindowSize) :
                                            getMinimizers(seq, numKmers, itr, minimizerWindowSize, graph, minCoverage);
                            
                            if (targetSketchesNullIndexes.isEmpty()) {
                                targetSketches.add(sketch);
                                targetSketchesAltIDs.add(-1);
                                bestTargetSketchID = targetSketches.size()-1;
                            }
                            else {
                                bestTargetSketchID = targetSketchesNullIndexes.poll();
                                targetSketches.set(bestTargetSketchID, sketch);
                                targetSketchesAltIDs.set(bestTargetSketchID, -1);
                            }
                        }
                        else {                            
                            if (c > 0 && !overlapSketchIdSet.isEmpty()) {
                                // combine overlapping clusters
                                                                
                                // combine FASTA files
//                                System.out.println("Combining " + overlapSketchIdSet.size() + " clusters into cluster " +  bestTargetSketchID);
                                int newTargetAltID = aggregateClusterFastas(bestTargetSketchID, overlapSketchIdSet, targetSketchesAltIDs, clusteredLongReadsDirectory);
                                targetSketchesAltIDs.set(bestTargetSketchID, newTargetAltID);
                                
                                
                                ArrayDeque<long[]> sketches = new ArrayDeque<>();
                                // add the sketch of the best target
                                sketches.add(targetSketches.get(bestTargetSketchID));
                                
                                // add the sketches of other targets
                                for (int i : overlapSketchIdSet) {
                                    sketches.add(targetSketches.get(i));
                                    targetSketches.set(i, null);
                                    targetSketchesNullIndexes.add(i);
                                    targetSketchesAltIDs.set(i, -1);
                                }
                                
                                long[] newSketch = combineSketches(sketches);
                                targetSketches.set(bestTargetSketchID, newSketch);
                            }
                        }
                    }
                    
                    /** add sequence to target cluster*/
                    String path = clusteredLongReadsDirectory + File.separator + bestTargetSketchID + FASTA_EXT;
                    FastaWriter writer = new FastaWriter(path, true);
                    writer.write(name, seq);
                    writer.close();
                }
                
                fr.close();
            }
        }
        
        //merge alt FASTA files
        System.out.println("Merging alternate cluster FASTAs ...");
        int targetSketchesSize = targetSketches.size();
        for (int i=0; i<targetSketchesSize; ++i) {
            int altID = targetSketchesAltIDs.get(i);
            if (altID >= 0) {
                mergeClusterFastas(i, altID, clusteredLongReadsDirectory);
            }
        }
        
        System.out.println("Reads are clustered in " + (targetSketches.size() - targetSketchesNullIndexes.size()) + " groups.");
    }
    
    public boolean assembleLongReads(String clusteredLongReadsDirectory, 
                                    String assembledLongReadsDirectory, 
                                    String assembledLongReadsCombined,
                                    int numThreads,
                                    boolean writeUracil,
                                    String minimapOptions,
                                    int minKmerCov,
                                    String txptNamePrefix,
                                    boolean stranded,
                                    int minTranscriptLength,
                                    boolean removeArtifacts) throws IOException {
        
        ArrayList<Integer> clusterIDs = new ArrayList<>();
        
        for (File f : new File(clusteredLongReadsDirectory).listFiles()) {
            String name = f.getName();
            if (name.endsWith(FASTA_EXT)) {
                clusterIDs.add(Integer.parseInt(name.substring(0, name.lastIndexOf(FASTA_EXT))));
            }
        }
        
        Collections.sort(clusterIDs);
        System.out.println("Total of " + clusterIDs.size() + " clusters to be assembled");
        
        ArrayList<Integer> errors = new ArrayList<>();
        for (int clusterID : clusterIDs) {
            String stampPath = assembledLongReadsDirectory + File.separator + clusterID + ".DONE";
            File stampFile = new File(stampPath);
            
            if (!stampFile.exists()) {
                String readsPath = clusteredLongReadsDirectory + File.separator + clusterID + FASTA_EXT;
                String tmpPrefix = assembledLongReadsDirectory + File.separator + clusterID;
                String concensusPath = assembledLongReadsDirectory + File.separator + clusterID + "_transcripts" + FASTA_EXT;

                System.out.println("Assembling cluster `" + clusterID + "`...");

                boolean ok = overlapLayoutConcensus(readsPath, tmpPrefix, concensusPath, numThreads, stranded, minimapOptions, maxIndelSize, removeArtifacts);
                if (!ok) {
                    System.out.println("*** Error assembling cluster `" + clusterID + "`!!! ***");
                    errors.add(clusterID);
                    //@TODO return false;
                }
                else {
                    touch(stampFile);
                }
            }
        }
        
        boolean ok = errors.isEmpty();
        String assembledLongReadsConcatenated = assembledLongReadsDirectory + File.separator + "all_transcripts" + FASTA_EXT;
        String tmpPrefix = assembledLongReadsDirectory + File.separator + "all_transcripts_overlap";
        
        if (ok) {
            Pattern raconRcPattern = Pattern.compile("RC:i:(\\d+)");
            
            // combine assembly files
            FastaWriter fout = new FastaWriter(assembledLongReadsConcatenated, false);
            FastaReader fin;
            for (int clusterID : clusterIDs) {
                String clusterAssemblyPath = assembledLongReadsDirectory + File.separator + clusterID + "_transcripts" + FASTA_EXT;
                fin = new FastaReader(clusterAssemblyPath);
                while(fin.hasNext()) {
                    String[] nameCommentSeq = fin.nextWithComment();
                    String comment = nameCommentSeq[1];
                    String seq = nameCommentSeq[2];
                                        
                    if (writeUracil) {
                        seq = seq.replace('T', 'U');
                    }
                    
                    String length = Integer.toString(seq.length());
                    
                    String coverage = "1";
                    if (!comment.isEmpty()) {
                        Matcher m = raconRcPattern.matcher(comment);
                        if (m.find()) {
                            coverage = m.group(1);
                        }
                    }
                    
                    fout.write(txptNamePrefix + clusterID + "_" + nameCommentSeq[0] +
                            " l=" + length + " c=" + coverage,
                            seq);
                }
                fin.close();
            }
            fout.close();
        }
        else {
            System.out.println("Cannot assemble the following clusters:");
            Collections.sort(errors);
            System.out.println(Arrays.toString(errors.toArray()));
            return false;
        }
        
        System.out.println("Inter-cluster assembly...");
        ok = overlapLayout(assembledLongReadsConcatenated, tmpPrefix, assembledLongReadsCombined, numThreads, stranded, "-r " + Integer.toString(2*maxIndelSize), k, percentIdentity, minTranscriptLength, maxIndelSize, removeArtifacts);
        
        return ok;
    }
    
    public static final int LENGTH_STRATUM_MIN_S_INDEX = 0;
    public static final int LENGTH_STRATUM_S_Q1_INDEX = 1;
    public static final int LENGTH_STRATUM_Q1_MED_INDEX = 2;
    public static final int LENGTH_STRATUM_MED_Q3_INDEX = 3;
    public static final int LENGTH_STRATUM_Q3_MAX_INDEX = 4;
    
    public static int getLongReadLengthStratumIndex(LengthStats stats, int minSeqLen, int testLen) {
        if (testLen < minSeqLen) {
            return LENGTH_STRATUM_MIN_S_INDEX;
        }
        else if (testLen < stats.q1) {
            return LENGTH_STRATUM_S_Q1_INDEX;
        }
        else if (testLen < stats.median) {
            return LENGTH_STRATUM_Q1_MED_INDEX;
        }
        else if (testLen < stats.q3) {
            return LENGTH_STRATUM_MED_Q3_INDEX;
        }
        else {
            return LENGTH_STRATUM_Q3_MAX_INDEX;
        }
    }
    
    private static final String[] LENGTH_STRATUM_NAMES = new String[]{"min_s", "s_q1", "q1_med", "med_q3", "q3_max"};
    
    public class CorrectedLongReadsWriterWorker implements Runnable {
        private final ArrayBlockingQueue<Sequence> inputQueue;
        private final int maxSampleSize;
        private final int minSeqLen;
        private final FastaWriter[][] writers;
        private LengthStats sampleLengthStats = null;
        private boolean terminateWhenInputExhausts = false;
        private long numCorrected = 0;
        private boolean successful = false;
        private Exception exception = null;
        private LongReadCorrectionWorker[] workersToWait = null;
        
        public CorrectedLongReadsWriterWorker(ArrayBlockingQueue<Sequence> inputQueue, FastaWriter[][] writers, int maxSampleSize, int minSeqLen) {
            this.inputQueue = inputQueue;
            this.writers = writers;
            this.maxSampleSize = maxSampleSize;
            this.minSeqLen = minSeqLen;
        }
                
        @Override
        public void run() {
            try {
                ArrayDeque<Sequence> sample = new ArrayDeque<>(maxSampleSize);

                while(!(terminateWhenInputExhausts && inputQueue.isEmpty())) {
                    Sequence seq = inputQueue.poll(1, TimeUnit.MILLISECONDS);
                    if (seq == null) {
                        continue;
                    }

                    sample.add(seq);

                    if (sample.size() == maxSampleSize) {
                        break;
                    }
                }

                int sampleSize = sample.size();
                int[] lengths = new int[sampleSize];
                int i = 0;
                for (Sequence seq : sample) {
                    lengths[i++] = seq.length;
                }
                                
                sampleLengthStats = getLengthStats(lengths);
                
                System.out.println("Corrected Read Lengths Distribution (n=" + sampleSize + ")");
                System.out.println("\tmin\tq1\tmed\tq3\tmax");
                System.out.println("\t" + sampleLengthStats.min + "\t" + 
                                            sampleLengthStats.q1 + "\t" +
                                            sampleLengthStats.median + "\t" +
                                            sampleLengthStats.q3 + "\t" +
                                            sampleLengthStats.max);

                /** write the sample sequences to file */
                for (Sequence seq : sample) {
                    
                    int lengthStratumIndex = getLongReadLengthStratumIndex(sampleLengthStats, minSeqLen, seq.length);
                    int covStatumIndex = getCoverageOrderOfMagnitude(seq.coverage);
                    
                    //String header = Long.toString(++numCorrected) + " l=" + Integer.toString(seq.length) + " c=" + Float.toString(seq.coverage);
                    ++numCorrected;
                    String header = seq.name + " l=" + Integer.toString(seq.length) + " c=" + Float.toString(seq.coverage);
                    
                    writers[covStatumIndex][lengthStratumIndex].write(header, seq.seq);
                }

                /** write the remaining sequences to file */
                while(!(terminateWhenInputExhausts && areDependentWorkersDone() && inputQueue.isEmpty())) {
                    Sequence seq = inputQueue.poll(1, TimeUnit.MILLISECONDS);
                    if (seq == null) {
                        continue;
                    }

                    int lengthStratumIndex = getLongReadLengthStratumIndex(sampleLengthStats, minSeqLen, seq.length);
                    int covStatumIndex = getCoverageOrderOfMagnitude(seq.coverage);

                    //String header = Long.toString(++numCorrected) + " l=" + Integer.toString(seq.length) + " c=" + Float.toString(seq.coverage);
                    ++numCorrected;
                    String header = seq.name + " l=" + Integer.toString(seq.length) + " c=" + Float.toString(seq.coverage);

                    writers[covStatumIndex][lengthStratumIndex].write(header, seq.seq);
                }
                
                successful = true;
            
            } catch (Exception ex) {
                System.out.println(ex.getMessage());
                exception = ex;
                successful = false;
            }
        }
        
        public void terminateWhenInputExhausts(LongReadCorrectionWorker[] workersToWait) {
            terminateWhenInputExhausts = true;
            this.workersToWait = workersToWait;
        }
        
        public boolean areDependentWorkersDone() {
            if (workersToWait == null) {
                return false;
            }
            else {
                int numRemaining = workersToWait.length;
                for (LongReadCorrectionWorker worker : workersToWait) {
                    if (worker.successful || worker.exception != null) {
                        --numRemaining;
                    }
                }
                
                return numRemaining <= 0;
            }
        }
        
        public boolean isSucessful() {
            return successful;
        }
        
        public Exception getExceptionCaught() {
            return exception;
        }
        
        public long getNumCorrected() {
            return numCorrected;
        }
        
        public LengthStats getSampleLengthStats() {
            return sampleLengthStats;
        }
    }
    
    public static class Sequence {
        String name;
        String seq;
        int length;
        float coverage;
        
        public Sequence(String name, String seq, int length, float coverage) {
            this.name = name;
            this.seq = seq;
            this.length = length;
            this.coverage = coverage;
        }
    }
    
    public class LongReadCorrectionWorker implements Runnable {
        private final ArrayBlockingQueue<String[]> inputQueue;
        private final ArrayBlockingQueue<Sequence> outputQueue;
        private boolean terminateWhenInputExhausts = false;
        private boolean successful = false;
        private Exception exception = null;
        private final int maxErrCorrItr;
        private final int minKmerCov;
        private final int minNumSolidKmers;
        private final boolean reverseComplement;
        private boolean trimArtifact = false;
        private long numArtifacts = 0;

        public LongReadCorrectionWorker(ArrayBlockingQueue<String[]> inputQueue,
                                        ArrayBlockingQueue<Sequence> outputQueue,
                                        int maxErrCorrItr, int minKmerCov, int minNumSolidKmers,
                                        boolean reverseComplement, boolean trimArtifact) {
            this.inputQueue = inputQueue;
            this.outputQueue = outputQueue;
            this.maxErrCorrItr = maxErrCorrItr;
            this.minKmerCov = minKmerCov;
            this.minNumSolidKmers = minNumSolidKmers;
            this.reverseComplement = reverseComplement;
            this.trimArtifact = trimArtifact;
        }
        
        @Override
        public void run() {
            while(!terminateWhenInputExhausts || !inputQueue.isEmpty()) {
                try {                    
                    String[] nameSeqPair = inputQueue.poll(1, TimeUnit.SECONDS);
                    if (nameSeqPair == null) {
                        continue;
                    }
                    
                    String seq = reverseComplement ? reverseComplement(nameSeqPair[1]) : nameSeqPair[1];
                    
                    ArrayList<Kmer> kmers = graph.getKmers(seq);
                    
                    ArrayList<Kmer> correctedKmers = correctLongSequence(kmers, 
                                                                        graph, 
                                                                        maxErrCorrItr, 
                                                                        maxCovGradient, 
                                                                        lookahead, 
                                                                        maxIndelSize, 
                                                                        percentIdentity, 
                                                                        minKmerCov,
                                                                        minNumSolidKmers,
                                                                        true);
                                        
                    if (correctedKmers != null && !correctedKmers.isEmpty()) {
                        if (trimArtifact) {
                            ArrayList<Kmer> trimmed = trimReverseComplementArtifact(correctedKmers,
                                                graph, strandSpecific, 150, maxIndelSize, percentIdentity, maxCovGradient);
                            if (trimmed.size() < correctedKmers.size()) {
                                ++numArtifacts;
                                correctedKmers = trimmed;
                            }
                        }
                        
                        float cov = getMinMedianKmerCoverage(correctedKmers, 200);
                        seq = graph.assemble(correctedKmers);
                        outputQueue.put(new Sequence(nameSeqPair[0], seq, seq.length(), cov));
                    }
                } catch (Exception ex) {
                    System.out.println(ex.getMessage());
                    ex.printStackTrace();
                    exception = ex;
                    successful = false;
                    return;
                }
            }
            
            successful = true;
        }
        
        public void terminateWhenInputExhausts() {
            terminateWhenInputExhausts = true;
        }
        
        public boolean isSucessful() {
            return successful;
        }
        
        public Exception getExceptionCaught() {
            return exception;
        }
    }
    
    public void correctLongReadsMultithreaded(String[] inputFastxPaths,
                                                FastaWriter[][] outFastaWriters,
                                                int minKmerCov,
                                                int maxErrCorrItr,
                                                int numThreads,
                                                int maxSampleSize,
                                                int minSeqLen,
                                                boolean reverseComplement,
                                                boolean trimArtifact) throws InterruptedException, IOException, Exception {
        long numReads = 0;
        
        MyExecutorService service = new MyExecutorService(numThreads+1, numThreads+1);
        
        int maxQueueSize = 100;
        ArrayBlockingQueue<String[]> inputQueue = new ArrayBlockingQueue<>(maxQueueSize);
        ArrayBlockingQueue<Sequence> outputQueue = new ArrayBlockingQueue<>(maxQueueSize);
                
        int numCorrectionWorkers = numThreads;
        LongReadCorrectionWorker[] correctionWorkers = new LongReadCorrectionWorker[numCorrectionWorkers];
        int minNumSolidKmers = 100;
        
        for (int i=0; i<numCorrectionWorkers; ++i) {
            LongReadCorrectionWorker worker = new LongReadCorrectionWorker(inputQueue, outputQueue, maxErrCorrItr, minKmerCov, minNumSolidKmers, reverseComplement, trimArtifact);
            correctionWorkers[i] = worker;
            service.submit(worker);
        }
        
        //FastaWriter writer = new FastaWriter(outFasta, true);
        CorrectedLongReadsWriterWorker writerWorker = new CorrectedLongReadsWriterWorker(outputQueue, outFastaWriters, maxSampleSize, minSeqLen);
        service.submit(writerWorker);
        
        for (FastxSequenceIterator itr = new FastxSequenceIterator(inputFastxPaths); itr.hasNext() ; ++numReads) {
            inputQueue.put(itr.nextWithName());
        }
        
        for (LongReadCorrectionWorker worker : correctionWorkers) {
            worker.terminateWhenInputExhausts();
        }
        
        writerWorker.terminateWhenInputExhausts(correctionWorkers);
        
        service.terminate();
        
        // check for errors
        if (!writerWorker.isSucessful()) {
            throw writerWorker.getExceptionCaught();
        }
        
        long numArtifacts = 0;
        for (LongReadCorrectionWorker worker : correctionWorkers) {
            if (!worker.isSucessful()) {
                throw worker.getExceptionCaught();
            }
            numArtifacts += worker.numArtifacts;
        }
        
        assert inputQueue.isEmpty() && outputQueue.isEmpty();
        
        System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " sequences.");
        long numCorrected = writerWorker.getNumCorrected();
        long numDiscarded = numReads - numCorrected;
        System.out.println("\tKept:      " + NumberFormat.getInstance().format(numCorrected) + "(" + numCorrected * 100f/numReads + "%)");
        System.out.println("\tDiscarded: " + NumberFormat.getInstance().format(numDiscarded) + "(" + numDiscarded * 100f/numReads + "%)");
        
        if (numArtifacts > 0) { 
            System.out.println("\tArtifacts: " + NumberFormat.getInstance().format(numArtifacts) + "(" + numArtifacts * 100f/numReads + "%)");
        }
    }
    
    public String polishSequence(String seq, int maxErrCorrItr, int minKmerCov, int minNumSolidKmers) {
        ArrayList<Kmer> kmers = graph.getKmers(seq);

        ArrayList<Kmer> correctedKmers = correctLongSequence(kmers, 
                                                            graph, 
                                                            maxErrCorrItr, 
                                                            maxCovGradient, 
                                                            lookahead, 
                                                            maxIndelSize, 
                                                            percentIdentity, 
                                                            minKmerCov,
                                                            minNumSolidKmers,
                                                            false);
        
        if (correctedKmers != null && !correctedKmers.isEmpty()) {
            return graph.assemble(correctedKmers);
        }
        
        return null;
    }

    private class FragmentWriters {
        boolean assemblePolyaTails = false;
        FastaWriter[] longFragmentsOut, shortFragmentsOut, unconnectedReadsOut, longPolyaFragmentsOut, shortPolyaFragmentsOut, unconnectedPolyaReadsOut;
        FastaWriter longSingletonsOut, shortSingletonsOut, unconnectedSingletonsOut, longPolyaSingletonsOut, shortPolyaSingletonsOut, unconnectedPolyaSingletonsOut;
        
        private FragmentWriters(FragmentPaths fragPaths, boolean assemblePolyaTails) throws IOException {
            this.assemblePolyaTails = assemblePolyaTails;
            longFragmentsOut = new FastaWriter[]{new FastaWriter(fragPaths.longFragmentsFastaPaths[0], true),
                                                new FastaWriter(fragPaths.longFragmentsFastaPaths[1], true),
                                                new FastaWriter(fragPaths.longFragmentsFastaPaths[2], true),
                                                new FastaWriter(fragPaths.longFragmentsFastaPaths[3], true),
                                                new FastaWriter(fragPaths.longFragmentsFastaPaths[4], true),
                                                new FastaWriter(fragPaths.longFragmentsFastaPaths[5], true)};

            shortFragmentsOut = new FastaWriter[]{new FastaWriter(fragPaths.shortFragmentsFastaPaths[0], true),
                                                new FastaWriter(fragPaths.shortFragmentsFastaPaths[1], true),
                                                new FastaWriter(fragPaths.shortFragmentsFastaPaths[2], true),
                                                new FastaWriter(fragPaths.shortFragmentsFastaPaths[3], true),
                                                new FastaWriter(fragPaths.shortFragmentsFastaPaths[4], true),
                                                new FastaWriter(fragPaths.shortFragmentsFastaPaths[5], true)};

            unconnectedReadsOut = new FastaWriter[]{new FastaWriter(fragPaths.unconnectedReadsFastaPaths[0], true),
                                                    new FastaWriter(fragPaths.unconnectedReadsFastaPaths[1], true),
                                                    new FastaWriter(fragPaths.unconnectedReadsFastaPaths[2], true),
                                                    new FastaWriter(fragPaths.unconnectedReadsFastaPaths[3], true),
                                                    new FastaWriter(fragPaths.unconnectedReadsFastaPaths[4], true),
                                                    new FastaWriter(fragPaths.unconnectedReadsFastaPaths[5], true)};

            longSingletonsOut = new FastaWriter(fragPaths.longSingletonsFasta, true);
            shortSingletonsOut = new FastaWriter(fragPaths.shortSingletonsFasta, true);
            unconnectedSingletonsOut = new FastaWriter(fragPaths.unconnectedSingletonsFasta, true);

            if (assemblePolyaTails) {
                longPolyaFragmentsOut = new FastaWriter[]{new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[0], true),
                                                            new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[1], true),
                                                            new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[2], true),
                                                            new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[3], true),
                                                            new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[4], true),
                                                            new FastaWriter(fragPaths.longPolyaFragmentsFastaPaths[5], true)};

                shortPolyaFragmentsOut = new FastaWriter[]{new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[0], true),
                                                            new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[1], true),
                                                            new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[2], true),
                                                            new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[3], true),
                                                            new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[4], true),
                                                            new FastaWriter(fragPaths.shortPolyaFragmentsFastaPaths[5], true)};

                unconnectedPolyaReadsOut = new FastaWriter[]{new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[0], true),
                                                            new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[1], true),
                                                            new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[2], true),
                                                            new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[3], true),
                                                            new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[4], true),
                                                            new FastaWriter(fragPaths.unconnectedPolyaReadsFastaPaths[5], true)};

                longPolyaSingletonsOut = new FastaWriter(fragPaths.longPolyaSingletonsFasta, true);
                shortPolyaSingletonsOut = new FastaWriter(fragPaths.shortPolyaSingletonsFasta, true);
                unconnectedPolyaSingletonsOut = new FastaWriter(fragPaths.unconnectedPolyaSingletonsFasta, true);
            }
        }
        
        public FastaWriter getWriter(Fragment frag, String seq) {
            if (frag.isUnconnectedRead) {
                boolean isPolya = assemblePolyaTails && polyATailPattern.matcher(seq).matches();

                if (frag.minCov == 1) {
                    if (isPolya) {
                        return unconnectedPolyaSingletonsOut;
                    }
                    else {
                        return unconnectedSingletonsOut;
                    }
                }
                else {
                    int m = getCoverageOrderOfMagnitude(frag.minCov);

                    if (isPolya) {
                        return unconnectedPolyaReadsOut[m];
                    }
                    else {
                        return unconnectedReadsOut[m];
                    }
                }
            }
            else {
                boolean isLong = frag.length >= longFragmentLengthThreshold;        
                boolean isPolya = assemblePolyaTails && polyATailPattern.matcher(seq).matches();

                if (frag.minCov == 1) {
                    if (isLong) {
                        if (isPolya) {
                            return longPolyaSingletonsOut;
                        }
                        else {
                            return longSingletonsOut;
                        }
                    }
                    else {
                        if (isPolya) {
                            return shortPolyaSingletonsOut;
                        }
                        else {
                            return shortSingletonsOut;
                        }
                    }
                }
                else {
                    int m = getCoverageOrderOfMagnitude(frag.minCov);

                    if (isLong) {
                        if (isPolya) {
                            return longPolyaFragmentsOut[m];
                        }
                        else {
                            return longFragmentsOut[m];
                        }
                    }
                    else {
                        if (isPolya) {
                            return shortPolyaFragmentsOut[m];
                        }
                        else {
                            return shortFragmentsOut[m];
                        }
                    }
                }
            }
        }
        
        public void closeAll() throws IOException {
            for (FastaWriter f : longFragmentsOut) {
                f.close();
            }
            
            for (FastaWriter f : shortFragmentsOut) {
                f.close();
            }
            
            for (FastaWriter f : unconnectedReadsOut) {
                f.close();
            }
            
            longSingletonsOut.close();
            shortSingletonsOut.close();
            unconnectedSingletonsOut.close();
            
            if (assemblePolyaTails) {
                for (FastaWriter f : longPolyaFragmentsOut) {
                    f.close();
                }

                for (FastaWriter f : shortPolyaFragmentsOut) {
                    f.close();
                }

                for (FastaWriter f : unconnectedPolyaReadsOut) {
                    f.close();
                }

                longPolyaSingletonsOut.close();
                shortPolyaSingletonsOut.close();
                unconnectedPolyaSingletonsOut.close();
            }
        }
    }
    
    private class FragmentWriterWorker implements Runnable {
        private int shortestFragmentLengthAllowed;
        private ArrayBlockingQueue<Fragment> fragments;
        private FragmentPaths outPaths;
        private FragmentWriters writers;
        private boolean keepGoing = true;
        private long readPairsConnected = 0;
        private long readPairsNotConnected = 0;
        
        private FragmentWriterWorker(ArrayBlockingQueue<Fragment> fragments, FragmentPaths outPaths, boolean assemblePolyaTails, int shortestFragmentLengthAllowed) throws IOException {
            this.fragments = fragments;
            this.outPaths = outPaths;
            this.writers = new FragmentWriters(outPaths, assemblePolyaTails);
            this.shortestFragmentLengthAllowed = shortestFragmentLengthAllowed;
        }
        
        @Override
        public void run() {
            long fragmentId = 0;
            long unconnectedReadId = 0;
            
            try {
                while (true) {
                    Fragment frag = fragments.poll(10, TimeUnit.MICROSECONDS);

                    if (frag == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
                        if (frag.isUnconnectedRead) {
                            ++readPairsNotConnected;
                            ++unconnectedReadId;

                            String seq = frag.left;
                            if (seq != null && seq.length() >= k) {
                                FastaWriter writer = writers.getWriter(frag, seq);
                                writer.write(Long.toString(unconnectedReadId) + "L", seq);
                            }

                            seq = frag.right;
                            if (seq != null && seq.length() >= k) {
                                FastaWriter writer = writers.getWriter(frag, seq);
                                writer.write(Long.toString(unconnectedReadId) + "R", seq);
                            }
                        }
                        else {
                            ++readPairsConnected;

                            int fragLen = frag.kmers.size()+ k - 1; // not using frag.length because it is the original length and frag may be extended
                            
                            if (fragLen >= shortestFragmentLengthAllowed && frag.minCov > 0) {
                                ArrayList<Kmer> fragKmers = frag.kmers;

                                if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                    graph.addPairedKmers(fragKmers);

                                    if (fragLen >= longFragmentLengthThreshold) {
                                        for (Kmer kmer : fragKmers) {
                                            screeningBf.add(kmer.getHash());
                                        }
                                    }

                                    String header = Long.toString(++fragmentId);
                                    if (debug) {
                                        header += " L=[" + frag.left + "] R=[" + frag.right + "]";
                                    }
                                    String seq = graph.assemble(fragKmers);

                                    FastaWriter writer = writers.getWriter(frag, seq);
                                    writer.write(header, seq);
                                }
                            }
                        }
                    }
                }
                
                writers.closeAll();
            }
            catch(Exception ex) {
                ex.printStackTrace();
                throw new RuntimeException(ex);
            }
        }
        
        public void stopWhenEmpty() {
            keepGoing = false;
        }
        
        public long getNumConnected() {
            return readPairsConnected;
        }
        
        public long getNumUnconnected() {
            return readPairsNotConnected;
        }
    }    
    
    private static class FragmentPaths {
        String[] longFragmentsFastaPaths;
        String[] shortFragmentsFastaPaths;
        String[] unconnectedReadsFastaPaths;
        String longSingletonsFasta;
        String shortSingletonsFasta;
        String unconnectedSingletonsFasta;
        String[] longPolyaFragmentsFastaPaths;
        String[] shortPolyaFragmentsFastaPaths;
        String[] unconnectedPolyaReadsFastaPaths;
        String longPolyaSingletonsFasta;
        String shortPolyaSingletonsFasta;
        String unconnectedPolyaSingletonsFasta;
        
        public FragmentPaths(String outdir, String name) {
            String longFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.long.";
            String shortFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.short.";
            String unconnectedReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.";
            String longPolyaFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.polya.long.";
            String shortPolyaFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.polya.short.";
            String unconnectedPolyaReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.polya.";
        
            longFragmentsFastaPaths = new String[]{longFragmentsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                        longFragmentsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                        longFragmentsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                        longFragmentsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                        longFragmentsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                        longFragmentsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};

            shortFragmentsFastaPaths = new String[]{shortFragmentsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};

            unconnectedReadsFastaPaths = new String[]{unconnectedReadsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};

            longPolyaFragmentsFastaPaths = new String[]{longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};

            shortPolyaFragmentsFastaPaths = new String[]{shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                                shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                                shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                                shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                                shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                                shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};

            unconnectedPolyaReadsFastaPaths = new String[]{unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[0] + FASTA_EXT,
                                                                unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[1] + FASTA_EXT,
                                                                unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[2] + FASTA_EXT,
                                                                unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[3] + FASTA_EXT,
                                                                unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[4] + FASTA_EXT,
                                                                unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[5] + FASTA_EXT};
        
            longSingletonsFasta = longFragmentsFastaPrefix + "01" + FASTA_EXT;
            shortSingletonsFasta = shortFragmentsFastaPrefix + "01" + FASTA_EXT;
            unconnectedSingletonsFasta = unconnectedReadsFastaPrefix + "01" + FASTA_EXT;
        
            longPolyaSingletonsFasta = longPolyaFragmentsFastaPrefix + "01" + FASTA_EXT;
            shortPolyaSingletonsFasta = shortPolyaFragmentsFastaPrefix + "01" + FASTA_EXT;
            unconnectedPolyaSingletonsFasta = unconnectedPolyaReadsFastaPrefix + "01" + FASTA_EXT;
        }
        
        public ArrayList asList(boolean assemblePolya) {
            ArrayList<String> paths = new ArrayList<>(longFragmentsFastaPaths.length + shortFragmentsFastaPaths.length + unconnectedReadsFastaPaths.length + 3);
            paths.addAll(Arrays.asList(longFragmentsFastaPaths));
            paths.addAll(Arrays.asList(shortFragmentsFastaPaths));
            paths.addAll(Arrays.asList(unconnectedReadsFastaPaths));
            paths.add(longSingletonsFasta);
            paths.add(shortSingletonsFasta);
            paths.add(unconnectedSingletonsFasta);

            if (assemblePolya) {
                paths.addAll(Arrays.asList(longPolyaFragmentsFastaPaths));
                paths.addAll(Arrays.asList(shortPolyaFragmentsFastaPaths));
                paths.addAll(Arrays.asList(unconnectedPolyaReadsFastaPaths));
                paths.add(longPolyaSingletonsFasta);
                paths.add(shortPolyaSingletonsFasta);
                paths.add(unconnectedPolyaSingletonsFasta);
            }
            
            return paths;
        }
        
        public void deleteAll() throws IOException {
            FileSystem sys = FileSystems.getDefault();
            
            for (String path : longFragmentsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            for (String path : shortFragmentsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            for (String path : unconnectedReadsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            for (String path : longPolyaFragmentsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            for (String path : shortPolyaFragmentsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            for (String path : unconnectedPolyaReadsFastaPaths) {
                Files.deleteIfExists(sys.getPath(path));
            }
            
            Files.deleteIfExists(sys.getPath(longSingletonsFasta));
            Files.deleteIfExists(sys.getPath(shortSingletonsFasta));
            Files.deleteIfExists(sys.getPath(unconnectedSingletonsFasta));
            Files.deleteIfExists(sys.getPath(longPolyaSingletonsFasta));
            Files.deleteIfExists(sys.getPath(shortPolyaSingletonsFasta));
            Files.deleteIfExists(sys.getPath(unconnectedPolyaSingletonsFasta));
        }
    }
        
    public int[] assembleFragmentsMultiThreaded(FastxFilePair[] fastqPairs, 
                                                FragmentPaths fragPaths,
                                                int bound,
                                                int minOverlap,
                                                int sampleSize, 
                                                int numThreads, 
                                                int maxErrCorrIterations,
                                                boolean extendFragments,
                                                int minKmerCov,
                                                boolean keepArtifact) throws FileFormatException, IOException, InterruptedException {
        
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }        
        
        long numParsed = 0;
        
        int shortestFragmentLengthAllowed = k;
        int leftReadLengthThreshold = k;
        int rightReadLengthThreshold = k;
        
        boolean assemblePolyaTails = this.minPolyATailLengthRequired > 0;
        
        ArrayBlockingQueue<PairedReadSegments> readPairs = new ArrayBlockingQueue<>(numThreads*2);
        ArrayBlockingQueue<Fragment> fragments = new ArrayBlockingQueue<>(sampleSize);
                
        FragmentAssembler[] workers = new FragmentAssembler[numThreads];
        Thread[] threads = new Thread[numThreads];
        for (int i=0; i<numThreads; ++i) {
            workers[i] = new FragmentAssembler(readPairs,
                                                fragments,
                                                bound,
                                                minOverlap,
                                                false, // do not store paired kmers
                                                maxErrCorrIterations, 
                                                leftReadLengthThreshold,
                                                rightReadLengthThreshold,
                                                extendFragments,
                                                minKmerCov,
                                                keepArtifact
                                               );
            threads[i] = new Thread(workers[i]);
            threads[i].start();
        }
        
        FastxPairSequenceIterator rin = new FastxPairSequenceIterator(fastqPairs, this.seqPattern, this.qualPatternFrag);
        while (rin.hasNext()) {
            if (fragments.remainingCapacity() == 0) {
                break;
            }
            
            readPairs.put(rin.next());
            ++numParsed;
        }
        
        // Calculate length stats
        ArrayList<Integer> fragLengths = new ArrayList<>(sampleSize);
        for (Fragment frag : fragments) {
            if (!frag.isUnconnectedRead) {
                fragLengths.add(frag.length);
            }
        }

        if (fragLengths.isEmpty()) {
            // use lengths of unconnected reads as fragment lengths
            for (Fragment frag : fragments) {
                fragLengths.add(frag.left.length());
                fragLengths.add(frag.right.length());
            }
        }
        
        int[] fragLengthsStats = getMinQ1MedianQ3Max(fragLengths);
        
        longFragmentLengthThreshold = fragLengthsStats[1];
        setPairedKmerDistance(longFragmentLengthThreshold);
        
        int newBound = getPairedReadsMaxDistance(fragLengthsStats);
        for (FragmentAssembler w : workers) {
            w.updateBound(newBound);
        }
        
        System.out.println("Fragment Lengths Distribution (n=" + fragLengths.size() + ")");
        System.out.println("\tmin\tQ1\tM\tQ3\tmax");
        System.out.println("\t" + fragLengthsStats[0] + "\t" + fragLengthsStats[1] + "\t" + fragLengthsStats[2] + "\t" + fragLengthsStats[3] + "\t" + fragLengthsStats[4]);
        System.out.println("Paired kmers distance:       " + (longFragmentLengthThreshold - k - minNumKmerPairs));
        System.out.println("Max graph traversal depth:   " + newBound);

        FragmentWriterWorker writerWorker = new FragmentWriterWorker(fragments, fragPaths, assemblePolyaTails, shortestFragmentLengthAllowed);
        Thread writerThread = new Thread(writerWorker);
        writerThread.start();
        
        while (rin.hasNext()) {
            readPairs.put(rin.next());
            ++numParsed;
        }
        
        for (FragmentAssembler w : workers) {
            w.stopWhenEmpty();
        }

        for (Thread t : threads) {
            t.join();
        }

        writerWorker.stopWhenEmpty();
        writerThread.join();
        
        long numConnected = writerWorker.getNumConnected();
        long numUnconnected = writerWorker.getNumUnconnected();
        long numDiscarded = numParsed - numConnected - numUnconnected;

        System.out.println("Parsed " + NumberFormat.getInstance().format(numParsed) + " read pairs.");
        System.out.println("\tconnected:\t" + NumberFormat.getInstance().format(numConnected) + "\t(" + numConnected*100f/numParsed + "%)");
        System.out.println("\tnot connected:\t" + NumberFormat.getInstance().format(numUnconnected) + "\t(" + numUnconnected*100f/numParsed + "%)");
        System.out.println("\tdiscarded:\t" + NumberFormat.getInstance().format(numDiscarded) + "\t(" + numDiscarded*100f/numParsed + "%)");
        System.out.println("Fragments paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
        System.out.println("Screening Bloom filter FPR:              " + screeningBf.getFPR() * 100 + " %");

        return fragLengthsStats;
    }

    public void updateGraphDesc(File graphFile) throws IOException {
        graph.saveDesc(graphFile);
    }
    
    public void savePairedKmersBloomFilter(File graphFile) throws IOException {
        graph.savePkbf(graphFile);
    }
    
    public void restorePairedKmersBloomFilter(File graphFile) throws IOException {
        graph.destroyPkbf();
        graph.restorePkbf(graphFile);
        graph.updateFragmentKmerDistance(graphFile);
    }
    
    private long assembleTranscriptsMultiThreadedHelper(String fragmentsFasta, 
                                                    TranscriptWriter writer, 
                                                    int sampleSize, 
                                                    int numThreads, 
                                                    boolean includeNaiveExtensions,
                                                    boolean extendBranchFreeFragmentsOnly,
                                                    boolean keepArtifact,
                                                    boolean keepChimera,
                                                    boolean reqFragKmersConsistency,
                                                    float minKmerCov) throws InterruptedException, IOException {
        
        long numFragmentsParsed = 0;
        FastaReader fin = new FastaReader(fragmentsFasta);

        ArrayBlockingQueue<String> fragmentsQueue = new ArrayBlockingQueue<>(numThreads*2, true);
        ArrayBlockingQueue<Transcript> transcriptsQueue = new ArrayBlockingQueue<>(numThreads*2, true);
        
        TranscriptAssemblyWorker[] workers = new TranscriptAssemblyWorker[numThreads];
        Thread[] threads = new Thread[numThreads];
        for (int i=0; i<numThreads; ++i) {
            workers[i] = new TranscriptAssemblyWorker(fragmentsQueue, transcriptsQueue, includeNaiveExtensions, extendBranchFreeFragmentsOnly, keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            threads[i] = new Thread(workers[i]);
            threads[i].start();
        }
        
        TranscriptWriterWorker writerWorker = new TranscriptWriterWorker(transcriptsQueue, writer);
        Thread writerThread = new Thread(writerWorker);
        writerThread.start();

        while (fin.hasNext()) {
            fragmentsQueue.put(fin.next());
            ++numFragmentsParsed;
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

    public void assembleSingleEndReads(String[] readPaths,
                                    boolean reverseComplement,
                                    String outFasta,
                                    String outFastaShort,
                                    int numThreads,
                                    int minTranscriptLength,
                                    boolean keepArtifact,
                                    boolean keepChimera,
                                    String txptNamePrefix,
                                    float minKmerCov,
                                    boolean writeUracil) throws IOException, InterruptedException {

        //boolean assemblePolya = minPolyATailLengthRequired > 0;
        /*@TODO support prioritized assembly of polya reads */
        
        boolean includeNaiveExtensions = true;
        boolean extendBranchFreeFragmentsOnly = false;

        FastaWriter fout = new FastaWriter(outFasta, false);
        FastaWriter foutShort = new FastaWriter(outFastaShort, false);
        TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength, maxTipLength, writeUracil);
        
        ArrayBlockingQueue<String> readsQueue = new ArrayBlockingQueue<>(numThreads*2);
        ArrayBlockingQueue<Transcript> transcriptsQueue = new ArrayBlockingQueue<>(numThreads*2);
        
        TranscriptAssemblyWorker[] workers = new TranscriptAssemblyWorker[numThreads];
        Thread[] threads = new Thread[numThreads];
        for (int i=0; i<numThreads; ++i) {
            workers[i] = new TranscriptAssemblyWorker(readsQueue, transcriptsQueue, includeNaiveExtensions, extendBranchFreeFragmentsOnly, keepArtifact, keepChimera, false, minKmerCov);
            threads[i] = new Thread(workers[i]);
            threads[i].start();
        }
        
        TranscriptWriterWorker writerWorker = new TranscriptWriterWorker(transcriptsQueue, writer);
        Thread writerThread = new Thread(writerWorker);
        writerThread.start();
        
        ArrayList<String> fastaPaths = new ArrayList<>();
        ArrayList<String> fastqPaths = new ArrayList<>();
        for (String p : readPaths) {
            if (FastaReader.isCorrectFormat(p)) {
                fastaPaths.add(p);
            }
            else if (FastqReader.isCorrectFormat(p)) {
                fastqPaths.add(p);
            }
        }
        
        if (!fastaPaths.isEmpty()) {
            String[] paths = new String[fastaPaths.size()];
            fastaPaths.toArray(paths);
            FastaFilteredSequenceIterator rin = new FastaFilteredSequenceIterator(paths, seqPattern, reverseComplement);

            while (rin.hasNext()) {
                readsQueue.put(rin.next());
            }
        }
        
        if (!fastqPaths.isEmpty()) {
            String[] paths = new String[fastqPaths.size()];
            fastqPaths.toArray(paths);
            FastqFilteredSequenceIterator rin = new FastqFilteredSequenceIterator(paths, seqPattern, qualPatternFrag, reverseComplement);

            while (rin.hasNext()) {
                readsQueue.put(rin.next());
            }
        }
        
        for (TranscriptAssemblyWorker w : workers) {
            w.stopWhenEmpty();
        }

        for (Thread t : threads) {
            t.join();
        }

        writerWorker.stopWhenEmpty();
        writerThread.join();
        
        fout.close();
        foutShort.close();
    }
    
    public void assembleTranscriptsMultiThreaded(FragmentPaths fragPaths,
                                                String outFasta,
                                                String outFastaShort,
                                                String graphFile,
                                                int numThreads,
                                                int sampleSize,
                                                int minTranscriptLength,
                                                boolean keepArtifact,
                                                boolean keepChimera,
                                                boolean reqFragKmersConsistency,
                                                String txptNamePrefix,
                                                float minKmerCov,
                                                String branchFreeExtensionThreshold,
                                                boolean writeUracil) throws IOException, InterruptedException {
        
        long numFragmentsParsed = 0;

        boolean assemblePolya = minPolyATailLengthRequired > 0;

        FastaWriter fout = new FastaWriter(outFasta, false);
        FastaWriter foutShort = new FastaWriter(outFastaShort, false);
        //TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength, sensitiveMode ? maxTipLength : Math.max(k, maxTipLength));
        TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength, maxTipLength, writeUracil);


        boolean allowNaiveExtension = true;
        boolean extendBranchFreeOnly;

        if (assemblePolya) {
            // extend LONG fragments
            for (int mag=fragPaths.longPolyaFragmentsFastaPaths.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
                String fragmentsFasta = fragPaths.longPolyaFragmentsFastaPaths[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly, 
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }

            // extend SHORT fragments
            for (int mag=fragPaths.shortPolyaFragmentsFastaPaths.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
                String fragmentsFasta = fragPaths.shortPolyaFragmentsFastaPaths[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }

            // extend UNCONNECTED reads
            for (int mag=fragPaths.unconnectedPolyaReadsFastaPaths.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
                String fragmentsFasta = fragPaths.unconnectedPolyaReadsFastaPaths[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }

            extendBranchFreeOnly = isLowerStratum(STRATUM_01, branchFreeExtensionThreshold);

            // extend LONG singleton fragments
            writer.setOutputPrefix(txptNamePrefix + "01.L.");
            System.out.println("Parsing `" + fragPaths.longPolyaSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.longPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

            // extend SHORT singleton fragments
            writer.setOutputPrefix(txptNamePrefix + "01.S.");
            System.out.println("Parsing `" + fragPaths.shortPolyaSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.shortPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

            // extend UNCONNECTED reads
            writer.setOutputPrefix(txptNamePrefix + "01.U.");
            System.out.println("Parsing `" + fragPaths.unconnectedPolyaSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.unconnectedPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
        }


        // extend LONG fragments
        for (int mag=fragPaths.longFragmentsFastaPaths.length-1; mag>=0; --mag) {
            writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
            String fragmentsFasta = fragPaths.longFragmentsFastaPaths[mag];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly, 
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
        }          

        // extend SHORT fragments
        for (int mag=fragPaths.shortFragmentsFastaPaths.length-1; mag>=0; --mag) {
            writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
            String fragmentsFasta = fragPaths.shortFragmentsFastaPaths[mag];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
        }

        // extend UNCONNECTED reads
        for (int mag=fragPaths.unconnectedReadsFastaPaths.length-1; mag>=0; --mag) {
            writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
            String fragmentsFasta = fragPaths.unconnectedReadsFastaPaths[mag];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
        }

        extendBranchFreeOnly = isLowerStratum(STRATUM_01, branchFreeExtensionThreshold);

        // extend LONG singleton fragments
        writer.setOutputPrefix(txptNamePrefix + "01.L.");
        System.out.println("Parsing `" + fragPaths.longSingletonsFasta + "`...");
        numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.longSingletonsFasta, writer, sampleSize, numThreads,
                                                                allowNaiveExtension, extendBranchFreeOnly,
                                                                keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

        // extend SHORT singleton fragments
        writer.setOutputPrefix(txptNamePrefix + "01.S.");
        System.out.println("Parsing `" + fragPaths.shortSingletonsFasta + "`...");
        numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.shortSingletonsFasta, writer, sampleSize, numThreads,
                                                                allowNaiveExtension, extendBranchFreeOnly,
                                                                keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

        // extend UNCONNECTED reads
        writer.setOutputPrefix(txptNamePrefix + "01.U.");
        System.out.println("Parsing `" + fragPaths.unconnectedSingletonsFasta + "`...");
        numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragPaths.unconnectedSingletonsFasta, writer, sampleSize, numThreads,
                                                                allowNaiveExtension, extendBranchFreeOnly,
                                                                keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

        fout.close();
        foutShort.close();

        System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
        System.out.println("Screening Bloom filter FPR:      " + screeningBf.getFPR() * 100 + " %");
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
        printVersionInfo(false);
        System.out.println();
        
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
    
    public static void printVersionInfo(boolean exit) {
        System.out.println(
                "RNA-Bloom v" + VERSION + "\n" +
                "Ka Ming Nip, Canada's Michael Smith Genome Sciences Centre, BC Cancer\n" +
                "Copyright 2018"
        );
        
        if (exit) {
            System.exit(0);
        }
    }

    private final static String FIELD_SEPARATOR = "\\s+"; // any white space character
    
    private static boolean getPooledReadPaths(String pooledReadPathsListFile,
            HashMap<String, ArrayList<String>> pooledLeftReadPaths,
            HashMap<String, ArrayList<String>> pooledRightReadPaths) throws FileNotFoundException, IOException {
        
        BufferedReader br = new BufferedReader(new FileReader(pooledReadPathsListFile));
        
        String line;
        int lineNumber = 0;
        while ((line = br.readLine()) != null) {
            ++lineNumber;
            line = line.trim();
            
            if (line.isEmpty()) {
                continue;
            }
            
            String[] entry = line.split(FIELD_SEPARATOR);
            
            //SAMPLE_ID LEFT_PATH RIGHT_PATH
            if (entry.length == 3) {
                String id = entry[0];

                ArrayList<String> paths = pooledLeftReadPaths.get(id);
                if (paths == null) {
                    paths = new ArrayList<>();
                    pooledLeftReadPaths.put(id, paths);
                }
                paths.add(entry[1]);

                paths = pooledRightReadPaths.get(id);
                if (paths == null) {
                    paths = new ArrayList<>();
                    pooledRightReadPaths.put(id, paths);
                }
                paths.add(entry[2]);
            }
            else {
                exitOnError("Pool reads path file has unexpected number of columns on line " + lineNumber + ":\n\t" + line);
                return false;
            }
        }
        
        br.close();
        return true;
    }
    
    private static void correctLongReads(RNABloom assembler, 
            String[] readFastxPaths, String[][] correctedLongReadFileNames, 
            int maxErrCorrItr, int minKmerCov, int numThreads, int sampleSize, int minSeqLen, boolean reverseComplement, boolean trimArtifact) throws InterruptedException, IOException, Exception {
        
        /* set up the file writers */
        final int numCovStrata = COVERAGE_ORDER.length;
        final int numLenStrata = LENGTH_STRATUM_NAMES.length;
        FastaWriter[][] writers = new FastaWriter[numCovStrata][numLenStrata];
        for (int c=0; c<numCovStrata; ++c) {
            for (int l=0; l<numLenStrata; ++l) {
                writers[c][l] = new FastaWriter(correctedLongReadFileNames[c][l], true);
            }
        }
        
        assembler.correctLongReadsMultithreaded(readFastxPaths, writers, minKmerCov, maxErrCorrItr, numThreads, sampleSize, minSeqLen, reverseComplement, trimArtifact);
        
        for (int i=0; i<writers.length; ++i) {
            for (int j=0; j<writers[i].length; ++j) {
                writers[i][j].close();
            }
        }
    }
    
    private static void clusterLongReads(RNABloom assembler, 
            String[][] correctedLongReadFileNames, String clusteredLongReadsDirectory,
            int sketchSize, int numThreads, float minKmerCov, boolean useCompressedMinimizers) throws IOException, InterruptedException {
        
        File outdir = new File(clusteredLongReadsDirectory);
        if (outdir.exists()) {
            for (File f : outdir.listFiles()) {
                f.delete();
            }
        }
        else {
            outdir.mkdirs();
        }
        
        assembler.clusterLongReads(correctedLongReadFileNames, clusteredLongReadsDirectory, sketchSize, numThreads, minKmerCov, useCompressedMinimizers);
    }
    
    private static boolean assembleLongReads(RNABloom assembler, 
            String clusteredLongReadsDirectory, String assembledLongReadsDirectory,
            String assembledLongReadsCombined,
            int numThreads, boolean forceOverwrite,
            boolean writeUracil, String minimapOptions, int minKmerCov, String txptNamePrefix, 
            boolean stranded, int minTranscriptLength, boolean removeArtifacts) throws IOException {
        
        File outdir = new File(assembledLongReadsDirectory);
        if (outdir.exists()) {
            if (forceOverwrite) {
                for (File f : outdir.listFiles()) {
                    f.delete();
                }
            }
        }
        else {
            outdir.mkdirs();
        }
        
        return assembler.assembleLongReads(clusteredLongReadsDirectory, assembledLongReadsDirectory, assembledLongReadsCombined,
                numThreads, writeUracil, minimapOptions, minKmerCov, txptNamePrefix, stranded, minTranscriptLength, removeArtifacts);
    }
    
    private static void assembleFragments(RNABloom assembler, boolean forceOverwrite,
            String outdir, String name, FastxFilePair[] fqPairs,
            long sbfSize, long pkbfSize, int sbfNumHash, int pkbfNumHash, int numThreads,
            int bound, int minOverlap, int sampleSize, int maxErrCorrItr, boolean extendFragments,
            int minKmerCoverage, boolean keepArtifact) throws FileFormatException, IOException, InterruptedException {

        final File fragsDoneStamp = new File(outdir + File.separator + STAMP_FRAGMENTS_DONE);
        
        if (forceOverwrite || !fragsDoneStamp.exists()) {
            FragmentPaths fragPaths = new FragmentPaths(outdir, name);
            fragPaths.deleteAll();

            assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
            assembler.setupFragmentPairedKmersBloomFilter(pkbfSize, pkbfNumHash);

            int[] fragStats = assembler.assembleFragmentsMultiThreaded(fqPairs, 
                                                                        fragPaths,
                                                                        bound, 
                                                                        minOverlap,
                                                                        sampleSize,
                                                                        numThreads,
                                                                        maxErrCorrItr,
                                                                        extendFragments,
                                                                        minKmerCoverage,
                                                                        keepArtifact);

            String fragStatsFile = outdir + File.separator + name + ".fragstats";
            String graphFile = outdir + File.separator + name + ".graph";
            
            assembler.updateGraphDesc(new File(graphFile));
            assembler.writeFragStatsToFile(fragStats, fragStatsFile);

            touch(fragsDoneStamp);
        }
        else {
            System.out.println("WARNING: Fragments were already assembled for \"" + name + "!");
        }
    }
    
    private static void splitFastaByLength(String inFasta, String outLongFasta, String outShortFasta, int lengthThreshold) throws IOException {
        FastaReader fin = new FastaReader(inFasta);
        FastaWriter foutLong = new FastaWriter(outLongFasta, false);
        FastaWriter foutShort = new FastaWriter(outShortFasta, false);
        while(fin.hasNext()) {
            String[] nameCommentSeq = fin.nextWithComment();
            
            String header = nameCommentSeq[0];
            if (!nameCommentSeq[1].isEmpty()) {
                header += " " + nameCommentSeq[1];
            }
            
            String seq = nameCommentSeq[2];
            
            if (seq.length() >= lengthThreshold) {
                foutLong.write(header, seq);
            }
            else {
                foutShort.write(header, seq);
            }
        }
        foutLong.close();
        foutShort.close();
        fin.close();
    }
    
    private static boolean mergePooledAssemblies(String outdir, String assemblyName, String[] sampleNames,
            String txptFileExt, String shortTxptFileExt, String txptNamePrefix,
            int k, int numThreads, boolean stranded, int maxIndelSize, int maxTipLength,
            float percentIdentity, boolean removeArtifacts, int txptLengthThreshold, boolean writeUracil) throws IOException {
        
        String concatenatedFasta = outdir + File.separator + "all_transcripts" + FASTA_EXT;
        String reducedFasta = outdir + File.separator + "all_transcripts_nr" + FASTA_EXT;
        String tmpPrefix = outdir + File.separator + "all_transcripts_overlap";
        String outLongFasta = outdir + File.separator + assemblyName + ".transcripts" + FASTA_EXT;
        String outShortFasta = outdir + File.separator + assemblyName + ".transcripts.short" + FASTA_EXT;

        // combine assembly files
        FastaWriter fout = new FastaWriter(concatenatedFasta, false);
        FastaReader fin;
        for (String sampleName : sampleNames) {
            String longTxptsPath = outdir + File.separator + sampleName + File.separator + sampleName + ".transcripts" + txptFileExt;
            String shortTxptsPath = outdir + File.separator + sampleName + File.separator + sampleName + ".transcripts" + shortTxptFileExt;
            
            fin = new FastaReader(longTxptsPath);
            while(fin.hasNext()) {
                String[] nameCommentSeq = fin.nextWithComment();
                //String comment = nameCommentSeq[1];
                String seq = nameCommentSeq[2];

                if (writeUracil) {
                    seq = seq.replace('T', 'U');
                }

                fout.write(txptNamePrefix + sampleName + "_" + nameCommentSeq[0], seq);
            }
            fin.close();
            
            fin = new FastaReader(shortTxptsPath);
            while(fin.hasNext()) {
                String[] nameCommentSeq = fin.nextWithComment();
                //String comment = nameCommentSeq[1];
                String seq = nameCommentSeq[2];

                if (writeUracil) {
                    seq = seq.replace('T', 'U');
                }

                fout.write(txptNamePrefix + sampleName + "_" + nameCommentSeq[0], seq);
            }
            fin.close();
        }
        fout.close();

        boolean ok = overlapLayout(concatenatedFasta, tmpPrefix, reducedFasta, numThreads,
                        stranded, "-r " + Integer.toString(maxIndelSize), maxTipLength, percentIdentity, 2*k,
                        maxIndelSize, removeArtifacts);
        
        splitFastaByLength(reducedFasta, outLongFasta, outShortFasta, txptLengthThreshold);
        
        Files.deleteIfExists(FileSystems.getDefault().getPath(concatenatedFasta));
        Files.deleteIfExists(FileSystems.getDefault().getPath(reducedFasta));
        
        return ok;
    }
    
    private static void assembleTranscripts(RNABloom assembler, boolean forceOverwrite,
            String outdir, String name, String txptNamePrefix, boolean strandSpecific,
            long sbfSize, int sbfNumHash, long pkbfSize, int pkbfNumHash, int numThreads, boolean noFragDBG,
            int sampleSize, int minTranscriptLength, boolean keepArtifact, boolean keepChimera, 
            boolean reqFragKmersConsistency, boolean restorePairedKmers,
            float minKmerCov, String branchFreeExtensionThreshold,
            boolean reduceRedundancy, boolean assemblePolya, boolean writeUracil,
            String[] refTranscriptFile) throws IOException, InterruptedException {
        
        final File txptsDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_DONE);
        final String transcriptsFasta = outdir + File.separator + name + ".transcripts" + FASTA_EXT;
        final String shortTranscriptsFasta = outdir + File.separator + name + ".transcripts.short" + FASTA_EXT;
                
        if (forceOverwrite || !txptsDoneStamp.exists()) {
            MyTimer timer = new MyTimer();

            FragmentPaths fragPaths = new FragmentPaths(outdir, name);
            
            final String graphFile = outdir + File.separator + name + ".graph";
            
            if (!noFragDBG) {
                if (assembler.isGraphInitialized()) {
                    assembler.clearDbgBf();
                }
                
                System.out.println("Rebuilding graph from assembled fragments...");
                timer.start();
                if (restorePairedKmers) {
                    assembler.setupFragmentPairedKmersBloomFilter(pkbfSize, pkbfNumHash);
                    assembler.updateFragmentKmerDistance(graphFile);
                }
                assembler.populateGraphFromFragments(fragPaths.asList(assemblePolya), strandSpecific, restorePairedKmers);
                System.out.println("Graph rebuilt in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                if (refTranscriptFile != null && refTranscriptFile.length > 0) {
                    System.out.println("Augmenting graph with reference transcripts...");
                    timer.start();
                    assembler.addPairedKmersFromSequences(refTranscriptFile);
                    System.out.println("Graph augmented in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                }
            }

            
            Files.deleteIfExists(FileSystems.getDefault().getPath(transcriptsFasta));
            Files.deleteIfExists(FileSystems.getDefault().getPath(shortTranscriptsFasta));

            System.out.println("Assembling transcripts...");
            timer.start();

            assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);

            assembler.assembleTranscriptsMultiThreaded(fragPaths,
                                                        transcriptsFasta, 
                                                        shortTranscriptsFasta,
                                                        graphFile,
                                                        numThreads,
                                                        sampleSize,
                                                        minTranscriptLength,
                                                        keepArtifact,
                                                        keepChimera,
                                                        reqFragKmersConsistency,
                                                        txptNamePrefix,
                                                        minKmerCov,
                                                        branchFreeExtensionThreshold,
                                                        writeUracil);
            
            System.out.println("Transcripts assembled in " + MyTimer.hmsFormat(timer.elapsedMillis()));

            touch(txptsDoneStamp);

            System.out.println("Assembled transcripts at `" + transcriptsFasta + "`");
        }
        else {
            System.out.println("WARNING: Transcripts were already assembled for \"" + name + "\"!");
        }
                
        if (reduceRedundancy) {
            final File nrTxptsDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_NR_DONE);
            
            if (forceOverwrite || !nrTxptsDoneStamp.exists()) {
                String tmpPrefix = outdir + File.separator + name + ".tmp";
                String nrTranscriptsFasta = outdir + File.separator + name + ".transcripts.nr" + FASTA_EXT;
                String shortNrTranscriptsFasta = outdir + File.separator + name + ".transcripts.nr.short" + FASTA_EXT;
                
                Files.deleteIfExists(FileSystems.getDefault().getPath(nrTranscriptsFasta));
                
                System.out.println("Reducing redundancy in assembled transcripts...");
                MyTimer timer = new MyTimer();
                timer.start();
                
                boolean ok = assembler.generateNonRedundantTranscripts(transcriptsFasta, shortTranscriptsFasta,
                        tmpPrefix, nrTranscriptsFasta, shortNrTranscriptsFasta, numThreads, !keepArtifact, minTranscriptLength);

                if (ok) {
                    System.out.println("Redundancy reduced in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                    touch(nrTxptsDoneStamp);
                }
                else {
                    exitOnError("Error during redundancy reduction!");
                }
            }
            else {
                System.out.println("WARNING: Redundancy reduction had completed previously for \"" + name + "\"!");
            }                
        }
    }
    
    private boolean generateNonRedundantTranscripts(String inLongFastas, String inShortFasta,
            String tmpPrefix, String outLongFasta, String outShortFasta,
            int numThreads, boolean removeArtifacts, int txptLengthThreshold) throws IOException {
        
        String concatenatedFasta = tmpPrefix + "_ava_cat" + FASTA_EXT;
        String reducedFasta = tmpPrefix + "_ava_cat_nr" + FASTA_EXT;

        // combine assembly files
        FastaWriter fout = new FastaWriter(concatenatedFasta, false);
        
        FastaReader fin = new FastaReader(inLongFastas);
        while(fin.hasNext()) {
            String[] nameCommentSeq = fin.nextWithComment();
            //String comment = nameCommentSeq[1];
            String seq = nameCommentSeq[2];
            fout.write(nameCommentSeq[0], seq);
        }
        fin.close();
        
        fin = new FastaReader(inShortFasta);
        while(fin.hasNext()) {
            String[] nameCommentSeq = fin.nextWithComment();
            //String comment = nameCommentSeq[1];
            String seq = nameCommentSeq[2];
            fout.write(nameCommentSeq[0], seq);
        }
        fin.close();
        
        fout.close();
        
        boolean ok = overlapLayout(concatenatedFasta, tmpPrefix, reducedFasta, 
                        numThreads, strandSpecific, "-r " + Integer.toString(maxIndelSize),
                        maxTipLength, percentIdentity, 2*k, maxIndelSize, removeArtifacts);
        
        splitFastaByLength(reducedFasta, outLongFasta, outShortFasta, txptLengthThreshold);
        
        Files.deleteIfExists(FileSystems.getDefault().getPath(concatenatedFasta));
        Files.deleteIfExists(FileSystems.getDefault().getPath(reducedFasta));
        
        return ok;
    }
    
    private static boolean hasNtcard() {
        try {
            String cmd = "ntcard --version";
            Runtime rt = Runtime.getRuntime();
            Process pr = rt.exec(cmd);
            int exitVal = pr.waitFor();
            return exitVal == 0;
        }
        catch (IOException | InterruptedException e) {
            return false;
        }
    }
        
    private static NTCardHistogram getNTCardHistogram(int threads, int k, String histogramPathPrefix, String readPathsFile, boolean forceOverwrite) throws IOException, InterruptedException {
        NTCardHistogram hist = null;
        String histogramPath = histogramPathPrefix + "_k" + k + ".hist";
        int exitVal = 0;
        
        if (forceOverwrite || !new File(histogramPath).isFile()) {
            String cmd = "ntcard -t " + threads + " -k " + k + " -c 65535 -p " + histogramPathPrefix + " @" + readPathsFile;
            Runtime rt = Runtime.getRuntime();
            System.out.println("Running command: `" + cmd + "`...");
            Process pr = rt.exec(cmd);
            exitVal = pr.waitFor();
        }
        
        if (exitVal == 0) {
            System.out.println("Parsing histogram file `" + histogramPath + "`...");
            hist = new NTCardHistogram(histogramPath);
            
        }
        
        return hist;
    }
    
    private static void printBloomFilterMemoryInfo(float dbgGB, float cbfGB, float pkbfGB, float sbfGB) {
        System.out.println(  "\nBloom filters          Memory (GB)");
        System.out.println(    "====================================");
        if (dbgGB > 0)
            System.out.println("de Bruijn graph:       " + dbgGB);
        if (cbfGB > 0)
            System.out.println("k-mer counting:        " + cbfGB);
        if (pkbfGB > 0) {
            System.out.println("paired k-mers (reads): " + pkbfGB);
            System.out.println("paired k-mers (frags): " + pkbfGB);
        }
        if (sbfGB > 0)
            System.out.println("screening:             " + sbfGB);
        System.out.println(    "====================================");
        System.out.println(    "Total:                 " + (dbgGB+cbfGB+2*pkbfGB+sbfGB));
    }
    
    private static final String STAMP_STARTED = "STARTED";
    private static final String STAMP_DBG_DONE = "DBG.DONE";
    private static final String STAMP_FRAGMENTS_DONE = "FRAGMENTS.DONE";
    private static final String STAMP_TRANSCRIPTS_DONE = "TRANSCRIPTS.DONE";
    private static final String STAMP_TRANSCRIPTS_NR_DONE = "TRANSCRIPTS_NR.DONE";
    private static final String STAMP_LONG_READS_CORRECTED = "LONGREADS.CORRECTED";
    private static final String STAMP_LONG_READS_CLUSTERED = "LONGREADS.CLUSTERED";
    private static final String STAMP_LONG_READS_ASSEMBLED = "LONGREADS.ASSEMBLED";
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {        
        MyTimer timer = new MyTimer();
                    
        // Based on: http://commons.apache.org/proper/commons-cli/usage.html
        CommandLineParser parser = new DefaultParser();

        Options options = new Options();

                
        Option optLeftReads = Option.builder("l")
                                    .longOpt("left")
                                    .desc("left reads file(s)")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optLeftReads);
        
        Option optRightReads = Option.builder("r")
                                    .longOpt("right")
                                    .desc("right reads file(s)")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optRightReads);
        
        Option optPooledAssembly = Option.builder("pool")
                                    .longOpt("pool")
                                    .desc("list of read files for pooled assembly")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optPooledAssembly);
        
        Option optLongReads = Option.builder("long")
                                    .desc("long reads file(s)\n(Requires `minimap2` and `racon` in PATH. Presets `-k 17 -c 3 -indel 10 -e 3 -p 0.8` unless each option is defined otherwise.)")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optLongReads);

        Option optRefTranscripts = Option.builder("ref")
                                    .desc("reference transcripts file(s) for guiding the assembly process")
                                    .hasArgs()
                                    .argName("FILE")
                                    .build();
        options.addOption(optRefTranscripts);
        
        Option optRevCompLeft = Option.builder("rcl")
                                    .longOpt("revcomp-left")
                                    .desc("reverse-complement left reads [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optRevCompLeft);

        Option optRevCompRight = Option.builder("rcr")
                                    .longOpt("revcomp-right")
                                    .desc("reverse-complement right reads [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optRevCompRight);

        Option optRevCompLong = Option.builder("rc")
                                    .longOpt("revcomp-long")
                                    .desc("reverse-complement long reads [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optRevCompLong);
        
        Option optStranded = Option.builder("ss")
                                    .longOpt("stranded")
                                    .desc("reads are strand specific [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optStranded);
        
        final String optNameDefault = "rnabloom";
        Option optName = Option.builder("n")
                                    .longOpt("name")
                                    .desc("assembly name [" + optNameDefault + "]")
                                    .hasArg(true)
                                    .argName("STR")
                                    .build();
        options.addOption(optName);
        
        final String optPrefixDefault = "";
        Option optPrefix = Option.builder("prefix")
                                    .desc("name prefix in FASTA header for assembled transcripts")
                                    .hasArg(true)
                                    .argName("STR")
                                    .build();
        options.addOption(optPrefix);
        
        Option optUracil = Option.builder("u")
                                    .longOpt("uracil")
                                    .desc("output uracils (U) in place of thymines (T) in assembled transcripts [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optUracil);
        
        final String optThreadsDefault = "2";
        Option optThreads = Option.builder("t")
                                    .longOpt("threads")
                                    .desc("number of threads to run [" + optThreadsDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optThreads);
        
        final String optOutdirDefault = System.getProperty("user.dir") + File.separator + "rnabloom_assembly";
        Option optOutdir = Option.builder("o")
                                    .longOpt("outdir")
                                    .desc("output directory [" + optOutdirDefault + "]")
                                    .hasArg(true)
                                    .argName("PATH")
                                    .build();
        options.addOption(optOutdir);
        
        Option optForce = Option.builder("f")
                                    .longOpt("force")
                                    .desc("force overwrite existing files [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optForce);
        
        final String optKmerSizeDefault = "25"; 
        Option optKmerSize = Option.builder("k")
                                    .longOpt("kmer")
                                    .desc("k-mer size [" + optKmerSizeDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optKmerSize);
        
        final String optStageDefault = "4";
        Option optStage = Option.builder("stage")
                                    .desc("assembly termination stage\n" +
                                            "short reads: [3]\n" +
                                            "1. construct graph\n" +
                                            "2. assemble fragments\n" +
                                            "3. assemble transcripts\n" +
                                            "long reads: [4]\n" + 
                                            "1. construct graph\n" +
                                            "2. correct reads\n" +
                                            "3. cluster reads\n" + 
                                            "4. assemble transcripts")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optStage);
        
        final String optBaseQualDbgDefault = "3";
        Option optBaseQualDbg = Option.builder("q")
                                    .longOpt("qual-dbg")
                                    .desc("minimum base quality in reads for constructing DBG [" + optBaseQualDbgDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBaseQualDbg);

        final String optBaseQualFragDefault = "3";
        Option optBaseQualFrag = Option.builder("Q")
                                    .longOpt("qual-frag")
                                    .desc("minimum base quality in reads for fragment reconstruction [" + optBaseQualFragDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBaseQualFrag);        
        
        final String optMinKmerCovDefault = "1"; 
        Option optMinKmerCov = Option.builder("c")
                                    .longOpt("mincov")
                                    .desc("minimum k-mer coverage [" + optMinKmerCovDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinKmerCov);
        
        final String optAllHashDefault = "2";
        Option optAllHash = Option.builder("hash")
                                    .desc("number of hash functions for all Bloom filters [" + optAllHashDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optAllHash); 
        
        Option optSbfHash = Option.builder("sh")
                                    .longOpt("sbf-hash")
                                    .desc("number of hash functions for screening Bloom filter [" + optAllHashDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optSbfHash); 
        
        Option optDbgbfHash = Option.builder("dh")
                                    .longOpt("dbgbf-hash")
                                    .desc("number of hash functions for de Bruijn graph Bloom filter [" + optAllHashDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optDbgbfHash);

        Option optCbfHash = Option.builder("ch")
                                    .longOpt("cbf-hash")
                                    .desc("number of hash functions for k-mer counting Bloom filter [" + optAllHashDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optCbfHash);
        
        Option optPkbfHash = Option.builder("ph")
                                    .longOpt("pkbf-hash")
                                    .desc("number of hash functions for paired k-mers Bloom filter [" + optAllHashDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optPkbfHash);        

        Option optNumKmers = Option.builder("nk")
                                    .longOpt("num-kmers")
                                    .desc("expected number of unique k-mers in input reads")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optNumKmers);
        
        Option optNtcard = Option.builder("ntcard")
                                    .desc("count unique k-mers in input reads with ntCard [false]\n(Requires `ntcard` in PATH. If this option is used along with `-long`, the value for `-c` is set automatically based on the ntCard histogram, unless `-c` is defined otherwise)")
                                    .hasArg(false)
                                    .build();
        options.addOption(optNtcard);        
        
        Option optAllMem = Option.builder("mem")
                                    .longOpt("memory")
                                    .desc("total amount of memory (GB) for all Bloom filters [auto]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optAllMem);
        
        Option optSbfMem = Option.builder("sm")
                                    .longOpt("sbf-mem")
                                    .desc("amount of memory (GB) for screening Bloom filter [auto]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optSbfMem);
        
        Option optDbgbfMem = Option.builder("dm")
                                    .longOpt("dbgbf-mem")
                                    .desc("amount of memory (GB) for de Bruijn graph Bloom filter [auto]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optDbgbfMem);

        Option optCbfMem = Option.builder("cm")
                                    .longOpt("cbf-mem")
                                    .desc("amount of memory (GB) for k-mer counting Bloom filter [auto]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optCbfMem);
        
        Option optPkbfMem = Option.builder("pm")
                                    .longOpt("pkbf-mem")
                                    .desc("amount of memory (GB) for paired kmers Bloom filter [auto]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optPkbfMem);

        final String optFprDefault = "0.01";
        Option optFpr = Option.builder("fpr")
                                    .longOpt("fpr")
                                    .desc("maximum allowable false-positive rate of Bloom filters [" + optFprDefault + "]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optFpr);
        
        Option optSaveBf = Option.builder("savebf")
                                    .desc("save graph (Bloom filters) from stage 1 to disk [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optSaveBf);  
        
        final String optTipLengthDefault = "5";
        Option optTipLength = Option.builder("tiplength")
                                    .desc("maximum branch length to be considered a tip [" + optTipLengthDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optTipLength);
        
        final String optLookaheadDefault = "3";
        Option optLookahead = Option.builder("lookahead")
                                    .desc("number of k-mers to look ahead during graph traversal [" + optLookaheadDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optLookahead);        
        
        final String optSampleDefault = "1000";
        Option optSample = Option.builder("sample")
                                    .desc("sample size for estimating read/fragment lengths [" + optSampleDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optSample);
        
        final String optErrCorrItrDefault = "1";
        Option optErrCorrItr = Option.builder("e")
                                    .longOpt("errcorritr")
                                    .desc("number of iterations of error-correction in reads [" + optErrCorrItrDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optErrCorrItr);
        
        final String optMaxCovGradDefault = "0.50";
        Option optMaxCovGrad = Option.builder("grad")
                                    .longOpt("maxcovgrad")
                                    .desc("maximum k-mer coverage gradient for error correction [" + optMaxCovGradDefault + "]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optMaxCovGrad);
        
        final String optIndelSizeDefault = "1";
        Option optIndelSize = Option.builder("indel")
                                    .desc("maximum size of indels to be collapsed [" + optIndelSizeDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optIndelSize);  

        final String optPercentIdentityDefault = "0.90";
        Option optPercentIdentity = Option.builder("p")
                                    .longOpt("percent")
                                    .desc("minimum percent identity of sequences to be collapsed [" + optPercentIdentityDefault + "]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optPercentIdentity);
                
        final String optMinLengthDefault = "200";
        Option optMinLength = Option.builder("length")
                                    .desc("minimum transcript length in output assembly [" + optMinLengthDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinLength);  
        
        Option optNoReduce = Option.builder("norr")
                                    .desc("skip redundancy reduction for assembled transcripts [false]\n(will not create 'transcripts.nr.fa')")
                                    .hasArg(false)
                                    .build();
        options.addOption(optNoReduce);
        
        final String optOverlapDefault = "10";
        Option optOverlap = Option.builder("overlap")
                                    .desc("minimum number of overlapping bases between read mates [" + optOverlapDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optOverlap);
        
        final String optBoundDefault = "500";
        Option optBound = Option.builder("bound")
                                    .desc("maximum distance between read mates [" + optBoundDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optBound);
        
        Option optExtend = Option.builder("extend")
                                    .desc("extend fragments outward during fragment reconstruction [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optExtend);

        Option optNoFragmentsConsistency = Option.builder("nofc")
                                    .desc("turn off assembly consistency with fragment paired k-mers [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optNoFragmentsConsistency);
        
        Option optSensitive = Option.builder("sensitive")
                                    .desc("assemble transcripts in sensitive mode [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optSensitive);
        
        Option optKeepArtifact = Option.builder("artifact")
                                    .desc("keep potential sequencing artifacts [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optKeepArtifact);

        Option optKeepChimera = Option.builder("chimera")
                                    .desc("keep potential chimeras [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optKeepChimera);
                
        final String optBranchFreeExtensionDefault = STRATUM_E0;
        final String optBranchFreeExtensionChoicesStr = String.join("|", STRATA);
        Option optBranchFreeExtensionThreshold = Option.builder("stratum")
                                    .desc("fragments lower than the specified stratum are extended only if they are branch-free in the graph [" + optBranchFreeExtensionDefault + "]")
                                    .hasArg(true)
                                    .argName(optBranchFreeExtensionChoicesStr)
                                    .build();
        options.addOption(optBranchFreeExtensionThreshold);        
        
        final String optMinKmerPairsDefault = "10";
        Option optMinKmerPairs = Option.builder("pair")
                                    .desc("minimum number of consecutive k-mer pairs for assembling transcripts [" + optMinKmerPairsDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinKmerPairs);  
                
        final String optPolyATailDefault = "0";
        Option optPolyATail = Option.builder("a")
                                    .longOpt("polya")
                                    .desc("prioritize assembly of transcripts with poly-A tails of the minimum length specified [" + optPolyATailDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optPolyATail);  
        
        final String optMinimapOptionsDefault = "-r 150";
        Option optMinimapOptions = Option.builder("mmopt")
                                    .desc("options for minimap2 [" + optMinimapOptionsDefault + "]\n(`-x` and `-t` are already in use)")
                                    .hasArg(true)
                                    .argName("OPTIONS")
                                    .build();
        options.addOption(optMinimapOptions);

        Option optDebug = Option.builder("debug")
                                    .desc("print debugging information [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optDebug);
        
//        Option optHomopolymerCompressed = Option.builder("hpcm")
//                                    .desc("use homopolymer-compressed minimizers in long-read clustering [false]\n(Requires `-long`)")
//                                    .hasArg(false)
//                                    .build();
//        options.addOption(optHomopolymerCompressed);
        
        Option optHelp = Option.builder("h")
                                    .longOpt("help")
                                    .desc("print this message and exits")
                                    .build();
        options.addOption(optHelp);
        
        Option optVersion = Option.builder("v")
                                    .longOpt("version")
                                    .desc("print version information and exits")
                                    .build();
        options.addOption(optVersion);
        

        try {
            CommandLine line = parser.parse(options, args);
            
            if (line.getOptions().length == 0 || line.hasOption(optHelp.getOpt())) {
                printHelp(options, false);
            }
            
            if (line.hasOption(optVersion.getOpt())) {
                printVersionInfo(true);
            }
            
            System.out.println("RNA-Bloom v" + VERSION + "\n" +
                               "args: " + Arrays.toString(args) + "\n");
            
            String branchFreeExtensionThreshold = line.getOptionValue(optBranchFreeExtensionThreshold.getOpt(), optBranchFreeExtensionDefault);
            if (!isValidStratumName(branchFreeExtensionThreshold)) {
                exitOnError("Unknown stratum name specified, \"" + branchFreeExtensionThreshold + "\"");
            }
            
            final boolean debug = line.hasOption(optDebug.getOpt());
            final int endstage = Integer.parseInt(line.getOptionValue(optStage.getOpt(), optStageDefault));
            final int numThreads = Integer.parseInt(line.getOptionValue(optThreads.getOpt(), optThreadsDefault));
            final boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            final String name = line.getOptionValue(optName.getOpt(), optNameDefault);
            final String outdir = line.getOptionValue(optOutdir.getOpt(), optOutdirDefault);
            
            System.out.println("name:   " + name);
            System.out.println("outdir: " + outdir);
            
            File f = new File(outdir);
            if (!f.exists()) {
                System.out.println("WARNING: Output directory does not exist!");
                f.mkdirs();
                System.out.println("Created output directory at `" + outdir + "`");
            }
            
            final String graphFile = outdir + File.separator + name + ".graph";
            
            File startedStamp = new File(outdir + File.separator + STAMP_STARTED);
            File dbgDoneStamp = new File(outdir + File.separator + STAMP_DBG_DONE);
            File fragsDoneStamp = new File(outdir + File.separator + STAMP_FRAGMENTS_DONE);
            File txptsDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_DONE);
            File txptsNrDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_NR_DONE);
            File longReadsCorrectedStamp = new File(outdir + File.separator + STAMP_LONG_READS_CORRECTED);
            File longReadsClusteredStamp = new File(outdir + File.separator + STAMP_LONG_READS_CLUSTERED);
            File longReadsAssembledStamp = new File(outdir + File.separator + STAMP_LONG_READS_ASSEMBLED);
            
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
                
                if (txptsNrDoneStamp.exists()) {
                    txptsNrDoneStamp.delete();
                }
                
                if (longReadsCorrectedStamp.exists()) {
                    longReadsCorrectedStamp.delete();
                }
                
                if (longReadsClusteredStamp.exists()) {
                    longReadsClusteredStamp.delete();
                }
                
                if (longReadsAssembledStamp.exists()) {
                    longReadsAssembledStamp.delete();
                }
            }
                        
            String[] leftReadPaths = line.getOptionValues(optLeftReads.getOpt());
            String[] rightReadPaths = line.getOptionValues(optRightReads.getOpt());
            String[] longReadPaths = line.getOptionValues(optLongReads.getOpt());
            String[] refTranscriptPaths = line.getOptionValues(optRefTranscripts.getOpt());
            
            final String pooledReadsListFile = line.getOptionValue(optPooledAssembly.getOpt());
            final boolean pooledGraphMode = pooledReadsListFile != null;
                        
            HashMap<String, ArrayList<String>> pooledLeftReadPaths = new HashMap<>();
            HashMap<String, ArrayList<String>> pooledRightReadPaths = new HashMap<>();
            
            float maxBfMem = 0;
            
            if (pooledGraphMode) {
                System.out.println("Pooled assembly mode is ON!");
                
                if (!new File(pooledReadsListFile).isFile()) {
                    exitOnError("Cannot find pooled read paths list `" + pooledReadsListFile + "`");
                }
                
                System.out.println("Parsing pool reads list file `" + pooledReadsListFile + "`...");
                boolean parseOK = getPooledReadPaths(pooledReadsListFile, pooledLeftReadPaths, pooledRightReadPaths);
                
                if (!parseOK) {
                    exitOnError("Incorrect format of pooled read paths list file!");
                }
                
                int numLeftIds = pooledLeftReadPaths.size();
                int numRightIds = pooledRightReadPaths.size();
                
                if (numLeftIds != numRightIds) {
                    exitOnError("Pooled read paths list file has disagreeing number of sample IDs for left (" + numLeftIds + ") and right (" + numRightIds + ") reads!");
                }
                
                if (numLeftIds == 0) {
                    exitOnError("Pooled read paths list file is empty!");
                }
                
                ArrayList<String> leftPathsQueue = new ArrayList<>();
                ArrayList<String> rightPathsQueue = new ArrayList<>();
                
                for (String id : pooledLeftReadPaths.keySet()) {
                    leftPathsQueue.addAll(pooledLeftReadPaths.get(id));
                    rightPathsQueue.addAll(pooledRightReadPaths.get(id));
                }
                
                leftReadPaths = new String[leftPathsQueue.size()];
                rightReadPaths = new String[rightPathsQueue.size()];
                
                leftPathsQueue.toArray(leftReadPaths);
                rightPathsQueue.toArray(rightReadPaths);
                
                checkInputFileFormat(leftReadPaths);
                checkInputFileFormat(rightReadPaths);
                
                double readFilesTotalBytes = 0;

                for (String fq : leftReadPaths) {
                    readFilesTotalBytes += new File(fq).length();
                }
                for (String fq : rightReadPaths) {
                    readFilesTotalBytes += new File(fq).length();
                }
                
                maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (Math.max(NUM_BYTES_1MB * 100, readFilesTotalBytes) / NUM_BYTES_1GB))));
            }
            else if (longReadPaths != null && longReadPaths.length > 0) {
                checkInputFileFormat(longReadPaths);
                
                double readFilesTotalBytes = 0;

                for (String fq : longReadPaths) {
                    readFilesTotalBytes += new File(fq).length();
                }
                
                maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (Math.max(NUM_BYTES_1MB * 100, readFilesTotalBytes) / NUM_BYTES_1GB))));
            }
            else {
                double readFilesTotalBytes = 0;
                
                if (leftReadPaths == null || leftReadPaths.length == 0) {
                    exitOnError("Please specify left read files!");
                }

//                if (rightReadPaths == null || rightReadPaths.length == 0) {
//                    exitOnError("Please specify right read files!");
//                }

                if (leftReadPaths != null && rightReadPaths != null &&
                        leftReadPaths.length > 0 && rightReadPaths.length > 0 &&
                        leftReadPaths.length != rightReadPaths.length) {
                    exitOnError("Read files are not paired properly!");
                }
                
                checkInputFileFormat(leftReadPaths);
                
                for (String p : leftReadPaths) {
                    readFilesTotalBytes += new File(p).length();
                }
                
                if (rightReadPaths != null && rightReadPaths.length > 0) {
                    checkInputFileFormat(rightReadPaths);
                    
                    for (String p : rightReadPaths) {
                        readFilesTotalBytes += new File(p).length();
                    }
                }
                                
                maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (Math.max(NUM_BYTES_1MB * 100, readFilesTotalBytes) / NUM_BYTES_1GB))));
            }
            
            final boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            final boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            final boolean revCompLong = line.hasOption(optRevCompLong.getOpt());
            final boolean strandSpecific = line.hasOption(optStranded.getOpt());
            final boolean writeUracil = line.hasOption(optUracil.getOpt());
            final boolean outputNrTxpts = !line.hasOption(optNoReduce.getOpt());
            final String minimapOptions = line.getOptionValue(optMinimapOptions.getOpt(), optMinimapOptionsDefault);
//            final boolean useCompressedMinimizers = line.hasOption(optHomopolymerCompressed.getOpt());
            final boolean useCompressedMinimizers = false;

            boolean hasLeftReadFiles = leftReadPaths != null && leftReadPaths.length > 0;
            boolean hasRightReadFiles = rightReadPaths != null && rightReadPaths.length > 0;
            boolean hasLongReadFiles = longReadPaths != null && longReadPaths.length > 0;
            boolean hasRefTranscriptFiles = refTranscriptPaths != null && refTranscriptPaths.length > 0;
            
            if ((hasLongReadFiles || outputNrTxpts) && !hasMinimap2()) {
                exitOnError("`minimap2` not found in PATH!");
            }

            if (hasLongReadFiles && !hasRacon()) {
                exitOnError("`racon` not found in PATH!");
            }
            
            String defaultK = hasLongReadFiles ? "17" : optKmerSizeDefault;
            final int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), defaultK));
            
            String defaultPercentIdentity = hasLongReadFiles ? "0.7" : optPercentIdentityDefault;
            final float percentIdentity = Float.parseFloat(line.getOptionValue(optPercentIdentity.getOpt(), defaultPercentIdentity));
            
            String defaultMaxIndelSize = hasLongReadFiles ? "10" : optIndelSizeDefault;
            final int maxIndelSize = Integer.parseInt(line.getOptionValue(optIndelSize.getOpt(), defaultMaxIndelSize));
            
            String defaultMaxErrCorrItr = hasLongReadFiles ? "3" : optErrCorrItrDefault;
            final int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), defaultMaxErrCorrItr));
            
            String defaultMinKmerCov = hasLongReadFiles ? "3" : optMinKmerCovDefault;
            int minKmerCov = Integer.parseInt(line.getOptionValue(optMinKmerCov.getOpt(), defaultMinKmerCov));
            
            final int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), optBaseQualDbgDefault));
            final int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), optBaseQualFragDefault));
                        
            float sbfGB = Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), Float.toString(maxBfMem * 1f / 8f)));
            float dbgGB = Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), Float.toString(maxBfMem * 1f / 8f)));
            float cbfGB = Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), Float.toString(maxBfMem * 4f / 8f)));
            float pkbfGB = Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), Float.toString(maxBfMem * 1f / 8f)));
                        
            long sbfSize = (long) (NUM_BITS_1GB * sbfGB);
            long dbgbfSize = (long) (NUM_BITS_1GB * dbgGB);
            long cbfSize = (long) (NUM_BYTES_1GB * cbfGB);
            long pkbfSize = (long) (NUM_BITS_1GB * pkbfGB);
            
            final int allNumHash = Integer.parseInt(line.getOptionValue(optAllHash.getOpt(), optAllHashDefault));
            final String allNumHashStr = Integer.toString(allNumHash);
            final int sbfNumHash = Integer.parseInt(line.getOptionValue(optSbfHash.getOpt(), allNumHashStr));
            final int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), allNumHashStr));
            final int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), allNumHashStr));
            final int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), allNumHashStr));
            
            final float maxFPR = Float.parseFloat(line.getOptionValue(optFpr.getOpt(), optFprDefault));
            final boolean saveGraph = line.hasOption(optSaveBf.getOpt());
            boolean storeReadPairedKmers = hasLeftReadFiles || hasRightReadFiles || hasRefTranscriptFiles;
                        
            long expNumKmers = -1L;
            NTCardHistogram hist = null;
            if (line.hasOption(optNtcard.getOpt())) {
                if (!hasNtcard()) {
                    exitOnError("`ntcard` not found in your PATH!");
                }
                
                System.out.println("\nK-mer counting with ntCard...");
                timer.start();
                 
                String ntcard_reads_list_file = outdir + File.separator + name + ".ntcard.readslist.txt";
                BufferedWriter writer = new BufferedWriter(new FileWriter(ntcard_reads_list_file, false));
                
                if (hasLeftReadFiles && leftReadPaths != null) {
                    for (String p : leftReadPaths) {
                        writer.write(p + '\n');
                    }
                }
                
                if (hasRightReadFiles && rightReadPaths != null) {
                    for (String p : rightReadPaths) {
                        writer.write(p + '\n');
                    }
                }
                
                if (hasLongReadFiles && longReadPaths != null) {
                    for (String p : longReadPaths) {
                        writer.write(p + '\n');
                    }
                }
                
                if (hasRefTranscriptFiles && refTranscriptPaths != null) {
                    for (String p : refTranscriptPaths) {
                        writer.write(p + '\n');
                    }
                }
                
                writer.close();
                
                String histogramPathPrefix = outdir + File.separator + name;
                
                hist = getNTCardHistogram(numThreads, k, histogramPathPrefix, ntcard_reads_list_file, forceOverwrite);
                if (hist == null) {
                    exitOnError("Error running ntCard!");
                }
                
                expNumKmers = hist.numKmers;
                System.out.println("Number of unique k-mers: " + NumberFormat.getInstance().format(expNumKmers));

                if (hasLongReadFiles) {
                    if (!line.hasOption(optMinKmerCov.getOpt())) {
                        minKmerCov = hist.covThreshold;
                    }
                    
                    System.out.println("Min k-mer coverage threshold: " + NumberFormat.getInstance().format(minKmerCov));
                }
                    
                if (expNumKmers <= 0) {
                    exitOnError("Cannot get number of unique k-mers from ntCard! (" + expNumKmers + ")");
                }
                
                System.out.println("K-mer counting completed in " + MyTimer.hmsFormat(timer.elapsedMillis()));                
            }
            else {
                expNumKmers = Long.parseLong(line.getOptionValue(optNumKmers.getOpt(), "-1"));
                System.out.println("Min k-mer coverage threshold: " + NumberFormat.getInstance().format(minKmerCov));
            }
            
            if (expNumKmers > 0) {
                sbfSize = BloomFilter.getExpectedSize(expNumKmers, maxFPR, sbfNumHash);
                sbfGB = sbfSize / (float) NUM_BITS_1GB;
                
                dbgbfSize = BloomFilter.getExpectedSize(expNumKmers, maxFPR, dbgbfNumHash);
                dbgGB = dbgbfSize / (float) NUM_BITS_1GB;
                
                if (hist != null) {
                    cbfSize = CountingBloomFilter.getExpectedSize(expNumKmers-hist.numSingletons, maxFPR, cbfNumHash);
                    cbfGB = cbfSize / (float) NUM_BYTES_1GB;
                }
                else {
                    cbfSize = CountingBloomFilter.getExpectedSize(expNumKmers, maxFPR, cbfNumHash);
                    cbfGB = cbfSize / (float) NUM_BYTES_1GB;
                }
                
                pkbfSize = PairedKeysBloomFilter.getExpectedSize(expNumKmers, maxFPR, pkbfNumHash);
                pkbfGB = pkbfSize / (float) NUM_BITS_1GB;
            }
            
            /**@TODO ensure that sbfNumHash and pkbfNumHash <= max(dbgbfNumHash, cbfNumHash) */
                        
            final int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), optOverlapDefault));
            final int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), optSampleDefault));
            final int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), optBoundDefault));
            final int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), optLookaheadDefault));
            final int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), optTipLengthDefault));
            final float maxCovGradient = Float.parseFloat(line.getOptionValue(optMaxCovGrad.getOpt(), optMaxCovGradDefault));
            final int minTranscriptLength = Integer.parseInt(line.getOptionValue(optMinLength.getOpt(), optMinLengthDefault));
            
            final int minPolyATail = Integer.parseInt(line.getOptionValue(optPolyATail.getOpt(), optPolyATailDefault));
//            if (minPolyATail > 0) {
//                maxErrCorrItr = 0;
//                branchFreeExtensionThreshold = STRATUM_01;
//            }
            
            boolean keepArtifact = line.hasOption(optKeepArtifact.getOpt());
            boolean keepChimera = line.hasOption(optKeepChimera.getOpt());
            final boolean sensitiveMode = line.hasOption(optSensitive.getOpt());
            if (sensitiveMode) {
                branchFreeExtensionThreshold = STRATUM_01;
                keepArtifact = true;
                keepChimera = true;
            }
            
            final boolean noFragDBG = false; //line.hasOption(optNoFragDBG.getOpt());
            final boolean reqFragKmersConsistency = !line.hasOption(optNoFragmentsConsistency.getOpt());
            final boolean extendFragments = line.hasOption(optExtend.getOpt());
            final int minNumKmerPairs = Integer.parseInt(line.getOptionValue(optMinKmerPairs.getOpt(), optMinKmerPairsDefault));
            final String txptNamePrefix = line.getOptionValue(optPrefix.getOpt(), optPrefixDefault);

            if (!storeReadPairedKmers) {
                pkbfGB = 0;
            }
            
            if (hasLongReadFiles) {
                sbfGB = 0;
            }
            
            printBloomFilterMemoryInfo(dbgGB, cbfGB, pkbfGB, sbfGB);
                        
            RNABloom assembler = new RNABloom(k, qDBG, qFrag, debug);
            assembler.setParams(strandSpecific, maxTipLen, lookahead, maxCovGradient, maxIndelSize, percentIdentity, minNumKmerPairs, minPolyATail);

            FileWriter writer = new FileWriter(startedStamp, false);
            writer.write(String.join(" ", args));
            writer.close();
                        
            if (!forceOverwrite && dbgDoneStamp.exists()) {
                System.out.println("WARNING: Graph was already constructed (k=" + k + ")!");
                
                if (endstage == 1) {
                    System.exit(0);
                }
                
                boolean fragmentsDone = fragsDoneStamp.exists();
                
                if (!fragmentsDone || (outputNrTxpts && !txptsNrDoneStamp.exists()) || !txptsDoneStamp.exists()) {
                    System.out.println("Loading graph from file `" + graphFile + "`...");
                    assembler.restoreGraph(new File(graphFile), noFragDBG || !fragmentsDone || (outputNrTxpts && !txptsNrDoneStamp.exists()));
                }
            }
            else {                
                ArrayList<String> forwardFilesList = new ArrayList<>();
                ArrayList<String> backwardFilesList = new ArrayList<>();
                ArrayList<String> longFilesList = new ArrayList<>();
                ArrayList<String> refFilesList = new ArrayList<>();
                
                if (hasLeftReadFiles) {
                    if (revCompLeft) {
                        backwardFilesList.addAll(Arrays.asList(leftReadPaths));
                    }
                    else {
                        forwardFilesList.addAll(Arrays.asList(leftReadPaths));
                    }
                }
                
                if (hasRightReadFiles) {
                    if (revCompRight) {
                        backwardFilesList.addAll(Arrays.asList(rightReadPaths));
                    }
                    else {
                        forwardFilesList.addAll(Arrays.asList(rightReadPaths));
                    }
                }
                
                if (hasLongReadFiles) {
                    longFilesList.addAll(Arrays.asList(longReadPaths));
                }
                
                if (hasRefTranscriptFiles) {
                    refFilesList.addAll(Arrays.asList(refTranscriptPaths));
                }
                       
                System.out.println("\n> Stage 1: Construct graph from reads (k=" + k + ")");
                timer.start();
                                
                assembler.initializeGraph(strandSpecific, 
                        dbgbfSize, cbfSize, pkbfSize, 
                        dbgbfNumHash, cbfNumHash, pkbfNumHash, false, storeReadPairedKmers);
                
                if (!hasLongReadFiles) {
                    assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                }
                
                assembler.populateGraph(forwardFilesList, backwardFilesList, longFilesList, refFilesList, strandSpecific, revCompLong, numThreads, false, storeReadPairedKmers);
                
                if (!assembler.withinMaxFPR(maxFPR)) {
                    System.out.println("WARNING: Bloom filter FPR is higher than the maximum allowed FPR (" + maxFPR*100 +"%)!");
                    
                    System.out.println("Adjusting Bloom filter sizes...");

                    long[] suggestedSizes = assembler.getOptimalBloomFilterSizes(maxFPR);

                    assembler.destroyAllBf();
                    
                    dbgbfSize = suggestedSizes[1];
                    cbfSize = suggestedSizes[2];
                    pkbfSize = suggestedSizes[3];
                    sbfSize = suggestedSizes[0];

                    dbgGB = dbgbfSize / (float) NUM_BITS_1GB;
                    cbfGB = cbfSize / (float) NUM_BYTES_1GB;
                    pkbfGB = pkbfSize / (float) NUM_BITS_1GB;
                    sbfGB = sbfSize / (float) NUM_BITS_1GB;
                    
                    if (!storeReadPairedKmers) {
                        pkbfGB = 0;
                    }

                    if (hasLongReadFiles) {
                        sbfGB = 0;
                    }

                    printBloomFilterMemoryInfo(dbgGB, cbfGB, pkbfGB, sbfGB);
                    
                    assembler.initializeGraph(strandSpecific, 
                            dbgbfSize, cbfSize, pkbfSize, 
                            dbgbfNumHash, cbfNumHash, pkbfNumHash, false, storeReadPairedKmers);
                    
                    if (!hasLongReadFiles) {
                        assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                    }
                    
                    System.out.println("Repopulate graph ...");
                    
                    assembler.populateGraph(forwardFilesList, backwardFilesList, longFilesList, refFilesList, strandSpecific, revCompLong, numThreads, false, storeReadPairedKmers);
                }    
                
                if (saveGraph) {
                    System.out.println("Saving graph to disk `" + graphFile + "`...");
                    assembler.saveGraph(new File(graphFile));
                    touch(dbgDoneStamp);
                }
                
                System.out.println("> Stage 1 completed in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                if (endstage <= 1) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
            }
                        

            if (pooledGraphMode) {
                // assemble fragments for each sample
                int numSamples = pooledLeftReadPaths.size();
                int sampleId = 0;
                
                System.out.println("\n> Stage 2: Assemble fragments for " + numSamples + " samples");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                String[] sampleNames = new String[numSamples];
                pooledLeftReadPaths.keySet().toArray(sampleNames);
                Arrays.sort(sampleNames);
                
                for (String sampleName : sampleNames) {
                    System.out.println(">> Working on \"" + sampleName + "\" (sample " + ++sampleId + " of " + numSamples + ")...");
                    
                    ArrayList<String> lefts = pooledLeftReadPaths.get(sampleName);
                    ArrayList<String> rights = pooledRightReadPaths.get(sampleName);
                    
                    FastxFilePair[] fqPairs = new FastxFilePair[lefts.size()];
                    for (int i=0; i<lefts.size(); ++i) {
                        fqPairs[i] = new FastxFilePair(lefts.get(i), rights.get(i), revCompLeft, revCompRight);
                    }

                    String sampleOutdir = outdir + File.separator + sampleName;
                    new File(sampleOutdir).mkdirs();
                    
                    
                    MyTimer sampleTimer = new MyTimer();
                    sampleTimer.start();
                    
                    assembleFragments(assembler, forceOverwrite,
                                    sampleOutdir, sampleName, fqPairs,
                                    sbfSize, pkbfSize, sbfNumHash, pkbfNumHash, numThreads,
                                    bound, minOverlap, sampleSize, maxErrCorrItr, extendFragments, minKmerCov, keepArtifact);
                    
                    System.out.println(">> Fragments assembled in " + MyTimer.hmsFormat(sampleTimer.elapsedMillis()) + "\n");
                }
                
                System.out.println("> Stage 2 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                touch(fragsDoneStamp);
                
                if (endstage <= 2) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
                
                // assemble transcripts for each sample
                sampleId = 0;
                System.out.println("\n> Stage 3: Assemble transcripts for " + numSamples + " samples");
                stageTimer.start();
                                
                for (String sampleName : sampleNames) {
                    System.out.println(">> Working on \"" + sampleName + "\" (sample " + ++sampleId + " of " + numSamples + ")...");
                    
                    String sampleOutdir = outdir + File.separator + sampleName;
                    
                    assembleTranscripts(assembler, forceOverwrite,
                                    sampleOutdir, sampleName, txptNamePrefix, strandSpecific,
                                    sbfSize, sbfNumHash, pkbfSize, pkbfNumHash, numThreads, noFragDBG,
                                    sampleSize, minTranscriptLength, keepArtifact, keepChimera,
                                    reqFragKmersConsistency, true, minKmerCov,
                                    branchFreeExtensionThreshold, outputNrTxpts, minPolyATail > 0, writeUracil,
                                    refTranscriptPaths);
                    
                    System.out.print("\n");
                }
                
                /*
                if (outputNrTxpts) {
                    // combine assembly files
                    
                    String txptFileExt = outputNrTxpts ? ".nr" + FASTA_EXT : FASTA_EXT;
                    String shortTxptFileExt = outputNrTxpts ? ".nr.short" + FASTA_EXT : FASTA_EXT;

                    System.out.println(">> Merging transcripts from all samples...");
                    MyTimer mergeTimer = new MyTimer();
                    mergeTimer.start();
                    
                    boolean ok = mergePooledAssemblies(outdir, name, sampleNames, 
                                    txptFileExt, shortTxptFileExt, txptNamePrefix,
                                    k, numThreads, strandSpecific, maxIndelSize,
                                    maxTipLen, percentIdentity, !keepArtifact, minTranscriptLength, writeUracil);
                                        
                    if (ok) {
                        System.out.println(">> Merging completed in " + MyTimer.hmsFormat(mergeTimer.elapsedMillis()));
                    }
                    else {
                        exitOnError("Error during assembly merging!");
                    }
                }
                */
                
                System.out.println("> Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));                
                
                touch(txptsDoneStamp);
            }
            else if (hasLongReadFiles) {
                int sketchSize = minTranscriptLength + k + 1;
                
                final String correctedLongReadFilePrefix = outdir + File.separator + name + ".longreads.corrected";

                final int numCovStrata = COVERAGE_ORDER.length;
                final int numLenStrata = LENGTH_STRATUM_NAMES.length;
                String[][] correctedLongReadFileNames = new String[numCovStrata][numLenStrata];
                for (int c=0; c<numCovStrata; ++c) {
                    String covStratumName = COVERAGE_ORDER[c];
                    for (int l=0; l<numLenStrata; ++l) {
                        String lengthStratumName = LENGTH_STRATUM_NAMES[l];

                        String correctedLongReadsFasta = correctedLongReadFilePrefix + "." + covStratumName + "." + lengthStratumName + FASTA_EXT;
                        correctedLongReadFileNames[c][l] = correctedLongReadsFasta;
                    }
                }

                System.out.println("\n> Stage 2: Correct long reads for \"" + name + "\"");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                if (forceOverwrite || !longReadsCorrectedStamp.exists()) {
                    /* set up the file writers */
                    for (String[] row : correctedLongReadFileNames) {
                        for (String fasta : row) {
                            Files.deleteIfExists(FileSystems.getDefault().getPath(fasta));
                        }
                    }
                    
                    correctLongReads(assembler,
                            longReadPaths, correctedLongReadFileNames,
                            maxErrCorrItr, minKmerCov, numThreads, sampleSize, minTranscriptLength,
                            revCompLong, !keepArtifact);
                    
                    touch(longReadsCorrectedStamp);
                }
                
                System.out.println("> Stage 2 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                if (endstage <= 2) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
                
                final String clusteredLongReadsDirectory = outdir + File.separator + name + ".longreads.clusters";
                
                System.out.println("\n> Stage 3: Cluster long reads for \"" + name + "\"");
                stageTimer.start();
                
                if (forceOverwrite || !longReadsClusteredStamp.exists()) {                    
                    clusterLongReads(assembler,
                            correctedLongReadFileNames, clusteredLongReadsDirectory,
                            sketchSize, numThreads, minKmerCov, useCompressedMinimizers);
                    
                    touch(longReadsClusteredStamp);
                }
                
                System.out.println("Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                if (endstage <= 3) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
                
                System.out.println("\n> Stage 4: Assemble long reads for \"" + name + "\"");
                stageTimer.start();
                
                final String assembledLongReadsDirectory = outdir + File.separator + name + ".longreads.assembly";
                final String assembledLongReadsCombinedFile = outdir + File.separator + name + ".transcripts" + FASTA_EXT;
                if (forceOverwrite || !longReadsAssembledStamp.exists()) {
                    assembler.destroyAllBf();
                                        
                    boolean ok = assembleLongReads(assembler,
                            clusteredLongReadsDirectory, assembledLongReadsDirectory, assembledLongReadsCombinedFile,
                            numThreads, forceOverwrite, writeUracil, minimapOptions, minKmerCov, txptNamePrefix, strandSpecific, minTranscriptLength, !keepArtifact);
                    
                    if (ok) {
                        touch(longReadsAssembledStamp);
                    }
                    else {
                        exitOnError("Error assembling long reads!");
                    }
                }
                
//                if (outputNrTxpts) {
//                    generateNonRedundantTranscripts(assembler, forceOverwrite, outdir, name, sbfSize, sbfNumHash);
//                }
                
                System.out.println("> Stage 4 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
            }
            else if (hasLeftReadFiles && !hasRightReadFiles) {
                // Note: no stage 2
                
                System.out.println("\n> Stage 3: Assemble transcripts for \"" + name + "\"");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                
                final String transcriptsFasta = outdir + File.separator + name + ".transcripts" + FASTA_EXT;
                final String shortTranscriptsFasta = outdir + File.separator + name + ".transcripts.short" + FASTA_EXT;
                
                assembler.assembleSingleEndReads(leftReadPaths,
                                    revCompLeft,
                                    transcriptsFasta,
                                    shortTranscriptsFasta,
                                    numThreads,
                                    minTranscriptLength,
                                    keepArtifact,
                                    keepChimera,
                                    txptNamePrefix,
                                    minKmerCov,
                                    writeUracil);
                
                System.out.println("> Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
            }
            else {
                FastxFilePair[] fqPairs = new FastxFilePair[leftReadPaths.length];
                for (int i=0; i<leftReadPaths.length; ++i) {
                    fqPairs[i] = new FastxFilePair(leftReadPaths[i], rightReadPaths[i], revCompLeft, revCompRight);
                }

                System.out.println("\n> Stage 2: Assemble fragments for \"" + name + "\"");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                assembleFragments(assembler, forceOverwrite,
                                    outdir, name, fqPairs,
                                    sbfSize, pkbfSize, sbfNumHash, pkbfNumHash, numThreads,
                                    bound, minOverlap, sampleSize, maxErrCorrItr, extendFragments, minKmerCov, keepArtifact);
                
                System.out.println("> Stage 2 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                if (endstage <= 2) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }

                System.out.println("\n> Stage 3: Assemble transcripts for \"" + name + "\"");
                stageTimer.start();
                                
                assembleTranscripts(assembler, forceOverwrite,
                                outdir, name, txptNamePrefix, strandSpecific,
                                sbfSize, sbfNumHash, pkbfSize, pkbfNumHash, numThreads, noFragDBG,
                                sampleSize, minTranscriptLength, keepArtifact, keepChimera,
                                reqFragKmersConsistency, true, minKmerCov, 
                                branchFreeExtensionThreshold, outputNrTxpts, minPolyATail > 0, writeUracil,
                                refTranscriptPaths);
                
                System.out.println("> Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
            }      
        }
        catch (Exception exp) {
            handleException(exp);
        }
        
        System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
    }
}
