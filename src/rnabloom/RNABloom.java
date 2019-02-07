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
import java.io.File;
import java.io.FileNotFoundException;
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;
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
import rnabloom.io.FastaPairReader;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.FastqPairReader;
import rnabloom.io.FastxFilePair;
import rnabloom.io.PairedReadSegments;
import rnabloom.io.FastqReader;
import rnabloom.io.FastqRecord;
import rnabloom.io.FastxPairReader;
import rnabloom.io.FileFormatException;
import rnabloom.util.GraphUtils;
import static rnabloom.util.GraphUtils.*;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author Ka Ming Nip
 */
public class RNABloom {
    public final static String VERSION = "1.0.0";
    
//    private final static long NUM_PARSED_INTERVAL = 100000;
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    public final static long NUM_BYTES_1MB = (long) pow(1024, 2);
    public final static long NUM_BYTES_1KB = (long) 1024;
    
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
    
    private void exitOnError(String msg) {
        System.out.println("ERROR: " + msg);
        System.exit(1);
    }
    
    private void handleException(Exception ex) {
        System.out.println("ERROR: " + ex.getMessage() );
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
    
    public void restoreGraph(File f, boolean loadDbgBits) {
        try {
            if (graph != null) {
                graph.destroy();
            }
            
            graph = new BloomFilterDeBruijnGraph(f, loadDbgBits);

            if (loadDbgBits) {
                dbgFPR = graph.getDbgbfFPR();
                System.out.println("DBG Bloom filter FPR:               " + dbgFPR * 100 + " %");
            }
            
            covFPR = graph.getCbfFPR();
            System.out.println("Counting Bloom filter FPR:           " + covFPR * 100 + " %");
            System.out.println("Read paired k-mers Bloom filter FPR: " + graph.getRpkbf().getFPR() * 100 + " %");
        } catch (Exception ex) {
            handleException(ex);
        }
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
    
    public class FastqToGraphWorker implements Runnable {
        
        private final int id;
        private long numReads = 0;
        private final FastqReader reader;
        private final NTHashIterator itr;
        private boolean successful = false;
        private boolean incrementIfPresent = false;
        
        public FastqToGraphWorker(int id, FastqReader reader, boolean stranded, boolean reverseComplement, int numHash, boolean incrementIfPresent) {            
            this.id = id;
            this.reader = reader;
            this.incrementIfPresent = incrementIfPresent;
            
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
            
            try {
                Matcher mQual = qualPatternDBG.matcher("");
                Matcher mSeq = seqPattern.matcher("");
                
                if (incrementIfPresent) {
                    try {
                        long[] hashVals = itr.hVals;
                        FastqRecord record = new FastqRecord();
                        
                        while (true) {
                            reader.nextWithoutName(record);
                            mQual.reset(record.qual);

                            while (mQual.find()) {
                                mSeq.reset(record.seq.substring(mQual.start(), mQual.end()));
                                while (mSeq.find()) {
                                    itr.start(mSeq.group());
                                    while (itr.hasNext()) {
                                        itr.next();
                                        graph.addCountIfPresent(hashVals); // Only insert if kmer is absent in the graph
                                    }
                                }
                            }
                            
                            ++numReads;
                        }
                    }
                    catch (NoSuchElementException e) {
                        //end of file
                        successful = true;
                    }
                }
                else {
                    try {
                        long[] hashVals = itr.hVals;
                        FastqRecord record = new FastqRecord();
                        
                        while (true) {
                            reader.nextWithoutName(record);
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
                            
                            ++numReads;
                        }
                    }
                    catch (NoSuchElementException e) {
                        //end of file
                        successful = true;
                    }
                }
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
        /* Support storing paired kmers */
        
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
                long[] phashVals = pitr.hVals3;
                
                if (FastqReader.isFastq(path)) {
                    FastqReader fr = new FastqReader(path);
                    Matcher mQual = qualPatternDBG.matcher("");
                    
                    try {
                        FastqRecord record = new FastqRecord();
                        
                        if (storeReadPairedKmers) {
                            while (true) {
                                ++numReads;

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
//                                            if (screeningBf.lookupThenAdd(hashVals)) {
                                                addFunction.accept(hashVals);
//                                            }
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
                            }
                        }
                        else {
                            while (true) {
                                ++numReads;

                                fr.nextWithoutName(record);
                                mQual.reset(record.qual);
                                mSeq.reset(record.seq);

                                while (mQual.find()) {
                                    mSeq.region(mQual.start(), mQual.end());
                                    while (mSeq.find()) {
                                        itr.start(record.seq, mSeq.start(), mSeq.end());
                                        while (itr.hasNext()) {
                                            itr.next();
//                                            if (screeningBf.lookupThenAdd(hashVals)) {
                                                addFunction.accept(hashVals);
//                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    catch (NoSuchElementException e) {
                        //end of file
                    }
                    
                    fr.close();
                }
                else if (FastaReader.isFasta(path)) {
                    FastaReader fr = new FastaReader(path);
                    
                    try {
                        String seq;
                        
                        if (storeReadPairedKmers) {
                            while (true) {
                                ++numReads;

                                seq = fr.next();
                                mSeq.reset(seq);
                                                                
                                while (mSeq.find()) {
                                    int start = mSeq.start();
                                    int end = mSeq.end();
                                    
                                    itr.start(seq, start, end);
                                    while (itr.hasNext()) {
                                        itr.next();
//                                        if (screeningBf.lookupThenAdd(hashVals)) {
                                            addFunction.accept(hashVals);
//                                        }
                                    }
                                    
                                    if (end - start - k + 1 >= kmerPairDistance) {
                                        pitr.start(seq, start, end);
                                        while (pitr.hasNext()) {
                                            pitr.next();
                                            graph.addSingleReadPairedKmer(phashVals);
                                        }
                                    }
                                }
                            }
                        }
                        else {
                            while (true) {
                                ++numReads;

                                seq = fr.next();
                                mSeq.reset(seq);

                                while (mSeq.find()) {
                                    itr.start(seq, mSeq.start(), mSeq.end());
                                    while (itr.hasNext()) {
                                        itr.next();
//                                        if (screeningBf.lookupThenAdd(hashVals)) {
                                            addFunction.accept(hashVals);
//                                        }
                                    }
                                }
                            }
                        }
                    }
                    catch (NoSuchElementException e) {
                        //end of file
                    }
                    
                    fr.close();
                }
                else {
                    System.out.println("Incompatible file format in `" + path + "`");
                    throw new RuntimeException();
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
        
    public int getMaxReadLength(String path) {
        int max = -1;
        
        try {
            if (FastqReader.isFastq(path)) {
                FastqRecord r = new FastqRecord();
                FastqReader fr = new FastqReader(path);
                for (int i=0; i< 100 && fr.hasNext(); ++i) {
                    fr.nextWithoutName(r);
                    max = Math.max(max, r.seq.length());
                }
                fr.close();
            }
            else if (FastaReader.isFasta(path)) {
                FastaReader fr = new FastaReader(path);
                for (int i=0; i< 100 && fr.hasNext(); ++i) {
                    max = Math.max(max, fr.next().length());
                }
                fr.close();                
            }
            else {
                System.out.println("Incompatible file format in `" + path + "`");
            }
        }
        catch (IOException ex) {
            System.out.println("ERROR: " + ex.getMessage());
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
    
    public void populateGraph2(Collection<String> forwardFastqs,
                            Collection<String> reverseFastqs,
                            boolean strandSpecific,
                            int numThreads,
                            boolean addCountsOnly) {
        
        long numReads = 0;
        int numHash = graph.getMaxNumHash();
        
        try {
            ArrayDeque<FastqToGraphWorker> threadPool = new ArrayDeque<>();
           
            for (String path : forwardFastqs) {
                ExecutorService service = Executors.newFixedThreadPool(numThreads);
                
                System.out.println("Parsing `" + path + "`...");
                FastqReader reader = new FastqReader(path);

                for (int threadId=1; threadId<=numThreads; ++threadId) {
                    FastqToGraphWorker t = new FastqToGraphWorker(++threadId, reader, strandSpecific, false, numHash, addCountsOnly);
                    service.submit(t);
                    threadPool.add(t);
                }

                service.shutdown();
                service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

                while (!threadPool.isEmpty()) {
                    numReads += threadPool.pop().getReadCount();
                }
            }

            for (String path : reverseFastqs) {
                ExecutorService service = Executors.newFixedThreadPool(numThreads);
                
                System.out.println("Parsing `" + path + "`...");
                FastqReader reader = new FastqReader(path);
                
                for (int threadId=1; threadId<=numThreads; ++threadId) {
                    FastqToGraphWorker t = new FastqToGraphWorker(++threadId, reader, strandSpecific, true, numHash, addCountsOnly);
                    service.submit(t);
                    threadPool.add(t);
                }
                
                service.shutdown();
                service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

                while (!threadPool.isEmpty()) {
                    numReads += threadPool.pop().getReadCount();
                }
            }

            System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");
            
            dbgFPR = graph.getDbgbfFPR();
            covFPR = graph.getCbfFPR();
            
        } catch (Exception ex) {
            handleException(ex);
        }
    }
        
    public void setReadKmerDistance(Collection<String> forwardReadPaths,
                                    Collection<String> reverseReadPaths) {
        
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
    
    public void populateGraph(Collection<String> forwardReadPaths,
                            Collection<String> reverseReadPaths,
                            boolean strandSpecific,
                            int numThreads,
                            boolean addCountsOnly,
                            boolean storeReadKmerPairs) {        
        
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

        service.shutdown();
        
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            
            for (SeqToGraphWorker t : threadPool) {
                numReads += t.getReadCount();
            }

            System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");            
        } catch (Exception ex) {
            handleException(ex);
        }
        
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
        return graph.getDbgbfFPR() <= fpr || graph.getCbfFPR() <= fpr || graph.getRpkbfFPR() <= fpr;
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
    
//    public void addCountsOnly(Collection<String> forwardFastqs,
//                            Collection<String> reverseFastqs,
//                            boolean strandSpecific,
//                            int numThreads) {        
//                
//        /** parse the reads */
//        
//        int numReads = 0;
//        int numHash = graph.getMaxNumHash();
//        
//        ExecutorService service = Executors.newFixedThreadPool(numThreads);
//        
//        ArrayList<SeqToGraphWorker> threadPool = new ArrayList<>();
//        int threadId = 0;
//           
//        for (String fastq : forwardFastqs) {
//            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, fastq, strandSpecific, false, numHash, true);
//            service.submit(t);
//            threadPool.add(t);
//        }
//
//        for (String fastq : reverseFastqs) {
//            SeqToGraphWorker t = new SeqToGraphWorker(++threadId, fastq, strandSpecific, true, numHash, true);
//            service.submit(t);
//            threadPool.add(t);
//        }
//
//        service.shutdown();
//        
//        try {
//            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            
//            for (SeqToGraphWorker t : threadPool) {
//                numReads += t.getReadCount();
//            }
//
//            System.out.println("Parsed " + NumberFormat.getInstance().format(numReads) + " reads in total.");
////            System.out.println("Screening Bloom filter FPR:  " + screeningBf.getFPR() * 100 + " %");
//            
////            screeningBf.destroy();
//            
//        } catch (Exception ex) {
//            handleException(ex);
//        }
//        
//        dbgFPR = graph.getDbgbfFPR();
//        covFPR = graph.getCbfFPR();
//    }
    
    public void populateGraphFromFragments(Collection<String> fastas, boolean strandSpecific, boolean newKmerSize) {
        /** insert into graph if absent */
        
        /** parse the fragments */
                          
        NTHashIterator itr = graph.getHashIterator(graph.getMaxNumHash());
        long[] hashVals = itr.hVals;

        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;

        try {
            for (String path : fastas) {
                FastaReader fin = new FastaReader(path);
                
                System.out.println("Parsing `" + path + "`...");
                
                String seq;
                try {
                    if (newKmerSize) {
                        while (true) {
                            seq = fin.next();

                            if (itr.start(seq)) {
                                while (itr.hasNext()) {
                                    itr.next();
                                    graph.addIfAbsent(hashVals);
//                                    screeningBf.add(hashVals);
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
                        // reuse existing k-mer counts and paired kmers
                        
                        while (true) {
                            seq = fin.next();

                            if (itr.start(seq)) {
                                while (itr.hasNext()) {
                                    itr.next();
                                    graph.addDbgOnly(hashVals);
//                                    screeningBf.add(hashVals);
                                }
                            }
                        }
                    }
                }
                catch (NoSuchElementException e) {
                    //end of file
                }

                fin.close();
            }
        } catch (Exception ex) {
            handleException(ex);
        }
                
        dbgFPR = graph.getDbgbfFPR();
        System.out.println("DBG Bloom filter FPR:      " + dbgFPR * 100 + " %");
//        covFPR = graph.getCbfFPR();
    }
    
    public void insertIntoDeBruijnGraph(String fasta) throws IOException, Exception { 
        
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;
        
        FastaReader fin = new FastaReader(fasta);

        try {
            while (true) {
                if (itr.start(fin.next())) {
                    while (itr.hasNext()) {
                        itr.next();
                        graph.addDbgOnly(hashVals);
                    }
                }
            }
        }
        catch (NoSuchElementException e) {
            //end of file
        }

        fin.close();
    }

    public void insertIntoDeBruijnGraphAndPairedKmers(String fasta) throws IOException, Exception {
                
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;

        PairedNTHashIterator pItr = graph.getPairedHashIterator();
        long[] hashVals1 = pItr.hVals1;
        long[] hashVals2 = pItr.hVals2;
        long[] hashVals3 = pItr.hVals3;
        
        FastaReader fin = new FastaReader(fasta);
        
        String seq;
        try {
            while (true) {
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
        catch(NoSuchElementException e) {
            // end of file
        }

        fin.close();
    }
    
    public void insertIntoDeBruijnGraph(String[] fastas) throws IOException, Exception { 
        
        NTHashIterator itr = graph.getHashIterator(graph.getDbgbfNumHash());
        long[] hashVals = itr.hVals;
        
        for (String fasta : fastas) {
            FastaReader fin = new FastaReader(fasta);
            
            try {
                while (true) {
                    if (itr.start(fin.next())) {
                        while (itr.hasNext()) {
                            itr.next();
                            graph.addDbgOnly(hashVals);
                        }
                    }
                }
            }
            catch(NoSuchElementException e) {
                // end of file
            }
            
            fin.close();
        }
    }
    
//    private void insertIntoDeBruijnGraphMultiThreaded(ArrayDeque<String> fragmentsFasta, int sampleSize, int numThreads) throws InterruptedException, IOException, Exception {
//        ArrayBlockingQueue<String> fragmentsQueue = new ArrayBlockingQueue<>(sampleSize, true);
//
//        FragmentDbgWorker[] workers = new FragmentDbgWorker[numThreads];
//        Thread[] threads = new Thread[numThreads];
//        for (int i=0; i<numThreads; ++i) {
//            workers[i] = new FragmentDbgWorker(fragmentsQueue);
//            threads[i] = new Thread(workers[i]);
//            threads[i].start();
//        }
//
//        for (String fa : fragmentsFasta) {
//            FastaReader fin = new FastaReader(fa);
//
//            while (fin.hasNext()) {
//                fragmentsQueue.put(fin.next());
//            }
//
//            fin.close();
//        }
//
//        for (FragmentDbgWorker w : workers) {
//            w.stopWhenEmpty();
//        }
//
//        for (Thread t : threads) {
//            t.join();
//        }
//    }
    
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
                        if (itr.start(fragment)) {
                            while (itr.hasNext()) {
                                itr.next();
                                graph.addDbgOnly(hashVals);
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
        private String prefix = "";
        private long cid = 0;
        
        public TranscriptWriter(FastaWriter fout, 
                                FastaWriter foutShort,
                                int minTranscriptLength,
                                int maxTipLength) {
            this.fout = fout;
            this.foutShort = foutShort;
            this.minTranscriptLength = minTranscriptLength;
            this.maxTipLength = maxTipLength;
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
                headerBuilder.append(" F=[");
                headerBuilder.append(fragment);
                headerBuilder.append("]");
                
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
                                transcriptSB.setCharAt(p+i, Character.toLowerCase(c));
                            }
                        }

                        transcript = transcriptSB.toString();
                        
                        headerBuilder.append("]");
                    }
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
    
    private class TranscriptAssemblyWorker implements Runnable {
        
        private final ArrayBlockingQueue<String> fragments;
        private final ArrayBlockingQueue<Transcript> transcripts;
        private boolean keepGoing = true;
        private boolean extendBranchFreeFragmentsOnly = false;
        private boolean keepArtifact = false;
        private boolean keepChimera = false;
        private boolean reqFragKmersConsistency = false;
        private final float minKmerCov;
        
        public TranscriptAssemblyWorker(ArrayBlockingQueue<String> fragments,
                                        ArrayBlockingQueue<Transcript> transcripts,
                                        boolean includeNaiveExtensions,
                                        boolean extendBranchFreeFragmentsOnly,
                                        boolean keepArtifact,
                                        boolean keepChimera,
                                        boolean reqFragKmersConsistency,
                                        float minKmerCov) {
            this.fragments = fragments;
            this.transcripts = transcripts;
            this.extendBranchFreeFragmentsOnly = extendBranchFreeFragmentsOnly;
            this.keepArtifact = keepArtifact;
            this.keepChimera = keepChimera;
            this.reqFragKmersConsistency = reqFragKmersConsistency;
            this.minKmerCov = minKmerCov;
        }

        public void stopWhenEmpty () {
            keepGoing = false;
        }

        private void storeConsistentReadSegments(String fragment, ArrayList<Kmer> txptKmers) throws InterruptedException {
            ArrayDeque<ArrayList<Kmer>> readSegments = breakWithReadPairedKmers(txptKmers, graph, lookahead);

            int numReadSegs = readSegments.size();

            if (numReadSegs == 1) {
                if (!keepArtifact) {
                    if (!isTemplateSwitch2(txptKmers, graph, screeningBf, lookahead, percentIdentity)) {
                        txptKmers = trimReverseComplementArtifact(txptKmers, k, 0.5f);
                        transcripts.put(new Transcript(fragment, txptKmers));
                        
//                        String seq = cutHairPinLoop(graph.assemble(txptKmers), k, minPercentIdentity);                        
//                        if (seq == null) {
//                            transcripts.put(new Transcript(fragment, txptKmers));
//                        }
//                        else {
//                            transcripts.put(new Transcript(fragment, graph.getKmers(seq)));
//                        }
                    }
                }
                else {
                    transcripts.put(new Transcript(fragment, txptKmers));
                }
            }
            else if (numReadSegs > 1) {
                int numFragKmers = getNumKmers(fragment, k);
                for (ArrayList<Kmer> r : readSegments) {
                    String rAssembled = graph.assemble(r);
                    if (r.size() >= numFragKmers && rAssembled.contains(fragment)) {
                        if (!keepArtifact) {
                            if (!isTemplateSwitch2(txptKmers, graph, screeningBf, lookahead, percentIdentity)) {
                                txptKmers = trimReverseComplementArtifact(txptKmers, k, 0.5f);
                                transcripts.put(new Transcript(fragment, txptKmers));
                                
//                                String seq = cutHairPinLoop(rAssembled, k, minPercentIdentity);
//                                if (seq == null) {
//                                    transcripts.put(new Transcript(fragment, r));
//                                }
//                                else {
//                                    transcripts.put(new Transcript(fragment, graph.getKmers(seq)));
//                                }
                            }
                        }
                        else {
                            transcripts.put(new Transcript(fragment, r));
                        }
                        
                        break;
                    }
                }
            }
        }
        
        @Override
        public void run() {
            try {
                int fragKmersDist = graph.getFragPairedKmerDistance();
                int maxEdgeClipLength = minPolyATailLengthRequired > 0 ? 0 : maxTipLength;
                boolean keepBluntEndArtifact = keepArtifact;
                
                while (true) {
                    String fragment = fragments.poll(10, TimeUnit.MICROSECONDS);
                    
                    if (fragment == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {                        
                        ArrayList<Kmer> kmers = graph.getKmers(fragment);
                        
                        if (!kmers.isEmpty()) {
                            ArrayList<Kmer> originalFragKmers = new ArrayList<>(kmers);

                            if ( (!extendBranchFreeFragmentsOnly || isBranchFree(kmers, graph, maxTipLength)) &&
                                 !represented(kmers,
                                                graph,
                                                screeningBf,
                                                lookahead,
                                                maxIndelSize,
                                                maxEdgeClipLength,
                                                percentIdentity) &&
                                 (keepChimera || !isFusion(kmers, graph, screeningBf, lookahead)) &&
                                 (keepBluntEndArtifact || !isBluntEndArtifact(kmers, graph, screeningBf, maxEdgeClipLength)) ) {
                                
                                extendPE(kmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs, maxCovGradient, minKmerCov);

                                if (kmers.size() > fragKmersDist) {
                                    if (reqFragKmersConsistency) {
                                        ArrayDeque<ArrayList<Kmer>> fragSegments = breakWithFragPairedKmers(kmers, graph);
                                        int numFragSegs = fragSegments.size();

                                        if (numFragSegs >= 1) {
                                            for (ArrayList<Kmer> seg : fragSegments) {
                                                if (numFragSegs == 1 || (seg.size() >= originalFragKmers.size() && new HashSet<>(seg).containsAll(originalFragKmers))) {
                                                    storeConsistentReadSegments(fragment, seg);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    else {
                                        storeConsistentReadSegments(fragment, kmers);
                                    }
                                }
                                else {
                                    storeConsistentReadSegments(fragment, kmers);
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

                    ArrayDeque<ArrayList<Kmer>> segments = breakWithReadPairedKmers(fragmentKmers, graph, lookahead);
                    
                    if (segments.size() != 1) {
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
        private PairedReadSegments p;
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
        
        public FragmentAssembler(PairedReadSegments p,
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
            
            this.p = p;
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
                ArrayList<Kmer> leftKmers = null;
                ArrayList<Kmer> rightKmers = null;
                                
                // connect segments of each read
                String left = getBestSegment(p.left, graph);
                
                if (left.length() >= this.leftReadLengthThreshold) {
                    if (!isLowComplexity2(left)) {
                        if (minKmerCov > 1) {
                            leftKmers = graph.getKmers(left, minKmerCov);
                        }
                        else {
                            leftKmers = graph.getKmers(left);
                        }
                    }
                } 
                                
//                if (leftKmers == null || leftKmers.isEmpty()) {
//                    return;
//                }
                
                String right = getBestSegment(p.right, graph);
                                
                if (right.length() >= this.rightReadLengthThreshold) {
                    if (!isLowComplexity2(right)) {
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
                    return;
                }
                
                /*
                if (leftBad) {
                    boolean hasComplexKmer = false;
                    float minCov = Float.MAX_VALUE;
                    if (rightKmers.size() >= lookahead) {
                        for (Kmer kmer : rightKmers) {
                            if (kmer.count < minCov) {
                                minCov = kmer.count;
                            }

                            if (!hasComplexKmer && !graph.isLowComplexity(kmer)) {
                                hasComplexKmer = true;
                            }
                        }
                    }

                    if (hasComplexKmer) {
                        outList.put(new Fragment("N", graph.assemble(rightKmers), null, 0, minCov, true));
                        return;
                    }
                }
                
                if (rightBad) {
                    boolean hasComplexKmer = false;
                    float minCov = Float.MAX_VALUE;
                    if (leftKmers.size() >= lookahead) {
                        for (Kmer kmer : leftKmers) {
                            if (kmer.count < minCov) {
                                minCov = kmer.count;
                            }

                            if (!hasComplexKmer && !graph.isLowComplexity(kmer)) {
                                hasComplexKmer = true;
                            }
                        }
                    }

                    if (hasComplexKmer) {
                        outList.put(new Fragment(graph.assemble(leftKmers), "N", null, 0, minCov, true));
                        return;
                    }
                }
                */
                
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
                if (fragmentKmers != null) {
                    if (extendFragments) {
                        fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength, minKmerCov);
                    }
                    
                    if (trimArtifact) {
                        fragmentKmers = trimReverseComplementArtifact(fragmentKmers, k, 0.5f);
                    }
                    
                    if (graph.getReadPairedKmerDistance() < fragmentKmers.size()) {
                        ArrayDeque<ArrayList<Kmer>> segments = breakWithReadPairedKmers(fragmentKmers, graph, lookahead);

                        if (segments.size() != 1) {
                            fragmentKmers = null;
                        }
                    }
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

                    if (!leftBad && leftKmers.size() >= lookahead) {
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

                    if (!rightBad && rightKmers.size() >= lookahead) {
                        for (Kmer kmer : rightKmers) {
                            if (kmer.count < minCov) {
                                minCov = kmer.count;
                            }

                            if (!hasComplexRightKmer && !graph.isLowComplexity(kmer)) {
                                hasComplexRightKmer = true;
                            }
                        }
                    }

                    if (hasComplexLeftKmer || hasComplexRightKmer) {
                        left = leftBad ? "" : graph.assemble(leftKmers);
                        right = rightBad ? "" : graph.assemble(rightKmers);
                        
                        outList.put(new Fragment(left, right, null, 0, minCov, true));
                    }
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
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
    
    public void setupFragmentPairedKmersBloomFilter(long numBits, int numHash) {
        graph.initializePairKmersBloomFilter(numBits, numHash);
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
                                                float minKmerCov) {
        
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
                try {
                    while (true) {
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
                                        m = getMinCoverageOrderOfMagnitude(frag.minCov);

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
                                                m = getMinCoverageOrderOfMagnitude(frag.minCov);

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
                }
                catch (NoSuchElementException e) {
                    // end of file
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
                                m = getMinCoverageOrderOfMagnitude(frag.minCov);

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
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
            System.out.println("Rescued " + NumberFormat.getInstance().format(rescuedReadPairs) + " read pairs.");
            System.out.println("Paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
            System.out.println("Screening Bloom filter FPR:    " + screeningBf.getFPR() * 100 + " %");
        }        
    }
    
    private final static String LABEL_SEPARATOR = ":";
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
            
            longFragmentLengthThreshold = fragStats[1];
        
        } catch (Exception ex) {
            handleException(ex);
        }
        
        return fragStats;
    }
    
    public void setPairedKmerDistance(int fragStatQ1) {
        graph.setFragPairedKmerDistance(fragStatQ1 - k - minNumKmerPairs);
    }
    
    public int getPairedReadsMaxDistance(int[] fragStats) {
        return fragStats[3] + ((fragStats[3] - fragStats[1]) * 3 / 2); // 1.5*IQR
    }
    
    private FastaWriter assignUnconnectedFastaWriter(Fragment frag,
                                        String seq,
                                        boolean assemblePolyaTails,
                                        FastaWriter[] unconnectedReadsOut,
                                        FastaWriter unconnectedSingeltonsOut,
                                        FastaWriter[] unconnectedPolyaReadsOut,
                                        FastaWriter unconnectedPolyaSingletonsOut) {
        
        boolean isPolya = assemblePolyaTails && polyATailPattern.matcher(seq).matches();
        
        if (frag.minCov == 1) {
            if (isPolya) {
                return unconnectedPolyaSingletonsOut;
            }
            else {
                return unconnectedSingeltonsOut;
            }
        }
        else {
            int m = getMinCoverageOrderOfMagnitude(frag.minCov);
            
            if (isPolya) {
                return unconnectedPolyaReadsOut[m];
            }
            else {
                return unconnectedReadsOut[m];
            }
        }
    }
    
    private FastaWriter assignFragmentFastaWriter(Fragment frag,
                                        String seq,
                                        boolean assemblePolyaTails,
                                        FastaWriter[] longFragmentsOut,
                                        FastaWriter[] shortFragmentsOut,
                                        FastaWriter longSingletonsOut,
                                        FastaWriter shortSingletonsOut,
                                        FastaWriter[] longPolyaFragmentsOut,
                                        FastaWriter[] shortPolyaFragmentsOut,
                                        FastaWriter longPolyaSingletonsOut,
                                        FastaWriter shortPolyaSingletonsOut) {
        
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
            int m = getMinCoverageOrderOfMagnitude(frag.minCov);
            
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
    
    public int[] assembleFragmentsMultiThreaded(FastxFilePair[] fastqPairs, 
                                                String[] longFragmentsFastaPaths,
                                                String[] shortFragmentsFastaPaths,
                                                String[] unconnectedReadsFastaPaths,
                                                String longSingletonsFasta,
                                                String shortSingletonsFasta,
                                                String unconnectedSingletonsFasta,
                                                String[] longPolyaFragmentsFastaPaths,
                                                String[] shortPolyaFragmentsFastaPaths,
                                                String[] unconnectedPolyaReadsFastaPaths,
                                                String longPolyaSingletonsFasta,
                                                String shortPolyaSingletonsFasta,
                                                String unconnectedPolyaSingletonsFasta,
                                                int bound,
                                                int minOverlap,
                                                int sampleSize, 
                                                int numThreads, 
                                                int maxErrCorrIterations,
                                                boolean extendFragments,
                                                int minKmerCov,
                                                boolean keepArtifact) {
        
        if (dbgFPR <= 0) {
            dbgFPR = graph.getDbgbf().getFPR();
        }
        
        if (covFPR <= 0) {
            covFPR = graph.getCbf().getFPR();
        }        
        
        long fragmentId = 0;
        long unconnectedReadId = 0;
        long readPairsParsed = 0;
        long readPairsConnected = 0;
        long readPairsNotConnected = 0;
        
        int maxTasksQueueSize = numThreads;
        int maxConcurrentSubmissions = numThreads + maxTasksQueueSize;
        
        int newBound = bound;
        int[] fragLengthsStats = null;
        boolean pairedKmerDistanceIsSet = false;
        int shortestFragmentLengthAllowed = k;
        
        boolean assemblePolyaTails = this.minPolyATailLengthRequired > 0;
                        
        try {
            
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
            
            FastaWriter[] longPolyaFragmentsOut = null;
            FastaWriter[]  shortPolyaFragmentsOut = null;
            FastaWriter[] unconnectedPolyaReadsOut = null;
            FastaWriter longPolyaSingletonsOut = null;
            FastaWriter shortPolyaSingletonsOut = null;
            FastaWriter unconnectedPolyaSingletonsOut = null;
            
            if (assemblePolyaTails) {
                longPolyaFragmentsOut = new FastaWriter[]{new FastaWriter(longPolyaFragmentsFastaPaths[0], true),
                                                            new FastaWriter(longPolyaFragmentsFastaPaths[1], true),
                                                            new FastaWriter(longPolyaFragmentsFastaPaths[2], true),
                                                            new FastaWriter(longPolyaFragmentsFastaPaths[3], true),
                                                            new FastaWriter(longPolyaFragmentsFastaPaths[4], true),
                                                            new FastaWriter(longPolyaFragmentsFastaPaths[5], true)};

                shortPolyaFragmentsOut = new FastaWriter[]{new FastaWriter(shortPolyaFragmentsFastaPaths[0], true),
                                                            new FastaWriter(shortPolyaFragmentsFastaPaths[1], true),
                                                            new FastaWriter(shortPolyaFragmentsFastaPaths[2], true),
                                                            new FastaWriter(shortPolyaFragmentsFastaPaths[3], true),
                                                            new FastaWriter(shortPolyaFragmentsFastaPaths[4], true),
                                                            new FastaWriter(shortPolyaFragmentsFastaPaths[5], true)};

                unconnectedPolyaReadsOut = new FastaWriter[]{new FastaWriter(unconnectedPolyaReadsFastaPaths[0], true),
                                                            new FastaWriter(unconnectedPolyaReadsFastaPaths[1], true),
                                                            new FastaWriter(unconnectedPolyaReadsFastaPaths[2], true),
                                                            new FastaWriter(unconnectedPolyaReadsFastaPaths[3], true),
                                                            new FastaWriter(unconnectedPolyaReadsFastaPaths[4], true),
                                                            new FastaWriter(unconnectedPolyaReadsFastaPaths[5], true)};
            
                longPolyaSingletonsOut = new FastaWriter(longPolyaSingletonsFasta, true);
                shortPolyaSingletonsOut = new FastaWriter(shortPolyaSingletonsFasta, true);
                unconnectedPolyaSingletonsOut = new FastaWriter(unconnectedPolyaSingletonsFasta, true);
            }
                        
            PairedReadSegments p;
            ArrayBlockingQueue<Fragment> fragments = new ArrayBlockingQueue<>(sampleSize);
            
            for (FastxFilePair fxPair: fastqPairs) {

                FastxPairReader fxpr = null;
                
                if (FastqReader.isFastq(fxPair.leftPath) && FastqReader.isFastq(fxPair.rightPath)) {
                    fxpr = new FastqPairReader(fxPair.leftPath, fxPair.rightPath, qualPatternFrag, seqPattern, fxPair.leftRevComp, fxPair.rightRevComp);
                }
                else if (FastaReader.isFasta(fxPair.leftPath) && FastaReader.isFasta(fxPair.rightPath)) {
                    fxpr = new FastaPairReader(fxPair.leftPath, fxPair.rightPath, seqPattern, fxPair.leftRevComp, fxPair.rightRevComp);
                }
                else {
                    throw new FileFormatException("Incompatible file format for `" + fxPair.leftPath + "` and `" + fxPair.rightPath + "`");
                }
                
                System.out.println("Parsing `" + fxPair.leftPath + "` and `" + fxPair.rightPath + "`...");

                int leftReadLengthThreshold = k;
                int rightReadLengthThreshold = k;
                
                // set up thread pool
                MyExecutorService service = new MyExecutorService(numThreads, maxTasksQueueSize);
                
                if (!pairedKmerDistanceIsSet) {
                                        
                    // Assembled initial sample of fragments
                    while (fxpr.hasNext()) {
                        p = fxpr.next();
                        ++readPairsParsed;
                                                
                        service.submit(new FragmentAssembler(p,
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
                        ));

                        if (fragments.remainingCapacity() == 0) {
                            break;
                        }
                    }
                    
                    if (!fxpr.hasNext() && fragments.isEmpty()) {
                        // entire file has been exhausted but no output fragments; wait for the workers to be done
                        service.terminate();
                    }
                    
                    if (fragments.isEmpty()) {
                        System.out.println("***************************************************************\n" +
                                           "* WARNING: Insufficient high-quality k-mers from input reads! *\n" + 
                                           "*          Output files may be empty!                         *\n" +
                                           "***************************************************************");
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
                    while (!fragments.isEmpty()) {
                        Fragment frag = fragments.poll();
                        
                        if (frag.isUnconnectedRead) {
                            ++readPairsNotConnected;
                            ++unconnectedReadId;
                            
                            String seq = frag.left;
                            if (seq != null && seq.length() >= k) {
                                FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                            unconnectedReadsOut, unconnectedSingletonsOut,
                                                            unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                                writer.write(Long.toString(unconnectedReadId) + "L", seq);
                            }
                            
                            seq = frag.right;
                            if (seq != null && seq.length() >= k) {
                                FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                            unconnectedReadsOut, unconnectedSingletonsOut,
                                                            unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                                writer.write(Long.toString(unconnectedReadId) + "R", seq);
                            }
                        }
                        else {
                            ++readPairsConnected;
                            
                            if (frag.length >= shortestFragmentLengthAllowed && frag.minCov > 0) {
                                ArrayList<Kmer> fragKmers = frag.kmers;

                                if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                    graph.addPairedKmers(fragKmers);
                                    
                                    if (frag.length >= longFragmentLengthThreshold) {
                                        for (Kmer kmer : fragKmers) {
                                            screeningBf.add(kmer.getHash());
                                        }
                                    }
                                    
                                    String header = Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]";
                                    String seq = graph.assemble(fragKmers);

                                    FastaWriter writer = assignFragmentFastaWriter(frag, seq, assemblePolyaTails,
                                                            longFragmentsOut, shortFragmentsOut, longSingletonsOut, shortSingletonsOut,
                                                            longPolyaFragmentsOut, shortPolyaFragmentsOut, longPolyaSingletonsOut, shortPolyaSingletonsOut);
                                    
                                    writer.write(header, seq);
                                }
                            }
                        }
                    }                    
                    
                    pairedKmerDistanceIsSet = true;
                }
                
                // assemble the remaining fragments in multi-threaded mode
                while (fxpr.hasNext()) {
                    p = fxpr.next();
                    ++readPairsParsed;
                        
                    service.submit(new FragmentAssembler(p,
                                                        fragments,
                                                        newBound, 
                                                        minOverlap, 
                                                        false, // do not store paired k-mers
                                                        maxErrCorrIterations, 
                                                        leftReadLengthThreshold, 
                                                        rightReadLengthThreshold,
                                                        extendFragments,
                                                        minKmerCov,
                                                        keepArtifact
                    ));

                    if (fragments.remainingCapacity() <= maxConcurrentSubmissions) {

                        // write fragments to file
                        for (int i=0; i<sampleSize; ++i) {
                            Fragment frag = fragments.poll();

                            if (frag == null) {
                                break;
                            }

                            if (frag.isUnconnectedRead) {
                                ++readPairsNotConnected;
                                ++unconnectedReadId;

                                String seq = frag.left;
                                if (seq != null && seq.length() >= k) {
                                    FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                                unconnectedReadsOut, unconnectedSingletonsOut,
                                                                unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                                    writer.write(Long.toString(unconnectedReadId) + "L", seq);
                                }

                                seq = frag.right;
                                if (seq != null && seq.length() >= k) {
                                    FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                                unconnectedReadsOut, unconnectedSingletonsOut,
                                                                unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                                    writer.write(Long.toString(unconnectedReadId) + "R", seq);
                                }
                            }
                            else {
                                ++readPairsConnected;

                                if (frag.length >= shortestFragmentLengthAllowed && frag.minCov > 0) {
                                    ArrayList<Kmer> fragKmers = frag.kmers;

                                    if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                        graph.addPairedKmers(fragKmers);
                                        
                                        if (frag.length >= longFragmentLengthThreshold) {
                                            for (Kmer kmer : fragKmers) {
                                                screeningBf.add(kmer.getHash());
                                            }
                                        }

                                        String header = Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]";
                                        String seq = graph.assemble(fragKmers);

                                        FastaWriter writer = assignFragmentFastaWriter(frag, seq, assemblePolyaTails,
                                                                longFragmentsOut, shortFragmentsOut, longSingletonsOut, shortSingletonsOut,
                                                                longPolyaFragmentsOut, shortPolyaFragmentsOut, longPolyaSingletonsOut, shortPolyaSingletonsOut);
                                        writer.write(header, seq);
                                    }
                                }
                            }
                        }
                    }
                }
                
                service.terminate();

                // write fragments to file
                while (!fragments.isEmpty()) {
                    Fragment frag = fragments.poll();
                    
                    if (frag.isUnconnectedRead) {
                        ++readPairsNotConnected;
                        ++unconnectedReadId;

                        String seq = frag.left;
                        if (seq != null && seq.length() >= k) {
                            FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                        unconnectedReadsOut, unconnectedSingletonsOut,
                                                        unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                            writer.write(Long.toString(unconnectedReadId) + "L", seq);
                        }

                        seq = frag.right;
                        if (seq != null && seq.length() >= k) {
                            FastaWriter writer = assignUnconnectedFastaWriter(frag, seq, assemblePolyaTails,
                                                        unconnectedReadsOut, unconnectedSingletonsOut,
                                                        unconnectedPolyaReadsOut, unconnectedPolyaSingletonsOut);
                            writer.write(Long.toString(unconnectedReadId) + "R", seq);
                        }
                    }
                    else {
                        ++readPairsConnected;
                        
                        if (frag.length >= shortestFragmentLengthAllowed && frag.minCov > 0) {
                            ArrayList<Kmer> fragKmers = frag.kmers;

                            if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                graph.addPairedKmers(fragKmers);
                                
                                if (frag.length >= longFragmentLengthThreshold) {
                                    for (Kmer kmer : fragKmers) {
                                        screeningBf.add(kmer.getHash());
                                    }
                                }

                                String header = Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]";
                                String seq = graph.assemble(fragKmers);

                                FastaWriter writer = assignFragmentFastaWriter(frag, seq, assemblePolyaTails,
                                                        longFragmentsOut, shortFragmentsOut, longSingletonsOut, shortSingletonsOut,
                                                        longPolyaFragmentsOut, shortPolyaFragmentsOut, longPolyaSingletonsOut, shortPolyaSingletonsOut);
                                writer.write(header, seq);
                            }
                        }
                    }
                }
                                
                fxpr.close();
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
            
            if (assemblePolyaTails) {
                unconnectedPolyaSingletonsOut.close();
                longPolyaSingletonsOut.close();
                shortPolyaSingletonsOut.close();

                for (FastaWriter out : longPolyaFragmentsOut) {
                    out.close();
                }

                for (FastaWriter out : shortPolyaFragmentsOut) {
                    out.close();
                }

                for (FastaWriter out : unconnectedPolyaReadsOut) {
                    out.close();
                }
            }
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(readPairsParsed) + " read pairs.");
            
            System.out.println("\tconnected:\t" + NumberFormat.getInstance().format(readPairsConnected) + "\t(" + readPairsConnected*100f/readPairsParsed + "%)");
            System.out.println("\tnot connected:\t" + NumberFormat.getInstance().format(readPairsNotConnected) + "\t(" + readPairsNotConnected*100f/readPairsParsed + "%)");
            
            long numDiscarded = readPairsParsed-readPairsConnected-readPairsNotConnected;
            System.out.println("\tdiscarded:\t" + NumberFormat.getInstance().format(numDiscarded) + "\t(" + numDiscarded*100f/readPairsParsed + "%)");
            
            System.out.println("Fragments paired kmers Bloom filter FPR: " + graph.getPkbfFPR() * 100   + " %");
            System.out.println("Screening Bloom filter FPR:              " + screeningBf.getFPR() * 100 + " %");
        }

        return fragLengthsStats;
    }

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
            graph.updateFragmentKmerDistance(graphFile);
        } catch (Exception ex) {
            handleException(ex);
        }
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
                                                    float minKmerCov) throws InterruptedException, IOException, Exception {
        
        long numFragmentsParsed = 0;
        FastaReader fin = new FastaReader(fragmentsFasta);

        ArrayBlockingQueue<String> fragmentsQueue = new ArrayBlockingQueue<>(sampleSize, true);
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

        try {
            while (true) {
                ++numFragmentsParsed;
                fragmentsQueue.put(fin.next());
            }
        }
        catch (NoSuchElementException e) {
            // end of file
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

    public void assembleTranscriptsMultiThreaded(String[] longFragmentsFastas, 
                                                String[] shortFragmentsFastas,
                                                String[] unconnectedReadsFastas,
                                                String longSingletonsFasta,
                                                String shortSingletonsFasta,
                                                String unconnectedSingletonsFasta,
                                                String[] longPolyaFragmentsFastas, 
                                                String[] shortPolyaFragmentsFastas,
                                                String[] unconnectedPolyaReadsFastas,
                                                String longPolyaSingletonsFasta,
                                                String shortPolyaSingletonsFasta,
                                                String unconnectedPolyaSingletonsFasta,
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
                                                String branchFreeExtensionThreshold) {
        
        long numFragmentsParsed = 0;

        try {
            boolean assemblePolya = minPolyATailLengthRequired > 0;

            FastaWriter fout = new FastaWriter(outFasta, false);
            FastaWriter foutShort = new FastaWriter(outFastaShort, false);
            //TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength, sensitiveMode ? maxTipLength : Math.max(k, maxTipLength));
            TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength, maxTipLength);

   
            boolean allowNaiveExtension = true;
            boolean extendBranchFreeOnly;
            
            if (assemblePolya) {
                // extend LONG fragments
                for (int mag=longPolyaFragmentsFastas.length-1; mag>=0; --mag) {
                    writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
                    String fragmentsFasta = longPolyaFragmentsFastas[mag];
                    System.out.println("Parsing `" + fragmentsFasta + "`...");
                    extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                    numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly, 
                                                                            keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
                }

                // extend SHORT fragments
                for (int mag=shortPolyaFragmentsFastas.length-1; mag>=0; --mag) {
                    writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
                    String fragmentsFasta = shortPolyaFragmentsFastas[mag];
                    System.out.println("Parsing `" + fragmentsFasta + "`...");
                    extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                    numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly,
                                                                            keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
                }

                // extend UNCONNECTED reads
                for (int mag=unconnectedPolyaReadsFastas.length-1; mag>=0; --mag) {
                    writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
                    String fragmentsFasta = unconnectedPolyaReadsFastas[mag];
                    System.out.println("Parsing `" + fragmentsFasta + "`...");
                    extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                    numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly,
                                                                            keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
                }

                extendBranchFreeOnly = isLowerStratum(STRATUM_01, branchFreeExtensionThreshold);

                // extend LONG singleton fragments
                writer.setOutputPrefix(txptNamePrefix + "01.L.");
                System.out.println("Parsing `" + longPolyaSingletonsFasta + "`...");
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(longPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

                // extend SHORT singleton fragments
                writer.setOutputPrefix(txptNamePrefix + "01.S.");
                System.out.println("Parsing `" + shortPolyaSingletonsFasta + "`...");
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(shortPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

                // extend UNCONNECTED reads
                writer.setOutputPrefix(txptNamePrefix + "01.U.");
                System.out.println("Parsing `" + unconnectedPolyaSingletonsFasta + "`...");
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(unconnectedPolyaSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }
            
            
            // extend LONG fragments
            for (int mag=longFragmentsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
                String fragmentsFasta = longFragmentsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly, 
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }          

            // extend SHORT fragments
            for (int mag=shortFragmentsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
                String fragmentsFasta = shortFragmentsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }
            
            // extend UNCONNECTED reads
            for (int mag=unconnectedReadsFastas.length-1; mag>=0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
                String fragmentsFasta = unconnectedReadsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                extendBranchFreeOnly = isLowerStratum(COVERAGE_ORDER[mag], branchFreeExtensionThreshold);
                numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly,
                                                                        keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            }

            extendBranchFreeOnly = isLowerStratum(STRATUM_01, branchFreeExtensionThreshold);

            // extend LONG singleton fragments
            writer.setOutputPrefix(txptNamePrefix + "01.L.");
            System.out.println("Parsing `" + longSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(longSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

            // extend SHORT singleton fragments
            writer.setOutputPrefix(txptNamePrefix + "01.S.");
            System.out.println("Parsing `" + shortSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(shortSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);

            // extend UNCONNECTED reads
            writer.setOutputPrefix(txptNamePrefix + "01.U.");
            System.out.println("Parsing `" + unconnectedSingletonsFasta + "`...");
            numFragmentsParsed += assembleTranscriptsMultiThreadedHelper(unconnectedSingletonsFasta, writer, sampleSize, numThreads,
                                                                    allowNaiveExtension, extendBranchFreeOnly,
                                                                    keepArtifact, keepChimera, reqFragKmersConsistency, minKmerCov);
            
            fout.close();
            foutShort.close();
            
        } catch (Exception ex) {
            handleException(ex);
        } finally {
            System.out.println("Parsed " + NumberFormat.getInstance().format(numFragmentsParsed) + " fragments.");
            System.out.println("Screening Bloom filter FPR:      " + screeningBf.getFPR() * 100 + " %");
        }
    }
    
    public void reduceSequenceRedundancy(String inFasta, String outFasta) {
        try {
            int numRemoved = reduceRedundancy(inFasta,
                                                outFasta,
                                                graph,
                                                screeningBf,
                                                lookahead,
                                                maxIndelSize,
                                                maxTipLength,
                                                percentIdentity);
            
            System.out.println("Removed " + NumberFormat.getInstance().format(numRemoved) + " redundant sequences.");
        }
        catch (Exception ex) {
            handleException(ex);
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

    public final static String FIELD_SEPARATOR = "\\s+"; // any white space character
    
    public static boolean getPooledReadPaths(String pooledReadPathsListFile,
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
                System.out.println("ERROR: Pool reads path file has unexpected number of columns on line " + lineNumber + ":\n\t" + line);
                return false;
            }
        }
        
        br.close();
        return true;
    }
    
    private static void assembleFragments(RNABloom assembler, boolean forceOverwrite,
            String outdir, String name, FastxFilePair[] fqPairs,
            long sbfSize, long pkbfSize, int sbfNumHash, int pkbfNumHash, int numThreads,
            int bound, int minOverlap, int sampleSize, int maxErrCorrItr, boolean extendFragments,
            int minKmerCoverage, boolean keepArtifact) {
        
        final String longFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.long.";
        final String shortFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.short.";
        final String unconnectedReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.";
        final String longPolyaFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.polya.long.";
        final String shortPolyaFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.polya.short.";
        final String unconnectedPolyaReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.polya.";
        
        final String[] longFragmentsFastaPaths = {longFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        longFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        longFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        longFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        longFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        longFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

        final String[] shortFragmentsFastaPaths = {shortFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        shortFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

        final String[] unconnectedReadsFastaPaths = {unconnectedReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        unconnectedReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        unconnectedReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        unconnectedReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        unconnectedReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        unconnectedReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

        final String[] longPolyaFragmentsFastaPaths = {longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

        final String[] shortPolyaFragmentsFastaPaths = {shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

        final String[] unconnectedPolyaReadsFastaPaths = {unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                        unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                        unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                        unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                        unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                        unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};
        
        final String longSingletonsFastaPath = longFragmentsFastaPrefix + "01.fa";
        final String shortSingletonsFastaPath = shortFragmentsFastaPrefix + "01.fa";
        final String unconnectedSingletonsFastaPath = unconnectedReadsFastaPrefix + "01.fa";
        
        final String longPolyaSingletonsFastaPath = longPolyaFragmentsFastaPrefix + "01.fa";
        final String shortPolyaSingletonsFastaPath = shortPolyaFragmentsFastaPrefix + "01.fa";
        final String unconnectedPolyaSingletonsFastaPath = unconnectedPolyaReadsFastaPrefix + "01.fa";

        final File fragsDoneStamp = new File(outdir + File.separator + STAMP_FRAGMENTS_DONE);
        
        final String fragStatsFile = outdir + File.separator + name + ".fragstats";
        final String graphFile = outdir + File.separator + name + ".graph";
        
        if (forceOverwrite || !fragsDoneStamp.exists()) {
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

            for (String fragmentsFasta : longPolyaFragmentsFastaPaths) {
                fragmentsFile = new File(fragmentsFasta);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
            }

            for (String fragmentsFasta : shortPolyaFragmentsFastaPaths) {
                fragmentsFile = new File(fragmentsFasta);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
            }

            for (String fragmentsFasta : unconnectedPolyaReadsFastaPaths) {
                fragmentsFile = new File(fragmentsFasta);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }
            }

            fragmentsFile = new File(longPolyaSingletonsFastaPath);
            if (fragmentsFile.exists()) {
                fragmentsFile.delete();
            }

            fragmentsFile = new File(shortPolyaSingletonsFastaPath);
            if (fragmentsFile.exists()) {
                fragmentsFile.delete();
            }

            fragmentsFile = new File(unconnectedPolyaSingletonsFastaPath);
            if (fragmentsFile.exists()) {
                fragmentsFile.delete();
            }

            assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
            assembler.setupFragmentPairedKmersBloomFilter(pkbfSize, pkbfNumHash);

            int[] fragStats = assembler.assembleFragmentsMultiThreaded(fqPairs, 
                                            longFragmentsFastaPaths, 
                                            shortFragmentsFastaPaths,
                                            unconnectedReadsFastaPaths,
                                            longSingletonsFastaPath,
                                            shortSingletonsFastaPath,
                                            unconnectedSingletonsFastaPath,
                                            longPolyaFragmentsFastaPaths, 
                                            shortPolyaFragmentsFastaPaths,
                                            unconnectedPolyaReadsFastaPaths,
                                            longPolyaSingletonsFastaPath,
                                            shortPolyaSingletonsFastaPath,
                                            unconnectedPolyaSingletonsFastaPath,
                                            bound, 
                                            minOverlap,
                                            sampleSize,
                                            numThreads,
                                            maxErrCorrItr,
                                            extendFragments,
                                            minKmerCoverage,
                                            keepArtifact);

            assembler.savePairedKmersBloomFilter(new File(graphFile));
            assembler.writeFragStatsToFile(fragStats, fragStatsFile);

            try {
                touch(fragsDoneStamp);
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(1);
            }
        }
        else {
            System.out.println("WARNING: Fragments were already assembled for \"" + name + "!");
        }
    }
    
    private static void assembleTranscripts(RNABloom assembler, boolean forceOverwrite,
            String outdir, String name, String txptNamePrefix, boolean strandSpecific,
            long sbfSize, int sbfNumHash, int numThreads, boolean noFragDBG,
            int sampleSize, int minTranscriptLength, boolean keepArtifact, boolean keepChimera, 
            boolean reqFragKmersConsistency, boolean restorePairedKmersBloomFilter,
            float minKmerCov, String branchFreeExtensionThreshold,
            boolean reduceRedundancy, boolean assemblePolya) {
        
        final File txptsDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_DONE);
        final File txptsNrDoneStamp = new File(outdir + File.separator + STAMP_TRANSCRIPTS_NR_DONE);
        final String transcriptsFasta = outdir + File.separator + name + ".transcripts.fa";
        final String shortTranscriptsFasta = outdir + File.separator + name + ".transcripts.short.fa";
                
        if (forceOverwrite || !txptsDoneStamp.exists()) {
            MyTimer timer = new MyTimer();
            
            if (restorePairedKmersBloomFilter) {
                assembler.restorePairedKmersBloomFilter(new File(outdir + File.separator + name + ".graph"));
            }
            
            final String longFragmentsFastaPrefix =         outdir + File.separator + name + ".fragments.long.";
            final String shortFragmentsFastaPrefix =        outdir + File.separator + name + ".fragments.short.";
            final String unconnectedReadsFastaPrefix =      outdir + File.separator + name + ".unconnected.";
            final String longPolyaFragmentsFastaPrefix =    outdir + File.separator + name + ".fragments.polya.long.";
            final String shortPolyaFragmentsFastaPrefix =   outdir + File.separator + name + ".fragments.polya.short.";
            final String unconnectedPolyaReadsFastaPrefix = outdir + File.separator + name + ".unconnected.polya.";

            final String[] longFragmentsFastaPaths = {longFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            longFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String[] shortFragmentsFastaPaths = {shortFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            shortFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String[] unconnectedReadsFastaPaths = {unconnectedReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            unconnectedReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String longSingletonsFastaPath = longFragmentsFastaPrefix + "01.fa";
            final String shortSingletonsFastaPath = shortFragmentsFastaPrefix + "01.fa";
            final String unconnectedSingletonsFastaPath = unconnectedReadsFastaPrefix + "01.fa";
            
            final String[] longPolyaFragmentsFastaPaths = {longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            longPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String[] shortPolyaFragmentsFastaPaths = {shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            shortPolyaFragmentsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String[] unconnectedPolyaReadsFastaPaths = {unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            unconnectedPolyaReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            final String longPolyaSingletonsFastaPath = longPolyaFragmentsFastaPrefix + "01.fa";
            final String shortPolyaSingletonsFastaPath = shortPolyaFragmentsFastaPrefix + "01.fa";
            final String unconnectedPolyaSingletonsFastaPath = unconnectedPolyaReadsFastaPrefix + "01.fa";
            
            final String graphFile = outdir + File.separator + name + ".graph";
            
            if (!noFragDBG) {
                if (assembler.isGraphInitialized()) {
                    assembler.clearDbgBf();
                }

                // repopulate with kmers from fragments
                ArrayList<String> fragmentPaths = new ArrayList<>(longFragmentsFastaPaths.length + shortFragmentsFastaPaths.length + unconnectedReadsFastaPaths.length + 3);
                fragmentPaths.addAll(Arrays.asList(longFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(shortFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(unconnectedReadsFastaPaths));
                fragmentPaths.add(longSingletonsFastaPath);
                fragmentPaths.add(shortSingletonsFastaPath);
                fragmentPaths.add(unconnectedSingletonsFastaPath);

                if (assemblePolya) {
                    fragmentPaths.addAll(Arrays.asList(longPolyaFragmentsFastaPaths));
                    fragmentPaths.addAll(Arrays.asList(shortPolyaFragmentsFastaPaths));
                    fragmentPaths.addAll(Arrays.asList(unconnectedPolyaReadsFastaPaths));
                    fragmentPaths.add(longPolyaSingletonsFastaPath);
                    fragmentPaths.add(shortPolyaSingletonsFastaPath);
                    fragmentPaths.add(unconnectedPolyaSingletonsFastaPath);
                }
                
                System.out.println("Rebuilding graph from assembled fragments...");
                timer.start();
                assembler.populateGraphFromFragments(fragmentPaths, strandSpecific, false);
                System.out.println("Graph rebuilt in " + MyTimer.hmsFormat(timer.elapsedMillis()));  
            }

            File transcriptsFile = new File(transcriptsFasta);
            if (transcriptsFile.exists()) {
                transcriptsFile.delete();
            }
            
            File shortTranscriptsFile = new File(shortTranscriptsFasta);
            if (shortTranscriptsFile.exists()) {
                shortTranscriptsFile.delete();
            }

            System.out.println("Assembling transcripts...");
            timer.start();

            assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);

            assembler.assembleTranscriptsMultiThreaded(longFragmentsFastaPaths, 
                                                        shortFragmentsFastaPaths,
                                                        unconnectedReadsFastaPaths,
                                                        longSingletonsFastaPath,
                                                        shortSingletonsFastaPath,
                                                        unconnectedSingletonsFastaPath,
                                                        longPolyaFragmentsFastaPaths, 
                                                        shortPolyaFragmentsFastaPaths,
                                                        unconnectedPolyaReadsFastaPaths,
                                                        longPolyaSingletonsFastaPath,
                                                        shortPolyaSingletonsFastaPath,
                                                        unconnectedPolyaSingletonsFastaPath,
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
                                                        branchFreeExtensionThreshold);
            
            System.out.println("Transcripts assembled in " + MyTimer.hmsFormat(timer.elapsedMillis()));

            try {
                touch(txptsDoneStamp);
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(1);
            }

            System.out.println("Assembled transcripts at `" + transcriptsFasta + "`");
        }
        else {
            System.out.println("WARNING: Transcripts were already assembled for \"" + name + "\"!");
        }
        
        if (reduceRedundancy) {
            final String transcriptsNrFasta = outdir + File.separator + name + ".transcripts.nr.fa";

            if (forceOverwrite || !txptsNrDoneStamp.exists()) {
                File transcriptsNrFile = new File(transcriptsNrFasta);
                if (transcriptsNrFile.exists()) {
                    transcriptsNrFile.delete();
                }
                
                MyTimer timer = new MyTimer();

                System.out.println("Reducing redundancy...");
                timer.start();
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                
                assembler.reduceSequenceRedundancy(transcriptsFasta, transcriptsNrFasta);

                System.out.println("Redundancy reduction completed in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                try {
                    touch(txptsNrDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
                
                System.out.println("Non-redundant transcripts at `" + transcriptsNrFasta + "`");
            }   
            else {
                System.out.println("WARNING: Redundancy reduction already completed for \"" + name + "\"!");
            }
        }
    }
    
    private static long getNumUniqueKmers(int threads, int k, String histogramPathPrefix, String[] leftReadPaths, String[] rightReadPaths, boolean forceOverwrite) throws IOException, InterruptedException {
        long numKmers = -1L;
        String histogramPath = histogramPathPrefix + "_k" + k + ".hist";
        int exitVal = 0;
        
        if (forceOverwrite || !new File(histogramPath).isFile()) {
            String cmd = "ntcard -t " + threads + " -k " + k + " -c 65535 -p " + histogramPathPrefix + " " + String.join(" ", leftReadPaths) + " " + String.join(" ", rightReadPaths);
            Runtime rt = Runtime.getRuntime();
            System.out.println("Running command: `" + cmd + "`...");
            Process pr = rt.exec(cmd);
            exitVal = pr.waitFor();
        }
        
        if (exitVal == 0) {
            System.out.println("Parsing histogram file `" + histogramPath + "`...");
            BufferedReader br = new BufferedReader(new FileReader(histogramPath));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.length() > 0) {
                    String[] cols = line.split("\t");
                    if (cols[0].equals("F0")) {
                        numKmers = Long.parseLong(cols[1]);
                        br.close();
                        break; 
                    }
                }
            }
            br.close();
        }
        
        return numKmers;
    }
    
    private static final String STAMP_STARTED = "STARTED";
    private static final String STAMP_DBG_DONE = "DBG.DONE";
    private static final String STAMP_FRAGMENTS_DONE = "FRAGMENTS.DONE";
    private static final String STAMP_TRANSCRIPTS_DONE = "TRANSCRIPTS.DONE";
    private static final String STAMP_TRANSCRIPTS_NR_DONE = "TRANSCRIPTS_NR.DONE";
    
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
                                    .desc("assembled transcript name prefix in FASTA header")
                                    .hasArg(true)
                                    .argName("STR")
                                    .build();
        options.addOption(optPrefix);
        
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
        
        final String optStageDefault = "3";
        Option optStage = Option.builder("stage")
                                    .desc("assembly end point, eg. '1' (construct graph), '2' (assemble fragments), '3' (assemble transcripts) [" + optStageDefault + "]")
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
                                    .desc("run `ntCard` to count the number of unique k-mers in input reads [false]")
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

        final String optFprDefault = "0.10";
        Option optFpr = Option.builder("fpr")
                                    .longOpt("fpr")
                                    .desc("maximum allowable false-positive rate of Bloom filters [" + optFprDefault + "]")
                                    .hasArg(true)
                                    .argName("DECIMAL")
                                    .build();
        options.addOption(optFpr);
        
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

        final String optSampleDefault = "1000";
        Option optSample = Option.builder("sample")
                                    .desc("sample size for estimating fragment lengths [" + optSampleDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optSample);
        
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
        
        Option optReduce = Option.builder("nr")
                                    .desc("output non-redundant transcripts in 'transcripts.nr.fa' [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optReduce);
        
        final String optErrCorrItrDefault = "1";
        Option optErrCorrItr = Option.builder("e")
                                    .longOpt("errcorritr")
                                    .desc("number of iterations of error-correction in reads [" + optErrCorrItrDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optErrCorrItr);        

        Option optExtend = Option.builder("extend")
                                    .desc("extend fragments outward during fragment reconstruction [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optExtend);

        Option optNoFragDBG = Option.builder("nofdbg")
                                    .desc("do not rebuild DBG from fragment k-mers [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optNoFragDBG);

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
                                    .desc("do not remove potential sequencing artifact [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optKeepArtifact);

        Option optKeepChimera = Option.builder("chimera")
                                    .desc("do not remove potential chimera [false]")
                                    .hasArg(false)
                                    .build();
        options.addOption(optKeepChimera);
        
        final String optBranchFreeExtensionDefault = STRATUM_E0;
        final String optBranchFreeExtensionChoicesStr = String.join("|", STRATA);
        Option optBranchFreeExtensionThreshold = Option.builder("stratum")
                                    .desc("fragments lower than the specified stratum are extended only if they are branch-free [" + optBranchFreeExtensionDefault + "]")
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
        
        final String optMinLengthDefault = "200";
        Option optMinLength = Option.builder("length")
                                    .desc("minimum transcript length in output assembly [" + optMinLengthDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optMinLength);  
        
        final String optPolyATailDefault = "0";
        Option optPolyATail = Option.builder("a")
                                    .longOpt("polya")
                                    .desc("prioritize assembly of transcripts with poly-A tails of the minimum length specified [" + optPolyATailDefault + "]")
                                    .hasArg(true)
                                    .argName("INT")
                                    .build();
        options.addOption(optPolyATail);  
        
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
                System.out.println("ERROR: Unknown stratum name specified, \"" + branchFreeExtensionThreshold + "\"");
                System.exit(1);
            }
            
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
            }
                        
            final String pooledReadsListFile = line.getOptionValue(optPooledAssembly.getOpt());
            String[] leftReadPaths = line.getOptionValues(optLeftReads.getOpt());
            String[] rightReadPaths = line.getOptionValues(optRightReads.getOpt());
            final boolean pooledGraphMode = pooledReadsListFile != null;
            
            HashMap<String, ArrayList<String>> pooledLeftReadPaths = new HashMap<>();
            HashMap<String, ArrayList<String>> pooledRightReadPaths = new HashMap<>();
            
            if (pooledGraphMode) {
                System.out.println("Pooled assembly mode is ON!");
                
                if (!new File(pooledReadsListFile).isFile()) {
                    System.out.println("ERROR: Cannot find pooled read paths list `" + pooledReadsListFile + "`");
                    System.exit(1);
                }
                
                System.out.println("Parsing pool reads list file `" + pooledReadsListFile + "`...");
                boolean parseOK = getPooledReadPaths(pooledReadsListFile, pooledLeftReadPaths, pooledRightReadPaths);
                
                if (!parseOK) {
                    System.out.println("ERROR: Incorrect format of pooled read paths list file!");
                    System.exit(1);
                }
                
                int numLeftIds = pooledLeftReadPaths.size();
                int numRightIds = pooledRightReadPaths.size();
                
                if (numLeftIds != numRightIds) {
                    System.out.println("ERROR: Pooled read paths list file has disagreeing number of sample IDs for left (" + numLeftIds + ") and right (" + numRightIds + ") reads!");
                    System.exit(1);
                }
                
                if (numLeftIds == 0) {
                    System.out.println("ERROR: Pooled read paths list file is empty!");
                    System.exit(1);                    
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
            }
            else {
                if (leftReadPaths == null || leftReadPaths.length == 0) {
                    System.out.println("ERROR: Please specify left read files!");
                    System.exit(1);
                }

                if (rightReadPaths == null || rightReadPaths.length == 0) {
                    System.out.println("ERROR: Please specify right read files!");
                    System.exit(1);
                }

                if (leftReadPaths.length != rightReadPaths.length) {
                    System.out.println("ERROR: Read files are not paired properly!");
                    System.exit(1);
                }
            }
            
            final boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            final boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            final boolean strandSpecific = line.hasOption(optStranded.getOpt());
            
            final int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), optKmerSizeDefault));
            
            final int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), optBaseQualDbgDefault));
            final int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), optBaseQualFragDefault));
            final int minKmerCov = Integer.parseInt(line.getOptionValue(optMinKmerCov.getOpt(), optMinKmerCovDefault));
            
            double leftReadFilesTotalBytes = 0;
            double rightReadFilesTotalBytes = 0;
            
            for (String fq : leftReadPaths) {
                leftReadFilesTotalBytes += new File(fq).length();
            }
            for (String fq : rightReadPaths) {
                rightReadFilesTotalBytes += new File(fq).length();
            }
            
            final float maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (Math.max(NUM_BYTES_1MB * 100, 0.80f * Math.min(leftReadFilesTotalBytes, rightReadFilesTotalBytes)) / NUM_BYTES_1GB))));
            float sbfGB = Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 8.5f)));
            float dbgGB = Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), Float.toString(maxBfMem * 1f / 8.5f)));
            float cbfGB = Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), Float.toString(maxBfMem * 6f / 8.5f)));
            float pkbfGB = Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 8.5f)));
                        
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
            
            long expNumKmers = -1L;
            if (line.hasOption(optNtcard.getOpt())) {
                System.out.println("\nK-mer counting with ntcard...");
                String histogramPathPrefix = outdir + File.separator + name;
                
                timer.start();
                expNumKmers = getNumUniqueKmers(numThreads, k, histogramPathPrefix, leftReadPaths, rightReadPaths, forceOverwrite);
                    
                if (expNumKmers < 0) {
                    System.out.println("ERROR: Cannot get number of unique k-mers from ntcard!");
                    System.exit(1);
                }
                
                System.out.println("Number of unique k-mers: " + NumberFormat.getInstance().format(expNumKmers));
                System.out.println("K-mer counting completed in " + MyTimer.hmsFormat(timer.elapsedMillis()));                
            }
            else {
                expNumKmers = Long.parseLong(line.getOptionValue(optNumKmers.getOpt(), "-1"));
            }
            
            if (expNumKmers > 0) {
                sbfSize = BloomFilter.getExpectedSize(expNumKmers, maxFPR, sbfNumHash);
                sbfGB = sbfSize / (float) NUM_BITS_1GB;
                
                dbgbfSize = BloomFilter.getExpectedSize(expNumKmers, maxFPR, dbgbfNumHash);
                dbgGB = dbgbfSize / (float) NUM_BITS_1GB;
                
                cbfSize = CountingBloomFilter.getExpectedSize(expNumKmers, maxFPR, cbfNumHash);
                cbfGB = cbfSize / (float) NUM_BYTES_1GB;
                
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
            final float percentIdentity = Float.parseFloat(line.getOptionValue(optPercentIdentity.getOpt(), optPercentIdentityDefault));
            final boolean outputNrTxpts = line.hasOption(optReduce.getOpt());
            final int maxIndelSize = Integer.parseInt(line.getOptionValue(optIndelSize.getOpt(), optIndelSizeDefault));
            int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), optErrCorrItrDefault));
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
            
            final boolean noFragDBG = line.hasOption(optNoFragDBG.getOpt());
            final boolean reqFragKmersConsistency = !line.hasOption(optNoFragmentsConsistency.getOpt());
            final boolean extendFragments = line.hasOption(optExtend.getOpt());
            final int minNumKmerPairs = Integer.parseInt(line.getOptionValue(optMinKmerPairs.getOpt(), optMinKmerPairsDefault));
            final String txptNamePrefix = line.getOptionValue(optPrefix.getOpt(), optPrefixDefault);
            

            System.out.println("\nBloom filters         Memory (GB)");
            System.out.println("====================================");
            System.out.println("de Bruijn graph:       " + dbgGB);
            System.out.println("k-mer counting:        " + cbfGB);
            System.out.println("paired k-mers (reads): " + pkbfGB);
            System.out.println("paired k-mers (frags): " + pkbfGB);
            System.out.println("screening:             " + sbfGB);
            System.out.println("====================================");
            System.out.println("Total:                 " + (dbgGB+cbfGB+2*pkbfGB+sbfGB));
            
            RNABloom assembler = new RNABloom(k, qDBG, qFrag);
            assembler.setParams(strandSpecific, maxTipLen, lookahead, maxCovGradient, maxIndelSize, percentIdentity, minNumKmerPairs, minPolyATail);

            try {
                FileWriter writer = new FileWriter(startedStamp, false);
                writer.write(String.join(" ", args));
                writer.close();
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(1);
            }
            
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
                
                if (revCompLeft) {
                    backwardFilesList.addAll(Arrays.asList(leftReadPaths));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(leftReadPaths));
                }
                
                if (revCompRight) {
                    backwardFilesList.addAll(Arrays.asList(rightReadPaths));
                }
                else {
                    forwardFilesList.addAll(Arrays.asList(rightReadPaths));
                }
                       
                System.out.println("\n* Stage 1: Construct graph from reads (k=" + k + ")");
                timer.start();
                
                assembler.initializeGraph(strandSpecific, 
                        dbgbfSize, cbfSize, pkbfSize, 
                        dbgbfNumHash, cbfNumHash, pkbfNumHash, false, true);
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads, false, true);
                
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
                    
                    System.out.println("Bloom filters          Memory (GB)");
                    System.out.println("====================================");
                    System.out.println("de Bruijn graph:       " + dbgGB);
                    System.out.println("k-mer counting:        " + cbfGB);
                    System.out.println("paired k-mers (reads): " + pkbfGB);
                    System.out.println("paired k-mers (frags): " + pkbfGB);
                    System.out.println("screening:             " + sbfGB);
                    System.out.println("====================================");
                    System.out.println("Total:                 " + (dbgGB+cbfGB+2*pkbfGB+sbfGB));
                    
                    assembler.initializeGraph(strandSpecific, 
                            dbgbfSize, cbfSize, pkbfSize, 
                            dbgbfNumHash, cbfNumHash, pkbfNumHash, false, true);
                    
                    assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                    
                    System.out.println("Repopulate graph ...");
                    
                    assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads, false, true);
                }    
                
                
                System.out.println("Saving graph to file `" + graphFile + "`...");
                assembler.saveGraph(new File(graphFile));
                
                System.out.println("* Stage 1 completed in " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                touch(dbgDoneStamp);
                
                if (endstage <= 1) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
            }
                        

            if (pooledGraphMode) {
                // assemble fragments for each sample
                int numSamples = pooledLeftReadPaths.size();
                int sampleId = 0;
                
                System.out.println("\n* Stage 2: Assemble fragments for " + numSamples + " samples");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                for (String sampleName : pooledLeftReadPaths.keySet()) {
                    System.out.println("** Working on \"" + sampleName + "\" (sample " + ++sampleId + " of " + numSamples + ")...");
                    
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
                                    bound, minOverlap, sampleSize, maxErrCorrItr, extendFragments, minKmerCov, !keepArtifact);
                    
                    System.out.println("** Fragments assembled in " + MyTimer.hmsFormat(sampleTimer.elapsedMillis()) + "\n");
                }
                
                System.out.println("* Stage 2 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                touch(fragsDoneStamp);
                
                if (endstage <= 2) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }
                
                // assemble transcripts for each sample
                sampleId = 0;
                System.out.println("\n* Stage 3: Assemble transcripts for " + numSamples + " samples");
                stageTimer.start();
                
                for (String sampleName : pooledLeftReadPaths.keySet()) {
                    System.out.println("** Working on \"" + sampleName + "\" (sample " + ++sampleId + " of " + numSamples + ")...");
                    
                    String sampleOutdir = outdir + File.separator + sampleName;
                    
                    assembleTranscripts(assembler, forceOverwrite,
                                    sampleOutdir, sampleName, txptNamePrefix, strandSpecific,
                                    sbfSize, sbfNumHash, numThreads, noFragDBG,
                                    sampleSize, minTranscriptLength, keepArtifact, keepChimera,
                                    reqFragKmersConsistency, true, minKmerCov,
                                    branchFreeExtensionThreshold, outputNrTxpts, minPolyATail > 0);
                    
                    System.out.print("\n");
                }
                
                System.out.println("* Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));                
                
                touch(txptsDoneStamp);
            }
            else {
                FastxFilePair[] fqPairs = new FastxFilePair[leftReadPaths.length];
                for (int i=0; i<leftReadPaths.length; ++i) {
                    fqPairs[i] = new FastxFilePair(leftReadPaths[i], rightReadPaths[i], revCompLeft, revCompRight);
                }

                System.out.println("\n* Stage 2: Assemble fragments for \"" + name + "\"");
                MyTimer stageTimer = new MyTimer();
                stageTimer.start();
                
                assembleFragments(assembler, forceOverwrite,
                                    outdir, name, fqPairs,
                                    sbfSize, pkbfSize, sbfNumHash, pkbfNumHash, numThreads,
                                    bound, minOverlap, sampleSize, maxErrCorrItr, extendFragments, minKmerCov, keepArtifact);
                
                System.out.println("* Stage 2 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
                
                if (endstage <= 2) {
                    System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
                    System.exit(0);
                }

                System.out.println("\n* Stage 3: Assemble transcripts for \"" + name + "\"");
                stageTimer.start();
                
                assembleTranscripts(assembler, forceOverwrite,
                                outdir, name, txptNamePrefix, strandSpecific,
                                sbfSize, sbfNumHash, numThreads, noFragDBG,
                                sampleSize, minTranscriptLength, keepArtifact, keepChimera,
                                reqFragKmersConsistency, false, minKmerCov, 
                                branchFreeExtensionThreshold, outputNrTxpts, minPolyATail > 0);
                
                System.out.println("* Stage 3 completed in " + MyTimer.hmsFormat(stageTimer.elapsedMillis()));
            }      
        }
        catch (Exception exp) {
            System.out.println("ERROR: " + exp.getMessage() );
            exp.printStackTrace();
            System.exit(1);
        }
        
        System.out.println("Total runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
    }
}
