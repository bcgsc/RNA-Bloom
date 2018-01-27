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
 * @author kmnip
 */
public class RNABloom {
    public final static String VERSION = "0.9.0";
    
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
//    private Pattern homoPolymerKmerPattern;
    private BloomFilterDeBruijnGraph graph = null;
    private BloomFilter screeningBf = null;

    private int maxTipLength;
    private int lookahead;
    private float maxCovGradient;
    private int maxIndelSize;
    private int minPolyATailLengthRequired;
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
        this.percentError = 1.0f - percentIdentity;
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
        graph.clearDbgbf();
        graph.clearCbf();
        graph.clearPkbf();
        graph.clearRpkbf();
        
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
        private int numReads = 0;
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
                kmerPairDistance = graph.getReadKmerDistance();
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
                for (int i=0; i< 10 && fr.hasNext(); ++i) {
                    fr.nextWithoutName(r);
                    max = Math.max(max, r.seq.length());
                }
                fr.close();
            }
            else if (FastaReader.isFasta(path)) {
                FastaReader fr = new FastaReader(path);
                for (int i=0; i< 10 && fr.hasNext(); ++i) {
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
        }
        
        if (readLength < 0) {
            exitOnError("Cannot determine read length from read files.");
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
        
        graph.setReadKmerDistance(readLength - k - minNumKmerPairs);
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
        
        int numReads = 0;
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
        
        dbgFPR = graph.getDbgbfFPR();
        covFPR = graph.getCbfFPR();
        
        System.out.println("Screening Bloom filter FPR:  " + screeningBf.getFPR() * 100 + " %");
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
                        if (!hasPolyATail) {
                            // skip this transcript because it does not contain poly A tail
                            return;
                        }
                        
                        pasPositions = getPolyASignalPositions(transcript, polyASignalPattern, polyATailOnlyMatchingPattern);
                    }
                    else {
                        hasPolyTHead = polyTHeadOnlyPattern.matcher(transcript).matches();
                        
                        if (!hasPolyATail && !hasPolyTHead) {
                            return;
                        }
                        
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
                    // add PAS and their positions to header
                    headerBuilder.append(" PAS=[");                    
                    
                    if (pasPositions != null && !pasPositions.isEmpty()) {
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
                    }
                    
                    headerBuilder.append("]");
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
                int fragKmersDist = graph.getPairedKmerDistance();
                
                while (true) {
                    String fragment = fragments.poll(10, TimeUnit.MICROSECONDS);
                    
//                    fragment = "AGGTGACCGGAATGTAGCGACCACTACTCGTGTAACTACTGGTAGACCTACAGTAGGACCTTGGCAAATCTGTCCAAATCAACAAACCGTATTCAAGGGGCAAGTCGGAGTAATTTATGAGATGCTTCCAGCCTCCACTCAACAAAATGCTATCAACGCATTTGTGGAGAATGTGCTATTGAATTCAAACTTCTACGGATTGGCTTTGGATATTCTAACCCATGACAATAGAACACTTGTCACTGCTATTCCATACCCAACAACTGACAAATACAATGTTCAAGGATATGGAAGTGCGAGAAGTGTAGATGAATTTAGAAATCAAGTTATTGCGTTCAATGGAATGATTCAACCTATTTCACAATCCAGCAGTATTTCAGATGCACTTTTGTACGTTACGACCAGTCTGCCAGCTGG";
                    
                    if (fragment == null) {
                        if (!keepGoing) {
                            break;
                        }
                    }
                    else {
//                        if (fragment.equals("CACAATCAACACATCACTGGACACATACAGCCACTTTCTCGCAGAAACATCTCATCTGAAAATTATGAACTTCAACATGAACTATAGCACCGAACCGTCTCAATACCAATCTTACGTCATCGCGAACGCTCCTTCCGCCACAGCGCCAATTTGTACAGAAATCACGTGCTGTGTTTGTTGAGCCTACGGGTCATGCATACGGACG")) {
//                            System.out.println(fragment);
//                        }
                        
                        if (minPolyATailLengthRequired > 0) {
                            if (!polyATailPattern.matcher(fragment).matches()) {
                                // skip this fragment because it does not contain poly A tail
                                continue;
                            }
                        }

                        ArrayList<Kmer> fragKmers = graph.getKmers(fragment);
                        
                        if (!fragKmers.isEmpty()) {
//                            ArrayList<Kmer2> fragKmers2 = correctErrorsSE(fragKmers, graph, lookahead, maxIndelSize, maxCovGradient, covFPR, percentIdentity);
//                            if (fragKmers2 != null) {
//                                fragKmers = fragKmers2;
//                            }
                            
                            ArrayList<Kmer> fragKmers2 = new ArrayList<>(fragKmers);

                            if ((!extendBranchFreeFragmentsOnly || isBranchFree(fragKmers, graph, maxTipLength)) &&
                                    !represented(fragKmers,
                                                        graph,
                                                        screeningBf,
                                                        lookahead,
                                                        maxIndelSize,
                                                        maxTipLength,
                                                        percentIdentity)) {
                                
                                if (includeNaiveExtensions) {
                                    extendWithPairedKmers(fragKmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs, maxCovGradient);
                                }
                                else {
                                    extendWithPairedKmersDFS(fragKmers, graph, lookahead, maxTipLength, screeningBf, maxIndelSize, percentIdentity, minNumKmerPairs, maxCovGradient);
                                }

                                if (fragKmers.size() > fragKmersDist) {
                                    ArrayDeque<ArrayList<Kmer>> fragSegments = breakWithPairedKmers(fragKmers, graph);
                                    int numFragSegs = fragSegments.size();

                                    if (numFragSegs >= 1) {
                                        for (ArrayList<Kmer> f : fragSegments) {
                                            if (numFragSegs == 1 || new HashSet<>(f).containsAll(fragKmers2)) {
                                                ArrayDeque<ArrayList<Kmer>> readSegments = breakWithReadPairedKmers(f, graph);

                                                int numReadSegs = readSegments.size();

                                                if (numReadSegs == 1) {
                                                    transcripts.put(new Transcript(fragment, f));
                                                }
                                                else if (numReadSegs > 1) {
                                                    for (ArrayList<Kmer> r : readSegments) {
                                                        if (new HashSet<>(r).containsAll(fragKmers2)) {
                                                            transcripts.put(new Transcript(fragment, r));
                                                            break;
                                                        }
                                                    }
                                                }

                                                break;
                                            }
                                        }
                                    }
                                }
                                else {
                                    ArrayDeque<ArrayList<Kmer>> readSegments = breakWithReadPairedKmers(fragKmers, graph);

                                    int numReadSegs = readSegments.size();

                                    if (numReadSegs == 1) {
                                        transcripts.put(new Transcript(fragment, fragKmers));
                                    }
                                    else if (numReadSegs > 1) {
                                        for (ArrayList<Kmer> r : readSegments) {
                                            if (new HashSet<>(r).containsAll(fragKmers2)) {
                                                transcripts.put(new Transcript(fragment, r));
                                                break;
                                            }
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
                                                        percentIdentity);

                    if (correctedReadPair.corrected) {
                        leftKmers = correctedReadPair.leftKmers;
                        rightKmers = correctedReadPair.rightKmers;
                    }
                }
                
                if (!leftKmers.isEmpty() && !rightKmers.isEmpty()) {

                    ArrayList<Kmer> fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound-k+1-leftKmers.size()-rightKmers.size(), lookahead, minOverlap, maxCovGradient, true);

                    ArrayDeque<ArrayList<Kmer>> segments = breakWithReadPairedKmers(fragmentKmers, graph);
                    
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
                                    fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength);
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
        private int polyXMinLen;
        private int polyXMaxMismatches;
        private boolean extendFragments;
        
        public FragmentAssembler(PairedReadSegments p,
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
                String left = getBestSegment(p.left, graph);
                
                if (left.length() < this.leftReadLengthThreshold) {
                    return;
                } 
                
                String right = getBestSegment(p.right, graph);
                
                if (right.length() < this.rightReadLengthThreshold) {
                    return;
                }

//                right = chompRightPolyX(right, polyXMinLen, polyXMaxMismatches);
//                
//                if (right.length() < this.rightReadLengthThreshold) {
//                    return;
//                }
                                
                ArrayList<Kmer> leftKmers = graph.getKmers(left);
                if (leftKmers.isEmpty()) {
                    return;
                }
                
                ArrayList<Kmer> rightKmers = graph.getKmers(right);
                if (rightKmers.isEmpty()) {
                    return;
                }
                                      
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

                ArrayList<Kmer> fragmentKmers = overlapAndConnect(leftKmers, rightKmers, graph, bound, lookahead, minOverlap, maxCovGradient, true);

                ArrayDeque<ArrayList<Kmer>> segments = breakWithReadPairedKmers(fragmentKmers, graph);

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
                            //if (extendFragments) {
                            if (extendFragments && minCov <= k*2 && isBranchFree(fragmentKmers, graph, maxTipLength)) {
                                fragmentKmers = naiveExtend(fragmentKmers, graph, maxTipLength);
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

                    if (hasComplexLeftKmer || hasComplexRightKmer) {
                        outList.put(new Fragment(graph.assemble(leftKmers), graph.assemble(rightKmers), null, 0, minCov, true));
                    }
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
                                                            extendFragments
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
            
            longFragmentLengthThreshold = fragStats[1];
        
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
    
    public int[] assembleFragmentsMultiThreaded(FastxFilePair[] fastqPairs, 
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
        
        if (graph.getReadKmerDistance() > 0) {
            System.out.println("Read paired-kmers Bloom filter FPR: " + graph.getRpkbf().getFPR() * 100 + " %");
        }
        
        System.out.println("Assembling fragments...");
        
        long fragmentId = 0;
        long unconnectedReadId = 0;
        long readPairsParsed = 0;
        
        int maxTasksQueueSize = numThreads;
        int maxConcurrentSubmissions = numThreads + maxTasksQueueSize;
        
        int newBound = bound;
        int[] fragLengthsStats = null;
        boolean pairedKmerDistanceIsSet = false;
        int shortestFragmentLengthAllowed = k;
                        
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
                int polyXMinLen = k;
                int polyXMaxMismatches = 1;
                
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
                                                            polyXMinLen,
                                                            polyXMaxMismatches,
                                                            extendFragments
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
                                ArrayList<Kmer> fragKmers = frag.kmers;

                                if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                    if (frag.minCov == 1) {
                                        graph.addPairedKmers(fragKmers);

                                        if (frag.length >= longFragmentLengthThreshold) {
                                            for (Kmer kmer : fragKmers) {
                                                screeningBf.add(kmer.getHash());
                                            }

                                            longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                        }
                                        else {
                                            shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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

                                                longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                            }
                                            else {
                                                shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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
                while (fxpr.hasNext()) {
                    p = fxpr.next();
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

                        if (fragments.remainingCapacity() <= maxConcurrentSubmissions) {

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
                                        ArrayList<Kmer> fragKmers = frag.kmers;

                                        if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                            if (frag.minCov == 1) {
                                                graph.addPairedKmers(fragKmers);

                                                if (frag.length >= longFragmentLengthThreshold) {
                                                    for (Kmer kmer : fragKmers) {
                                                        screeningBf.add(kmer.getHash());
                                                    } 

                                                    longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                                }
                                                else {
                                                    shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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

                                                        longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                                    }
                                                    else {
                                                        shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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
                            ArrayList<Kmer> fragKmers = frag.kmers;

                            if (!containsAllKmers(screeningBf, fragKmers) || !graph.containsAllPairedKmers(fragKmers)) {
                                if (frag.minCov == 1) {
                                    graph.addPairedKmers(fragKmers);

                                    if (frag.length >= longFragmentLengthThreshold) {
                                        for (Kmer kmer : fragKmers) {
                                            screeningBf.add(kmer.getHash());
                                        }

                                        longSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                    }
                                    else {
                                        shortSingletonsOut.write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
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

                                            longFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                        }
                                        else {
                                            shortFragmentsOut[m].write(Long.toString(++fragmentId) + " L=[" + frag.left + "] R=[" + frag.right + "]", graph.assemble(frag.kmers));
                                        }
                                    }
                                }
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
                                                    boolean extendBranchFreeFragmentsOnly) throws InterruptedException, IOException, Exception {
        
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
                                                boolean sensitiveMode,
                                                String txptNamePrefix) {
        
        long numFragmentsParsed = 0;

        try {

            if (minPolyATailLengthRequired > 0) {
                System.out.println("Assembling polyadenylated transcripts only...");
            }
            else {
                System.out.println("Assembling transcripts...");
            }
        

            FastaWriter fout = new FastaWriter(outFasta, false);
            FastaWriter foutShort = new FastaWriter(outFastaShort, false);
            TranscriptWriter writer = new TranscriptWriter(fout, foutShort, minTranscriptLength);

   
            boolean allowNaiveExtension = true;
            boolean extendBranchFreeOnly = false;
            
            // extend LONG fragments
            
            for (int mag=longFragmentsFastas.length-1; mag>0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".L.");
                String fragmentsFasta = longFragmentsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly);
            }          

            // extend SHORT fragments
            
            for (int mag=shortFragmentsFastas.length-1; mag>0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".S.");
                String fragmentsFasta = shortFragmentsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly);
            }
            
            // extend UNCONNECTED reads
            
            for (int mag=unconnectedReadsFastas.length-1; mag>0; --mag) {
                writer.setOutputPrefix(txptNamePrefix + "E" + mag + ".U.");
                String fragmentsFasta = unconnectedReadsFastas[mag];
                System.out.println("Parsing `" + fragmentsFasta + "`...");
                numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                            allowNaiveExtension, extendBranchFreeOnly);
            }
            
            if (sensitiveMode) {
                System.out.println("Sensitive assembly mode is ON...");
            }
            else {
                // be extra careful with extending low coverage fragments (ie. 01, E0)
                allowNaiveExtension = false;
                extendBranchFreeOnly = true;
            }
            
            
            // extend LONG fragments
            
            writer.setOutputPrefix(txptNamePrefix + "E0.L.");
            String fragmentsFasta = longFragmentsFastas[0];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);

            // extend SHORT fragments
            
            writer.setOutputPrefix(txptNamePrefix + "E0.S.");
            fragmentsFasta = shortFragmentsFastas[0];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);
            
            // extend UNCONNECTED reads

            writer.setOutputPrefix(txptNamePrefix + "E0.U.");
            fragmentsFasta = unconnectedReadsFastas[0];
            System.out.println("Parsing `" + fragmentsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(fragmentsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);
            
            
            // extend LONG singleton fragments

            writer.setOutputPrefix(txptNamePrefix + "01.L.");
            System.out.println("Parsing `" + longSingletonsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(longSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);

            // extend SHORT singleton fragments

            writer.setOutputPrefix(txptNamePrefix + "01.S.");
            System.out.println("Parsing `" + shortSingletonsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(shortSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);

            // extend UNCONNECTED reads

            writer.setOutputPrefix(txptNamePrefix + "01.U.");
            System.out.println("Parsing `" + unconnectedSingletonsFasta + "`...");
            numFragmentsParsed += extendFragmentsMultiThreadedHelper(unconnectedSingletonsFasta, writer, sampleSize, numThreads,
                                                                        allowNaiveExtension, extendBranchFreeOnly);
            
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
                "RNA-Bloom version " + VERSION + "\n" +
                "Ka Ming Nip, Canada's Michael Smith Genome Sciences Centre\n" +
                "Copyright 2017"
        );
        
        if (exit) {
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
//        final String FRAGMENTS_K2_DONE = "FRAGMENTS.K2.DONE";
//        final String TRANSFRAGS_DONE = "TRANSFRAGS.DONE";
        final String TRANSCRIPTS_DONE = "TRANSCRIPTS.DONE";
        
        MyTimer timer = new MyTimer();
        
        System.out.println("args: " + Arrays.toString(args));
        
        // -left /home/gengar/test_data/GAPDH/GAPDH_2.fq.gz -right /home/gengar/test_data/GAPDH/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /home/gengar/test_assemblies/GAPDH
        // -left /home/gengar/test_data/SRR1360926/SRR1360926_2.fastq.gz -right /home/gengar/test_data/SRR1360926/SRR1360926_1.fastq.gz -revcomp-right -stranded -name SRR1360926 -outdir /home/gengar/test_assemblies/SRR1360926
        // -left /home/gengar/test_data/SRR1360926/SRR1360926.RNF213.2.fq.gz -right /home/gengar/test_data/SRR1360926/SRR1360926.RNF213.1.fq.gz -revcomp-right -stranded -name RNF213 -outdir /home/gengar/test_assemblies/RNF213
        
        // -t 1 -pair 10 -mem 2 -left /projects/btl2/kmnip/rna-bloom/example/SRX983106/SRR1957705_2.fastq.gz.trim.fq.gz -right /projects/btl2/kmnip/rna-bloom/example/SRX983106/SRR1957705_1.fastq.gz.trim.fq.gz -revcomp-right -stranded -name SRR953119 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/c.elegans/SRR953119/apr07.debug
        // -left /projects/btl2/kmnip/rna-bloom/tests/GAPDH_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/tests/GAPDH_1.fq.gz -revcomp-right -stranded -name gapdh -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/gapdh
        // -left /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_2.fq.gz -right /projects/btl2/kmnip/rna-bloom/example/SRP043027/trimmed_mod_1.fq.gz -revcomp-right -stranded -name SRR1360926 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/SRR1360926
        // -left /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_2.fq.gz -right /projects/btl2/kmnip/ENCODE/MCF-7_nucleus_all_1.fq.gz -revcomp-right -stranded -name mcf7 -outdir /projects/btl2/kmnip/rna-bloom/tests/java_assemblies/mcf7        
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
        
//        Option optKmerSize2 = Option.builder("k2")
//                                    .longOpt("kmer2")
//                                    .desc("2nd kmer size")
//                                    .hasArg(true)
//                                    .argName("INT")
//                                    .build();
//        options.addOption(optKmerSize2);
        
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

        Option optSensitive = Option.builder("sensitive")
                                    .longOpt("sensitive")
                                    .desc("assemble transcripts in sensitive mode")
                                    .hasArg(false)
                                    .build();
        options.addOption(optSensitive);
        
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
        
        Option optPolyATail = Option.builder("a")
                                    .longOpt("polya")
                                    .desc("only assemble transcripts with poly-A tails of the minimum length specified")
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
            
            final int numThreads = Integer.parseInt(line.getOptionValue(optThreads.getOpt(), "2"));
            final boolean forceOverwrite = line.hasOption(optForce.getOpt());
            
            final String name = line.getOptionValue(optName.getOpt(), "rnabloom");
            final String outdir = line.getOptionValue(optOutdir.getOpt(), System.getProperty("user.dir") + File.separator + name + "_assembly");
            
            final String longFragmentsFastaPrefix =      outdir + File.separator + name + ".fragments.long.";
            final String shortFragmentsFastaPrefix =     outdir + File.separator + name + ".fragments.short.";
            final String unconnectedReadsFastaPrefix =   outdir + File.separator + name + ".unconnected.";
//            final String unconnectedK2ReadsFastaPrefix = outdir + File.separator + name + ".unconnected.k2.";
            final String transcriptsFasta =              outdir + File.separator + name + ".transcripts.fa";
            final String shortTranscriptsFasta =         outdir + File.separator + name + ".transcripts.short.fa";
            final String graphFile = outdir + File.separator + name + ".graph";
//            final String graphK2File = outdir + File.separator + name + ".k2.graph";
            final String fragStatsFile = outdir + File.separator + name + ".fragstats";
            
            File startedStamp = new File(outdir + File.separator + STARTED);
            File dbgDoneStamp = new File(outdir + File.separator + DBG_DONE);
            File fragsDoneStamp = new File(outdir + File.separator + FRAGMENTS_DONE);
//            File fragsK2DoneStamp = new File(outdir + File.separator + FRAGMENTS_K2_DONE);
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
                
//                if (fragsK2DoneStamp.exists()) {
//                    fragsK2DoneStamp.delete();
//                }
                
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
            
            final boolean revCompLeft = line.hasOption(optRevCompLeft.getOpt());
            final boolean revCompRight = line.hasOption(optRevCompRight.getOpt());
            final boolean strandSpecific = line.hasOption(optStranded.getOpt());
            
            final int k = Integer.parseInt(line.getOptionValue(optKmerSize.getOpt(), "25"));
//            final int k2 = Integer.parseInt(line.getOptionValue(optKmerSize2.getOpt(), Integer.toString(k)));
            
            final int qDBG = Integer.parseInt(line.getOptionValue(optBaseQualDbg.getOpt(), "3"));
            final int qFrag = Integer.parseInt(line.getOptionValue(optBaseQualFrag.getOpt(), "3"));
            
            double leftReadFilesTotalBytes = 0;
            for (String fq : fastqsLeft) {
                leftReadFilesTotalBytes += new File(fq).length();
            }
            double rightReadFilesTotalBytes = 0;
            for (String fq : fastqsRight) {
                rightReadFilesTotalBytes += new File(fq).length();
            }
            
            final float maxBfMem = (float) Float.parseFloat(line.getOptionValue(optAllMem.getOpt(), Float.toString((float) (Math.max(NUM_BYTES_1MB * 100, 0.75f * Math.min(leftReadFilesTotalBytes, rightReadFilesTotalBytes)) / NUM_BYTES_1GB))));
            final float sbfGB = Float.parseFloat(line.getOptionValue(optSbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 8f)));
            final float dbgGB = Float.parseFloat(line.getOptionValue(optDbgbfMem.getOpt(), Float.toString(maxBfMem * 1f / 8f)));
            final float cbfGB = Float.parseFloat(line.getOptionValue(optCbfMem.getOpt(), Float.toString(maxBfMem * 6f / 8f)));
            final float pkbfGB = Float.parseFloat(line.getOptionValue(optPkbfMem.getOpt(), Float.toString(maxBfMem * 0.5f / 8f)));
            
            final long sbfSize = (long) (NUM_BITS_1GB * sbfGB);
            final long dbgbfSize = (long) (NUM_BITS_1GB * dbgGB);
            final long cbfSize = (long) (NUM_BYTES_1GB * cbfGB);
            final long pkbfSize = (long) (NUM_BITS_1GB * pkbfGB);
            
            final int allNumHash = Integer.parseInt(line.getOptionValue(optAllHash.getOpt(), "2"));
            final String allNumHashStr = Integer.toString(allNumHash);
            final int sbfNumHash = Integer.parseInt(line.getOptionValue(optSbfHash.getOpt(), allNumHashStr));
            final int dbgbfNumHash = Integer.parseInt(line.getOptionValue(optDbgbfHash.getOpt(), allNumHashStr));
            final int cbfNumHash = Integer.parseInt(line.getOptionValue(optCbfHash.getOpt(), allNumHashStr));
            final int pkbfNumHash = Integer.parseInt(line.getOptionValue(optPkbfHash.getOpt(), allNumHashStr));
            
            /**@TODO ensure that sbfNumHash and pkbfNumHash <= max(dbgbfNumHash, cbfNumHash) */
                        
            final int minOverlap = Integer.parseInt(line.getOptionValue(optOverlap.getOpt(), "10"));
            final int sampleSize = Integer.parseInt(line.getOptionValue(optSample.getOpt(), "1000"));
            final int bound = Integer.parseInt(line.getOptionValue(optBound.getOpt(), "500"));
            final int lookahead = Integer.parseInt(line.getOptionValue(optLookahead.getOpt(), "3"));
            final int maxTipLen = Integer.parseInt(line.getOptionValue(optTipLength.getOpt(), Integer.toString(k)));
            final float maxCovGradient = Float.parseFloat(line.getOptionValue(optMaxCovGrad.getOpt(), "0.5"));
            final float percentIdentity = Float.parseFloat(line.getOptionValue(optPercentIdentity.getOpt(), "0.95"));
            final int maxIndelSize = Integer.parseInt(line.getOptionValue(optIndelSize.getOpt(), "1"));
            final int maxErrCorrItr = Integer.parseInt(line.getOptionValue(optErrCorrItr.getOpt(), "1"));
            final int minTranscriptLength = Integer.parseInt(line.getOptionValue(optMinLength.getOpt(), "200"));
            final int minPolyATail = Integer.parseInt(line.getOptionValue(optPolyATail.getOpt(), "0"));
            final boolean sensitiveMode = line.hasOption(optSensitive.getOpt());
            final boolean extendFragments = line.hasOption(optExtend.getOpt());
            final int minNumKmerPairs = Integer.parseInt(line.getOptionValue(optMinKmerPairs.getOpt(), "10"));
            final String txptNamePrefix = line.getOptionValue(optPrefix.getOpt(), "");
            

            System.out.println("Bloom filters      Memory (GB)");
            System.out.println("==============================");
            System.out.println("de Bruijn graph:   " + dbgGB);
            System.out.println("kmer counting:     " + cbfGB);
            System.out.println("paired kmers (SE): " + pkbfGB);
            System.out.println("paired kmers (PE): " + pkbfGB);
            System.out.println("screening:         " + sbfGB);
            System.out.println("==============================");
            System.out.println("Total:             " + (dbgGB+cbfGB+2*pkbfGB+sbfGB));
            
            System.out.println("name:    " + name);
            System.out.println("outdir:  " + outdir);
            
            File f = new File(outdir);
            if (!f.exists()) {
                System.out.println("WARNING: Output directory does not exist!");
                f.mkdirs();
                System.out.println("Created output directory at `" + outdir + "`");
            }

            RNABloom assembler = new RNABloom(k, qDBG, qFrag);
            assembler.setParams(strandSpecific, maxTipLen, lookahead, maxCovGradient, maxIndelSize, percentIdentity, minNumKmerPairs, minPolyATail);

            try {
                touch(startedStamp);
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(1);
            }
            
            if (!forceOverwrite && dbgDoneStamp.exists()) {
                System.out.println("WARNING: Graph was already generated (k=" + k + ")!");
                
                if (!fragsDoneStamp.exists() || !txptsDoneStamp.exists()) {
                    System.out.println("Loading graph from file `" + graphFile + "`...");
                    assembler.restoreGraph(new File(graphFile));
                }
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
                       
                System.out.println("Building graph from reads (k=" + k + ")...");
                timer.start();
                
                assembler.initializeGraph(strandSpecific, 
                        dbgbfSize, cbfSize, pkbfSize, 
                        dbgbfNumHash, cbfNumHash, pkbfNumHash, false, true);
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads, false, true);
                
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
                        

            FastxFilePair[] fqPairs = new FastxFilePair[fastqsLeft.length];
            for (int i=0; i<fastqsLeft.length; ++i) {
                fqPairs[i] = new FastxFilePair(fastqsLeft[i], fastqsRight[i], revCompLeft, revCompRight);
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
                
                timer.start();
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                assembler.setupFragmentPairedKmersBloomFilter(pkbfSize, pkbfNumHash);
                
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
                
                assembler.savePairedKmersBloomFilter(new File(graphFile));
                assembler.writeFragStatsToFile(fragStats, fragStatsFile);
                
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                                
                try {
                    touch(fragsDoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
            else {
                System.out.println("WARNING: Fragments were already assembled (k=" + k + ")!");
            }

            // Connect unconnected reads by using the 2nd kmer size
            /*
            assembler.setK(k2);
            
            String[] unconnectedK2ReadsFastaPaths = {unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[0] + ".fa",
                                            unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[1] + ".fa",
                                            unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[2] + ".fa",
                                            unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[3] + ".fa",
                                            unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[4] + ".fa",
                                            unconnectedK2ReadsFastaPrefix + COVERAGE_ORDER[5] + ".fa"};

            String unconnectedK2SingletonsFastaPath = unconnectedK2ReadsFastaPrefix + "01.fa";

            if (!forceOverwrite && fragsK2DoneStamp.exists()) {
                System.out.println("WARNING: Read pairs were already rescued (k=" + k2 + ")!");
                
                if (!txptsDoneStamp.exists()) {
                    System.out.println("Loading graph from file `" + graphK2File + "`...");
                    File tmp = new File(graphK2File);
                    assembler.restoreGraph(tmp);
                    assembler.restorePairedKmersBloomFilter(tmp);
                }
            }
            else {
                boolean newKmerSize = k!=k2;
                
                if (assembler.isGraphInitialized()) {
                    assembler.clearDbgBf();
                    if (newKmerSize) {
                        assembler.clearCBf();
                        assembler.clearPkBf();
                        assembler.clearRpkBf();
                    }
                }
                else {
                    assembler.initializeGraph(strandSpecific, 
                            dbgbfSize, cbfSize, pkbfSize, 
                            dbgbfNumHash, cbfNumHash, pkbfNumHash, true, true);
                }
                
                assembler.setupKmerScreeningBloomFilter(sbfSize, sbfNumHash);
                
                // adjust paired kmer distance
                int[] fragStats = assembler.restoreFragStatsFromFile(fragStatsFile);
                assembler.setPairedKmerDistance(fragStats[1]);
                int fragSizeBound = fragStats[3] + 3*(fragStats[3] - fragStats[1])/2;

                // repopulate with NEW kmers from fragments
                ArrayList<String> fragmentPaths = new ArrayList<>(longFragmentsFastaPaths.length + shortFragmentsFastaPaths.length + unconnectedReadsFastaPaths.length + 3);
                fragmentPaths.addAll(Arrays.asList(longFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(shortFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(unconnectedReadsFastaPaths));
                fragmentPaths.add(longSingletonsFastaPath);
                fragmentPaths.add(shortSingletonsFastaPath);
                fragmentPaths.add(unconnectedSingletonsFastaPath);

                System.out.println("Rebuilding graph from assembled fragments (k=" + k2 + ")...");
                timer.start();
                assembler.repopulateGraph(fragmentPaths, strandSpecific, newKmerSize);
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));    
                
                if (newKmerSize) {
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

                    System.out.println("Counting kmers in reads (k=" + k2 + ")...");
                    timer.start();
                    assembler.populateGraph(forwardFilesList, backwardFilesList, strandSpecific, numThreads, true, true);
                    System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                }
                
                // Remove existing output files
                
                File fragmentsFile = new File(unconnectedK2SingletonsFastaPath);
                if (fragmentsFile.exists()) {
                    fragmentsFile.delete();
                }

                for (String fragmentsFasta : unconnectedK2ReadsFastaPaths) {
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
                
                timer.start();
                assembler.rescueUnconnectedMultiThreaded(allUnconnectedReads, 
                                                    longFragmentsFastaPaths,
                                                    shortFragmentsFastaPaths,
                                                    unconnectedK2ReadsFastaPaths,
                                                    longSingletonsFastaPath,
                                                    shortSingletonsFastaPath,
                                                    unconnectedK2SingletonsFastaPath,
                                                    fragSizeBound,
                                                    minOverlap,
                                                    sampleSize, 
                                                    numThreads, 
                                                    maxErrCorrItr,
                                                    extendFragments);
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));
                
                // Save DBG-Bf and pk-Bf to disk
                System.out.println("Saving graph to file `" + graphK2File + "`...");
                assembler.saveGraph(new File(graphK2File));
                assembler.savePairedKmersBloomFilter(new File(graphK2File));
                
                // Touch stamp
                try {
                    touch(fragsK2DoneStamp);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
            */
            
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
                if (assembler.isGraphInitialized()) {
                    assembler.clearDbgBf();
                }
//                else {
//                    assembler.initializeGraph(strandSpecific, 
//                            dbgbfSize, cbfSize, pkbfSize, 
//                            dbgbfNumHash, cbfNumHash, pkbfNumHash, true, true);
//                }
                
                // adjust paired kmer distance
                int[] fragStats = assembler.restoreFragStatsFromFile(fragStatsFile);

                // repopulate with kmers from fragments
                ArrayList<String> fragmentPaths = new ArrayList<>(longFragmentsFastaPaths.length + shortFragmentsFastaPaths.length + unconnectedReadsFastaPaths.length + 3);
                fragmentPaths.addAll(Arrays.asList(longFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(shortFragmentsFastaPaths));
                fragmentPaths.addAll(Arrays.asList(unconnectedReadsFastaPaths));
                fragmentPaths.add(longSingletonsFastaPath);
                fragmentPaths.add(shortSingletonsFastaPath);
                fragmentPaths.add(unconnectedSingletonsFastaPath);

                System.out.println("Rebuilding graph from assembled fragments (k=" + k + ")...");
                timer.start();
                assembler.populateGraphFromFragments(fragmentPaths, strandSpecific, false);
                System.out.println("Time elapsed: " + MyTimer.hmsFormat(timer.elapsedMillis()));  
                
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
                                                            unconnectedReadsFastaPaths,
                                                            longSingletonsFastaPath,
                                                            shortSingletonsFastaPath,
                                                            unconnectedSingletonsFastaPath,
                                                            transcriptsFasta, 
                                                            shortTranscriptsFasta,
                                                            graphFile,
                                                            numThreads,
                                                            sampleSize,
                                                            minTranscriptLength,
                                                            sensitiveMode,
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
            else {
                System.out.println("WARNING: Transcripts were already assembled!");
            }
            
        }
        catch (ParseException exp) {
            System.out.println("ERROR:" + exp.getMessage() );
            System.exit(1);
        }
        
        System.out.println("Total Runtime: " + MyTimer.hmsFormat(timer.totalElapsedMillis()));
    }
}
