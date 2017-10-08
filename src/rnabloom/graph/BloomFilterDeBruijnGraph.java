/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.graph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.PairedKeysPartitionedBloomFilter;
import rnabloom.bloom.hash.*;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author kmnip
 */
public class BloomFilterDeBruijnGraph {
    
    private BloomFilter dbgbf;
    private CountingBloomFilter cbf;
    private PairedKeysPartitionedBloomFilter pkbf = null;
    
    private int dbgbfNumHash;
    private int cbfNumHash;
    private int dbgbfCbfMaxNumHash;
    private final HashFunction2 hashFunction;
    private int k;
    private int kMinus1;
    private boolean stranded;
    private int pairedKmersDistance = -1;
    private long pkbfNumBits;
    private int pkbfNumHash;
//    private KmerBitsUtils2 bitsUtils;
    
    private final static String FILE_DESC_EXTENSION = ".desc";
    private final static String FILE_DBGBF_EXTENSION = ".dbgbf";
    private final static String FILE_CBF_EXTENSION = ".cbf";
    private final static String FILE_PKBF_LEFT_EXTENSION = ".lkbf";
    private final static String FILE_PKBF_RIGHT_EXTENSION = ".rkbf";
    private final static String FILE_PKBF_PAIR_EXTENSION = ".pkbf";

    
    private final static String LABEL_SEPARATOR = ":";
    private final static String LABEL_DBGBF_CBF_NUM_HASH = "dbgbfCbfMaxNumHash";
    private final static String LABEL_K = "k";
    private final static String LABEL_STRANDED = "stranded";
//    private final static String LABEL_SEED = "seed";
    private final static String LABEL_PAIRED_KMER_DIST = "pairedKmersDistance";
    private final static String LABEL_PKBF_NUM_BITS = "pkbfNumBits";
    private final static String LABEL_PKBF_NUM_HASH = "pkbfNumHash";
    
    public BloomFilterDeBruijnGraph(long dbgbfNumBits,
                                    long cbfNumBytes,
                                    long pkbfNumBits,
                                    int dbgbfNumHash,
                                    int cbfNumHash,
                                    int pkbfNumHash,
                                    int k,
                                    boolean stranded) {
        this.k = k;
        this.kMinus1 = k-1;
//        this.bitsUtils = new KmerBitsUtils2(k);
        this.stranded = stranded;
        this.dbgbfNumHash = dbgbfNumHash;
        this.cbfNumHash = cbfNumHash;
        this.dbgbfCbfMaxNumHash = Math.max(dbgbfNumHash, cbfNumHash);
        if (stranded) {
            this.hashFunction = new HashFunction2(k);
        }
        else {
            this.hashFunction = new CanonicalHashFunction2(k);
        }
        this.dbgbf = new BloomFilter(dbgbfNumBits, dbgbfNumHash, this.hashFunction);
        this.cbf = new CountingBloomFilter(cbfNumBytes, cbfNumHash, this.hashFunction);
        this.pkbfNumBits = pkbfNumBits;
        this.pkbfNumHash = pkbfNumHash;
    }
    
    public BloomFilterDeBruijnGraph(File graphFile) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(graphFile));
        String line;
        while ((line = br.readLine()) != null) {
            String[] entry = line.split(LABEL_SEPARATOR);
            String key = entry[0];
            String val = entry[1];
            switch(key) {
                case LABEL_DBGBF_CBF_NUM_HASH:
                    dbgbfCbfMaxNumHash = Integer.parseInt(val);
                    break;
                case LABEL_K:
                    k = Integer.parseInt(val);
                    kMinus1 = k-1;
//                    bitsUtils = new KmerBitsUtils2(k);
                    break;
                case LABEL_STRANDED:
                    stranded = Boolean.parseBoolean(val);
                    break;
                case LABEL_PAIRED_KMER_DIST:
                    pairedKmersDistance = Integer.parseInt(val);
                    break;
                case LABEL_PKBF_NUM_BITS:
                    pkbfNumBits = Long.parseLong(val);
                    break;
                case LABEL_PKBF_NUM_HASH:
                    pkbfNumHash = Integer.parseInt(val);
                    break;
            }
        }
        br.close();
        
        if (stranded) {
            this.hashFunction = new HashFunction2(k);
        }
        else {
            this.hashFunction = new CanonicalHashFunction2(k);
        }
        
        String dbgbfBitsPath = graphFile.getPath() + FILE_DBGBF_EXTENSION;
        String dbgbfDescPath = dbgbfBitsPath + FILE_DESC_EXTENSION;
        dbgbf = new BloomFilter(new File(dbgbfDescPath), new File(dbgbfBitsPath), hashFunction);
        
        String cbfBitsPath = graphFile.getPath() + FILE_CBF_EXTENSION;
        String cbfDescPath = cbfBitsPath + FILE_DESC_EXTENSION;
        cbf = new CountingBloomFilter(new File(cbfDescPath), new File(cbfBitsPath), hashFunction);
        
        dbgbfNumHash = dbgbf.getNumHash();
        cbfNumHash = cbf.getNumHash();
        
        String leftBitsPath = graphFile.getPath() + FILE_PKBF_LEFT_EXTENSION;
        String rightBitsPath = graphFile.getPath() + FILE_PKBF_RIGHT_EXTENSION;
        String pairBitsPath = graphFile.getPath() + FILE_PKBF_PAIR_EXTENSION;
        String pkbfDescPath = pairBitsPath + FILE_DESC_EXTENSION;
        
        File leftBitsFile = new File(leftBitsPath);
        File rightBitsFile = new File(rightBitsPath);
        File pairBitsFile = new File(pairBitsPath);
        File pkbfDescFile = new File(pkbfDescPath);
        
        if (pkbfDescFile.isFile() && leftBitsFile.isFile() && rightBitsFile.isFile() && pairBitsFile.isFile()) {
            pkbf = new PairedKeysPartitionedBloomFilter(pkbfDescFile, leftBitsFile, rightBitsFile, pairBitsFile, hashFunction);
        }
    }

    public HashFunction2 getHashFunction() {
        return this.hashFunction;
    }
    
    public int getDbgbfNumHash () {
        return dbgbfNumHash;
    }

    public int getCbfNumHash () {
        return cbfNumHash;
    }
    
    public int getPkbfNumHash () {
        return pkbfNumHash;
    }
    
    public int getMaxNumHash() {
        return dbgbfCbfMaxNumHash;
    }
    
    public void destroy() {
        dbgbf.destroy();
        cbf.destroy();
        if (pkbf != null) {
            pkbf.destroy();
        }        
    }
    
    public void clearAllBf() {
        dbgbf.empty();
        cbf.empty();
        pkbf.empty();
    }
    
    public void clearDbgbf() {
        dbgbf.empty();
    }
    
    public void clearCbf() {
        cbf.empty();
    }
    
    public void clearPkbf() {
        pkbf.empty();
    }

    public void destroyDbgbf() {
        if (dbgbf != null) {
            dbgbf.destroy();
        }
    }
    
    public void destroyCbf() {
        if (cbf != null) {
            cbf.destroy();
        }
    }
    
    public void destroyPkbf() {
        if (pkbf != null) {
            pkbf.destroy();
        }
    }
    
    public BloomFilter getDbgbf() {
        return dbgbf;
    }

    public CountingBloomFilter getCbf() {
        return cbf;
    }

    public PairedKeysPartitionedBloomFilter getPkbf() {
        return pkbf;
    }
    
    public boolean isStranded() {
        return stranded;
    }    
    
    private void saveDesc(File graphFile) throws IOException {
        FileWriter writer = new FileWriter(graphFile);
        writer.write(LABEL_DBGBF_CBF_NUM_HASH + LABEL_SEPARATOR + dbgbfCbfMaxNumHash + "\n" +
                    LABEL_STRANDED + LABEL_SEPARATOR + stranded + "\n" +
                    LABEL_K + LABEL_SEPARATOR + k + "\n" +
                    LABEL_PAIRED_KMER_DIST + LABEL_SEPARATOR + pairedKmersDistance + "\n" +
                    LABEL_PKBF_NUM_BITS + LABEL_SEPARATOR + pkbfNumBits + "\n" +
                    LABEL_PKBF_NUM_HASH + LABEL_SEPARATOR + pkbfNumHash + "\n");
        writer.close();
    }
    
    public void save(File graphFile) throws IOException {
        /** write graph desc*/
        saveDesc(graphFile);
        
        /** write Bloom filters*/
        String dbgbfBitsPath = graphFile.getPath() + FILE_DBGBF_EXTENSION;
        String dbgbfDescPath = dbgbfBitsPath + FILE_DESC_EXTENSION;
        dbgbf.save(new File(dbgbfDescPath), new File(dbgbfBitsPath));
        
        String cbfBitsPath = graphFile.getPath() + FILE_CBF_EXTENSION;
        String cbfDescPath = cbfBitsPath + FILE_DESC_EXTENSION;
        cbf.save(new File(cbfDescPath), new File(cbfBitsPath));
        
        if (pkbf != null) {
            savePkbf(graphFile);
        }
    }
    
    public void savePkbf(File graphFile) throws IOException {
        /** update the graph desc because kmer pair distance is updated*/
        saveDesc(graphFile);
        
        String leftBitsPath = graphFile.getPath() + FILE_PKBF_LEFT_EXTENSION;
        String rightBitsPath = graphFile.getPath() + FILE_PKBF_RIGHT_EXTENSION;
        String pairBitsPath = graphFile.getPath() + FILE_PKBF_PAIR_EXTENSION;
        String pkbfDescPath = pairBitsPath + FILE_DESC_EXTENSION;
        
        pkbf.save(new File(pkbfDescPath), new File(leftBitsPath), new File(rightBitsPath), new File(pairBitsPath));
    }
    
    public void restorePkbf(File graphFile) throws IOException {
        String leftBitsPath = graphFile.getPath() + FILE_PKBF_LEFT_EXTENSION;
        String rightBitsPath = graphFile.getPath() + FILE_PKBF_RIGHT_EXTENSION;
        String pairBitsPath = graphFile.getPath() + FILE_PKBF_PAIR_EXTENSION;
        String pkbfDescPath = pairBitsPath + FILE_DESC_EXTENSION;
        
        pkbf = new PairedKeysPartitionedBloomFilter(new File(pkbfDescPath), new File(leftBitsPath), new File(rightBitsPath), new File(pairBitsPath), hashFunction);
    }
    
    public void initializePairKmersBloomFilter() {
        this.pkbf = new PairedKeysPartitionedBloomFilter(pkbfNumBits, pkbfNumHash, this.hashFunction);
    }
    
    public void setPairedKmerDistance(int d) {
        this.pairedKmersDistance = d;
    }
    
    public int getPairedKmerDistance() {
        return this.pairedKmersDistance;
    }
    
    public int getK() {
        return k;
    }
    
    public void setK(int k) {
        this.k = k;
        this.kMinus1 = k-1;
        this.hashFunction.setK(k);
    }
    
    public int getKMinus1() {
        return kMinus1;
    }
        
    public boolean isLowComplexity(Kmer2 kmer) {
        return isLowComplexity2(kmer.bytes);
    }
    
    public void add(String kmer) {
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        hashFunction.getHashValues(kmer, dbgbfCbfMaxNumHash, hashVals);
        add(hashVals);
    }
    
    public void add(final long[] hashVals) {
        dbgbf.add(hashVals);
        cbf.increment(hashVals);        
    }
    
    public void addIfAbsent(final long[] hashVals) {
        if (!dbgbf.lookup(hashVals) || cbf.getCount(hashVals) == 0) {
            add(hashVals);
        }
    }
    
    public void addDbgOnly(final long[] hashVals) {
        dbgbf.add(hashVals);
    }
    
    public void addCountOnly(final long[] hashVals) {
        cbf.increment(hashVals);
    }
    
//    public void addFragmentKmersFromSeq(String seq) {
//        final int numKmers = getNumKmers(seq, k);
//        for (int i=0; i<numKmers; ++i) {
//            pkbf.add(seq.substring(i, i+k));
//        }
//    }
//    
//    public void addFragmentKmers(ArrayList<Kmer> kmers) {
//        for (Kmer kmer : kmers) {
//            pkbf.add(kmer.hashVals);
//        }
//    }
    
    public void addPairedKmers(long[] leftHashVals, long[] rightHashVals, long[] pairingHashVals) {
        pkbf.add(leftHashVals, rightHashVals, pairingHashVals);
    }
    
    public void addPairedKmersFromSeq(String seq) {
        // add paired kmers
        final int upperBound = getNumKmers(seq, k) - pairedKmersDistance;
        
        if (upperBound > 0) {
            for (int i=0; i<upperBound; ++i) {
                pkbf.add(seq.substring(i, i+k), seq.substring(i+pairedKmersDistance, i+k+pairedKmersDistance));
            }
        }
    }
        
    public void addPairedKmers(ArrayList<Kmer2> kmers) {
        Kmer2 left, right;
        
        // add paired kmers
        final int upperBound = kmers.size() - pairedKmersDistance;
        
        if (upperBound > 0) {
            for (int i=0; i<upperBound; ++i) {
                left = kmers.get(i);
                right = kmers.get(i+pairedKmersDistance);

                pkbf.add(left.getHash(),
                        right.getHash(),
                        left.getKmerPairHashValue(right));
            }
        }
    }
     
    public boolean containsAllPairedKmers(ArrayList<Kmer2> kmers) {
        Kmer2 left, right;
        
        // add paired kmers
        final int upperBound = kmers.size() - pairedKmersDistance;
        
        if (upperBound <= 0) {
            return false;
        }
        
        for (int i=0; i<upperBound; ++i) {
            left = kmers.get(i);
            right = kmers.get(i+pairedKmersDistance);
            
            if (!pkbf.lookup(left.getHash(),
                            right.getHash(),
                            left.getKmerPairHashValue(right))) {
                return false;
            }
        }
        
        return true;
    }
    
//    public void addFragmentKmersAndPairedKmersFromSeq(String seq) {
//        final int numKmers = getNumKmers(seq, k);
//        
//        // add kmers
//        for (int i=0; i<numKmers; ++i) {
//            pkbf.add(seq.substring(i, i+k));
//        }
//
//        // add paired kmers
//        final int upperBound = numKmers-pairedKmersDistance;
//        for (int i=0; i<upperBound; ++i) {
//            pkbf.addPair(seq.substring(i, i+k), seq.substring(i+pairedKmersDistance, i+k+pairedKmersDistance));
//        }
//    }
    
//    public void addFragmentKmersAndPairedKmers(ArrayList<Kmer> kmers) {
//        addFragmentKmers(kmers);
//        addPairedKmers(kmers);
//    }
    
//    public boolean lookupFragmentKmer(final long[] hashVals) {
//        return pkbf.lookup(hashVals);
//    }
//    
//    public boolean lookupPairedKmers(String kmer1, String kmer2) {
//        return pkbf.lookup(kmer1, kmer2);
//    }
    
    public boolean lookupLeftKmer(final long hashVal) {
        return pkbf.lookupLeft(hashVal);
    }

    public boolean lookupRightKmer(final long hashVal) {
        return pkbf.lookupRight(hashVal);
    }
    
    public boolean lookupLeftKmer(final long[] hashVals) {
        return pkbf.lookupLeft(hashVals);
    }
    
    public boolean lookupRightKmer(final long[] hashVals) {
        return pkbf.lookupRight(hashVals);
    }
    
    public boolean lookupKmerPair(Kmer2 left, Kmer2 right) {
        return pkbf.lookup(left.getHash(),
                            right.getHash(),
                            left.getKmerPairHashValue(right));
    }
    
    public boolean lookupKmerPairing(Kmer2 left, Kmer2 right) {
        return pkbf.lookupPair(left.getKmerPairHashValue(right));
    }
    
    public boolean contains(String kmer) {
        return dbgbf.lookup(kmer);
    }
    
    public boolean contains(final long[] hashVals) {
        return dbgbf.lookup(hashVals);
    }

    public void increment(String kmer) {
        cbf.increment(kmer);
    }
    
    public float getCount(String kmer) {
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        hashFunction.getHashValues(kmer, dbgbfCbfMaxNumHash, hashVals);
        return getCount(hashVals);
    }
    
    public float getCount(final long[] hashVals) {
        if (dbgbf.lookup(hashVals)) {
            return cbf.getCount(hashVals);
        }
        else {
            return 0;
        }        
    }

    public float getDbgbfFPR() {
        return dbgbf.getFPR();
    }

    public float getCbfFPR() {
        return cbf.getFPR();
    }
    
    public float getPkbfFPR() {
        return pkbf.getFPR();
    }
    
    public float getFPR() {
        return dbgbf.getFPR() * cbf.getFPR();
    }
    
    public Kmer2 getKmer(String kmer) {
        return hashFunction.getKmer(kmer, dbgbfCbfMaxNumHash, this);
    }
    
//    public class Kmer {
//        public byte[] bytes;
//        public float count;
//        public long fHashVal;
//        public long rHashVal;
//        
////        public long[] hashVals;
////        public ArrayDeque<Kmer> predecessors = null;
////        public ArrayDeque<Kmer> successors = null;
//                
//        public Kmer(String seq, float count, long fHashVal, long rHashVal) {
//            this.bytes = stringToBytes(seq, k);
//            this.count = count;
//            this.fHashVal = fHashVal;
//            this.rHashVal = rHashVal;
//        }
//        
//        public Kmer(byte[] bytes, float count, long fHashVal, long rHashVal) {
//            this.bytes = bytes;
//            this.count = count;
//            this.fHashVal = fHashVal;
//            this.rHashVal = rHashVal;
//        }
//        
//        public boolean equals(Kmer other) {
//            return Arrays.equals(bytes, other.bytes);
//        }
//        
//        @Override
//        public String toString() {
//            return bytesToString(bytes, bytes.length);
//        }
//        
//        @Override
//        public int hashCode(){
//            return (int) Math.min(fHashVal, rHashVal); 
//        }
//
//        @Override
//        public boolean equals(Object obj) {
//            if (this == obj) {
//                return true;
//            }
//            if (obj == null) {
//                return false;
//            }
//            if (getClass() != obj.getClass()) {
//                return false;
//            }
//            return Arrays.equals(bytes, ((Kmer)obj).bytes);
//        }
//    }
    
    public String getPrefix(String kmer) {
        return kmer.substring(0, kMinus1);
    }
    
    public String getSuffix(String kmer) {
        return kmer.substring(1, k);
    }
    
    public CharSequence getPrefixCharSeq(String kmer) {
        return kmer.subSequence(0, kMinus1);
    }
    
    public CharSequence getSuffixCharSeq(String kmer) {
        return kmer.subSequence(1, k);
    }
    
//    public boolean hasPredecessors(Kmer kmer) {
//        PredecessorsNTHashIterator itr = hashFunction.getPredecessorsNTHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.fHashVal, (char) kmer.bytes[kMinus1]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                return true;
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getPredecessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[kMinus1]);
////        for (int i=0; i<4; ++i) {
////            if (dbgbf.lookup(allHashVals[i]) && cbf.getCount(allHashVals[i]) > 0) {
////                return true;
////            }
////        }
//        return false;
//    }
//
//    public boolean hasAtLeastXPredecessors(Kmer kmer, int x) {
//        int result = 0;
//        
//        PredecessorsNTHashIterator itr = hashFunction.getPredecessorsNTHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.fHashVal, (char) kmer.bytes[kMinus1]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                if (++result >= x) {
//                    return true;
//                }
//            }
//        }
//        
//        return false;
//    }
//    
//    public int getNumPredecessors(Kmer kmer) {
//        int result = 0;
//        
//        PredecessorsNTHashIterator itr = hashFunction.getPredecessorsNTHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.fHashVal, (char) kmer.bytes[kMinus1]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                ++result;
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getPredecessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[kMinus1]);
////        
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (dbgbf.lookup(hashVals) && cbf.getCount(hashVals) > 0) {
////                ++result;
////            }
////        }
//        
//        return result;
//    }
//    
//    public ArrayDeque<Kmer> getPredecessors(Kmer kmer) {
//        
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//  
//        PredecessorsNTHashIterator itr = hashFunction.getPredecessorsNTHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.fHashVal, (char) kmer.bytes[kMinus1]);
//        long[] hVals = itr.hVals;
//        
//        float count;
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals)) {
//                count = cbf.getCount(hVals);
//                if (count > 0) {
//                    byte[] bytes = shiftRight(kmer.bytes, k);
//                    bytes[0] = (byte) c;
//                    result.add(new Kmer(bytes, count, hVals));
//                }
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getPredecessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[kMinus1]);
////        
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    byte[] bytes = shiftRight(kmer.bytes, k);
////                    bytes[0] = NUCLEOTIDE_BYTES[i];
////                    result.add(new Kmer(bytes, count, hashVals));
////                }
////            }
////        }
//        
////        String prefix = getPrefix(kmer.seq);
////        long[] hashVals = new long[dbgbfCbfMaxNumHash];
////        for (char c : NUCLEOTIDES) {
////            String v = c + prefix;
////            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
////            if (dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    result.add(new Kmer(v, count, hashVals));
////                }
////            }
////        }
//        
////        kmer.predecessors = result;
//        
//        return result;
//    }
//    
//    public ArrayDeque<Kmer> getPredecessors(Kmer kmer, BloomFilter bf) {        
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//                       
//        PredecessorsNTHashIterator itr = hashFunction.getPredecessorsNTHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[kMinus1]);
//        long[] hVals = itr.hVals;
//        
//        float count;
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (bf.lookup(hVals) && dbgbf.lookup(hVals)) {
//                count = cbf.getCount(hVals);
//                if (count > 0) {
//                    byte[] bytes = shiftRight(kmer.bytes, k);
//                    bytes[0] = (byte) c;
//                    result.add(new Kmer(bytes, count, hVals));
//                }
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getPredecessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[kMinus1]);
////        
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (bf.lookup(hashVals) && dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    byte[] bytes = shiftRight(kmer.bytes, k);
////                    bytes[0] = NUCLEOTIDE_BYTES[i];
////                    result.add(new Kmer(bytes, count, hashVals));
////                }
////            }
////        }
//        
//        return result;
//    }
//    
//    public ArrayDeque<String> getPredecessors(String kmer) {
//        ArrayDeque<String> result = new ArrayDeque<>(4);
//        final String prefix = getPrefix(kmer);
//        long[] hashVals = new long[dbgbfCbfMaxNumHash];
//        for (char c : NUCLEOTIDES) {
//            String v = c + prefix;
//            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
//            if (dbgbf.lookup(hashVals)) {
//                float count = cbf.getCount(hashVals);
//                if (count > 0) {
//                    result.add(v);
//                }
//            }
//        }
//        return result;
//    }
//
//    public boolean hasSuccessors(Kmer kmer) {
//        SuccessorsNTHashIterator itr = hashFunction.getSuccessorsHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[0]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                return true;
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getSuccessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[0]);
////        for (int i=0; i<4; ++i) {
////            if (dbgbf.lookup(allHashVals[i]) && cbf.getCount(allHashVals[i]) > 0) {
////                return true;
////            }
////        }
//
//        return false;
//    }
//    
//    public boolean hasAtLeastXSuccessors(Kmer kmer, int x) {
//        int result = 0;
//
//        SuccessorsNTHashIterator itr = hashFunction.getSuccessorsHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[0]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                if (++result >= x) {
//                    return true;
//                }
//            }
//        }
//        
//        return false;
//    }
//    
//    public int getNumSuccessors(Kmer kmer) {
//        int result = 0;
//
//        SuccessorsNTHashIterator itr = hashFunction.getSuccessorsHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[0]);
//        long[] hVals = itr.hVals;
//        
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals) && cbf.getCount(hVals) > 0) {
//                ++result;
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getSuccessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[0]);
////                
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (dbgbf.lookup(hashVals) && cbf.getCount(hashVals) > 0) {
////                ++result;
////            }
////        }
//        
//        return result;
//    }
//    
//    public ArrayDeque<Kmer> getSuccessors(Kmer kmer) {        
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//        
//        SuccessorsNTHashIterator itr = hashFunction.getSuccessorsHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[0]);
//        long[] hVals = itr.hVals;
//        
//        float count;
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (dbgbf.lookup(hVals)) {
//                count = cbf.getCount(hVals);
//                if (count > 0) {
//                    byte[] bytes = shiftLeft(kmer.bytes, k);
//                    bytes[kMinus1] = (byte) c;
//                    result.add(new Kmer(bytes, count, hVals));
//                }
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getSuccessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[0]);
////                
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    byte[] bytes = shiftLeft(kmer.bytes, k);
////                    bytes[kMinus1] = NUCLEOTIDE_BYTES[i];
////                    result.add(new Kmer(bytes, count, hashVals));
////                }
////            }
////        }
//        
////        final String suffix = getSuffix(kmer.seq);
////        long[] hashVals = new long[dbgbfCbfMaxNumHash];
////        for (char c : NUCLEOTIDES) {
////            String v = suffix + c;
////            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
////            if (dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    result.add(new Kmer(v, count, hashVals));
////                }
////            }
////        }
//        
////        kmer.successors = result;
//        
//        return result;
//    }
//    
//    public ArrayDeque<Kmer> getSuccessors(Kmer kmer, BloomFilter bf) {
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//
//        SuccessorsNTHashIterator itr = hashFunction.getSuccessorsHashIterator(dbgbfCbfMaxNumHash);
//        itr.start(kmer.hashVals, (char) kmer.bytes[0]);
//        long[] hVals = itr.hVals;
//        
//        float count;
//        for (char c : NUCLEOTIDES) {
//            itr.next(c);
//            if (bf.lookup(hVals) && dbgbf.lookup(hVals)) {
//                count = cbf.getCount(hVals);
//                if (count > 0) {
//                    byte[] bytes = shiftLeft(kmer.bytes, k);
//                    bytes[kMinus1] = (byte) c;
//                    result.add(new Kmer(bytes, count, hVals));
//                }
//            }
//        }
//        
////        long[][] allHashVals = hashFunction.getSuccessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, (char) kmer.bytes[0]);
////        
////        long[] hashVals;
////        for (int i=0; i<4; ++i) {
////            hashVals = allHashVals[i];
////            if (bf.lookup(hashVals) && dbgbf.lookup(hashVals)) {
////                float count = cbf.getCount(hashVals);
////                if (count > 0) {
////                    byte[] bytes = shiftLeft(kmer.bytes, k);
////                    bytes[kMinus1] = NUCLEOTIDE_BYTES[i];
////                    result.add(new Kmer(bytes, count, hashVals));
////                }
////            }
////        }
//        
//        return result;
//    }
//    
//    public ArrayDeque<String> getSuccessors(String kmer) {
//        ArrayDeque<String> result = new ArrayDeque<>(4);
//        final String suffix = getSuffix(kmer);
//        long[] hashVals = new long[dbgbfCbfMaxNumHash];
//        for (char c : NUCLEOTIDES) {
//            String v = suffix + c;
//            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
//            if (dbgbf.lookup(hashVals)) {
//                float count = cbf.getCount(hashVals);
//                if (count > 0) {
//                    result.add(v);
//                }
//            }
//        }
//        return result;
//    }
//
//    public ArrayDeque<Kmer> getLeftVariants(Kmer kmer) {
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//        
//        LeftVariantsNTHashIterator itr = hashFunction.getLeftVariantsNTHashIterator(dbgbfCbfMaxNumHash);
//        char charOut = (char) kmer.bytes[0];
//        itr.start(kmer.hashVals, charOut);
//        long[] hVals = itr.hVals;
//        
//        byte[] bytes;
//        float count;
//        for (char charIn : getAltNucleotides(charOut)) {
//            itr.next(charIn);
//            count = getCount(hVals);
//            if (count > 0) {
//                bytes = Arrays.copyOf(kmer.bytes, k);
//                bytes[0] = (byte) charIn;
//                result.add(new Kmer(bytes, count, hVals));
//            }
//        }
//        
////        final String suffix = getSuffix(kmer.toString());
////        String v;
////        float count;
////        
////        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
////        for (char c : getAltNucleotides(kmer.bytes[0])) {
////            v = c + suffix;
////
////            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
////
////            count = getCount(hashVals);
////            if (count > 0) {
////                result.add(new Kmer(v, count, hashVals, k));
////            }
////        }
//                
//        return result;
//    }
    
    public ArrayDeque<String> getLeftVariants(String kmer) {
        ArrayDeque<String> result = new ArrayDeque<>(4);
        final String suffix = getSuffix(kmer);
        String v;
        for (char c : getAltNucleotides(kmer.charAt(0))) {
            v = c + suffix;

            if (contains(v)) {
                result.add(v);
            }
        }
        return result;
    }

//    public ArrayDeque<Kmer> getRightVariants(Kmer kmer) {
//        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
//        
//        RightVariantsNTHashIterator itr = hashFunction.getRightVariantsNTHashIterator(dbgbfCbfMaxNumHash);
//        char charOut = (char) kmer.bytes[kMinus1];
//        itr.start(kmer.hashVals, charOut);
//        long[] hVals = itr.hVals;
//        
//        byte[] bytes;
//        float count;
//        for (char charIn : getAltNucleotides(charOut)) {
//            itr.next(charIn);
//            count = getCount(hVals);
//            if (count > 0) {
//                bytes = Arrays.copyOf(kmer.bytes, k);
//                bytes[kMinus1] = (byte) charIn;
//                result.add(new Kmer(bytes, count, hVals));
//            }
//        }
//        
////        final String prefix = getPrefix(kmer.toString());
////        String v;
////        float count;
////        
////        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
////        for (char c : getAltNucleotides(kmer.bytes[kMinus1])) {
////            v = prefix + c;
////
////            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
////
////            count = getCount(hashVals);
////            if (count > 0) {
////                result.add(new Kmer(v, count, hashVals, k));
////            }
////        }
//                
//        return result;
//    }
    
    public ArrayDeque<String> getRightVariants(String kmer) {
        ArrayDeque<String> result = new ArrayDeque<>(4);
        final String prefix = getPrefix(kmer);
        String v;
        for (char c : getAltNucleotides(kmer.charAt(kMinus1))) {
            v = prefix + c;
            if (contains(v)) {
                result.add(v);
            }
        }
        return result;
    }
    
    public float[] getCounts(String[] kmers){
        final int numKmers = kmers.length;
        float[] counts = new float[numKmers];
        for (int i=0; i<numKmers; ++i) {
            counts[i] = getCount(kmers[i]);
        }
        return counts;
    }
    
//    public float[] getCounts(String seq) {
//        final int numKmers = getNumKmers(seq, k);
//        
//        float[] counts = new float[numKmers];
//        
//        NTHashIterator itr = getHashIterator();
//        itr.start(seq);
//        long[] hVals = itr.hVals;
//        
//        for (int i=0; i<numKmers; ++i) {
//            itr.next();
//            counts[i] = getCount(hVals);
//        }
//        
//        /*
//        String[] kmers = kmerize(seq, k);
//        for (int i=0; i<numKmers; ++i) {
//            counts[i] = getCount(kmers[i]);
//        }
//        */
//        
//        return counts;
//    }
    
//    public float getMinCount(String seq) {
//        final int numKmers = getNumKmers(seq, k);
//        
//        NTHashIterator itr = getHashIterator();
//        itr.start(seq);
//        long[] hVals = itr.hVals;
//        
//        float minCount = 0;
//        
//        if (numKmers > 0) {
//            itr.next();
//            minCount = getCount(hVals);
//
//            float c;
//            for (int i=1; i<numKmers; ++i) {
//                itr.next();
//                c = getCount(hVals);
//                if (c < minCount) {
//                    minCount = c;
//                }
//            }
//        }
//        
//        return minCount;
//    }
        
//    public boolean isValidSeq(String seq) {
//        NTHashIterator itr = getHashIterator();
//        itr.start(seq);
//        long[] hVals = itr.hVals;
//        
//        while (itr.hasNext()) {
//            itr.next();
//            if (!dbgbf.lookup(hVals)) {
//                return false;
//            }
//        }
//                
//        return true;
//    }
    
    public NTHashIterator getHashIterator() {
        return hashFunction.getHashIterator(this.dbgbfCbfMaxNumHash);
    }
    
    public NTHashIterator getHashIterator(int numHash) {
        return hashFunction.getHashIterator(numHash);
    }
    
    public PairedNTHashIterator getPairedHashIterator() {
        return hashFunction.getPairedHashIterator(this.pkbfNumHash, this.pairedKmersDistance);
    }
    
    public ArrayList<Kmer2> getKmers(String seq) {                
        return hashFunction.getKmers(seq, this.dbgbfCbfMaxNumHash, this);
}
    
    public String assemble(ArrayDeque<Kmer2> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size() + kMinus1);
        
        for (byte b : kmers.getFirst().bytes) {
            sb.append((char) b);
        }
        
        for (Kmer2 e : kmers) {
            sb.append((char) e.bytes[kMinus1]);
        }
        
        return sb.toString();
    }
    
    public String assemble(ArrayList<Kmer2> kmers, int start, int end) {
        StringBuilder sb = new StringBuilder(end - start + kMinus1);
        
        for (byte b : kmers.get(start).bytes) {
            sb.append((char) b);
        }
        
        for (int i=start+1; i<end; ++i) {
            sb.append((char) kmers.get(i).bytes[kMinus1]);
        }
        
        return sb.toString();
    }
    
    public String assembleReverseOrder(ArrayDeque<Kmer2> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size() + kMinus1);
        
        Iterator<Kmer2> itr = kmers.descendingIterator();
        
        for (byte b : itr.next().bytes) {
            sb.append((char) b);
        }
        
        while (itr.hasNext()) {
            sb.append((char) itr.next().bytes[kMinus1]);
        }
        
        return sb.toString();
    }
    
    public String assemble(Collection<Kmer2> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size() + kMinus1);
        
        Iterator<Kmer2> itr = kmers.iterator();
        
        for (byte b : itr.next().bytes) {
            sb.append((char) b);
        }
        
        while (itr.hasNext()) {
            sb.append((char) itr.next().bytes[kMinus1]);
        }
                
        return sb.toString();
    }
  
    public String assembleFirstBase(Collection<Kmer2> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size() + kMinus1);
        
        for (Kmer2 kmer : kmers) {
            sb.append((char) kmer.bytes[0]);
        }
        
        return sb.toString();
    }

    public String assembleLastBase(Collection<Kmer2> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size() + kMinus1);
        
        for (Kmer2 kmer : kmers) {
            sb.append((char) kmer.bytes[kMinus1]);
        }
        
        return sb.toString();
    }
}
