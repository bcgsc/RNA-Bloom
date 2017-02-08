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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.PairedKeysBloomFilter;
import rnabloom.bloom.hash.CanonicalHashFunction2;
import rnabloom.bloom.hash.HashFunction2;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.util.SeqUtils.KmerSeqIterator;
import static rnabloom.util.SeqUtils.NUCLEOTIDES;
import static rnabloom.util.SeqUtils.getNumKmers;

/**
 *
 * @author kmnip
 */
public class BloomFilterDeBruijnGraph {
    
    private final BloomFilter dbgbf;
    private final CountingBloomFilter cbf;
    private PairedKeysBloomFilter pkbf = null;
    
    private int dbgbfCbfMaxNumHash;
    private final HashFunction2 hashFunction;
    private int k;
    private int overlap;
    private boolean stranded;
    private int pairedKmersDistance = -1;
    private long pkbfNumBits;
    private int pkbfNumHash;
    
    private final static String FILE_DESC_EXTENSION = ".desc";
    private final static String FILE_DBGBF_EXTENSION = ".dbgbf";
    private final static String FILE_CBF_EXTENSION = ".cbf";
    private final static String FILE_PKBF_EXTENSION = ".pkbf";
    
    private final static String LABEL_SEPARATOR = ":";
    private final static String LABEL_DBGBF_CBF_NUM_HASH = "dbgbfCbfMaxNumHash";
    private final static String LABEL_K = "k";
    private final static String LABEL_STRANDED = "stranded";
    private final static String LABEL_SEED = "seed";
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
        this.overlap = k-1;
        this.stranded = stranded;
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
                    overlap = k-1;
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
        
        String pkbfBitsPath = graphFile.getPath() + FILE_PKBF_EXTENSION;
        String pkbfDescPath = pkbfBitsPath + FILE_DESC_EXTENSION;
        File pkbfBitsFile = new File(pkbfBitsPath);
        File pkbfDescFile = new File(pkbfDescPath);
        if (pkbfDescFile.isFile() && pkbfBitsFile.isFile()) {
            pkbf = new PairedKeysBloomFilter(pkbfDescFile, pkbfBitsFile, hashFunction);
        }
    }

    public HashFunction2 getHashFunction() {
        return this.hashFunction;
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

    public PairedKeysBloomFilter getPkbf() {
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
        
        String pkbfBitsPath = graphFile.getPath() + FILE_PKBF_EXTENSION;
        String pkbfDescPath = pkbfBitsPath + FILE_DESC_EXTENSION;
        pkbf.save(new File(pkbfDescPath), new File(pkbfBitsPath));
    }
    
    public void restorePkbf(File graphFile) throws IOException {
        String pkbfBitsPath = graphFile.getPath() + FILE_PKBF_EXTENSION;
        String pkbfDescPath = pkbfBitsPath + FILE_DESC_EXTENSION;
        pkbf = new PairedKeysBloomFilter(new File(pkbfDescPath), new File(pkbfBitsPath), hashFunction);
    }
    
    public void initializePairKmersBloomFilter() {
        this.pkbf = new PairedKeysBloomFilter(pkbfNumBits, pkbfNumHash, this.hashFunction);
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
    
    public void add(String kmer) {
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        hashFunction.getHashValues(kmer, dbgbfCbfMaxNumHash, hashVals);
        dbgbf.add(hashVals);
        cbf.increment(hashVals);
    }
    
    public void add(final long[] hashVals) {
        dbgbf.add(hashVals);
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
    
    public void addPairedKmersFromSeq(String seq) {
        // add paired kmers
        final int upperBound = getNumKmers(seq, k) - pairedKmersDistance;
        for (int i=0; i<upperBound; ++i) {
            pkbf.addSingleAndPair(seq.substring(i, i+k), seq.substring(i+pairedKmersDistance, i+k+pairedKmersDistance));
        }
    }
    
    public void addPairedKmers(ArrayList<Kmer> kmers) {
        // add paired kmers
        final int upperBound = kmers.size() - pairedKmersDistance;
        for (int i=0; i<upperBound; ++i) {
            pkbf.addSingleAndPair(kmers.get(i).hashVals, kmers.get(i+pairedKmersDistance).hashVals);
        }
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
    
    public boolean lookupFragmentKmer(String kmer) {
        return pkbf.lookup(kmer);
    }

    public boolean lookupFragmentKmer(final long[] hashVals) {
        return pkbf.lookup(hashVals);
    }
    
    public boolean lookupPairedKmers(String kmer1, String kmer2) {
        return pkbf.lookupSingleAndPair(kmer1, kmer2);
    }
    
    public boolean lookupPairedKmers(long[] hashVals1, long[] hashVals2) {
        return pkbf.lookupSingleAndPair(hashVals1, hashVals2);
    }
    
    public boolean lookupKmerPairing(long[] hashVals1, long[] hashVals2) {
        return pkbf.lookupPair(hashVals1, hashVals2);
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
    
    public Kmer getKmer(String kmer) {
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        hashFunction.getHashValues(kmer, dbgbfCbfMaxNumHash, hashVals);
        return new Kmer(kmer, this.getCount(hashVals), hashVals, false);
    }
    
    public static class Kmer {
        public String seq;
        public float count;
        public long[] hashVals;
        public LinkedList<Kmer> predecessors = null;
        public LinkedList<Kmer> successors = null;
                
        public Kmer(String seq, float count, long[] hashVals) {
            this.seq = seq;
            this.count = count;
            this.hashVals = Arrays.copyOf(hashVals, hashVals.length);
        }
        
        public Kmer(String seq, float count, long[] hashVals, boolean copyHashVals) {
            this.seq = seq;
            this.count = count;
            if (copyHashVals) {
                this.hashVals = Arrays.copyOf(hashVals, hashVals.length);
            }
            else {
                this.hashVals = hashVals;
            }
        }
        
        public boolean equals(Kmer other) {
            return this.seq.equals(other.seq);
        }
    }
    
    public String getPrefix(String kmer) {
        return kmer.substring(0, overlap);
    }
    
    public String getSuffix(String kmer) {
        return kmer.substring(1, k);
    }
    
    public CharSequence getPrefixCharSeq(String kmer) {
        return kmer.subSequence(0, overlap);
    }
    
    public CharSequence getSuffixCharSeq(String kmer) {
        return kmer.subSequence(1, k);
    }
    
    public LinkedList<Kmer> getPredecessors(Kmer kmer) {
        if (kmer.predecessors != null) {
            return kmer.predecessors;
        }
        
        LinkedList<Kmer> result = new LinkedList<>();
        
        final StringBuilder buffer = new StringBuilder(k);
        buffer.append('A');
        buffer.append(getPrefix(kmer.seq));
               
        long[][] allHashVals = hashFunction.getPredecessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, kmer.seq.charAt(k-1));
        
        for (int i=0; i<4; ++i) {
            if (dbgbf.lookup(allHashVals[i])) {
                float count = cbf.getCount(allHashVals[i]);
                if (count > 0) {
                    buffer.setCharAt(0, NUCLEOTIDES[i]);
                    result.add(new Kmer(buffer.toString(), count, allHashVals[i], false));
                }
            }
        }
        
//        String prefix = getPrefix(kmer.seq);
//        long[] hashVals = new long[dbgbfCbfMaxNumHash];
//        for (char c : NUCLEOTIDES) {
//            String v = c + prefix;
//            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
//            if (dbgbf.lookup(hashVals)) {
//                float count = cbf.getCount(hashVals);
//                if (count > 0) {
//                    result.add(new Kmer(v, count, hashVals));
//                }
//            }
//        }
        
        kmer.predecessors = result;
        
        return result;
    }
    
    public LinkedList<String> getPredecessors(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        final String prefix = getPrefix(kmer);
        long[] hashVals = new long[dbgbfCbfMaxNumHash];
        for (char c : NUCLEOTIDES) {
            String v = c + prefix;
            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
            if (dbgbf.lookup(hashVals)) {
                float count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(v);
                }
            }
        }
        return result;
    }

    public LinkedList<Kmer> getSuccessors(Kmer kmer) {
        if (kmer.successors != null) {
            return kmer.successors;
        }
        
        LinkedList<Kmer> result = new LinkedList<>();
        
        final StringBuilder buffer = new StringBuilder(k);
        buffer.append(getSuffix(kmer.seq));
        buffer.append('A');
               
        long[][] allHashVals = hashFunction.getSuccessorsHashValues(dbgbfCbfMaxNumHash, kmer.hashVals, kmer.seq.charAt(0));
        
        int lastIndex = k-1;
        for (int i=0; i<4; ++i) {
            if (dbgbf.lookup(allHashVals[i])) {
                float count = cbf.getCount(allHashVals[i]);
                if (count > 0) {
                    buffer.setCharAt(lastIndex, NUCLEOTIDES[i]);
                    result.add(new Kmer(buffer.toString(), count, allHashVals[i], false));
                }
            }
        }
        
//        final String suffix = getSuffix(kmer.seq);
//        long[] hashVals = new long[dbgbfCbfMaxNumHash];
//        for (char c : NUCLEOTIDES) {
//            String v = suffix + c;
//            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
//            if (dbgbf.lookup(hashVals)) {
//                float count = cbf.getCount(hashVals);
//                if (count > 0) {
//                    result.add(new Kmer(v, count, hashVals));
//                }
//            }
//        }
        
        kmer.successors = result;
        
        return result;
    }
    
    public LinkedList<String> getSuccessors(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        final String suffix = getSuffix(kmer);
        long[] hashVals = new long[dbgbfCbfMaxNumHash];
        for (char c : NUCLEOTIDES) {
            String v = suffix + c;
            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
            if (dbgbf.lookup(hashVals)) {
                float count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(v);
                }
            }
        }
        return result;
    }

    public LinkedList<Kmer> getLeftVariants(Kmer kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
        
        final String suffix = getSuffix(kmer.seq);
        String v;
        float count;
        
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        for (char c : NUCLEOTIDES) {
            v = c + suffix;
            
            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
            
            count = getCount(hashVals);
            if (count > 0) {
                result.add(new Kmer(v, count, hashVals));
            }
        }
        
        return result;
    }
    
    public LinkedList<String> getLeftVariants(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        final String suffix = getSuffix(kmer);
        String v;
        for (char c : NUCLEOTIDES) {
            v = c + suffix;
            
            if (contains(v)) {
                result.add(v);
            }
        }
        return result;
    }

    public LinkedList<Kmer> getRightVariants(Kmer kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
        
        final String prefix = getPrefix(kmer.seq);
        String v;
        float count;
        
        final long[] hashVals = new long[dbgbfCbfMaxNumHash];
        for (char c : NUCLEOTIDES) {
            v = prefix + c;
            
            hashFunction.getHashValues(v, dbgbfCbfMaxNumHash, hashVals);
            
            count = getCount(hashVals);
            if (count > 0) {
                result.add(new Kmer(v, count, hashVals));
            }
        }
        
        return result;
    }
    
    public LinkedList<String> getRightVariants(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        final String prefix = getPrefix(kmer);
        String v;
        for (char c : NUCLEOTIDES) {
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
    
    public float[] getCounts(String seq) {
        final int numKmers = getNumKmers(seq, k);
        
        float[] counts = new float[numKmers];
        
        NTHashIterator itr = getHashIterator();
        itr.start(seq);
        long[] hVals = itr.hVals;
        
        for (int i=0; i<numKmers; ++i) {
            itr.next();
            counts[i] = getCount(hVals);
        }
        
        /*
        String[] kmers = kmerize(seq, k);
        for (int i=0; i<numKmers; ++i) {
            counts[i] = getCount(kmers[i]);
        }
        */
        
        return counts;
    }
    
    public float getMinCount(String seq) {
        final int numKmers = getNumKmers(seq, k);
        
        NTHashIterator itr = getHashIterator();
        itr.start(seq);
        long[] hVals = itr.hVals;
        
        float minCount = 0;
        
        if (numKmers > 0) {
            itr.next();
            minCount = getCount(hVals);

            float c;
            for (int i=1; i<numKmers; ++i) {
                itr.next();
                c = getCount(hVals);
                if (c < minCount) {
                    minCount = c;
                }
            }
        }
        
        return minCount;
    }
        
    public float[] getMinMedianMaxKmerCoverage(String seq) {
        float[] minMedianMax = new float[3];
        
        KmerSeqIterator itr = new KmerSeqIterator(seq, k);
        final int numKmers = itr.numKmers;
        final int halfNumKmers = numKmers/2;
        ArrayList<Float> counts = new ArrayList<>(numKmers);
        
        while (itr.hasNext()) {
            counts.add(cbf.getCount(itr.next()));
        }
        
        Collections.sort(counts);
        
        minMedianMax[0] = counts.get(0);
        minMedianMax[2] = counts.get(numKmers-1);
        
        if (numKmers % 2 == 0) {
            minMedianMax[1] = (counts.get(halfNumKmers) + counts.get(halfNumKmers -1))/2.0f;
        }
        else {
            minMedianMax[1] = counts.get(halfNumKmers);
        }
        
        return minMedianMax;
    }
        
    public boolean isValidSeq(String seq) {
        NTHashIterator itr = getHashIterator();
        itr.start(seq);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            if (!dbgbf.lookup(hVals)) {
                return false;
            }
        }
                
        return true;
    }
    
    public NTHashIterator getHashIterator() {
        return hashFunction.getHashIterator(this.dbgbfCbfMaxNumHash);
    }
    
    public NTHashIterator getHashIterator(int numHash) {
        return hashFunction.getHashIterator(numHash);
    }
    
    public ArrayList<Kmer> getKmers(String seq) {        
        ArrayList<Kmer> result = new ArrayList<>();
        
        NTHashIterator itr = getHashIterator();
        itr.start(seq);
        long[] hVals = itr.hVals;
        int i;
        while (itr.hasNext()) {
            itr.next();
            i = itr.getPos();
            if (dbgbf.lookup(hVals)) {
                result.add(new Kmer(seq.substring(i, i+k), cbf.getCount(hVals), hVals));
            }
            else {
                result.add(new Kmer(seq.substring(i, i+k), 0f, hVals));
            }
        }
        
        return result;
    }
    
    public class KmerIterator implements Iterator<Kmer> {
        private String seq;
        private int i;
        public int numKmers;
        private NTHashIterator itr;
        private long[] hVals;

        public KmerIterator(String seq) {
            this.seq = seq;
            this.numKmers = seq.length() - k + 1;
            this.itr = hashFunction.getHashIterator(dbgbfCbfMaxNumHash);
            itr.start(seq);
            hVals = itr.hVals;
        }

        @Override
        public boolean hasNext() {
            return itr.hasNext();
        }

        @Override
        public Kmer next() {
            itr.next();
            
            i = itr.getPos();
            
            if (dbgbf.lookup(hVals)) {
                return new Kmer(seq.substring(i, i+k), cbf.getCount(hVals), hVals);
            }
            
            return new Kmer(seq.substring(i, i+k), 0, hVals);
        }
    }
}
