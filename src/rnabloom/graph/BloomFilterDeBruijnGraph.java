/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.PairedKeysBloomFilter;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.hash.SmallestStrandHashFunction;
import static rnabloom.util.SeqUtils.kmerize;

/**
 *
 * @author kmnip
 */
public class BloomFilterDeBruijnGraph {
    
    private final static char[] NUCLEOTIDES = new char[] {'A','C','G','T'};
    private final BloomFilter dbgbf;
    private final CountingBloomFilter cbf;
    private final PairedKeysBloomFilter pkbf;
    private final int maxHash;
    private final HashFunction hashFunction;
    private final int k;
    private final int overlap;
    private int pairedKmersDistance;
    
    public BloomFilterDeBruijnGraph(long dbgbfNumBits,
                                    long cbfNumBytes,
                                    long pkbfNumBits,
                                    int dbgbfNumHash,
                                    int cbfNumHash,
                                    int pkbfNumHash,
                                    int seed,
                                    int k,
                                    boolean stranded) {
        this.maxHash = Math.max(Math.max(dbgbfNumHash, cbfNumHash), pkbfNumHash);
        this.k = k;
        this.overlap = k-1;
        if (stranded) {
            this.hashFunction = new HashFunction(maxHash, seed, k);
        }
        else {
            this.hashFunction = new SmallestStrandHashFunction(maxHash, seed, k);
        }
        this.dbgbf = new BloomFilter(dbgbfNumBits, dbgbfNumHash, this.hashFunction);
        this.cbf = new CountingBloomFilter(cbfNumBytes, cbfNumHash, this.hashFunction);
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
        final long[] hashVals = hashFunction.getHashValues(kmer);
        dbgbf.add(hashVals);
        cbf.increment(hashVals);
    }
    
    public void addKmersFromSeq(String seq) {
        final int numKmers = seq.length()-k+1;
        
        for (int i=0; i<numKmers; ++i) {
            long[] hashVals = hashFunction.getHashValues(seq.substring(i, i+k));
            dbgbf.add(hashVals);
            cbf.increment(hashVals);
        }
    }
    
    public void addFragmentKmersFromSeq(String seq) {
        int numKmers = seq.length()-k+1;
        
        for (int i=0; i<numKmers; ++i) {
            pkbf.add(hashFunction.getHashValues(seq.substring(i, i+k)));
        }
    }
    
    public void addPairedKmersFromSeq(String seq) {
        int numKmers = seq.length()-k+1;
        long[][] hashValues = new long[numKmers][hashFunction.getNumHash()];
        
        // add kmers
        for (int i=0; i<numKmers; ++i) {
            long[] v = hashFunction.getHashValues(seq.substring(i, i+k));
            pkbf.add(v);
            hashValues[i] = v;
        }
        
        // add paired kmers
        int upperBound = numKmers-pairedKmersDistance;
        for (int i=0; i<upperBound; ++i) {
            pkbf.addPair(hashValues[i], hashValues[i+pairedKmersDistance]);
        }
    }
    
    public boolean lookupFragmentKmer(String kmer) {
        return pkbf.lookup(kmer);
    }

    public boolean lookupPairedKmers(String kmer1, String kmer2) {
        return pkbf.lookupSingleAndPair(kmer1, kmer2);
    }
    
    public boolean contains(String kmer) {
        return dbgbf.lookup(kmer);
    }

    public void increment(String kmer) {
        cbf.increment(kmer);
    }
    
    public float getCount(String kmer) {
        final long[] hashVals = hashFunction.getHashValues(kmer);
        if (dbgbf.lookup(hashVals)) {
            return cbf.getCount(hashVals);
        }
        else {
            return 0;
        }
    }

    public float getFPR() {
        return dbgbf.getFPR() * cbf.getFPR();
    }
    
    public Kmer getKmer(String kmer) {
        return new Kmer(kmer, this.getCount(kmer));
    }
    
    public static class Kmer {
        public String seq;
        public float count;
        
        public Kmer(String seq, float count) {
            this.seq = seq;
            this.count = count;
        }
        
        public boolean equals(Kmer other) {
            return this.seq.equals(other.seq);
        }
    }
    
    private String getPrefix(String kmer) {
        return kmer.substring(0, overlap);
    }
    
    private String getSuffix(String kmer) {
        return kmer.substring(1, k);
    }
    
    public LinkedList<Kmer> getPredecessors(Kmer kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
        String prefix = getPrefix(kmer.seq);
        String v;
        long[] hashVals;
        float count;
        for (char c : NUCLEOTIDES) {
            v = c + prefix;
            hashVals = this.hashFunction.getHashValues(v);
            if (dbgbf.lookup(hashVals)) {
                count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(new Kmer(v, count));
                }
            }
        }
        return result;
    }
    
    public LinkedList<String> getPredecessors(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        String prefix = getPrefix(kmer);
        String v;
        long[] hashVals;
        float count;
        for (char c : NUCLEOTIDES) {
            v = c + prefix;
            hashVals = this.hashFunction.getHashValues(v);
            if (dbgbf.lookup(hashVals)) {
                count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(v);
                }
            }
        }
        return result;
    }

    public LinkedList<Kmer> getSuccessors(Kmer kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
        String suffix = getSuffix(kmer.seq);
        String v;
        long[] hashVals;
        float count;
        for (char c : NUCLEOTIDES) {
            v = suffix + c;
            hashVals = this.hashFunction.getHashValues(v);
            if (dbgbf.lookup(hashVals)) {
                count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(new Kmer(v, count));
                }
            }
        }
        return result;
    }
    
    public LinkedList<String> getSuccessors(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        String suffix = getSuffix(kmer);
        String v;
        long[] hashVals;
        float count;
        for (char c : NUCLEOTIDES) {
            v = suffix + c;
            hashVals = this.hashFunction.getHashValues(v);
            if (dbgbf.lookup(hashVals)) {
                count = cbf.getCount(hashVals);
                if (count > 0) {
                    result.add(v);
                }
            }
        }
        return result;
    }
    
    public LinkedList<String> getLeftVariants(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        String suffix = getSuffix(kmer);
        String v;
        for (char c : NUCLEOTIDES) {
            v = c + suffix;
            if (contains(v)) {
                result.add(v);
            }
        }
        return result;
    }

    public LinkedList<String> getRightVariants(String kmer) {
        LinkedList<String> result = new LinkedList<>();
        
        String prefix = getPrefix(kmer);
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
        int numKmers = kmers.length;
        float[] counts = new float[numKmers];
        for (int i=0; i<numKmers; ++i) {
            counts[i] = getCount(kmers[i]);
        }
        return counts;
    }
    
    public float getMedianKmerCoverage(String seq){
        return getMedianKmerCoverage(kmerize(seq, k));
    }
    
    public float getMedianKmerCoverage(String[] kmers) {
        int numKmers = kmers.length;
        int halfNumKmers = numKmers/2;
        
        ArrayList<Float> counts = new ArrayList<>(numKmers);
        for (String kmer : kmers) {
            counts.add(cbf.getCount(kmer));
        }
        
        Collections.sort(counts);
        
        if (numKmers % 2 == 0) {
            return (counts.get(halfNumKmers) + counts.get(halfNumKmers -1))/2.0f;
        }
        
        return counts.get(halfNumKmers);
    }
        
    public boolean isValidSeq(String seq) {
        return containsAll(kmerize(seq, k));
    }
    
    public boolean isValidSeq(Collection<Kmer> kmers) {
        for (Kmer kmer : kmers) {
            if (kmer.count == 0) {
                return false;
            }
        }
        return true;
    }
    
    public boolean containsAll(String[] kmers) {
        for (String kmer : kmers) {
            if (!contains(kmer)) {
                return false;
            }
        }
        return true;
    }
    
    public ArrayList<Kmer> getKmers(String seq) {
        final int numKmers = seq.length()-k+1;
        ArrayList<Kmer> result = new ArrayList<>(numKmers);
        
        for (int i=0; i<numKmers; ++i) {
            result.add(getKmer(seq.substring(i, i+k)));
        }
        
        return result;
    }
}
