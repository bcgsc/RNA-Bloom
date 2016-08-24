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
import java.util.List;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.hash.SmallestStrandHashFunction;
import static rnabloom.util.SeqUtils.kmerize;

/**
 *
 * @author kmnip
 */
public class BloomFilterDeBruijnGraph {
    
    private final BloomFilter dbgbf;
    private final CountingBloomFilter cbf;
    private final int maxHash;
    private final HashFunction hashFunction;
    private final int k;
    private final int overlap;
    private final static char[] NUCLEOTIDES = new char[] {'A','C','G','T'};
        
    public BloomFilterDeBruijnGraph(long dbgbfSize,
                                    int cbfSize,
                                    int dbgbfNumHash,
                                    int cbfNumHash,
                                    int seed,
                                    int k,
                                    boolean stranded) {
        this.maxHash = Math.max(dbgbfNumHash, cbfNumHash);
        this.k = k;
        this.overlap = k-1;
        if (stranded) {
            this.hashFunction = new HashFunction(maxHash, seed, k);
        }
        else {
            this.hashFunction = new SmallestStrandHashFunction(maxHash, seed, k);
        }
        this.dbgbf = new BloomFilter(dbgbfSize, dbgbfNumHash, this.hashFunction);
        this.cbf = new CountingBloomFilter(cbfSize, cbfNumHash, this.hashFunction);
    }
    
    public int getK() {
        return k;
    }
    
    public void add(String kmer) {
        final long[] hashVals = hashFunction.getHashValues(kmer);
        dbgbf.add(hashVals);
        cbf.increment(hashVals);
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
    
    public String getPrefix(String kmer) {
        return kmer.substring(0, overlap);
    }
    
    public String getSuffix(String kmer) {
        return kmer.substring(1, k);
    }

    public char getFirstBase(String kmer) {
        return kmer.charAt(0);
    }
    
    public char getLastBase(String kmer) {
        return kmer.charAt(overlap);
    }
    
    public LinkedList<Kmer> getPredecessors(Kmer kmer) {
        return getPredecessors(kmer.seq);
    }
    
    public LinkedList<Kmer> getPredecessors(String kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
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
                    result.add(new Kmer(v, count));
                }
            }
        }
        return result;
    }

    public LinkedList<Kmer> getSuccessors(Kmer kmer) {
        return getSuccessors(kmer.seq);
    }
    
    public LinkedList<Kmer> getSuccessors(String kmer) {
        LinkedList<Kmer> result = new LinkedList<>();
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
                    result.add(new Kmer(v, count));
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
        final int numKmers = seq.length() - k + 1;
        ArrayList<Kmer> result = new ArrayList<>(numKmers);
        
        for (int i=0; i<numKmers; ++i) {
            result.add(getKmer(seq.substring(i, i+k)));
        }
        
        return result;
    }
}
