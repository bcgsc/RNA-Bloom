/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.graph;

import java.util.ArrayList;
import java.util.List;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.hash.SmallestStrandHashFunction;

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
    
    public void add(String kmer) {
        final long[] hashVals = hashFunction.getHashValues(kmer);
        dbgbf.add(hashVals);
        cbf.increment(hashVals);
    }

    public boolean lookup(String kmer) {
        return dbgbf.lookup(kmer);
    }

    public void increment(String key) {
        cbf.increment(key);
    }
    
    public float getCount(String kmer) {
        return cbf.getCount(kmer);
    }

    public float getFPR() {
        return dbgbf.getFPR() * cbf.getFPR();
    }
    
    public List<String> getPredecessors(String kmer) {
        List<String> result = new ArrayList<>(4);
        String prefix = kmer.substring(0, overlap);
        String v;
        for (char c : NUCLEOTIDES) {
            v = c + prefix;
            if (lookup(v)) {
                result.add(v);
            }
        }
        return result;
    }

    public List<String> getSuccessors(String kmer) {
        List<String> result = new ArrayList<>(4);
        String suffix = kmer.substring(1, k);
        String v;
        for (char c : NUCLEOTIDES) {
            v = suffix + c;
            if (lookup(v)) {
                result.add(v);
            }
        }
        return result;
    }
    
    public List<String> getLeftVariants(String kmer) {
        List<String> result = new ArrayList<>(4);
        String suffix = kmer.substring(1, k);
        String v;
        for (char c : NUCLEOTIDES) {
            v = c + suffix;
            if (lookup(v)) {
                result.add(v);
            }
        }
        return result;
    }

    public List<String> getRightVariants(String kmer) {
        List<String> result = new ArrayList<>(4);
        
        String prefix = kmer.substring(0, overlap);
        String v;
        for (char c : NUCLEOTIDES) {
            v = prefix + c;
            if (lookup(v)) {
                result.add(v);
            }
        }
        return result;
    }
    
    
}
