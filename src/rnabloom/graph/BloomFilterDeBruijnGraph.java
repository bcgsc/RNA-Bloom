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
        return cbf.getCount(kmer);
    }

    public float getFPR() {
        return dbgbf.getFPR() * cbf.getFPR();
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
    
    public ArrayList<Kmer> getPredecessors(Kmer kmer) {
        return getPredecessors(kmer.seq);
    }
    
    public ArrayList<Kmer> getPredecessors(String kmer) {
        ArrayList<Kmer> result = new ArrayList<>(4);
        String prefix = kmer.substring(0, overlap);
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

    public ArrayList<Kmer> getSuccessors(Kmer kmer) {
        return getSuccessors(kmer.seq);
    }
    
    public ArrayList<Kmer> getSuccessors(String kmer) {
        ArrayList<Kmer> result = new ArrayList<>(4);
        String suffix = kmer.substring(1, k);
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
    
    public ArrayList<String> getLeftVariants(String kmer) {
        ArrayList<String> result = new ArrayList<>(4);
        String suffix = kmer.substring(1, k);
        String v;
        for (char c : NUCLEOTIDES) {
            v = c + suffix;
            if (contains(v)) {
                result.add(v);
            }
        }
        return result;
    }

    public ArrayList<String> getRightVariants(String kmer) {
        ArrayList<String> result = new ArrayList<>(4);
        
        String prefix = kmer.substring(0, overlap);
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
    
    public float getMeanKmerCoverage(String seq) {
        return getMeanKmerCoverage(kmerize(seq, k));
    }
    
    public float getMeanKmerCoverage(String[] kmers) {
        float count = 0;
        for (String kmer : kmers) {
            count += cbf.getCount(kmer);
        }
        return count/kmers.length;
    }
    
    public float getMeanKmerCoverage(ArrayList<Kmer> kmers) {
        float count = 0;
        for (Kmer kmer : kmers) {
            count += cbf.getCount(kmer.seq);
        }
        return count/kmers.size();
    }
    
    public boolean isValidSeq(String seq) {
        return containsAll(kmerize(seq, k));
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
        String[] kmers = kmerize(seq, this.k);
        ArrayList<Kmer> result = new ArrayList<>(kmers.length);
        
        for (String kmer : kmers) {
            if (this.contains(kmer)) {
                result.add(new Kmer(kmer, this.getCount(kmer)));
            }
            else {
                result.add(new Kmer(kmer, 0));
            }
        }
        
        return result;
    }
}
