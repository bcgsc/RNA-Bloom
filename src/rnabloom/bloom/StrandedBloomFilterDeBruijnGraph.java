/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import util.hash.MurmurHash3;
import static util.hash.MurmurHash3.murmurhash3_x64_128;

/**
 *
 * @author kmnip
 */
public class StrandedBloomFilterDeBruijnGraph implements DeBruijnGraphBloomFilterInterface {
    
    private BloomFilter dbgbf;
    private CountingBloomFilter cbf;
    private int seed;
    private int maxHash;
    private int k;
    
    public StrandedBloomFilterDeBruijnGraph(long dbgbfSize,
                                            int cbfSize,
                                            int dbgbfNumHash,
                                            int cbfNumHash,
                                            int seed,
                                            int k) {
        this.seed = seed;
        this.maxHash = Math.max(dbgbfNumHash, cbfNumHash);
        this.k = k;
        this.dbgbf = new BloomFilter(dbgbfSize, dbgbfNumHash, seed, k);
        this.cbf = new CountingBloomFilter(cbfSize, cbfNumHash, seed, k);
    }
    
    public long[] getHashValues(final String kmer) {
        final byte[] b = kmer.getBytes();
        final long[] hashVals = new long[maxHash];
        murmurhash3_x64_128(b, 0, k, seed, maxHash, hashVals);
                
        return hashVals;
    }
    
    @Override
    public void add(String kmer) {
        final long[] hashVals = getHashValues(kmer);
        dbgbf.add(hashVals);
        cbf.increment(hashVals);
    }

    @Override
    public boolean lookup(String kmer) {
        return dbgbf.lookup(kmer);
    }

    @Override
    public void increment(String kmer) {
        cbf.increment(kmer);
    }

    @Override
    public float getCount(String kmer) {
        return cbf.getCount(kmer);
    }

    @Override
    public float getFPR() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String[] getPredecessors(String kmer) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String[] getSuccessors(String kmer) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
