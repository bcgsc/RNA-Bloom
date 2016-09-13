/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.MurmurHash3.murmurhash3_x64_128;

/**
 *
 * @author kmnip
 */
public class HashFunction {
    protected final int numHash;
    protected final int seed;
    protected final int k;
    
    public HashFunction(int numHash, int seed, int k) {
        this.numHash = numHash;
        this.seed = seed;
        this.k = k;
    }
    
    public long[] getHashValues(final String kmer) {
        return getHashValues(kmer, numHash);
    }
    
    public long[] getHashValues(final String kmer, final int numHash) {
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(kmer.getBytes(), 0, k, seed, numHash, hashVals);
                
        return hashVals;
    }
    
    public long[] getHashValues(final String kmer1, final String kmer2) {
        return getHashValues(kmer1, kmer2, numHash);
    }
    
    public long[] getHashValues(final String kmer1, final String kmer2, final int numHash) {
        final long[] hashVals1 = new long[numHash];
        murmurhash3_x64_128(kmer1.getBytes(), 0, k, seed, numHash, hashVals1);
        
        final long[] hashVals2 = new long[numHash];
        murmurhash3_x64_128(kmer2.getBytes(), 0, k, seed, numHash, hashVals2);
        
        final long[] hashVal = new long[numHash];
        for (int i=0; i<numHash; ++i) {
            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
        }
        
        return hashVal;
    }
    
    public int getNumHash() {
        return numHash;
    }
    
    public int getSeed() {
        return seed;
    }
    
    public int getK() {
        return k;
    }
    
    public static long combineHashValues(long a, long b) {
        // See: http://stackoverflow.com/a/27952689
        a ^= b + 0x9e3779b9 + (a << 6) + (b >> 2);
        return a < 0 ? -a : a;
    }
}
