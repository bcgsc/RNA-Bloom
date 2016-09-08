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
    private final int numHash;
    private final int seed;
    private final int k;
    
    public HashFunction(int numHash, int seed, int k) {
        this.numHash = numHash;
        this.seed = seed;
        this.k = k;
    }
    
    public long[] getHashValues(final String kmer) {
        final byte[] b = kmer.getBytes();
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(b, 0, k, seed, numHash, hashVals);
                
        return hashVals;
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
