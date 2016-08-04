/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static util.hash.MurmurHash3.murmurhash3_x64_128;

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
}
