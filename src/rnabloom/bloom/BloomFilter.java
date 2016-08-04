/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.pow;
import static java.lang.Math.exp;
import static util.hash.MurmurHash3.murmurhash3_x64_128;

/**
 *
 * @author kmnip
 */
public class BloomFilter implements BloomFilterInterface {    
    protected final LargeBitBuffer bitArray;
    protected final int numHash;
    protected final int seed;
    protected final long size;
    protected final int keyLength;
        
    public BloomFilter(long size, int numHash, int seed, int keyLength) {
        
        this.size = size;
        this.bitArray = new LargeBitBuffer(size);
        this.numHash = numHash;
        this.seed = seed;
        this.keyLength = keyLength;
    }
        
    @Override
    public synchronized void add(final String key) {
        final byte[] b = key.getBytes();
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(b, 0, keyLength, seed, numHash, hashVals);
        
        add(hashVals);
    }
    
    public synchronized void add(final long[] hashVals){
        for (int h=0; h<numHash; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }

    @Override
    public boolean lookup(final String key) {
        final byte[] b = key.getBytes();
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(b, 0, keyLength, seed, numHash, hashVals);
        
        return lookup(hashVals);
    }

    public boolean lookup(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(hashVals[h] % size)) {
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        return (float) pow(1 - exp(-numHash * bitArray.popCount() / size), numHash);
    }
    
}
