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
    protected final long[] longArray;
    protected final int numHash;
    protected final int seed;
    protected final long size;
    protected final int keyLength;
    
    protected static final long MAX_SIZE = (long) Integer.MAX_VALUE * (long) Long.SIZE;
    
    public BloomFilter(long size, int numHash, int seed, int keyLength) {
        if (size > MAX_SIZE) {
            throw new UnsupportedOperationException("Size is too large.");
        }
        
        this.size = size;
        this.longArray = new long[(int) (size/Long.SIZE) + 1];
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
        long i;
        for (int h=0; h<numHash; ++h) {
            i = hashVals[h] % size;
            longArray[(int) (i/Long.SIZE)] |= (1 << (int) (i % Long.SIZE));
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
        long i;
        for (int h=0; h<numHash; ++h) {
            i = hashVals[h] % size;
            if ((longArray[(int) (i/Long.SIZE)] & (1 << (int) (i % Long.SIZE))) == 0) {
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
        
        long n = 0;
        for (long l : longArray) {
            n += Long.bitCount(l);
        }
        
        return (float) pow(1 - exp(-numHash * n / size), numHash);
    }
    
}
