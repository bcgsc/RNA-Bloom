/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.random;
import static java.lang.Math.scalb;
import static util.hash.MurmurHash3.murmurhash3_x64_128;

/**
 *
 * @author kmnip
 */
public class CountingBloomFilter implements CountingBloomFilterInterface {
    protected final byte[] counts;
    protected final int numHash;
    protected final int seed;
    protected final int size;
    protected final int keyLength;
    
    protected static final long MAX_SIZE = (long) Integer.MAX_VALUE;
    
    private static final byte MANTISSA = 3;
    private static final byte MANTI_MASK = 0xFF >> (8 - MANTISSA);
    private static final byte ADD_MASK = 0x80 >> (7 - MANTISSA);
    
    public CountingBloomFilter(int size, int numHash, int seed, int keyLength) {
        if (size > MAX_SIZE) {
            throw new UnsupportedOperationException("Size is too large.");
        }
        
        this.size = size;
        this.counts = new byte[size];
        this.numHash = numHash;
        this.seed = seed;
        this.keyLength = keyLength;
    }
        
    @Override
    public synchronized void increment(String key) {
        // generate hash values
        final byte[] b = key.getBytes();
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(b, 0, keyLength, seed, numHash, hashVals);
        
        increment(hashVals);
    }
    
    public synchronized void increment(long[] hashVals) {
        // find the smallest count at all hash positions
        byte min = counts[(int) (hashVals[0] % size)];
        byte c;
        int h;
        for (h=1; h<numHash; ++h) {
            c = counts[(int) (hashVals[h] % size)];
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
            }
        }
        
        // increment the smallest count
        byte updated = min;
        if (min <= MANTI_MASK) {
            ++updated;
        }
        else {
            int shiftVal = (1 << ((updated >> MANTISSA) - 1));
            if ((int) (random() * Integer.MAX_VALUE) % shiftVal == 0) {
                ++updated;
            }
        }
        
        // update the smallest count only
        int index;
        for (h=1; h<numHash; ++h) {
            index = (int) (hashVals[h] % size);

            if (counts[index] == min) {
                counts[index] = updated;
            }
        }
    }

    @Override
    public float getCount(String key) {
        // generate hash values
        final byte[] b = key.getBytes();
        final long[] hashVals = new long[numHash];
        murmurhash3_x64_128(b, 0, keyLength, seed, numHash, hashVals);
        
        return getCount(hashVals);
    }
    
    public float getCount(long[] hashVals) {
        // find the smallest count
        byte min = counts[(int) (hashVals[0] % size)];
        byte c;
        for (int h=1; h<numHash; ++h) {
            c = counts[(int) (hashVals[h] % size)];
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
            }
        }
        
        // return the float value
        if (min <= MANTI_MASK) {
            return (float) min;
        }
        
        return scalb((min & MANTI_MASK) | ADD_MASK, (min >> MANTISSA) - 1);
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        int n = 0;
        for (byte b : counts) {
            if (b > 0) {
                ++n;
            }
        }
        
        return (float) pow(1 - exp(-numHash * n / size), numHash);
    }
    
}
