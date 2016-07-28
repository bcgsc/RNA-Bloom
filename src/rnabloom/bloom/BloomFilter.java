/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.pow;
import static java.lang.Math.exp;
import static util.hash.MurmurHash3.murmurhash3_x64_128;
import static util.hash.MurmurHash3.LongPair;

/**
 *
 * @author kmnip
 */
public class BloomFilter extends AbstractBloomFilter {
    protected final long[] long_array;
    protected final int num_hash;
    protected final int seed;
    protected final long size;
    protected final int key_length;
    
    protected static final long max_size = (long) Integer.MAX_VALUE * (long) Long.SIZE;
    
    public BloomFilter(long size, int num_hash, int seed, int key_length) {
        if (size > max_size) {
            throw new UnsupportedOperationException("Size is too large.");
        }
        
        this.size = size;
        this.long_array = new long[(int) size/Long.SIZE + 1];
        this.num_hash = num_hash;
        this.seed = seed;
        this.key_length = key_length;
    }
        
    @Override
    public void add(final String key) {
        final byte[] b = key.getBytes();
        
        long bit_index;
        final LongPair p = new LongPair();
        
        for (int i=0; i<num_hash; ++i) {
            murmurhash3_x64_128(b, 0, key_length, seed+i, p);
            
            bit_index = (p.val1 >>> 1) % size;
            long_array[(int) bit_index/Long.SIZE] |= (1 << (int) bit_index % Long.SIZE);
            
            if (++i<num_hash) {
                bit_index = (p.val2 >>> 1) % size;
                long_array[(int) bit_index/Long.SIZE] |= (1 << (int) bit_index % Long.SIZE);
            }
        }
    }

    @Override
    public boolean lookup(final String key) {
        final byte[] b = key.getBytes();
        
        long bit_index;
        final LongPair p = new LongPair();
        
        for (int i=0; i<num_hash; ++i) {
            murmurhash3_x64_128(b, 0, key_length, seed+i, p);
            
            bit_index = (p.val1 >>> 1) % size;
            if ((long_array[(int) bit_index/Long.SIZE] & (1 << (int) bit_index % Long.SIZE)) == 0) {
                return false;
            }
            
            if (++i<num_hash) {
                bit_index = (p.val2 >>> 1) % size;
                if ((long_array[(int) bit_index/Long.SIZE] & (1 << (int) bit_index % Long.SIZE)) == 0) {
                    return false;
                }
            }
        }
        
        return true;
    }

    @Override
    public void save(String path) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void restore(String path) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        long n = 0;
        for (long l : long_array) {
            n += Long.bitCount(l);
        }
        
        return (float) pow(1 - exp(-num_hash * n / size), num_hash);
    }
    
}
