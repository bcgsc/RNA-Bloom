/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static util.hash.MurmurHash3.murmurhash3_x86_32;
import static java.lang.Math.random;
import static java.lang.Math.scalb;
import static util.hash.MurmurHash3.murmurhash3_x86_32;
import static java.lang.Math.scalb;

/**
 *
 * @author kmnip
 */
public class CountingBloomFilter extends AbstractCountingBloomFilter {
    protected final byte[] counts;
    protected final int num_hash;
    protected final int seed;
    protected final int[] seeds;
    protected final int size;
    protected final int key_length;
    
    private static final byte MANTISSA = 3;
    private static final byte MANTI_MASK = 0xFF >> (8 - MANTISSA);
    private static final byte ADD_MASK = 0x80 >> (7 - MANTISSA);
    
    public CountingBloomFilter(int size, int num_hash, int seed, int key_length) {        
        this.size = size;
        this.counts = new byte[size];
        this.num_hash = num_hash;
        this.seed = seed;
        this.seeds = new int[num_hash];
        for (int i=0; i<num_hash; ++i) {
            this.seeds[i] = seed + i;
        }
        this.key_length = key_length;
    }
    
    private byte min_val(String key) {
        byte min = counts[(murmurhash3_x86_32(key, 0, key_length, seeds[0]) >>> 1) % size];
        
        byte c;
        for (int i=1; i<num_hash; ++i) {
            c = counts[(murmurhash3_x86_32(key, 0, key_length, seeds[i]) >>> 1) % size];
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
            }
        }

        return min;        
    }
    
    @Override
    public void increment(String key) {
        final byte min = min_val(key);
        
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
        
        int index;
        for (int s : seeds) {
            index = (murmurhash3_x86_32(key, 0, key_length, s) >>> 1) % size;

            if (counts[index] == min) {
                counts[index] = updated;
            }
        }
    }

    @Override
    public float lookup(String key) {
        final byte min = min_val(key);
        
        if (min <= MANTI_MASK) {
            return (float) min;
        }
        
        return scalb((min & MANTI_MASK) | ADD_MASK, (min >> MANTISSA) - 1);
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
        
        int n = 0;
        for (byte b : counts) {
            if (b > 0) {
                ++n;
            }
        }
        
        return (float) pow(1 - exp(-num_hash * n / size), num_hash);
    }
    
}
