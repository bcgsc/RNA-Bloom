/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import static java.lang.Math.pow;
import static java.lang.Math.exp;

/**
 *
 * @author kmnip
 */
public class BloomFilter implements BloomFilterInterface {    
    protected AbstractLargeBitBuffer bitArray;
    protected final int numHash;
    protected final long size;
    protected final HashFunction hashFunction;
        
    public BloomFilter(long size, int numHash, HashFunction hashFunction) {
        
        this.size = size;
        try {
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
        this.numHash = numHash;
        this.hashFunction = hashFunction;
    }
        
    @Override
    public void add(final String key) {
        add(hashFunction.getHashValues(key));
    }
    
    public void add(final long[] hashVals){
        for (int h=0; h<numHash; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }

    @Override
    public boolean lookup(final String key) {        
        return lookup(hashFunction.getHashValues(key));
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
