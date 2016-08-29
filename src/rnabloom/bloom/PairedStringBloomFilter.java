/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.exp;
import static java.lang.Math.pow;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author kmnip
 */
public class PairedStringBloomFilter implements BloomFilterInterface {    
    
    protected AbstractLargeBitBuffer bitArray;
    protected final long size;
    protected final int numHash1;    
    protected final HashFunction hashFunction1;
    protected final int numHash2;
    protected final HashFunction hashFunction2;
        
    public PairedStringBloomFilter(long size, int numHash1, HashFunction hashFunction1, int numHash2, HashFunction hashFunction2) {
        
        this.size = size;
        try {
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
        this.numHash1 = numHash1;
        this.hashFunction1 = hashFunction1;
        this.numHash2 = numHash2;
        this.hashFunction2 = hashFunction2;
    }
        
    @Override
    public void add(final String key) {
        add(hashFunction1.getHashValues(key));
    }
    
    public void add(final long[] hashVals){
        for (int h=0; h<numHash1; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }

    public void addPair(final String key1, final String key2) {
        long[] hashVals = hashFunction2.getHashValues(key1 + key2);
        for (int h=0; h<numHash2; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }
    
    @Override
    public boolean lookup(final String key) {        
        return lookup(hashFunction1.getHashValues(key));
    }

    public boolean lookup(final long[] hashVals) {
        for (int h=0; h<numHash1; ++h) {
            if (!bitArray.get(hashVals[h] % size)) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean lookupPair(final String key1, final String key2) {
        long[] hashVals = hashFunction2.getHashValues(key1 + key2);
        for (int h=0; h<numHash2; ++h) {
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
        
        return (float) pow(1 - exp(-numHash1 * bitArray.popCount() / size), numHash1);
    }

    public float getPairedFPR() {
        return (float) pow(1 - exp(-numHash2 * bitArray.popCount() / size), numHash2);
    }
    
}
