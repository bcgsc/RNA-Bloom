/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import static java.lang.Math.exp;
import static java.lang.Math.pow;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author kmnip
 */
public class PairedStringBloomFilter extends BloomFilter {    
    
    protected final int numHash2;
    protected final HashFunction hashFunction2;
        
    public PairedStringBloomFilter(long size, int numHash1, HashFunction hashFunction1, int numHash2, HashFunction hashFunction2) {
        
        super(size, numHash1, hashFunction1);
        this.numHash2 = numHash2;
        this.hashFunction2 = hashFunction2;
    }

    public void addPair(final String key1, final String key2) {
        long[] hashVals = hashFunction2.getHashValues(key1 + key2);
        for (int h=0; h<numHash2; ++h) {
            bitArray.set(hashVals[h] % size);
        }
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

    public float getPairedFPR() {
        return (float) pow(1 - exp(-numHash2 * bitArray.popCount() / size), numHash2);
    }
    
}
