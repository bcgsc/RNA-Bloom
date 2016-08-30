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
public class PairedKeysBloomFilter extends BloomFilter {
    
    protected final int numHash2;
    protected final HashFunction hashFunction2;
        
    public PairedKeysBloomFilter(long size, int singleKeyNumHash, HashFunction singleKeyHashFunction, int pairedKeysNumHash, HashFunction pairedKeysHashFunction) {
        super(size, singleKeyNumHash, singleKeyHashFunction);
        this.numHash2 = pairedKeysNumHash;
        this.hashFunction2 = pairedKeysHashFunction;
    }
    
    public void addPair(final String key1, final String key2) {
        long[] hashVals = hashFunction2.getHashValues(key1 + key2);
        for (int h=0; h<numHash2; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }

    public void addSingleAndPair(final String key1, final String key2) {
        this.add(key1);
        this.add(key2);
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

    public boolean lookupSingleAndPair(final String key1, final String key2) {
        return this.lookup(key1) && this.lookup(key2) && this.lookupPair(key1, key2);
    }
    
    public float getPairedFPR() {
        return (float) pow(1 - exp(-numHash2 * bitArray.popCount() / size), numHash2);
    }
    
}
