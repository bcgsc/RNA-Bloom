/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author kmnip
 */
public class PairedKeysBloomFilter extends BloomFilter {
     
    public PairedKeysBloomFilter(long size, int numHash, HashFunction hashFunction) {
        super(size, numHash, hashFunction);
    }
    
    public void addPair(final String key1, final String key2) {
        addPair(super.hashFunction.getHashValues(key1), super.hashFunction.getHashValues(key2));
    }

    public void addPair(final long[] hash1, final long[] hash2) {
        for (int h=0; h<numHash; ++h) {
            bitArray.set(HashFunction.combineHashValues(hash1[h], hash2[h]) % size);
        }
    }
    
    public void addSingleAndPair(final String key1, final String key2) {
        long[] hash1 = super.hashFunction.getHashValues(key1);
        long[] hash2 = super.hashFunction.getHashValues(key2);
        
        for (int h=0; h<numHash; ++h) {
            bitArray.set(hash1[h] % size);
            bitArray.set(hash2[h] % size);
            bitArray.set(HashFunction.combineHashValues(hash1[h], hash2[h]) % size);
        }
    }
        
    public boolean lookupPair(final String key1, final String key2) {
        long[] hash1 = super.hashFunction.getHashValues(key1);
        long[] hash2 = super.hashFunction.getHashValues(key2);
        
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(HashFunction.combineHashValues(hash1[h], hash2[h]) % size)) {
                return false;
            }
        }
        
        return true;
    }

    public boolean lookupSingleAndPair(final String key1, final String key2) {
        long[] hash1 = super.hashFunction.getHashValues(key1);
        long[] hash2 = super.hashFunction.getHashValues(key2);
        
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(hash1[h] % size) ||
                    !bitArray.get(hash2[h] % size) ||
                    !bitArray.get(HashFunction.combineHashValues(hash1[h], hash2[h]) % size)) {
                return false;
            }
        }
        
        return true;
    }
        
}
