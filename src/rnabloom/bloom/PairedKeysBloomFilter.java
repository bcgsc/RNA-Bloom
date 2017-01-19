/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import rnabloom.bloom.hash.HashFunction2;

/**
 *
 * @author kmnip
 */
public class PairedKeysBloomFilter extends BloomFilter {
     
    public PairedKeysBloomFilter(long size, int numHash, HashFunction2 hashFunction) {
        super(size, numHash, hashFunction);
    }
    
    public PairedKeysBloomFilter(File desc, File bits, HashFunction2 hashFunction) throws FileNotFoundException, IOException {
        super(desc, bits, hashFunction);
    }
    
    public void addPair(String key1, String key2) {
        long[] hash3 = hashFunction.getHashValues(key1, key2, numHash);
        for (int h=0; h<numHash; ++h) {
            bitArray.set(getIndex(hash3[h]));
        }        
    }
    
    public void addPair(final long[] hvals1, final long[] hvals2) {
        long[] hash3 = hashFunction.getHashValues(hvals1, hvals2, numHash);
        for (int h=0; h<numHash; ++h) {
            bitArray.set(getIndex(hash3[h]));
        }   
    }
    
    public void addSingleAndPair(String key1, String key2) {
        long[] hash1 = new long[numHash];
        hashFunction.getHashValues(key1, numHash, hash1);
        long[] hash2 = new long[numHash];
        hashFunction.getHashValues(key2, numHash, hash2);
        long[] hash3 = hashFunction.getHashValues(key1, key2, numHash);
        
        for (int h=0; h<numHash; ++h) {
            bitArray.set(getIndex(hash1[h]));
            bitArray.set(getIndex(hash2[h]));
            bitArray.set(getIndex(hash3[h]));
        }
    }

    public void addSingleAndPair(long[] hvals1, long[] hvals2) {
        long[] hash3 = hashFunction.getHashValues(hvals1, hvals2, numHash);
        
        for (int h=0; h<numHash; ++h) {
            bitArray.set(getIndex(hvals1[h]));
            bitArray.set(getIndex(hvals2[h]));
            bitArray.set(getIndex(hash3[h]));
        }
    }
        
    public boolean lookupPair(String key1, String key2) {
        long[] hash = hashFunction.getHashValues(key1, key2, numHash);
        
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(getIndex(hash[h]))) {
                return false;
            }
        }
        
        return true;
    }

    public boolean lookupSingleAndPair(String key1, String key2) {
        long[] hash1 = new long[numHash];
        hashFunction.getHashValues(key1, numHash, hash1);
        long[] hash2 = new long[numHash];
        hashFunction.getHashValues(key2, numHash, hash2);
        long[] hash3 = hashFunction.getHashValues(key1, key2, numHash);
        
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(getIndex(hash1[h])) ||
                    !bitArray.get(getIndex(hash2[h])) ||
                    !bitArray.get(getIndex(hash3[h]))) {
                return false;
            }
        }
        
        return true;
    }
        
}
