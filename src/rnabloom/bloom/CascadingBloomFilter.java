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
public class CascadingBloomFilter implements BloomFilterInterface {
    protected final BloomFilter[] bfs;
    protected final BloomFilter topLevelBf;
    protected final int numLevels;
    protected final long size; // number of bits
    protected final long partitionSize;
    protected final int numHash;
    protected final HashFunction hashFunction;
    
    public CascadingBloomFilter(long size, int numHash, HashFunction hashFunction, int numLevels) {
        this.numLevels = numLevels;
        bfs = new BloomFilter[numLevels];
        this.size = size;
        this.partitionSize = size/numLevels;
        this.numHash = numHash;
        this.hashFunction = hashFunction;
        for (int i=0; i<numLevels; ++i) {
            bfs[i] = new BloomFilter(partitionSize, numHash, hashFunction);
        }
        topLevelBf = bfs[numLevels-1];
    }
    
    public int getNumLevels() {
        return numLevels;
    }
    
    public BloomFilter getBloomFilter(int level) {
        return bfs[level];
    }
    
    @Override
    public void add(String key) {
        long[] hashVals = new long[numHash];
        this.hashFunction.getHashValues(key, numHash, hashVals);
        add(hashVals);
    }

    public void add(long[] hashVals) {
        for (BloomFilter bf : bfs) {
            if (!bf.lookupThenAdd(hashVals)) {
                break;
            }
        }
    }
    
    @Override
    public boolean lookup(String key) {
        long[] hashVals = new long[numHash];
        this.hashFunction.getHashValues(key, numHash, hashVals);
        return lookup(hashVals);
    }

    public boolean lookup(long[] hashVals) {
        return topLevelBf.lookup(hashVals);
    }

    public boolean lookup(long[] hashVals, int level) {
        return bfs[level].lookup(hashVals);
    }
    
    public boolean lookupThenAdd(long[] hashVals) {
        for (BloomFilter bf : bfs) {
            if (!bf.lookupThenAdd(hashVals)) {
                return false;
            }
        }
        return true;
    }
    
    @Override
    public float getFPR() {
        return topLevelBf.getFPR();
    }
    
    public float getFPR(int level) {
        return bfs[level].getFPR();
    }
    
    public void empty() {
        for (BloomFilter bf : bfs) {
            bf.empty();
        }
    }
    
    public void destroy() {
        for (BloomFilter bf : bfs) {
            bf.destroy();
        }
    }
    
    public boolean equivalent(CascadingBloomFilter bf) {
        boolean paramsEqual = this.numLevels != bf.numLevels ||
                                this.size != bf.size || 
                                this.numHash != bf.numHash ||
                                this.partitionSize != bf.partitionSize;
        
        if (paramsEqual) {
            for (int i=0; i<this.numLevels; ++i) {
                if (!this.bfs[i].equivalent(bf.getBloomFilter(i))) {
                    return false;
                }
            }
        }
        else {
            return false;
        }
        
        return true;
    }
    
}
