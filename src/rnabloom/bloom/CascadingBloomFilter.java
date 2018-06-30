/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom.bloom;

import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author Ka Ming Nip
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
