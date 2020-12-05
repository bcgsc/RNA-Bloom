/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;

/**
 *
 * @author Ka Ming Nip
 */
public class SimpleBloomFilter implements BloomFilterInterface {
    protected AbstractLargeBitBuffer bitArray;
    protected long size;
    protected long popcount = -1;
    
    public SimpleBloomFilter(long size) {
        this.size = size;
        try {
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
    }
    
    protected long getIndex(String key) {
        long hashCode = (long) key.hashCode();
        return (hashCode + (long) Integer.MAX_VALUE + 1L) % size;
    }
    
    @Override
    public void add(String key) {
        bitArray.set(getIndex(key));
    }

    @Override
    public boolean lookup(String key) {
        return bitArray.get(getIndex(key));
    }
    
    public boolean lookupAndAdd(String key) {
        return bitArray.getAndSet(getIndex(key));
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        popcount = bitArray.popCount();
        return (float) ((double) popcount / (double) size);
    }
    
    public void empty() {
        if (this.bitArray != null) {
            this.bitArray.empty();
        }
    }
    
    public void destroy() {
        if (this.bitArray != null) {
            this.bitArray.destroy();
            this.bitArray = null;
        }
    }
}
