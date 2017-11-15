/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;

/**
 *
 * @author kmnip
 */
public class NTHashIterator {
    protected CharSequence seq;
    protected final int k;
    protected final int kMod64;
    protected final int h;
    protected int start = -1;
    protected int pos = -1;
    protected int max = -2;
    public long[] hVals = null;

    public NTHashIterator(int k, int h) {
        this.k = k;
        this.kMod64 = k%64;
        this.h = h;
        this.hVals = new long[h];
    }
    
    public boolean start(CharSequence seq) {
        return start(seq, 0, seq.length());
    }

    public boolean start(CharSequence seq, int start, int end) {
        this.seq = seq;
        this.start = start-1;
        this.pos = this.start;
        this.max = end - k;
        
        return max >= 0;
    }
    
    public void next() {
        if (pos == start) {
            // hash first kmer
            NTM64(seq, k, h, start+1, hVals);
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals, kMod64);
        }
        else {
            hVals = null;
        }
        ++pos;
    }

    public boolean hasNext() {
        return pos < max;
    }
    
    public int getPos() {
        return pos;
    }
    
    public int getMax() {
        return max;
    }
}
