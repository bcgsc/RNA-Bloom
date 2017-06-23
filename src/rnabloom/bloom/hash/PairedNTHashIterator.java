/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;

/**
 *
 * @author gengar
 */
public class PairedNTHashIterator {
    protected CharSequence seq;
    protected final int k;
    protected final int kMod64;
    protected final int h;
    protected int pos = -1;
    protected int max = -2;
    protected int d = 0;
    public long[] hVals1 = null;
    public long[] hVals2 = null;

    public PairedNTHashIterator(int k, int h, int d) {
        this.k = k;
        this.kMod64 = k%64;
        this.h = h;
        
        this.hVals1 = new long[h];
        this.hVals2 = new long[h];
        this.d = d;
    }
    
    public void start(CharSequence seq) {
        this.seq = seq;
        this.pos = -1;
        this.max = seq.length() - k - d;
    }

    public void next() {
        if (pos == -1) {
            NTM64(seq.subSequence(0, k), k, h, hVals1);
            NTM64(seq.subSequence(d, d+k), k, h, hVals2);
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals1, kMod64);
            NTM64(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, hVals2, kMod64);
        }
        else {
            hVals1 = null;
            hVals2 = null;
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
