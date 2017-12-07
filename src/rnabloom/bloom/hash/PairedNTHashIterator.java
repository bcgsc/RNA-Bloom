/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import java.util.Arrays;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;
import static rnabloom.bloom.hash.NTHash.NTM64;
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
    protected int start = -1;
    protected int pos = -1;
    protected int max = -2;
    protected int d = 0;
    public long[] hVals1 = null;
    public long[] hVals2 = null;
    public long[] hVals3 = null;

    public PairedNTHashIterator(int k, int h, int d) {
        this.k = k;
        this.kMod64 = k%64;
        this.h = h;
        
        this.hVals1 = new long[h];
        this.hVals2 = new long[h];
        this.hVals3 = new long[h];
        this.d = d;
    }
    
    public boolean start(CharSequence seq) {
        return start(seq, 0, seq.length());
    }
    
    public boolean start(CharSequence seq, int start, int end) {
        this.seq = seq;
        this.start = start-1;
        this.pos = this.start;
        this.max = end - k - d;
        
        return max >= 0;
    }

    public void next() {
        if (pos == -1) {
            NTM64(seq, k, h, 0, hVals1);
            NTM64(seq, k, h, d, hVals2);
            NTM64(combineHashValues(hVals1[0], hVals2[0]), hVals3, k, h);
            ++pos;
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals1, kMod64);
            NTM64(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, hVals2, kMod64);
            NTM64(combineHashValues(hVals1[0], hVals2[0]), hVals3, k, h);
            ++pos;
        }
        else {
            Arrays.fill(hVals1, 0);
            Arrays.fill(hVals2, 0);
            Arrays.fill(hVals3, 0);
            hVals1 = null;
            hVals2 = null;
            hVals3 = null;
        }
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
