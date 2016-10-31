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
    protected final int h;
    protected int pos;
    protected int max;
    public long[] hVals = null;

    public NTHashIterator(CharSequence seq, int k, int h) {
        this.seq = seq;
        this.k = k;
        this.h = h;
        this.hVals = new long[h];
        this.pos = -1;
        this.max = seq.length() - k;
    }

    public void setSeq(CharSequence seq) {
        this.seq = seq;
        this.pos = -1;
        this.max = seq.length() - k;
    }

    public void next() {
        if (pos == -1) {
            NTM64(seq.subSequence(0, k), k, h, hVals);
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals);
        }
        else {
            hVals = null;
        }
        ++pos;
    }

    public boolean hasNext() {
        return pos < max;
    }
}
