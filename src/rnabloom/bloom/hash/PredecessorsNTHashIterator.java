/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.msTab;
import static rnabloom.bloom.hash.NTHash.multiSeed;
import static rnabloom.bloom.hash.NTHash.multiShift;

/**
 *
 * @author gengar
 */
public class PredecessorsNTHashIterator {
    protected int k;
    protected int kMinus1Mod64;
    protected int h;
    public long[] hVals = null;
    protected long tmpVal;
    
    public PredecessorsNTHashIterator(final int k, final int h) {
        this.k = k;
        this.kMinus1Mod64 = (k-1)%64;
        this.h = h;
        this.hVals = new long[h];
    }
    
    public void start(final long[] hVals, char charOut) {
        tmpVal = Long.rotateRight(hVals[0], 1) ^ msTab[charOut][63];
    }
    
    public void next(char charIn) {
        long bVal = tmpVal ^ msTab[charIn][kMinus1Mod64];
        hVals[0] = bVal;
        long tVal;
        for(int i=1; i<h; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVals[i] = tVal;
        }
    }
}
