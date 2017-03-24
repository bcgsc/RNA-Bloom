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
public class SuccessorsNTHashIterator {
    protected int k;
    protected int h;
    public long[] hVals = null;
    protected long tmpVal;
    
    public SuccessorsNTHashIterator(final int k, final int h) {
        this.k = k;
        this.h = h;
        this.hVals = new long[h];
    }
    
    public void start(final long[] hVals, char charOut) {
        tmpVal = Long.rotateLeft(hVals[0], 1) ^ msTab[charOut][k%64];
    }
    
    public void next(char charIn) {
        long bVal = tmpVal ^ msTab[charIn][0];
        hVals[0] = bVal;
        long tVal;
        for(int i=1; i<h; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVals[i] = tVal;
        }
    }
}
