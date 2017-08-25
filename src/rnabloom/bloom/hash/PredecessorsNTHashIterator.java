/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.bloom.hash.NTHash.msTab;

/**
 *
 * @author Ka Ming Nip
 */
public class PredecessorsNTHashIterator {
    protected int k;
    protected int kMinus1Mod64;
    protected long tmpVal;
    protected int numHash;
    public long[] hVals;
    
    public PredecessorsNTHashIterator(final int k, final int numHash) {
        this.k = k;
        this.kMinus1Mod64 = (k-1)%64;
        this.numHash = numHash;
        this.hVals = new long[numHash];
    }
    
    public void start(final long fHashVal, final char charOut) {
        tmpVal = Long.rotateRight(fHashVal, 1) ^ msTab[charOut][63];
    }
    
    public void next(final char charIn) {
        NTM64(tmpVal ^ msTab[charIn][kMinus1Mod64], hVals, k, numHash);
    }
}
