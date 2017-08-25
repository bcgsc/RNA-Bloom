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
public class SuccessorsNTHashIterator {
    protected int k;
    protected long tmpVal;
    protected int numHash;
    public long[] hVals;
    
    public SuccessorsNTHashIterator(final int k, final int numHash) {
        this.k = k;
        this.numHash = numHash;
        this.hVals = new long[numHash];
    }
    
    public void start(final long fHashVal, final char charOut) {
        tmpVal = Long.rotateLeft(fHashVal, 1) ^ msTab[charOut][k%64];
    }
    
    public void next(final char charIn) {
        NTM64(tmpVal ^ msTab[charIn][0], hVals, k, numHash);
    }
}
