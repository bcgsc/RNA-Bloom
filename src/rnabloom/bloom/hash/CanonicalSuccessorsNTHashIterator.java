/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.bloom.hash.NTHash.cpOff;
import static rnabloom.bloom.hash.NTHash.msTab;
import static rnabloom.util.SeqUtils.NUCLEOTIDES_BYTES;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalSuccessorsNTHashIterator {
    protected int k;
    protected int kMinus1Mod64;
    protected int i = -1;
    protected long tmpValF, tmpValR;
    public long fHashVal, rHashVal;
    protected int numHash;
    public long[] hVals;
    
    public CanonicalSuccessorsNTHashIterator(final int k, final int numHash) {
        this.k = k;
        this.kMinus1Mod64 = (k-1)%64;
        this.numHash = numHash;
        this.hVals = new long[numHash];
    }
    
    public boolean hasNext() {
        return i < 3;
    }
    
    public void start(final long fHashVal, final long rHashVal, final byte charOut) {
        tmpValF = Long.rotateLeft(fHashVal, 1) ^ msTab[charOut][k%64];
        tmpValR = Long.rotateRight(rHashVal, 1) ^ msTab[charOut&cpOff][63];
        i = -1;
    }
    
    public void next() {
        byte charIn = NUCLEOTIDES_BYTES[++i];
        
        fHashVal = tmpValF ^ msTab[charIn][0];
        rHashVal = tmpValR ^ msTab[charIn&cpOff][kMinus1Mod64];
        
        NTM64(Math.min(fHashVal, rHashVal), hVals, k, numHash);
    }
    
    public byte currentChar() {
        return NUCLEOTIDES_BYTES[i];
    }
}
