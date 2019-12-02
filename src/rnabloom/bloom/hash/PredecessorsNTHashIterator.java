/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
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
import static rnabloom.bloom.hash.NTHash.msTab;
import static rnabloom.util.SeqUtils.NUCLEOTIDES_BYTES;

/**
 *
 * @author Ka Ming Nip
 */
public class PredecessorsNTHashIterator {
    protected int k;
    protected int i = -1;
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
    
    public boolean hasNext() {
        return i < 3;
    }
    
    public void start(final long fHashVal, final byte charOut) {
        tmpVal = Long.rotateRight(fHashVal, 1) ^ msTab[charOut][63];
        i = -1;
    }
    
    public void next() {
        NTM64(tmpVal ^ msTab[NUCLEOTIDES_BYTES[++i]][kMinus1Mod64], hVals, k, numHash);
    }
    
    public byte currentChar() {
        return NUCLEOTIDES_BYTES[i];
    }
}
