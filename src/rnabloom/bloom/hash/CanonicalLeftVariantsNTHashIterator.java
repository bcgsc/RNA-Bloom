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
import static rnabloom.bloom.hash.NTHash.cpOff;
import static rnabloom.bloom.hash.NTHash.msTab;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalLeftVariantsNTHashIterator {
    protected int k;
    protected int kMinus1Mod64;
    protected long tmpValF, tmpValR;
    public long fHashVal, rHashVal;
    protected int numHash;
    public long[] hVals;
    
    public CanonicalLeftVariantsNTHashIterator(final int k, final int numHash) {
        this.k = k;
        this.kMinus1Mod64 = (k-1)%64;
        this.numHash = numHash;
        this.hVals = new long[numHash];
    }
    
    public void start(final long fHashVal, final long rHashVal, final byte charOut) {
        tmpValF = fHashVal ^ msTab[charOut][kMinus1Mod64];
        tmpValR = rHashVal ^ msTab[charOut&cpOff][0];
    }
    
    public void next(final byte charIn) {        
        fHashVal = tmpValF ^ msTab[charIn][kMinus1Mod64];
        rHashVal = tmpValR ^ msTab[charIn&cpOff][0];
        
        NTM64(Math.min(fHashVal, rHashVal), hVals, k, numHash);
    }
}
