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

import java.util.Arrays;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;
import static rnabloom.bloom.hash.NTHash.NTM64;

/**
 *
 * @author Ka Ming Nip
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
    public long[] hValsL = null;
    public long[] hValsR = null;
    public long[] hValsP = null;

    public PairedNTHashIterator(int k, int h, int d) {
        this.k = k;
        this.kMod64 = k%64;
        this.h = h;
        
        this.hValsL = new long[h];
        this.hValsR = new long[h];
        this.hValsP = new long[h];
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
        if (pos == start) {
            ++pos;
            NTM64(seq, k, h, pos, hValsL);
            NTM64(seq, k, h, pos+d, hValsR);
            NTM64(combineHashValues(hValsL[0], hValsR[0]), hValsP, k, h);
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hValsL, kMod64);
            NTM64(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, hValsR, kMod64);
            NTM64(combineHashValues(hValsL[0], hValsR[0]), hValsP, k, h);
            ++pos;
        }
        else {
            Arrays.fill(hValsL, 0);
            Arrays.fill(hValsR, 0);
            Arrays.fill(hValsP, 0);
            hValsL = null;
            hValsR = null;
            hValsP = null;
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
