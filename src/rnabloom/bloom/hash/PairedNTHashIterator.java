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
        if (pos == start) {
            ++pos;
            NTM64(seq, k, h, pos, hVals1);
            NTM64(seq, k, h, pos+d, hVals2);
            NTM64(combineHashValues(hVals1[0], hVals2[0]), hVals3, k, h);
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
