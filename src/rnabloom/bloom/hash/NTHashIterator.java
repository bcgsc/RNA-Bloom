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
import static rnabloom.bloom.hash.NTHash.NTM64;

/**
 *
 * @author Ka Ming Nip
 */
public class NTHashIterator {
    protected CharSequence seq;
    protected final int k;
    protected final int kMod64;
    protected final int h;
    protected int start = -1;
    protected int pos = -1;
    protected int max = -2;
    public long[] hVals = null;

    public NTHashIterator(int k, int h) {
        this.k = k;
        this.kMod64 = k%64;
        this.h = h;
        this.hVals = new long[h];
    }
    
    public boolean start(CharSequence seq) {
        return start(seq, 0, seq.length());
    }

    public boolean start(CharSequence seq, int start, int end) {
        this.seq = seq;
        this.start = start-1;
        this.pos = this.start;
        this.max = end - k;
        
        return max >= 0;
    }
    
    public void next() {
        if (pos == start) {
            // hash first kmer
            NTM64(seq, k, h, ++pos, hVals);
        }
        else if (pos < max) {
            NTM64(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals, kMod64);
            ++pos;
        }
        else {
            Arrays.fill(hVals, 0);
            hVals = null;
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
