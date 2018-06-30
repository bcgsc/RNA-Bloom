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
import static rnabloom.bloom.hash.NTHash.NTMC64;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalPairedNTHashIterator extends PairedNTHashIterator {
    public long[] frVal1 = new long[2];
    public long[] frVal2 = new long[2];
    
    public CanonicalPairedNTHashIterator(int k, int h, int d) {
        super(k, h, d);
    }
    
    @Override
    public void next() {
        if (pos == start) {
            ++pos;
            NTMC64(seq, k, h, pos, frVal1, hVals1);
            NTMC64(seq, k, h, pos+d, frVal2, hVals2);
            NTM64(Math.min(combineHashValues(frVal1[0], frVal2[0]), combineHashValues(frVal2[1], frVal1[1])), hVals3, k, h);
        }
        else if (pos < max) {
            NTMC64(seq.charAt(pos), seq.charAt(pos+k), k, h, frVal1, hVals1);
            NTMC64(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, frVal2, hVals2);
            NTM64(Math.min(combineHashValues(frVal1[0], frVal2[0]), combineHashValues(frVal2[1], frVal1[1])), hVals3, k, h);
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
}
