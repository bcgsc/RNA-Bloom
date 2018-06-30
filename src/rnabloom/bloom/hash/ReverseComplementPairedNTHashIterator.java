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
import static rnabloom.bloom.hash.NTHash.NTM64RC;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author Ka Ming Nip
 */
public class ReverseComplementPairedNTHashIterator extends PairedNTHashIterator {

    public ReverseComplementPairedNTHashIterator(int k, int h, int d) {
        super(k, h, d);
    }
    
    @Override
    public void next() {
        if (pos == start) {
            ++pos;
            NTM64RC(seq, k, h, pos, hVals1);
            NTM64RC(seq, k, h, pos+d, hVals2);
            NTM64(combineHashValues(hVals2[0], hVals1[0]), hVals3, k, h);
        }
        else if (pos < max) {
            NTM64RC(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals1);
            NTM64RC(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, hVals2);
            NTM64(combineHashValues(hVals2[0], hVals1[0]), hVals3, k, h);
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
