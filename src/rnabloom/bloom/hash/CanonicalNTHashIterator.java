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
import static rnabloom.bloom.hash.NTHash.NTMC64;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalNTHashIterator extends NTHashIterator {
    public final long[] frhval = new long[2];
    private final int kMinus1Mod64;

    public CanonicalNTHashIterator(int k, int h) {
        super(k, h);
        kMinus1Mod64 = (k-1)%64;
    }

    @Override
    public void next() {
        if (pos == start) {
            NTMC64(seq, k, h, ++pos, frhval, hVals);
        }
        else if (pos < max) {
            NTMC64(seq.charAt(pos), seq.charAt(pos+k), k, h, frhval, hVals, kMod64, kMinus1Mod64);
            ++pos;
        }
        else {
            Arrays.fill(hVals, 0);
            hVals = null;
        }
    }
}
