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

import static rnabloom.bloom.hash.NTHash.NTM64RC;

/**
 *
 * @author Ka Ming Nip
 */
public class ReverseComplementNTHashIterator extends NTHashIterator {
    public ReverseComplementNTHashIterator(int k, int h) {
        super(k, h);
    }

    @Override
    public void next() {
        if (pos == start) {
            NTM64RC(seq, k, h, ++pos, hVals);
        }
        else if (pos < max) {
            NTM64RC(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals);
            ++pos;
        }
        else {
            hVals = null;
        }
    }
}
