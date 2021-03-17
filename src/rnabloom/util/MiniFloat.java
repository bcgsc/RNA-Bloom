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
package rnabloom.util;

import static java.lang.Math.random;
import static java.lang.Math.scalb;

/**
 * Based on https://github.com/bcgsc/abyss/blob/master/LogKmerCount/plc.h
 * @author Ka Ming Nip
 */
public class MiniFloat {
    private static final byte MANTISSA = 3;
    private static final byte MANTI_MASK = 0xFF >> (8 - MANTISSA);
    private static final byte ADD_MASK = 0x80 >> (7 - MANTISSA);
    
    public static byte increment(byte b) {
        if (b <= MANTI_MASK ||
                (b < Byte.MAX_VALUE &&
                    (int) (random() * Integer.MAX_VALUE) % (1 << ((b >> MANTISSA) - 1)) == 0)) {
            return (byte) (b + 1);
        }
        return b;
    }
    
    public static float toFloat(byte b) {
        if (b <= MANTI_MASK) {
            return (float) b;
        }
        return scalb((b & MANTI_MASK) | ADD_MASK, (b >> MANTISSA) - 1);
    }
}
