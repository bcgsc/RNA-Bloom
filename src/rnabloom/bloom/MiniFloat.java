/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

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
