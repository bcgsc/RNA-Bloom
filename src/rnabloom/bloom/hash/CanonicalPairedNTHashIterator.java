/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.HashFunction2.combineHashValues;
import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.bloom.hash.NTHash.NTMC64;

/**
 *
 * @author gengar
 */
public class CanonicalPairedNTHashIterator extends PairedNTHashIterator {
    public long[] frVal1 = new long[2];
    public long[] frVal2 = new long[2];
    
    public CanonicalPairedNTHashIterator(int k, int h, int d) {
        super(k, h, d);
    }
    
    @Override
    public void next() {
        if (pos == -1) {
            NTMC64(seq, k, h, 0, frVal1, hVals1);
            NTMC64(seq, k, h, d, frVal2, hVals2);
            NTM64(Math.min(combineHashValues(frVal1[0], frVal2[0]), combineHashValues(frVal2[1], frVal1[1])), hVals3, k, h);
        }
        else if (pos < max) {
            NTMC64(seq.charAt(pos), seq.charAt(pos+k), k, h, frVal1, hVals1);
            NTMC64(seq.charAt(pos+d), seq.charAt(pos+k+d), k, h, frVal2, hVals2);
            NTM64(Math.min(combineHashValues(frVal1[0], frVal2[0]), combineHashValues(frVal2[1], frVal1[1])), hVals3, k, h);
        }
        else {
            hVals1 = null;
            hVals2 = null;
            hVals3 = null;
        }
        ++pos;
    }
}
