/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import java.util.Arrays;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;
import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.bloom.hash.NTHash.NTM64RC;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author kmnip
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
