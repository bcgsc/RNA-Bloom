/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTMC64;

/**
 *
 * @author kmnip
 */
public class CanonicalNTHashIterator extends NTHashIterator {
    private final long[] frhval = new long[2];

    public CanonicalNTHashIterator(CharSequence seq, int k, int h) {
        super(seq, k, h);
    }

    @Override
    public void next() {
        if (pos == -1) {
            NTMC64(seq.subSequence(0, k), k, h, frhval, hVals);
        }
        else if (pos < max) {
            NTMC64(seq.charAt(pos), seq.charAt(pos+k), k, h, frhval, hVals);
        }
        else {
            hVals = null;
        }
        ++pos;
    }
}
