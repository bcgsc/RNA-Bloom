/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64RC;
import static rnabloom.bloom.hash.NTHash.NTM64RC;

/**
 *
 * @author kmnip
 */
public class ReverseComplementNTHashIterator extends NTHashIterator {
    public ReverseComplementNTHashIterator(int k, int h) {
        super(k, h);
    }

    @Override
    public void next() {
        if (pos == -1) {
            NTM64RC(seq.subSequence(0, k), k, h, hVals);
        }
        else if (pos >= 0) {
            NTM64RC(seq.charAt(pos), seq.charAt(pos+k), k, h, hVals);
        }
        else {
            hVals = null;
        }
        ++pos;
    }
}
