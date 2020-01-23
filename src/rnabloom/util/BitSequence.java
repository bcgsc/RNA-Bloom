/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.BitSet;
import static rnabloom.util.KmerBitsUtils.bitsToSeq;
import static rnabloom.util.KmerBitsUtils.seqToBits;

/**
 *
 * @author kmnip
 */
public class BitSequence implements Comparable<Object> {
    int length;
    BitSet bits;

    public BitSequence(String seq) {
        this.length = seq.length();
        this.bits = seqToBits(seq);
    }

    @Override
    public String toString() {
        return bitsToSeq(this.bits, this.length, false);
    }

    public String toString(boolean useUracil) {
        return bitsToSeq(this.bits, this.length, useUracil);
    }
    
    public boolean equals(BitSequence other) {
        return bits.equals(other.bits);
    }
    
    @Override
    public int compareTo(Object other) {
        return ((BitSequence) other).length - this.length;
    }
}
