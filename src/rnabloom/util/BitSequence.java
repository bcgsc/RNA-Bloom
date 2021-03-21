/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.Arrays;
import static rnabloom.util.SeqBitsUtils.bitsToRevCompSeq;
import static rnabloom.util.SeqBitsUtils.bitsToSeq;
import static rnabloom.util.SeqBitsUtils.seqToBits;
import static rnabloom.util.SeqBitsUtils.seqToBitsParallelized;

/**
 *
 * @author Ka Ming Nip
 */
public class BitSequence implements Comparable<Object> {
    public int length;
    private byte[] bits;

    public BitSequence(String seq) {
        this.length = seq.length();
        this.bits = seqToBitsParallelized(seq);
    }

    @Override
    public String toString() {
        return bitsToSeq(this.bits, this.length);
    }
    
    public String subString(int start, int end) {
        return bitsToSeq(this.bits, this.length, start, end);
    }

    public String subStringRevComp(int start, int end) {
        return bitsToRevCompSeq(this.bits, this.length, start, end);
    }
    
    public boolean equals(BitSequence other) {
        return Arrays.equals(bits, other.bits);
    }
    
    @Override
    public int compareTo(Object other) {
        return ((BitSequence) other).length - this.length;
    }
}
