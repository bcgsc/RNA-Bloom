/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.BitSet;

/**
 *
 * @author kmnip
 */
public class KmerBitsUtils {
    private static final float LOW_COMPLEXITY_THRESHOLD = 0.87f;
    
    private int k;
    private int numBits;
    private int lastBitIndex, secondLastBitIndex;
    private int t1, t2, t3;
    
    
    public KmerBitsUtils(int k) {
        this.k = k;
        this.numBits = k * 2;
        
        this.lastBitIndex = this.numBits - 1;
        this.secondLastBitIndex = this.numBits - 2;
        
        this.t1 = Math.round(k * LOW_COMPLEXITY_THRESHOLD);
        this.t2 = Math.round(k/2 * LOW_COMPLEXITY_THRESHOLD);
        this.t3 = Math.round(k/3 * LOW_COMPLEXITY_THRESHOLD);
    }
    
    public BitSet seqToBits(String seq) {
        BitSet bits = new BitSet(numBits);
        
        for (int i=0; i<k; ++i) {
            switch (seq.charAt(i)) {
                case 'A':
                    // 00
                    break;
                case 'C':
                    // 01
                    bits.set(2*i + 1);
                    break;
                case 'G':
                    // 10
                    bits.set(2*i);
                    break;
                case 'T':
                    // 11
                    bits.set(2*i, 2*i + 2);
                    break;
            }
        }
        
        return bits;
    }
    
    public String bitsToSeq(BitSet bits) {
        StringBuilder sb = new StringBuilder(k);
        
        for (int i=0; i<numBits; ++i) {
            if (bits.get(i)) {
                if (bits.get(++i)) {
                    // 11
                    sb.append('T');
                }
                else {
                    // 10
                    sb.append('G');
                }
            }
            else {
                if (bits.get(++i)) {
                    // 01
                    sb.append('C');
                }
                else {
                    // 00
                    sb.append('A');
                }
            }
        }
        
        return sb.toString();
    }
    
    private int getBaseRank(BitSet bits, int index) {
        index *= 2; // bit index
        if (bits.get(index)) {
            if (bits.get(index+1)) {
                // 11
                return 3;
            }
            else {
                // 10
                return 2;
            }
        }
        else {
            if (bits.get(index+1)) {
                // 01
                return 1;
            }
            else {
                // 00
                return 0;
            }
        }
    }
    
    public char getBase(BitSet bits, int index) {
        index *= 2; // bit index
        if (bits.get(index)) {
            if (bits.get(index+1)) {
                // 11
                return 'T';
            }
            else {
                // 10
                return 'G';
            }
        }
        else {
            if (bits.get(index+1)) {
                // 01
                return 'C';
            }
            else {
                // 00
                return 'A';
            }
        }
    }
    
    public char getFirstBase(BitSet bits) {
        if (bits.get(0)) {
            if (bits.get(1)) {
                // 11
                return 'T';
            }
            else {
                // 10
                return 'G';
            }
        }
        else {
            if (bits.get(1)) {
                // 01
                return 'C';
            }
            else {
                // 00
                return 'A';
            }
        }
    }
    
    public char getLastBase(BitSet bits) {
        if (bits.get(secondLastBitIndex)) {
            if (bits.get(lastBitIndex)) {
                // 11
                return 'T';
            }
            else {
                // 10
                return 'G';
            }
        }
        else {
            if (bits.get(lastBitIndex)) {
                // 01
                return 'C';
            }
            else {
                // 00
                return 'A';
            }
        }
    }
    
    public void setFirstBase(BitSet bits, char c) {
        switch (c) {
            case 'A':
                // 00
                bits.set(0, 2, false);
                break;
            case 'C':
                // 01
                bits.set(0, false);
                bits.set(1);
                break;
            case 'G':
                // 10
                bits.set(0);
                bits.set(1, false);
                break;
            case 'T':
                // 11
                bits.set(0, 2);
                break;
        }
    }
    
    public void setLastBase(BitSet bits, char c) {
        switch (c) {
            case 'A':
                // 00
                bits.set(secondLastBitIndex, numBits, false);
                break;
            case 'C':
                // 01
                bits.set(secondLastBitIndex, false);
                bits.set(lastBitIndex);
                break;
            case 'G':
                // 10
                bits.set(secondLastBitIndex);
                bits.set(lastBitIndex, false);
                break;
            case 'T':
                // 11
                bits.set(secondLastBitIndex, numBits);
                break;
        }
    }
    
    public BitSet shiftLeft(BitSet bits) {
        BitSet newBits = new BitSet(numBits);
        
        for (int i=2; i<numBits; ++i) {
            if (bits.get(i)) {
                newBits.set(i-2);
            }
        }
        
        return newBits;
    }
    
    public BitSet shiftRight(BitSet bits) {
        BitSet newBits = new BitSet(numBits);
        
        for (int i=0; i<secondLastBitIndex; ++i) {
            if (bits.get(i)) {
                newBits.set(i+2);
            }
        }
        
        return newBits;
    }
        
    public boolean isLowComplexity(BitSet bits) {
        byte nf1[]     = new byte[4];
        byte nf2[][]   = new byte[4][4];
        byte nf3[][][] = new byte[4][4][4];
        
        int c3 = getBaseRank(bits, 0);
        int c2 = getBaseRank(bits, 1);
        int c1 = getBaseRank(bits, 2);
        
        ++nf1[c3];
        ++nf1[c2];
        ++nf1[c1];
        
        ++nf2[c3][c2];
        ++nf2[c2][c1];
        
        ++nf3[c3][c2][c1];
                
        for (int i=3; i<k; ++i) {
            c3 = c2;
            c2 = c1;
            c1 = getBaseRank(bits, i);
            
            ++nf1[c1];
            ++nf2[c2][c1];
            ++nf3[c3][c2][c1];
        }
        
        // homopolymer runs
        for (byte n : nf1) {
            if (n >= t1) {
                return true;
            }
        }
        
        // di-nucleotide content
        if (nf1[0]+nf1[1]>t1 || nf1[0]+nf1[2]>t1 || nf1[0]+nf1[3]>t1 || 
                nf1[1]+nf1[2]>t1 || nf1[1]+nf1[3]>t1 || nf1[2]+nf1[3]>t1) {
            return true;
        }
        
        // di-nucleotide repeat
        for (byte[] n1 : nf2) {
            for (byte n : n1) {
                if (n >= t2) {
                    return true;
                }
            }
        }
        
        // tri-nucleotide repeat
        for (byte[][] n2 : nf3) {
            for (byte[] n1 : n2) {
                for (byte n : n1) {
                    if (n >= t3) {
                        return true;
                    }
                }
            }
        }
        
        return false;
    }
}
