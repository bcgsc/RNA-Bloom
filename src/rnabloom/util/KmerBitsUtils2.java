/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
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

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author Ka Ming Nip
 */
public class KmerBitsUtils2 {
    private static final float LOW_COMPLEXITY_THRESHOLD = 0.87f;
    private static final int NUM_BITS_PER_BASE = 2;
    private static final int NUM_BASES_PER_BYTE = Byte.SIZE / NUM_BITS_PER_BASE;
    
    private final int k;
    private final int lastByteIndex, numBytes, numFullBytes, numRemainingBases, numRemainingBits;
    private final int lastBaseBitIndex1, lastBaseBitIndex2;
    private final int t1, t2, t3;
    
    private static final byte[] BIT_MASKS = new byte[]{(byte) 1, (byte) 2, (byte)(1<<2), (byte)(1<<3), (byte)(1<<4), (byte)(1<<5), (byte)(1<<6), (byte)(1<<7)};
    private static final byte[] NOT_BIT_MASKS = new byte[]{~(byte) 1, ~(byte) 2, ~(byte)(1<<2), ~(byte)(1<<3), ~(byte)(1<<4), ~(byte)(1<<5), ~(byte)(1<<6), ~(byte)(1<<7)};

    private static final byte[][] BYTE_BASE_MASKS = getByteBaseMasks();
    
    private final byte lastBaseBitMask1, lastBaseBitMask2, lastBaseNotBitMask1, lastBaseNotBitMask2;

    
    public KmerBitsUtils2(int k) {
        this.k = k;
        this.numFullBytes = k / NUM_BASES_PER_BYTE;
        this.numRemainingBases = k % NUM_BASES_PER_BYTE;
        this.numRemainingBits = NUM_BITS_PER_BASE * numRemainingBases;
        
        if (numRemainingBases > 0) {
            this.numBytes = numFullBytes + 1;
            this.lastByteIndex = numFullBytes;
            this.lastBaseBitIndex2 = numRemainingBits - 1;
            this.lastBaseBitIndex1 = lastBaseBitIndex2 - 1;
        }
        else {
            this.numBytes = numFullBytes;
            this.lastByteIndex = numFullBytes - 1;
            this.lastBaseBitIndex1 = 6;
            this.lastBaseBitIndex2 = 7;
        }
        
        this.lastBaseBitMask1 = BIT_MASKS[lastBaseBitIndex1];
        this.lastBaseBitMask2 = BIT_MASKS[lastBaseBitIndex2];
        
        this.lastBaseNotBitMask1 = NOT_BIT_MASKS[lastBaseBitIndex1];
        this.lastBaseNotBitMask2 = NOT_BIT_MASKS[lastBaseBitIndex2];
        
        this.t1 = Math.round(k * LOW_COMPLEXITY_THRESHOLD);
        this.t2 = Math.round(k/2 * LOW_COMPLEXITY_THRESHOLD);
        this.t3 = Math.round(k/3 * LOW_COMPLEXITY_THRESHOLD);
    }
    
    public static byte[][] getByteBaseMasks() {
        byte[][] baseMasks = new byte[NUM_BASES_PER_BYTE][4];
        
        for (int i=0; i<NUM_BASES_PER_BYTE; ++i) {
            // A
            // do nothing
            
            // C
            baseMasks[i][1] |= BIT_MASKS[i*2+1];
            
            // G
            baseMasks[i][2] |= BIT_MASKS[i*2];

            // T
            baseMasks[i][3] |= BIT_MASKS[i*2];
            baseMasks[i][3] |= BIT_MASKS[i*2+1];
        }
        
        return baseMasks;
    }
    
    public byte[] kmerToBits(String kmer) {
        
        byte[] bits = new byte[numBytes];
        
        int byteIndex, baseIndex;
        for (int i=0; i<k; ++i) {
            byteIndex = i * NUM_BITS_PER_BASE / Byte.SIZE;
            baseIndex = i % NUM_BASES_PER_BYTE;
            
            switch (kmer.charAt(i)) {
                // A: do nothing
                case 'C':
                    // 01
                    bits[byteIndex] |= BYTE_BASE_MASKS[baseIndex][1];
                    break;
                case 'G':
                    // 10
                    bits[byteIndex] |= BYTE_BASE_MASKS[baseIndex][2];
                    break;
                case 'T':
                    // 11
                    bits[byteIndex] |= BYTE_BASE_MASKS[baseIndex][3];
                    break;
            }
        }
        
        return bits;
    }
    
    public byte[] seqToBits(String seq, int start) {
        // get bits for kmer at start
        byte[] bits = new byte[numBytes];
        
        for (int i=0; i<numFullBytes; ++i) {
            byte b = 0;
            for (int j=0; j<NUM_BASES_PER_BYTE; ++j) {
                
                switch (seq.charAt(start++)) {
                    // A: do nothing
                    case 'C':
                        // 01
                        b |= BYTE_BASE_MASKS[j][1];
                        break;
                    case 'G':
                        // 10
                        b |= BYTE_BASE_MASKS[j][2];
                        break;
                    case 'T':
                        // 11
                        b |= BYTE_BASE_MASKS[j][3];
                        break;
                }
            }
            bits[i] = b;
        }
        
        if (numRemainingBases > 0) {
            byte b = 0;
            for (int j=0; j<numRemainingBases; ++j) {
                
                switch (seq.charAt(start++)) {
                    // A: do nothing
                    case 'C':
                        // 01
                        b |= BYTE_BASE_MASKS[j][1];
                        break;
                    case 'G':
                        // 10
                        b |= BYTE_BASE_MASKS[j][2];
                        break;
                    case 'T':
                        // 11
                        b |= BYTE_BASE_MASKS[j][3];
                        break;
                }
            }
            bits[numFullBytes] = b;
        }
        
        return bits;
    }
        
    public String bitsToKmer(byte[] bits) {
        StringBuilder sb = new StringBuilder(k);
        
        for (int i=0; i<numFullBytes; ++i) {
            byte b = bits[i];
            
            for (int j=0; j<Byte.SIZE; ++j) {
                if ((b & BIT_MASKS[j]) == 0) {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 00
                        sb.append('A');
                    }
                    else {
                        // 01
                        sb.append('C');
                    }
                }
                else {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 10
                        sb.append('G');
                    }
                    else {
                        // 11
                        sb.append('T');
                    }
                }
            }
        }
        
        if (numRemainingBases > 0) {
            byte b = bits[numFullBytes]; // get last byte
            for (int j=0; j<numRemainingBits; ++j) {
                if ((b & BIT_MASKS[j]) == 0) {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 00
                        sb.append('A');
                    }
                    else {
                        // 01
                        sb.append('C');
                    }
                }
                else {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 10
                        sb.append('G');
                    }
                    else {
                        // 11
                        sb.append('T');
                    }
                }
            }
        }
        
        return sb.toString();
    }
    
    public int[] bitsToRanks(byte[] bits) {
        int[] ranks = new int[k];
        
        int r = 0;
        
        for (int i=0; i<numFullBytes; ++i) {
            byte b = bits[i];
            
            for (int j=0; j<Byte.SIZE; ++j) {
                if ((b & BIT_MASKS[j]) == 0) {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 00
                        ranks[r++] = 0;
                    }
                    else {
                        // 01
                        ranks[r++] = 1;
                    }
                }
                else {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 10
                        ranks[r++] = 2;
                    }
                    else {
                        // 11
                        ranks[r++] = 3;
                    }
                }
            }
        }
        
        if (numRemainingBases > 0) {
            byte b = bits[numFullBytes]; // get last byte
            for (int j=0; j<numRemainingBits; ++j) {
                if ((b & BIT_MASKS[j]) == 0) {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 00
                        ranks[r++] = 0;
                    }
                    else {
                        // 01
                        ranks[r++] = 1;
                    }
                }
                else {
                    if ((b & BIT_MASKS[++j]) == 0) {
                        // 10
                        ranks[r++] = 2;
                    }
                    else {
                        // 11
                        ranks[r++] = 3;
                    }
                }
            }
        }
        
        return ranks;
    }
    
    public char getBase(byte[] bits, int index) {
        byte b = bits[index / NUM_BASES_PER_BYTE];
        int bitIndex = (index * NUM_BITS_PER_BASE) % Byte.SIZE;
        if ((b & BIT_MASKS[bitIndex]) == 0) {
            if ((b & BIT_MASKS[bitIndex+1]) == 0) {
                // 00
                return 'A';
            }
            else {
                // 01
                return 'C';
            }
        }
        else {
            if ((b & BIT_MASKS[bitIndex+1]) == 0) {
                // 10
                return 'G';
            }
            else {
                // 11
                return 'T';
            }
        }
    }
    
    public char getFirstBase(byte[] bits) {
        byte b = bits[0];
        if ((b & 1) == 0) {
            if ((b & 2) == 0) {
                // 00
                return 'A';
            }
            else {
                // 01
                return 'C';
            }
        }
        else {
            if ((b & 2) == 0) {
                // 10
                return 'G';
            }
            else {
                // 11
                return 'T';
            }
        }
    }
    
    public char getLastBase(byte[] bits) {
        byte b = bits[lastByteIndex];
        if ((b & lastBaseBitMask1) == 0) {
            if ((b & lastBaseBitMask2) == 0) {
                // 00
                return 'A';
            }
            else {
                // 01
                return 'C';
            }
        }
        else {
            if ((b & lastBaseBitMask2) == 0) {
                // 10
                return 'G';
            }
            else {
                // 11
                return 'T';
            }
        }
    }

    public void setBase(byte[] bits, int index, char c) {
        int byteIndex = index / NUM_BASES_PER_BYTE;
        byte b = bits[byteIndex];
        int bitIndex = (index * NUM_BITS_PER_BASE) % Byte.SIZE;
        
        switch (c) {
            case 'A':
                // 00
                b &= NOT_BIT_MASKS[bitIndex];
                b &= NOT_BIT_MASKS[bitIndex+1];
                break;
            case 'C':
                // 01
                b &= NOT_BIT_MASKS[bitIndex];
                b |= BIT_MASKS[bitIndex+1];
                break;
            case 'G':
                // 10
                b |= BIT_MASKS[bitIndex];
                b &= NOT_BIT_MASKS[bitIndex+1];
                break;
            case 'T':
                // 11
                b |= BIT_MASKS[bitIndex];
                b |= BIT_MASKS[bitIndex+1];
                break;
        }
        
        bits[byteIndex] = b;
    }
    
    public void setFirstBase(byte[] bits, char c) {
        byte b = bits[0];
        
        switch (c) {
            case 'A':
                // 00
                b &= NOT_BIT_MASKS[0];
                b &= NOT_BIT_MASKS[1];
                break;
            case 'C':
                // 01
                b &= NOT_BIT_MASKS[0];
                b |= BIT_MASKS[1];
                break;
            case 'G':
                // 10
                b |= BIT_MASKS[0];
                b &= NOT_BIT_MASKS[1];
                break;
            case 'T':
                // 11
                b |= BIT_MASKS[0];
                b |= BIT_MASKS[1];
                break;
        }
        
        bits[0] = b;
    }
    
    public void setLastBase(byte[] bits, char c) {
        byte b = bits[lastByteIndex];
        
        switch (c) {
            case 'A':
                // 00
                b &= lastBaseNotBitMask1;
                b &= lastBaseNotBitMask2;
                break;
            case 'C':
                // 01
                b &= lastBaseNotBitMask1;
                b |= lastBaseBitMask2;
                break;
            case 'G':
                // 10
                b |= lastBaseBitMask1;
                b &= lastBaseNotBitMask2;
                break;
            case 'T':
                // 11
                b |= lastBaseBitMask1;
                b |= lastBaseBitMask2;
                break;
        }
        
        bits[lastByteIndex] = b;
    }
    
    public byte[] shiftLeft(byte[] bits) {
        byte[] newBits = new byte[numBytes];
        
        for (int i=0; i<lastByteIndex; ++i) {
            byte b = bits[i];
            byte b2 = bits[i+1];
            
            b >>= 2;
            
            if ((b2 & 1) == 0) {
                if ((b2 & 2) == 0) {
                    // 00
                    b &= NOT_BIT_MASKS[6];
                    b &= NOT_BIT_MASKS[7];
                }
                else {
                    // 01
                    b &= NOT_BIT_MASKS[6];
                    b |= BIT_MASKS[7];
                }
            }
            else {
                if ((b2 & 2) == 0) {
                    // 10
                    b |= BIT_MASKS[6];
                    b &= NOT_BIT_MASKS[7];
                }
                else {
                    // 11
                    b |= BIT_MASKS[6];
                    b |= BIT_MASKS[7];
                }
            }
            
            newBits[i] = b;
        }
        
        newBits[lastByteIndex] = (byte) (bits[lastByteIndex] >> 2);
        
        return newBits;
    }
    
    public byte[] shiftRight(byte[] bits) {
        byte[] newBits = new byte[numBytes];
        
        for (int i=lastByteIndex; i>0; --i) {
            byte b = bits[i];
            byte b2 = bits[i-1];
            
            b <<= 2;
            
            if ((b2 & BIT_MASKS[6]) == 0) {
                if ((b2 & BIT_MASKS[7]) == 0) {
                    // 00
                    b &= NOT_BIT_MASKS[0];
                    b &= NOT_BIT_MASKS[1];
                }
                else {
                    // 01
                    b &= NOT_BIT_MASKS[0];
                    b |= BIT_MASKS[1];
                }
            }
            else {
                if ((b2 & BIT_MASKS[7]) == 0) {
                    // 10
                    b |= BIT_MASKS[0];
                    b &= NOT_BIT_MASKS[1];
                }
                else {
                    // 11
                    b |= BIT_MASKS[0];
                    b |= BIT_MASKS[1];
                }
            }
            
            newBits[i] = b;
        }
        
        newBits[0] = (byte) (bits[0] << 2);
        
        if (numRemainingBits > 0) {
            // empty the overflow bits
            byte b = newBits[lastByteIndex];
            b &= NOT_BIT_MASKS[numRemainingBits];
            b &= NOT_BIT_MASKS[numRemainingBits+1];
            newBits[lastByteIndex] = b;
        }
        
        return newBits;
    }
        
    public boolean isLowComplexity(byte[] bits) {
        byte nf1[]     = new byte[4];
        byte nf2[][]   = new byte[4][4];
        byte nf3[][][] = new byte[4][4][4];
        
        int[] ranks = bitsToRanks(bits);
        
        int c3 = ranks[0];
        int c2 = ranks[1];
        int c1 = ranks[2];
        
        ++nf1[c3];
        ++nf1[c2];
        ++nf1[c1];
        
        ++nf2[c3][c2];
        ++nf2[c2][c1];
        
        ++nf3[c3][c2][c1];
                
        for (int i=3; i<k; ++i) {
            c3 = c2;
            c2 = c1;
            c1 = ranks[i];
            
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
    
    public byte[] copyOf(byte[] bits) {
        return Arrays.copyOf(bits, numBytes);
    }
    
    public static void main(String[] args) {
        Random r = new Random();
        char[] nucleotides = new char[] {'A', 'C', 'G', 'T'};
        
        int seqLen = r.nextInt(50);
        StringBuilder sb = new StringBuilder(seqLen);

        for (int i=0; i<seqLen; ++i) {
            sb.append(nucleotides[r.nextInt(4)]);
        }
        String seq = sb.toString();
        
//        String seq = "CCCAGTAAACAAGCAAATCATGCAC";
        System.out.println(seq);
        
        KmerBitsUtils2 utils = new KmerBitsUtils2(seq.length());
        
        byte[] bits = utils.kmerToBits(seq);
        
//        byte[] newBits = utils.shiftRight(bits);
        System.out.println(utils.bitsToKmer(bits));
        
        
    }
}
