/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.ThreadLocalRandom;

/**
 *
 * @author kmnip
 */
public class NucleotideBitsUtils {
    
    private static final byte[] BIT_MASKS = new byte[]{(byte) 1, (byte) 2, (byte)(1<<2), (byte)(1<<3), (byte)(1<<4), (byte)(1<<5), (byte)(1<<6), (byte)(1<<7)};
    private static final char[][] BYTE_TETRAMER_ARRAY = getByteToTetramerArray();
    
    public static char[] byteToTetramer(byte b) {
        return BYTE_TETRAMER_ARRAY[((int)b)+128];
    }
    
    private static char[][] getByteToTetramerArray() {
        char[][] result = new char[128*2][4];
        
        char[] chars = new char[]{'A', 'C', 'G', 'T'};
        
        char[] tetramer = new char[4];
        for (char c0 : chars) {
            tetramer[0] = c0;
            for (char c1 : chars) {
                tetramer[1] = c1;
                for (char c2 : chars) {
                    tetramer[2] = c2;
                    for (char c3 : chars) {
                        tetramer[3] = c3;
                        byte b = tetramerToByte(tetramer);
                        result[((int)b)+128] = Arrays.copyOf(tetramer, 4);
                    }
                }
            }
        }
        
        return result;
    }
    
    private static byte tetramerToByte(char[] tetramer) {
        byte bits = 0;

        for (int i=0; i<4; ++i) {
            switch (tetramer[3-i]) {
                // A: do nothing
                case 'C':
                    // 01
                    bits |= BIT_MASKS[i*2];
                    break;
                case 'G':
                    // 10
                    bits |= BIT_MASKS[i*2+1];
                    break;
                case 'T':
                    // 11
                    bits |= BIT_MASKS[i*2];
                    bits |= BIT_MASKS[i*2+1];
                    break;
            }
        }
        
        return bits;
    }
    
    public static byte[] seqToByteArray(String seq) {
        int seqLen = seq.length();
        int remainder = seqLen % 4;
        int numBytesForSeq = remainder > 0 ? seqLen/4 +1 : seqLen/4;
        byte[] result = new byte[4 + numBytesForSeq]; // 4 bytes for seqLen
        
        byte[] seqLenAsByteArray = ByteBuffer.allocate(4).putInt(seqLen).array();
        for (int i=0; i<4; ++i) {
            result[i] = seqLenAsByteArray[i];
        }
        
        char[] tetramer = new char[4];
        int lastSeqIndex = seqLen - remainder;
        for (int i=0; i<lastSeqIndex; i+=4) {
            tetramer[0] = seq.charAt(i);
            tetramer[1] = seq.charAt(i+1);
            tetramer[2] = seq.charAt(i+2);
            tetramer[3] = seq.charAt(i+3);
            result[i/4+4] = tetramerToByte(tetramer);
        }
        
        if (remainder > 0) {
            for (int i=0; i<remainder; ++i) {
                tetramer[i] = seq.charAt(lastSeqIndex+i);
            }
            for (int i=remainder; i<4; ++i) {
                tetramer[i] = 'A';
            }
            result[seqLen/4+4] = tetramerToByte(tetramer);
        }
        
        return result;
    }
    
    public static String byteArrayToSeq(int seqLen, byte[] b) {
        //int seqLen = ByteBuffer.wrap(Arrays.copyOfRange(b, 0, 4)).getInt();
        StringBuilder sb = new StringBuilder(seqLen);
        int byteArrLastIndex = b.length - 1;
        for (int i=0; i<byteArrLastIndex; ++i) {
            sb.append(byteToTetramer(b[i]));
        }
        
        int remainder = seqLen % 4;
        char[] lastTetramer = byteToTetramer(b[byteArrLastIndex]);
        if (remainder == 0) {
            sb.append(lastTetramer);
        }
        else {
            for (int i=0; i<remainder; ++i) {
                sb.append(lastTetramer[i]);
            }    
        }
                
        return sb.toString();
    }
    
    public static BitSet seqToBitset(String seq) {
        int len = seq.length();
        BitSet bits = new BitSet(len);
        
        for (int i=0; i<len; ++i) {
            switch (seq.charAt(i)) {
                case 'A': case 'a':
                    // 00
                    break;
                case 'C': case 'c':
                    // 01
                    bits.set(2*i + 1);
                    break;
                case 'G': case 'g':
                    // 10
                    bits.set(2*i);
                    break;
                case 'T': case 't': case 'U': case 'u':
                    // 11
                    bits.set(2*i, 2*i + 2);
                    break;
                default:
                    // use a random nucleotide
                    switch (ThreadLocalRandom.current().nextInt(0, 4)) {
                        case 0:
                            break;
                        case 1:
                            bits.set(2*i + 1);
                            break;
                        case 2:
                            bits.set(2*i);
                            break;
                        case 3:
                            bits.set(2*i, 2*i + 2);
                            break;
                    }
                    break;
            }
        }
        
        return bits;
    }
    
    public static String bitsetToSeq(BitSet bits, int strLen, boolean useUracil) {
        int numBits = strLen*2;
        StringBuilder sb = new StringBuilder(strLen);
        
        char utBase = useUracil ? 'U' : 'T';
        
        for (int i=0; i<numBits; ++i) {
            if (bits.get(i)) {
                if (bits.get(++i)) {
                    // 11
                    sb.append(utBase);
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
    
    public static void main(String[] args) {
        String seq = "CACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGCAGCTCCCAGATGGGTCTTGTCCCAGGTGCGGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCTCTTATACAAATCTCCGAGCCCACGAGAGCTCTCAGGATATCGTATGGATGAGACAGCTTGAACACACAA";
        
        int len = seq.length();
        BitSet bits = seqToBitset(seq);
        String seq2 = bitsetToSeq(bits, len, false);
        
        System.out.println(seq + "\n" + seq2);
        System.out.println(seq.equals(seq2));
        
        byte[] byteArr = seqToByteArray(seq);
        int seqLen = ByteBuffer.wrap(Arrays.copyOfRange(byteArr, 0, 4)).getInt();
        byte[] seqByteArr = Arrays.copyOfRange(byteArr, 4, byteArr.length);
        String seq3 = byteArrayToSeq(seqLen, seqByteArr);
        
        System.out.println(seq + "\n" + seq3);
        System.out.println(seq.equals(seq3));
    }
}
