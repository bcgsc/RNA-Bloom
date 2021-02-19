/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.concurrent.ThreadLocalRandom;

/**
 *
 * @author Ka Ming Nip
 */
public class SeqBitsUtils {
    private static final byte[][][][] BYTE_LOOKUP_TABLE = getByteLookupTable();
    private static final String[] TETRAMER_LOOKUP_TABLE = getTetramerLookupTable();
    
    private static byte[][][][] getByteLookupTable() {
        byte[][][][] lookupTable = new byte[4][4][4][4];
        for (int h=0; h<4; ++h) {
            // XX 00 00 00
            byte b0 = (byte) (h << 6);
            for (int i=0; i<4; ++i) {
                // 00 XX 00 00
                byte b1 = (byte) (i << 4);
                for (int j=0; j<4; ++j) {
                    // 00 00 XX 00
                    byte b2 = (byte) (j << 2);
                    for (int k=0; k<4; ++k) {
                        byte b = (byte) k;
                        b |= b2;
                        b |= b1;
                        b |= b0;
                        lookupTable[h][i][j][k] = b;
                    }
                }
            }
        }
        return lookupTable;
    }
    
    private static char getBase(int i) {
        switch(i) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                return getBase(ThreadLocalRandom.current().nextInt(0, 4));
        }
    }
    
    private static int getIndex(char c) {
        switch(c) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 1;
            case 'G':
            case 'g':
                return 2;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                return 3;
            default:
                return ThreadLocalRandom.current().nextInt(0, 4);
        }
    }
    
    private static String[] getTetramerLookupTable() {
        String[] lookupTable = new String[256];
        char[] tetramer = new char[4];
        for (int h=0; h<4; ++h) {
            // XX 00 00 00
            byte b0 = (byte) (h << 6);
            tetramer[0] = getBase(h);
            for (int i=0; i<4; ++i) {
                // 00 XX 00 00
                byte b1 = (byte) (i << 4);
                tetramer[1] = getBase(i);
                for (int j=0; j<4; ++j) {
                    // 00 00 XX 00
                    byte b2 = (byte) (j << 2);
                    tetramer[2] = getBase(j);
                    for (int k=0; k<4; ++k) {
                        byte b = (byte) k;
                        b |= b2;
                        b |= b1;
                        b |= b0;
                        tetramer[3] = getBase(k);
                        lookupTable[b + 128] = new String(tetramer);
                    }
                }
            }
        }
        return lookupTable;
    }
    
    public static byte[] seqToBits(String seq) {
        final int[] bitShifts = new int[]{6, 4, 2, 0};
        int seqLen = seq.length();
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        boolean hasExtra = numExtraBases > 0;
        int bLen = hasExtra ? numFullBytes + 1 : numFullBytes;
        byte[] bytes = new byte[bLen]; 
        
        int bIndex = 0;
        for (int i=0; i+3<seqLen; i+=4) {
            bytes[bIndex] = BYTE_LOOKUP_TABLE[getIndex(seq.charAt(i))][getIndex(seq.charAt(i+1))][getIndex(seq.charAt(i+2))][getIndex(seq.charAt(i+3))];
            ++bIndex;
        }
        
        if (hasExtra) {
            byte b = 0;
            int offset = seqLen - numExtraBases;
            for (int i=0; i<numExtraBases; ++i) {
                b |= (byte) (getIndex(seq.charAt(i + offset)) << bitShifts[i]);
            }
            bytes[numFullBytes] = b;
        }
        
        return bytes;
    }
    
    public static String bitsToSeq(byte[] bytes, int seqLen) {
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        
        StringBuilder sb = new StringBuilder(seqLen);
        for (int i=0; i<numFullBytes; ++i) {
            sb.append(TETRAMER_LOOKUP_TABLE[bytes[i] + 128]);
        }
        
        if (numExtraBases > 0) {
            sb.append(TETRAMER_LOOKUP_TABLE[bytes[numFullBytes] + 128].substring(0, numExtraBases));
        }
        
        return sb.toString();
    }
    
    public static void main(String[] args) {
        //debug
    }
}
