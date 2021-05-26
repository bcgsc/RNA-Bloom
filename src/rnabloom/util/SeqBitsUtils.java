/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;

/**
 *
 * @author Ka Ming Nip
 */
public class SeqBitsUtils {
    private static final byte[][][][] BYTE_LOOKUP_TABLE = getByteLookupTable();
    private static final String[] TETRAMER_LOOKUP_TABLE = getTetramerLookupTable();
    private static final byte[] REVERSE_COMPLEMENT_BYTES = getRevCompBytes();
    private static final int BYTE_OFFSET = 128;

    private static byte[][][][] getByteLookupTable() {
        byte[][][][] lookupTable = new byte[4][4][4][4];
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                for (int k=0; k<4; ++k) {
                    for (int l=0; l<4; ++l) {
                        lookupTable[i][j][k][l] = getByte(i, j, k, l);
                    }
                }
            }
        }
        return lookupTable;
    }
    
//    private static byte[][][][] getByteLookupTable() {
//        byte[][][][] lookupTable = new byte[4][4][4][4];
//        for (int h=0; h<4; ++h) {
//            // XX 00 00 00
//            byte b0 = (byte) (h << 6);
//            for (int i=0; i<4; ++i) {
//                // 00 XX 00 00
//                byte b1 = (byte) (i << 4);
//                for (int j=0; j<4; ++j) {
//                    // 00 00 XX 00
//                    byte b2 = (byte) (j << 2);
//                    for (int k=0; k<4; ++k) {
//                        byte b = (byte) k;
//                        b |= b2;
//                        b |= b1;
//                        b |= b0;
//                        lookupTable[h][i][j][k] = b;
//                    }
//                }
//            }
//        }
//        return lookupTable;
//    }
    
//    private static byte[] getRevCompBytes() {
//        byte[] lookupTable = new byte[256];
//        for (int h=0; h<4; ++h) {
//            // XX 00 00 00
//            byte b0 = (byte) (h << 6);
//            byte b0RC = (byte) (3-h);
//            for (int i=0; i<4; ++i) {
//                // 00 XX 00 00
//                byte b1 = (byte) (i << 4);
//                byte b1RC = (byte) ((3-i) << 2);
//                for (int j=0; j<4; ++j) {
//                    // 00 00 XX 00
//                    byte b2 = (byte) (j << 2);
//                    byte b2RC = (byte) ((3-j) << 4);
//                    for (int k=0; k<4; ++k) {
//                        byte b = (byte) k;
//                        b |= b2;
//                        b |= b1;
//                        b |= b0;
//                        byte bRC = (byte) ((3-k) << 6);
//                        bRC |= b0RC;
//                        bRC |= b1RC;
//                        bRC |= b2RC;
//                        lookupTable[b + 128] = bRC;
//                    }
//                }
//            }
//        }
//        return lookupTable;
//    }
    
    private static byte[] getRevCompBytes() {
        byte[] lookupTable = new byte[256];
        
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                for (int k=0; k<4; ++k) {
                    for (int l=0; l<4; ++l) {
                        byte b = getByte(i, j, k, l);
                        lookupTable[b + BYTE_OFFSET] = getRevCompByte(i, j, k, l);
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
    
    private static byte getByte(int i, int j, int k, int l) {
        return (byte) (i*64 + j*16 + k*4 + l - BYTE_OFFSET);
    }

    private static byte getRevCompByte(int i, int j, int k, int l) {
        return (byte) ((3-l)*64 + (3-k)*16 + (3-j)*4 + (3-i) - BYTE_OFFSET);
    }
    
//    private static String[] getTetramerLookupTable() {
//        String[] lookupTable = new String[256];
//        char[] tetramer = new char[4];
//        for (int h=0; h<4; ++h) {
//            // XX 00 00 00
//            byte b0 = (byte) (h << 6);
//            tetramer[0] = getBase(h);
//            for (int i=0; i<4; ++i) {
//                // 00 XX 00 00
//                byte b1 = (byte) (i << 4);
//                tetramer[1] = getBase(i);
//                for (int j=0; j<4; ++j) {
//                    // 00 00 XX 00
//                    byte b2 = (byte) (j << 2);
//                    tetramer[2] = getBase(j);
//                    for (int k=0; k<4; ++k) {
//                        byte b = (byte) k;
//                        b |= b2;
//                        b |= b1;
//                        b |= b0;
//                        tetramer[3] = getBase(k);
//                        lookupTable[b + 128] = new String(tetramer);
//                    }
//                }
//            }
//        }
//        return lookupTable;
//    }
    
   private static String[] getTetramerLookupTable() {
        String[] lookupTable = new String[256];
        char[] tetramer = new char[4];
        for (int i=0; i<4; ++i) {
            tetramer[0] = getBase(i);
            for (int j=0; j<4; ++j) {
                tetramer[1] = getBase(j);
                for (int k=0; k<4; ++k) {
                    tetramer[2] = getBase(k);
                    for (int l=0; l<4; ++l) {
                        tetramer[3] = getBase(l);
                        byte b = getByte(i, j, k, l);
                        lookupTable[b + BYTE_OFFSET] = new String(tetramer);
                    }
                }
            }
        }
        return lookupTable;
    }
    
    public static String getTetramer(byte b) {
        return TETRAMER_LOOKUP_TABLE[b + BYTE_OFFSET];
    }
    
    public static byte reverseComplement(byte b) {
        return REVERSE_COMPLEMENT_BYTES[b + BYTE_OFFSET];
    }
    
    public static int seqLenToNumBytes(int seqLen) {
        return seqLen % 4 > 0 ? seqLen/4+1 : seqLen/4;
    }
    
    public static byte[] intToFourBytes(int i) {
        return ByteBuffer.allocate(4).putInt(i).array();
    }
    
    public static int fourBytesToInt(byte[] b) {
        return ByteBuffer.wrap(b).getInt();
    }
        
    public static byte[] seqToBits(String seq) {
        int seqLen = seq.length();
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        boolean hasExtra = numExtraBases > 0;
        int bLen = hasExtra ? numFullBytes + 1 : numFullBytes;
        byte[] bytes = new byte[bLen]; 
        
        int bIndex = 0;
        for (int i=0; i+3<seqLen; i+=4) {
            bytes[bIndex] = BYTE_LOOKUP_TABLE[getIndex(seq.charAt(i))][getIndex(seq.charAt(i+1))][getIndex(seq.charAt(i+2))][getIndex(seq.charAt(i+3))];
//            bytes[bIndex] = getByte(getIndex(seq.charAt(i)), getIndex(seq.charAt(i+1)),
//                    getIndex(seq.charAt(i+2)), getIndex(seq.charAt(i+3)));
            ++bIndex;
        }
        
        if (hasExtra) {
            int offset = seqLen - numExtraBases;
            int i = getIndex(seq.charAt(offset));
            int j = 1 < numExtraBases ? getIndex(seq.charAt(offset+1)) : 0;
            int k = 2 < numExtraBases ? getIndex(seq.charAt(offset+2)) : 0;
            int l = 3 < numExtraBases ? getIndex(seq.charAt(offset+3)) : 0;
            bytes[numFullBytes] = BYTE_LOOKUP_TABLE[i][j][k][l];
//            bytes[numFullBytes] = getByte(i, j, k, l);
        }
        
        return bytes;
    }
    
    public static byte[] seqToBitsParallelized(String seq) {
        int seqLen = seq.length();
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        boolean hasExtra = numExtraBases > 0;
        int bLen = hasExtra ? numFullBytes + 1 : numFullBytes;
        byte[] bytes = new byte[bLen]; 
        
        IntStream.range(0, numFullBytes).parallel().forEach(e -> {
            int i = e * 4;
            bytes[e] = BYTE_LOOKUP_TABLE[getIndex(seq.charAt(i))][getIndex(seq.charAt(i+1))][getIndex(seq.charAt(i+2))][getIndex(seq.charAt(i+3))];
//            bytes[e] = getByte(getIndex(seq.charAt(i)), getIndex(seq.charAt(i+1)),
//                    getIndex(seq.charAt(i+2)), getIndex(seq.charAt(i+3)));
        });
        
        if (hasExtra) {
            int offset = seqLen - numExtraBases;
            int i = getIndex(seq.charAt(offset));
            int j = 1 < numExtraBases ? getIndex(seq.charAt(offset+1)) : 0;
            int k = 2 < numExtraBases ? getIndex(seq.charAt(offset+2)) : 0;
            int l = 3 < numExtraBases ? getIndex(seq.charAt(offset+3)) : 0;
            bytes[numFullBytes] = BYTE_LOOKUP_TABLE[i][j][k][l];
//            bytes[numFullBytes] = getByte(i, j, k, l);
        }
        
        return bytes;
    }
    
    public static String bitsToSeq(byte[] bytes, int seqLen) {
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        
        StringBuilder sb = new StringBuilder(seqLen);
        for (int i=0; i<numFullBytes; ++i) {
            sb.append(getTetramer(bytes[i]));
        }
        
        if (numExtraBases > 0) {
            sb.append(getTetramer(bytes[numFullBytes]).substring(0, numExtraBases));
        }
        
        return sb.toString();
    }
    
    public static String bitsToSeqParallelized(byte[] bytes, int seqLen) {
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        
        char[] sb = new char[seqLen];
        IntStream.range(0, numFullBytes).parallel().forEach(e -> {
            String t = getTetramer(bytes[e]);
            int i = e*4;
            sb[i] = t.charAt(0);
            sb[i+1] = t.charAt(1);
            sb[i+2] = t.charAt(2);
            sb[i+3] = t.charAt(3);
        });
        
        if (numExtraBases > 0) {
            String t = getTetramer(bytes[numFullBytes]);
            int offset = numFullBytes*4;
            for (int i=0; i<numExtraBases; ++i) {
                sb[offset + i] = t.charAt(i);
            }
        }
        
        return new String(sb);
    }
    
    public static String bitsToSeq(byte[] bytes, int seqLen, int seqStart, int seqEnd) {
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        
        int startIndex = seqStart / 4;
        int startExtraBases = seqStart % 4;

        int endIndex = seqEnd / 4;
        int endExtraBases = seqEnd % 4;
        if (endIndex == numFullBytes) {
            endExtraBases = Math.min(endExtraBases, numExtraBases);
        }

        StringBuilder sb = new StringBuilder(seqEnd - seqStart);
        if (startExtraBases == 0) {
            sb.append(getTetramer(bytes[startIndex]));
        }
        else {
            sb.append(getTetramer(bytes[startIndex]).substring(startExtraBases));
        }
        
        for (int i=startIndex+1; i<endIndex; ++i) {
            sb.append(getTetramer(bytes[i]));
        }
        
        if (endExtraBases > 0) {
            sb.append(getTetramer(bytes[endIndex]).substring(0, endExtraBases));
        }
        
        return sb.toString();
    }
    
    public static String bitsToRevCompSeq(byte[] bytes, int seqLen, int seqStart, int seqEnd) {
        int numFullBytes = seqLen / 4;
        int numExtraBases = seqLen % 4;
        
        int startIndex = seqStart / 4;
        int startExtraBases = seqStart % 4;

        int endIndex = seqEnd / 4;
        int endExtraBases = seqEnd % 4;
        if (endIndex == numFullBytes) {
            endExtraBases = Math.min(endExtraBases, numExtraBases);
        }

        StringBuilder sb = new StringBuilder(seqEnd - seqStart);
        if (endExtraBases > 0) {
            sb.append(getTetramer(reverseComplement(bytes[endIndex])).substring(4-endExtraBases));
        }

        for (int i=endIndex-1; i>startIndex; --i) {
            sb.append(getTetramer(reverseComplement(bytes[i])));
        }
        
        if (startExtraBases == 0) {
            sb.append(getTetramer(reverseComplement(bytes[startIndex])));
        }
        else {
            sb.append(getTetramer(reverseComplement(bytes[startIndex])).substring(0, 4-startExtraBases));
        }
                
        return sb.toString();
    }
    
    public static void writeToFile(ArrayList<? extends BitSequence> list, String outFasta) throws IOException {
        FastaWriter writer = new FastaWriter(outFasta, false);
        
        int cid = 0;
        for (BitSequence b : list) {
            writer.write("r" + ++cid, b.toString());
        }
        
        writer.close();
    }
    
    public static ArrayList<BitSequence> getListFromFile(String inFasta) throws IOException {
        ArrayList<BitSequence> bitSeqs = new ArrayList<>();
        FastaReader reader = new FastaReader(inFasta);
        while (reader.hasNext()) {
            bitSeqs.add(new BitSequence(reader.next()));
        }
        reader.close();
        
        return bitSeqs;
    }
    
    public static void main(String[] args) {
        //debug
        
    }
}
