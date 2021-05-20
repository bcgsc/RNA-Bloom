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

import java.util.AbstractCollection;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.PrimitiveIterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.io.FastqRecord;

/**
 *
 * @author Ka Ming Nip
 */
public final class SeqUtils {
    private final static int CHAR_A_INT = (int) 'A';
    private final static int CHAR_C_INT = (int) 'C';
    private final static int CHAR_G_INT = (int) 'G';
    private final static int CHAR_T_INT = (int) 'T';
    private final static int CHAR_U_INT = (int) 'U';
    
    private final static byte CHAR_A_BYTE = (byte) 'A';
    private final static byte CHAR_C_BYTE = (byte) 'C';
    private final static byte CHAR_G_BYTE = (byte) 'G';
    private final static byte CHAR_T_BYTE = (byte) 'T';
    private final static byte CHAR_U_BYTE = (byte) 'U';
    public final static byte[] NUCLEOTIDES_BYTES = new byte[] {CHAR_A_BYTE, CHAR_C_BYTE, CHAR_G_BYTE, CHAR_T_BYTE};
    
    public static final char GAP_CHAR = 'N';
    
    public final static char[] NUCLEOTIDES = new char[] {'A','C','G','T'};
    public final static char[] A_ALT_NUCLEOTIDES = new char[] {'C','G','T'};
    public final static char[] C_ALT_NUCLEOTIDES = new char[] {'A','G','T'};
    public final static char[] G_ALT_NUCLEOTIDES = new char[] {'A','C','T'};
    public final static char[] T_ALT_NUCLEOTIDES = new char[] {'A','C','G'};
    public final static byte[] A_ALT_NUCLEOTIDES_BYTES = new byte[] {CHAR_C_BYTE,CHAR_G_BYTE,CHAR_T_BYTE};
    public final static byte[] C_ALT_NUCLEOTIDES_BYTES = new byte[] {CHAR_A_BYTE,CHAR_G_BYTE,CHAR_T_BYTE};
    public final static byte[] G_ALT_NUCLEOTIDES_BYTES = new byte[] {CHAR_A_BYTE,CHAR_C_BYTE,CHAR_T_BYTE};
    public final static byte[] T_ALT_NUCLEOTIDES_BYTES = new byte[] {CHAR_A_BYTE,CHAR_C_BYTE,CHAR_G_BYTE};
    
    private static final float LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ = 0.95f;
    private static final float LOW_COMPLEXITY_THRESHOLD_LONG_SEQ = 0.89f;
    
    public static final byte[] stringToBytes(String seq, int len) {
        byte[] arr = new byte[len];
        
        byte b;
        for (int i=0; i<len; ++i) {
            b = (byte) seq.charAt(i);
            if (b == CHAR_U_INT) {
                arr[i] = CHAR_T_INT;
            }
            else {
                arr[i] = b;
            }
        }
        
        return arr;
    }
    
    public static final String bytesToString(byte[] bytes, int len) {
        StringBuilder sb = new StringBuilder(len);
        
        for (byte b : bytes) {
            sb.append((char) b);
        }
        
        return sb.toString();
    }
    
    public static final String bytesToString(byte[] bytes, int start, int end) {
        StringBuilder sb = new StringBuilder(end-start);
        
        for (int i=start; i<end; ++i) {
            sb.append((char) bytes[i]);
        }
        
        return sb.toString();
    }
    
    public static final String getHeadTailSummary(String seq, int k) {
        int fragLen = seq.length();
        if (fragLen > 2*k) {
            return seq.substring(0, k) + "{" + Integer.toString(fragLen-2*k) + "}" + seq.substring(fragLen-k, fragLen);
        }
        
        return seq;
    }
    
    public static final byte[] shiftRight(byte[] bytes, int len) {
        byte[] bytes2 = new byte[len];
        System.arraycopy(bytes, 0, bytes2, 1, len-1);
        return bytes2;
    }
    
    public static final byte[] shiftLeft(byte[] bytes, int len) {
        byte[] bytes2 = new byte[len];
        System.arraycopy(bytes, 1, bytes2, 0, len-1);
        return bytes2;
    }
    
    public static int countBits(byte b) {
        b = (byte) ((byte) (b & 0b01010101) + (byte)(((byte) (b >>> 1)) & 0b01010101));
        b = (byte) ((byte) (b & 0b00110011) + (byte)(((byte) (b >>> 2)) & 0b00110011));
        b = (byte) ((byte) (b & 0b00001111) + (byte)(((byte) (b >>> 4)) & 0b00001111));
        
        return b;
    }
    
    public static final byte[] getAltNucleotides(byte c) {
        switch (c) {
            case CHAR_A_BYTE:
                return A_ALT_NUCLEOTIDES_BYTES;
            case CHAR_C_BYTE:
                return C_ALT_NUCLEOTIDES_BYTES;
            case CHAR_G_BYTE:
                return G_ALT_NUCLEOTIDES_BYTES;
            case CHAR_T_BYTE:
                return T_ALT_NUCLEOTIDES_BYTES;
            case CHAR_U_BYTE:
                return T_ALT_NUCLEOTIDES_BYTES;
            default:
                return NUCLEOTIDES_BYTES;
        }
    }
    
    public static final char[] getAltNucleotides(char c) {
        switch (c) {
            case 'A':
                return A_ALT_NUCLEOTIDES;
            case 'C':
                return C_ALT_NUCLEOTIDES;
            case 'G':
                return G_ALT_NUCLEOTIDES;
            case 'T':
                return T_ALT_NUCLEOTIDES;
            case 'U':
                return T_ALT_NUCLEOTIDES;
            default:
                return NUCLEOTIDES;
        }
    }
    
    public static float getPercentIdentity(String a, String b) {
        int aLen = a.length();
        int bLen = b.length();

        int d = getDistance(a, b, aLen, bLen);
        
        if (aLen <= bLen) {
            return ((float) (bLen - d))/(float)bLen;
        }
        
        return ((float) (aLen - d))/(float)aLen;
    }
    
    public static float getPercentIdentity(byte[] a, byte[] b) {
        int aLen = a.length;
        int bLen = b.length;

        int d = getDistance(a, b, aLen, bLen);
        
        if (aLen <= bLen) {
            return ((float) (bLen - d))/(float)bLen;
        }
        
        return ((float) (aLen - d))/(float)aLen;
    }
    
    private static int getDistance(String s, String t, int sLen, int tLen) {
        // compute the Levenshtein Distance
        // https://en.wikipedia.org/wiki/Levenshtein_distance
        
        // degenerate cases
        if (s.equals(t)) return 0;
        if (sLen == 0) return tLen;
        if (tLen == 0) return sLen;

        // create two work vectors of integer distances
        int[] v0 = new int[tLen+1];
        int[] v1 = new int[tLen+1];
        
        // initialize v0 (the previous row of distances)
        // this row is A[0][i]: edit distance for an empty s
        // the distance is just the number of characters to delete from t
        for (int i=0; i<=tLen; ++i) {
            v0[i] = i;
        }

        for (int i=0; i<sLen; ++i) {
            // calculate v1 (current row distances) from the previous row v0

            // first element of v1 is A[i+1][0]
            //   edit distance is delete (i+1) chars from s to match empty t
            
            v1[0] = i+1;
            
            for (int j=0; j<tLen; ++j) {
                v1[j+1] = min3(v1[j  ]+1,
                               v0[j+1]+1,
                               v0[j  ]+(s.charAt(i) == t.charAt(j) ? 0 : 1));
            }
            
            // copy v1 (current row) to v0 (previous row) for next iteration
            System.arraycopy(v1, 0, v0, 0, tLen);
        }

        return v1[tLen];
    }
    
    private static int getDistance(byte[] s, byte[] t, int sLen, int tLen) {
        // degenerate cases
        if (Arrays.equals(s, t)) return 0;
        if (sLen == 0) return tLen;
        if (tLen == 0) return sLen;
        
        // create two work vectors of integer distances
        int[] v0 = new int[tLen+1];
        int[] v1 = new int[tLen+1];
        
        // initialize v0 (the previous row of distances)
        // this row is A[0][i]: edit distance for an empty s
        // the distance is just the number of characters to delete from t
        for (int i=0; i<=tLen; ++i) {
            v0[i] = i;
        }

        for (int i=0; i<sLen; ++i) {
            // calculate v1 (current row distances) from the previous row v0

            // first element of v1 is A[i+1][0]
            //   edit distance is delete (i+1) chars from s to match empty t
            
            v1[0] = i+1;
            
            for (int j=0; j<tLen; ++j) {
                v1[j+1] = min3(v1[j  ]+1,
                               v0[j+1]+1,
                               v0[j  ]+(s[i] == t[j] ? 0 : 1));
            }
            
            // copy v1 (current row) to v0 (previous row) for next iteration
            System.arraycopy(v1, 0, v0, 0, tLen);
        }

        return v1[tLen];
    }
    
    private static int min3(int a, int b, int c) {
        return Math.min(a, Math.min(b, c));
    }
    
    public static final int getNumGC(String seq) {
        int numGC = 0;
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int c;
        
        while (itr.hasNext()) {
            c = itr.nextInt();
            switch(c) {
                case CHAR_C_INT:
                    ++numGC;
                    break;
                case CHAR_G_INT:
                    ++numGC;
                    break;
            }
        }
        
        return numGC;
    }
    
    public static final float getGCContent(String seq) {
        return (float) getNumGC(seq) / seq.length();
    }
    
    private static int nucleotideArrayIndex(byte b) {
        switch(b) {
            case CHAR_A_BYTE:
                return 0;
            case CHAR_C_BYTE:
                return 1;
            case CHAR_G_BYTE:
                return 2;
            case CHAR_T_BYTE:
                return 3;
            case CHAR_U_BYTE:
                return 3;
            default:
                return -1;
        }
    }
    
    private static int nucleotideArrayIndex(int c) {
        switch(c) {
            case CHAR_A_INT:
                return 0;
            case CHAR_C_INT:
                return 1;
            case CHAR_G_INT:
                return 2;
            case CHAR_T_INT:
                return 3;
            case CHAR_U_INT:
                return 3;
            default:
                return -1;
        }
    }
    
    public static int[] countN(String seq) {
        int[] result = new int[4];
        countN(seq, result);
        return result;
    }
    
    public static void countN(String seq, int[] result) {
        result[0] = 0;
        result[1] = 0;
        result[2] = 0;
        result[3] = 0;
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int i;
        while (itr.hasNext()) {
            i = nucleotideArrayIndex(itr.nextInt());
            if (i >= 0) {
                ++result[i];
            }
        }
    }

    public static final boolean isHomopolymer(byte[] seq) {
        int len = seq.length;
        if (len == 0) {
            return false;
        }
        
        byte b = seq[0];
        for (int i=1; i<len; ++i) {
            if (b != seq[i]) {
                return false;
            }
        }
        
        return true;
    }

    public static final boolean isLowComplexity2(byte[] bytes) {
        int length = bytes.length;
        int t1 = Math.min(Byte.MAX_VALUE, Math.round(length * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ));
        int t2 = Math.min(Byte.MAX_VALUE, Math.round(length * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ / 2));
        int t3 = Math.min(Byte.MAX_VALUE, Math.round(length * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ / 3));
        
        byte nf1[]     = new byte[4];
        byte nf2[][]   = new byte[4][4];
        byte nf3[][][] = new byte[4][4][4];
        
        int c3 = nucleotideArrayIndex(bytes[0]);
        int c2 = nucleotideArrayIndex(bytes[1]);
        int c1 = nucleotideArrayIndex(bytes[2]);
        
        ++nf1[c3];
        ++nf1[c2];
        ++nf1[c1];
        
        if (c3 != c2) {
            ++nf2[c3][c2];
        }
        if (c2 != c1) {
            ++nf2[c2][c1];
        }
        
        if (c3 != c2 || c2 != c1 || c3 != c1) {
            ++nf3[c3][c2][c1];
        }
        
        for (int i=3; i<length; ++i) {
            c3 = c2;
            c2 = c1;
            c1 = nucleotideArrayIndex(bytes[i]);
            
            if (++nf1[c1] >= t1)
                return true; // homopolymer runs
            if (c2 != c1 && ++nf2[c2][c1] >= t2)
                return true; // di-nucleotide repeat
            if ((c3 != c2 || c2 != c1 || c3 != c1) && ++nf3[c3][c2][c1] >= t3)
                return true; // tri-nucleotide repeat
        }
        
        // check bias in di/tri-nucleotide content
        return (nf1[0]+nf1[1]>=t1 || nf1[0]+nf1[2]>=t1 || nf1[0]+nf1[3]>=t1 || 
                nf1[1]+nf1[2]>=t1 || nf1[1]+nf1[3]>=t1 || nf1[2]+nf1[3]>=t1);
    }
    
    public static final boolean isRepeat(String seq) {
        final float repeatThreshold = 0.9f;
        final int length = seq.length();
        
         // homopolymer runs
        final int t1 = Math.round(length * repeatThreshold);
        byte[] nf1 = new byte[4];
        for (int i=0; i<length; ++i) {
            int n = nucleotideArrayIndex(seq.charAt(i));
            if (++nf1[n] >= t1)
                return true;
        }
        
        // di-nucleotide repeat
        final int t2 = Math.round(length/2 * repeatThreshold);
        for (int start=0; start<2; ++start) {
            byte[][] nf2 = new byte[4][4];
            for (int i=start; i<length-1; i+=2) {
                int n1 = nucleotideArrayIndex(seq.charAt(i));
                int n2 = nucleotideArrayIndex(seq.charAt(i+1));
                if (++nf2[n1][n2] >= t2)
                    return true;
            }
        }
        
        // tri-nucleotide repeat
        final int t3 = Math.round(length/3 * repeatThreshold);
        for (int start=0; start<3; ++start) {
            byte[][][] nf3 = new byte[4][4][4];
            for (int i=start; i<length-2; i+=3) {
                int n1 = nucleotideArrayIndex(seq.charAt(i));
                int n2 = nucleotideArrayIndex(seq.charAt(i+1));
                int n3 = nucleotideArrayIndex(seq.charAt(i+2));
                if (++nf3[n1][n2][n3] >= t3)
                    return true;
            }
        }
        
        return false;
    }
    
    public static final boolean isRepeat(byte[] bytes) {
        final float repeatThreshold = 0.9f;
        final int length = bytes.length;
        
         // homopolymer runs
        final int t1 = Math.round(length * repeatThreshold);
        byte[] nf1 = new byte[4];
        for (int i=0; i<length; ++i) {
            int n = nucleotideArrayIndex(bytes[i]);
            if (++nf1[n] >= t1)
                return true;
        }
        
        // di-nucleotide repeat
        final int t2 = Math.round(length/2 * repeatThreshold);
        for (int start=0; start<2; ++start) {
            byte[][] nf2 = new byte[4][4];
            for (int i=start; i<length-1; i+=2) {
                int n1 = nucleotideArrayIndex(bytes[i]);
                int n2 = nucleotideArrayIndex(bytes[i+1]);
                if (++nf2[n1][n2] >= t2)
                    return true;
            }
        }
        
        // tri-nucleotide repeat
        final int t3 = Math.round(length/3 * repeatThreshold);
        for (int start=0; start<3; ++start) {
            byte[][][] nf3 = new byte[4][4][4];
            for (int i=start; i<length-2; i+=3) {
                int n1 = nucleotideArrayIndex(bytes[i]);
                int n2 = nucleotideArrayIndex(bytes[i+1]);
                int n3 = nucleotideArrayIndex(bytes[i+2]);
                if (++nf3[n1][n2][n3] >= t3)
                    return true;
            }
        }
        
        return false;
    }

    public static final boolean isLowComplexityShort(String seq) {
        int length = seq.length();
        int t1 = Math.min(Short.MAX_VALUE, Math.round(length * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ));
        int t2 = Math.min(Short.MAX_VALUE, Math.round(length/2 * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ));
        int t3 = Math.min(Short.MAX_VALUE, Math.round(length/3 * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ));
        
        short nf1[]     = new short[4];
        short nf2[][]   = new short[4][4];
        short nf3[][][] = new short[4][4][4];
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int c3 = nucleotideArrayIndex(itr.nextInt());
        int c2 = nucleotideArrayIndex(itr.nextInt());
        int c1 = nucleotideArrayIndex(itr.nextInt());
        
        ++nf1[c3];
        ++nf1[c2];
        ++nf1[c1];
        
        ++nf2[c3][c2];
        ++nf2[c2][c1];
        
        ++nf3[c3][c2][c1];
        
        while (itr.hasNext()) {
            c3 = c2;
            c2 = c1;
            c1 = nucleotideArrayIndex(itr.nextInt());
            
            if (++nf1[c1] >= t1)
                return true; // homopolymer runs
            if (++nf2[c2][c1] >= t2)
                return true; // di-nucleotide repeat
            if (++nf3[c3][c2][c1] >= t3)
                return true; // tri-nucleotide repeat
        }
        
        // di-nucleotide content
        if (nf1[0]+nf1[1]>=t1 || nf1[0]+nf1[2]>=t1 || nf1[0]+nf1[3]>=t1 || 
                nf1[1]+nf1[2]>=t1 || nf1[1]+nf1[3]>=t1 || nf1[2]+nf1[3]>=t1) {
            return true;
        }
        
        return false;
    }

    public static final boolean isLowComplexity2N(byte[] seq) {
        int length = seq.length;
        if (length <= 2) {
            return false;
        }
        
        int t1 = Math.round(length * LOW_COMPLEXITY_THRESHOLD_SHORT_SEQ);
        
        byte nf1[] = new byte[4];
        
        for (int i=0; i<length; ++i) {
            if (++nf1[nucleotideArrayIndex(seq[i])] >= t1)
                return true;
        }
        
        // di-nucleotide content
        return (nf1[0]+nf1[1]>=t1 || nf1[0]+nf1[2]>=t1 || nf1[0]+nf1[3]>=t1 || 
                nf1[1]+nf1[2]>=t1 || nf1[1]+nf1[3]>=t1 || nf1[2]+nf1[3]>=t1);
    }
    
    public static final boolean isLowComplexity2N(String seq) {
        int length = seq.length();
        if (length <= 2) {
            return false;
        }
        
        int t1 = Math.round(length * LOW_COMPLEXITY_THRESHOLD_LONG_SEQ);
        
        int nf1[] = new int[4];
        
        for (PrimitiveIterator.OfInt itr = seq.chars().iterator(); itr.hasNext();) {
            if (++nf1[nucleotideArrayIndex(itr.nextInt())] >= t1)
                return true;
        }
        
        // di-nucleotide content
        return (nf1[0]+nf1[1]>=t1 || nf1[0]+nf1[2]>=t1 || nf1[0]+nf1[3]>=t1 || 
                nf1[1]+nf1[2]>=t1 || nf1[1]+nf1[3]>=t1 || nf1[2]+nf1[3]>=t1);
    }
    
    public static final boolean isLowComplexityLong(String seq) {
        int length = seq.length();
        if (length <= 6) {
            return false;
        }
        
        int t1 = Math.round(length * LOW_COMPLEXITY_THRESHOLD_LONG_SEQ);
        int t2 = Math.round(length * LOW_COMPLEXITY_THRESHOLD_LONG_SEQ / 2.0f);
        int t3 = Math.round(length * LOW_COMPLEXITY_THRESHOLD_LONG_SEQ / 3.0f);
        
        int nf1[]     = new int[4];
        int nf2[][]   = new int[4][4];
        int nf3[][][] = new int[4][4][4];
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int c3 = nucleotideArrayIndex(itr.nextInt());
        int c2 = nucleotideArrayIndex(itr.nextInt());
        int c1 = nucleotideArrayIndex(itr.nextInt());

        ++nf1[c3];
        ++nf1[c2];
        ++nf1[c1];
        
        if (c3 != c2) {
            ++nf2[c3][c2];
        }
        if (c2 != c1) {
            ++nf2[c2][c1];
        }
        
        if (c3 != c2 || c2 != c1 || c3 != c1) {
            ++nf3[c3][c2][c1];
        }
        
        while (itr.hasNext()) {
            c3 = c2;
            c2 = c1;
            c1 = nucleotideArrayIndex(itr.nextInt());
            
            if (++nf1[c1] >= t1)
                return true; // homopolymer runs
            if (c2 != c1 && ++nf2[c2][c1] >= t2)
                return true; // di-nucleotide repeat
            if ((c3 != c2 || c2 != c1 || c3 != c1) && ++nf3[c3][c2][c1] >= t3)
                return true; // tri-nucleotide repeat
        }
        
        // check bias in di-nucleotide content
        return (nf1[0]+nf1[1]>=t1 || nf1[0]+nf1[2]>=t1 || nf1[0]+nf1[3]>=t1 || 
                nf1[1]+nf1[2]>=t1 || nf1[1]+nf1[3]>=t1 || nf1[2]+nf1[3]>=t1);
    }
    
    public static final boolean isLowComplexityLongWindowed(String seq) {
        int seqLen = seq.length();
        int windowSize = 50;
        int numWindows = seqLen/windowSize;
        if (numWindows >= 4) {
            int offset = (seqLen % windowSize) / 2;
            int numLowComplexityWindows = 0;

            for (int i=0; i<numWindows; ++i) {
                int start = i*windowSize + offset;
                int end = start + windowSize;
                String window = seq.substring(start, end);
                if (isLowComplexityLong(window)) {
                    ++numLowComplexityWindows;
                }
            }

            return numLowComplexityWindows >= Math.floor(0.75 * numWindows);
        }
        
        return isLowComplexityLong(seq);
    }
    
    public static boolean isARich(String seq, int start, int end, float minFraction) {
        return isNRich(seq, start, end, minFraction, 'A');
    }

    public static boolean isTRich(String seq, int start, int end, float minFraction) {
        return isNRich(seq, start, end, minFraction, 'T');
    }
    
    public static boolean isNRich(String seq, int start, int end, float minFraction, char n) {
        int threshold = Math.round(minFraction * (end - start));
        int count = 0;
        
        for (int i=start; i<end; ++i) {
            if (seq.charAt(i) == n) {
                if (++count >= threshold) {
                    return true; 
                }
            }    
        }
        
        return false;
    }
    
    public static ArrayList<String> trimLowComplexityRegions(String seq, int minLowComplexityLength) {
        ArrayList<String> segments = new ArrayList<>();
        
        int seqLen = seq.length();
        int windowSize = 100;
        int numWindows = seqLen/windowSize;
        if (numWindows >= 5) {
            boolean[] complexityStatus = new boolean[numWindows];
            
            int offset = (seqLen % windowSize) / 2;
            int numLowComplexityWindows = 0;
            
            for (int i=0; i<numWindows; ++i) {
                int start = i*windowSize + offset;
                int end = start + windowSize;
                String window = seq.substring(start, end);
                if (isLowComplexityLong(window)) {
                    ++numLowComplexityWindows;
                    complexityStatus[i] = false;
                }
                else {
                    complexityStatus[i] = true;
                }
            }

            if (numLowComplexityWindows * windowSize >= minLowComplexityLength) {
                int start = -1;
                int end = -1;
                for (int i=0; i<numWindows; ++i) {
                    if (complexityStatus[i]) {
                        if (start < 0) {
                            start = i;
                        }
                        end = i;
                    }
                    else if (start >= 0) {
                        if (start == 0) {
                            segments.add(seq.substring(0, offset + (end+1) * windowSize));
                        }
                        else {
                            segments.add(seq.substring(offset + start * windowSize, offset + (end+1) * windowSize));
                        }
                        
                        start = -1;
                        end = -1;
                    }
                }

                if (start == 0) {
                    segments.add(seq.substring(0, seqLen));
                }
                else if (start > 0) {
                    segments.add(seq.substring(offset + start * windowSize, seqLen));
                }
            }
            else {
                segments.add(seq);
            }
        }
        else if (!isLowComplexityLongWindowed(seq)) {
            segments.add(seq);
        }
        
        return segments;
    }
    
    public static String trimLowComplexityEdges(String seq) {
        int seqLen = seq.length();
        int windowSize = 50;
        int numWindows = seqLen/windowSize;
        if (numWindows >= 4) {
            int offset = (seqLen % windowSize) / 2;

            int numHeadWindows = 0;
            int skip = 0;
            for (int i=0; i<numWindows; ++i) {
                int start = i*windowSize + offset;
                int end = start + windowSize;
                String window = seq.substring(start, end);
                if (isLowComplexityLong(window)) {
                    if (skip > 0) {
                        numHeadWindows += 2;
                        skip = 0;
                    }
                    else {
                        ++numHeadWindows;
                    }
                }
                else if (skip == 0) {
                    ++skip;
                }
                else {
                    break;
                }
            }
            
            int numTailWindows = 0;
            skip = 0;
            for (int i=numWindows-1; i>numHeadWindows; --i) {
                int start = i*windowSize + offset;
                int end = start + windowSize;
                String window = seq.substring(start, end);
                if (isLowComplexityLong(window)) {
                    if (skip > 0) {
                        numTailWindows += 2;
                        skip = 0;
                    }
                    else {
                        ++numTailWindows;
                    }
                }
                else if (skip == 0) {
                    ++skip;
                }
                else {
                    break;
                }
            }
            
            if (numHeadWindows > 0 || numTailWindows > 0) {
                int start = 0;
                int end = seqLen;
                if (numHeadWindows > 0) {
                    start = offset + numHeadWindows*windowSize;
                    if (isTRich(seq, 0, start, 0.8f)) {
                        start = 0;
                    }
                }
                
                if (numTailWindows > 0) {
                    end = seqLen - offset - numTailWindows*windowSize;
                    if (isARich(seq, end, seqLen, 0.8f)) {
                        end = seqLen;
                    }
                }
                
                if (start > 0 || end < seqLen) {
                    return seq.substring(start, end);
                }
            }
        }
        
        return seq;
    }
    
    public static int getHomoPolymerCompressedLength(String seq) {
        int length = 0;
        if (!seq.isEmpty()) {
            PrimitiveIterator.OfInt itr = seq.chars().iterator();
            int prev = itr.nextInt();
            ++length;

            while (itr.hasNext()) {
                int curr = itr.nextInt();
                if (curr != prev) {
                    ++length;
                    prev = curr;
                }
            }
        }
        
        return length;
    }
    
    public static String chompRightPolyX(String seq, int minLen, int tolerance) {
        int seqLen = seq.length();
        
        int numA = 0;
        int numC = 0;
        int numG = 0;
        int numT = 0;
        
        int i, count;
        for (i=seqLen-1; i>=0; --i) {
            switch (seq.charAt(i)) {
                case 'A':
                    count = ++numA;
                    break;
                case 'C':
                    count = ++numC;
                    break;
                case 'G':
                    count = ++numG;
                    break;
                case 'T':
                    count = ++numT;
                    break;
                case 'U':
                    count = ++numT;
                    break;
                default:
                    count = 0;
            }
        
            if (count > tolerance && count < seqLen - i - tolerance) {
                break;
            }
        }
        
        if (++i < seqLen-minLen) {
            char maxChar = 'A';
            int maxCount = numA;
            
            if (numC > maxCount) {
                maxChar = 'C';
                maxCount = numC;
            }
            
            if (numG > maxCount) {
                maxChar = 'G';
                maxCount = numG;
            }
            
            if (numT > maxCount) {
                maxChar = 'T';
                maxCount = numT;
            }
            
            for (int j=i+tolerance; j>i; --j) {
                if (seq.charAt(j) != maxChar) {
                    i = j+1;
                }
            }
            
            return seq.substring(0, i+minLen/2); // retain half of specified mininum length of the polymer
        }
        
        return seq;
    }
    
    public static final int getNumKmers(String seq, int k) {
        return seq.length() - k + 1;
    }
    
    public static final int getSeqLength(int numKmers, int k) {
        return k + numKmers - 1;
    }
    
    public static final String getFirstKmer(String seq, int k) {
        return seq.substring(0, k);
    }
    
    public static final String getLastKmer(String seq, int k) {
        int seqLen = seq.length();
        return seq.substring(seqLen-k, seqLen);
    }
    
    public static class KmerSeqIterator implements Iterator<String> {
        private String seq;
        private int k;
        private int i = 0;
        public int numKmers;
        
        public KmerSeqIterator(int k) {
            this.k = k;
        }

        public KmerSeqIterator(String seq, int k) {
            this.k = k;
            initialize(seq);
        }
        
        public final void initialize(String seq) {
            this.seq = seq;
            this.numKmers = seq.length() - k + 1;
            i = 0;
        }

        @Override
        public boolean hasNext() {
            return i < numKmers;
        }

        @Override
        public String next() {
            int j = i++;
            return seq.substring(j, j+k);
        }
        
        public void reset() {
            i = 0;
        }
    }
    
    public static final String[] kmerize(String seq, int k) {
        final int numKmers = seq.length() - k +1;
        final String[] kmers = new String[numKmers];
        
        for (int i=0; i<numKmers; ++i) {
            kmers[i] = seq.substring(i, i+k);
        }
        
        return kmers;
    }
    
    public static final void kmerizeToCollection(String seq, int k, AbstractCollection<String> kmers) {
        final int numKmers = seq.length()-k+1;
        
        for (int i=0; i<numKmers; ++i) {
            kmers.add(seq.substring(i, i+k));
        }
    }
        
    public static final String reverseComplement(String seq) {
        int seqLen = seq.length();
        char[] rc = new char[seqLen];
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int i = seqLen;
        int c;
        
        while (itr.hasNext()) {
            c = itr.nextInt();
            switch(c) {
                case CHAR_A_INT:
                    rc[--i] = 'T';
                    break;
                case CHAR_C_INT:
                    rc[--i] = 'G';
                    break;
                case CHAR_G_INT:
                    rc[--i] = 'C';
                    break;
                case CHAR_T_INT:
                    rc[--i] = 'A';
                    break;
                case CHAR_U_INT:
                    rc[--i] = 'A';
                    break;
                default:
                    rc[--i] = 'N';
            }
        }
        
        return new String(rc);
    }
    
    public static final String cutHairPinLoop(final String seq, final int seedSize, final float minPercentIdentity) {
        int seqLen = seq.length();
        int halfLen = seqLen/2;
        
        final int maxSeedSearchDepth = Math.min(halfLen, 200);
        final int maxLoopLength = Math.max(200, halfLen);
        final int maxLoopDiameter = maxLoopLength/2;
        
        int lastSeedIndex = maxSeedSearchDepth - seedSize;
        
        for (int i=0; i<lastSeedIndex; i+=seedSize) {
            String seed = seq.substring(i, i+seedSize);
            int seedRevCompIndex = seq.indexOf(reverseComplement(seed), i+seedSize);

            if (seedRevCompIndex >= 0) {
                int endIndex = seedRevCompIndex + seedSize;
                int halfIndex = (i + endIndex)/2;

                if (i+seedSize >= seedRevCompIndex-maxLoopLength) {
                   if (halfIndex < halfLen) {
                        return seq.substring(halfIndex);
                    }
                    else {
                        return seq.substring(0, halfIndex);
                    }
                }
                else {
                    String left = seq.substring(i, halfIndex-maxLoopDiameter);
                    String right = seq.substring(halfIndex+maxLoopDiameter, endIndex);

                    float pid = getPercentIdentity(left, reverseComplement(right));
                    
                    if (pid >= minPercentIdentity) {
                        if (halfIndex < halfLen) {
                            return seq.substring(halfIndex);
                        }
                        else {
                            return seq.substring(0, halfIndex);
                        }
                    }
                }
                
                break;
            }
        }
        
        lastSeedIndex = seqLen-maxSeedSearchDepth;
        
        for (int i=seqLen-seedSize; i>=lastSeedIndex; i-=seedSize) {
            String seed = seq.substring(i, i+seedSize);
            int seedRevCompIndex = seq.lastIndexOf(reverseComplement(seed));

            if (seedRevCompIndex >= 0 && seedRevCompIndex < i-seedSize) {
                int halfIndex = (seedRevCompIndex + i + seedSize)/2;
                
                if (seedRevCompIndex+seedSize >= i-maxLoopLength) {
                    if (halfIndex < halfLen) {
                        return seq.substring(halfIndex);
                    }
                    else {
                        return seq.substring(0, halfIndex);
                    }
                }
                else {
                    String left = seq.substring(seedRevCompIndex, halfIndex-maxLoopDiameter);
                    String right = seq.substring(halfIndex+maxLoopDiameter, i+seedSize);

                    float pid = getPercentIdentity(left, reverseComplement(right));
                    
                    if (pid >= minPercentIdentity) {
                        if (halfIndex < halfLen) {
                            return seq.substring(halfIndex);
                        }
                        else {
                            return seq.substring(0, halfIndex);
                        }
                    }
                }
                
                break;
            }
        }
        
        return null;
    }
    
    public static final byte complement(byte b) {
        switch(b) {
            case CHAR_A_BYTE:
                return CHAR_T_BYTE;
            case CHAR_C_BYTE:
                return CHAR_G_BYTE;
            case CHAR_G_BYTE:
                return CHAR_C_BYTE;
            case CHAR_T_BYTE:
                return CHAR_A_BYTE;
            case CHAR_U_BYTE:
                return CHAR_A_BYTE;
        }
        
        return b;
    }
    
    public static final boolean isReverseComplement(byte[] seq1, byte[] seq2) {
        int len = seq1.length;
        
        if (len == 0 || seq2.length != len) {
            return false;
        }
        
        for (int i=0; i<len; ++i) {
            if (complement(seq1[i]) != seq2[len-1-i]) {
                return false;
            }
        }
        
        return true;
    }
    
    public static final byte[] reverseComplement(byte[] seq) {
        int seqLen = seq.length;
        byte[] rc = new byte[seqLen];
        
        int start = seqLen;
        for (byte b : seq) {
            rc[--start] = complement(b);
        }
        
        return rc;
    }
    
    public static final String reverseComplement(byte[] seq, int start, int end) {
        StringBuilder sb = new StringBuilder(end-start);
        
        for (int i=end-1; i>= start; --i) {
            sb.append((char) complement(seq[i]));
        }
        
        return sb.toString();
    }
        
    public static final String smallestStrand(String seq) {
        String rc = reverseComplement(seq);
        
        if (seq.compareTo(rc) > 0) {
            return rc;
        }
        else {
            return seq;
        }
    }
    
    public static final String[] smallestStrand(String seq1, String seq2) {
        String rc2 = reverseComplement(seq2);
        int seq1Vrc2 = seq1.compareTo(rc2);
        
        if (seq1Vrc2 > 0) {
            return new String[]{rc2, reverseComplement(seq1)};
        }
        else if (seq1Vrc2 == 0) {
            String rc1 = reverseComplement(seq1);
            if (seq2.compareTo(rc1) > 0) {
                return new String[]{rc2, rc1};
            }
        }
        
        return new String[]{seq1, seq2};
        
        /*
        String rc1 = reverseComplement(seq1);
        String rc2 = reverseComplement(seq2);
        
        if ((seq1 + seq2).compareTo((rc2 + rc1)) > 0) {
            return new String[]{rc2, rc1};
        }
        else {
            return new String[]{seq1, seq2};
        }
        */
    }
    
    public static String overlapMaximally(String left, String right, int minOverlap) {
        String prefix = right.substring(0, minOverlap);
        int leftLength = left.length();
        int rightLength = right.length();
        
        int lowerL = 0;
        int maxLi = leftLength - minOverlap;
        
        while (lowerL >= 0 && lowerL <= maxLi) {
            lowerL = left.indexOf(prefix, lowerL);
            
            if (lowerL >= 0) {
                int upperL = lowerL+rightLength;
                if (upperL < leftLength) {
                    // could `right` be contained in `left`?
                    String leftTail = left.substring(lowerL+minOverlap, upperL);
                    String rightTail = right.substring(minOverlap);
                    if (leftTail.equals(rightTail)) {
                        // overlap: rightLength
                        return left;
                    }
                }
                else {
                    upperL = leftLength;
                    String leftTail = left.substring(lowerL+minOverlap, upperL);
                    String rightTail = right.substring(minOverlap, upperL-lowerL);
                    if (leftTail.equals(rightTail)) {
                        // overlap: leftLength - lowerL
                        return left + right.substring(upperL-lowerL);
                    }
                }
                ++lowerL;
            }
        }

        if (leftLength >= rightLength && left.contains(right)) {
            return left;
        }
        
        if (leftLength < rightLength && right.contains(left)) {
            return right;
        }
        
        return null;
    }
        
    public static String overlapMinimally(String left, String right, int minOverlap) {
        int li = left.length() - minOverlap;
        String suffix = left.substring(li);
        int ri = 0;
        int maxRi = right.length() - minOverlap;
        
        while (ri >= 0 && ri <= maxRi) {
            ri = right.indexOf(suffix, ri);
            
            if (ri == 0 || (ri > 0 && right.substring(0, ri).equals(left.substring(li-ri, li)))) {
                return left + right.substring(ri+minOverlap);
            }
            else if (ri >= 0) {
                ++ri;
            }
        }
        
        return null;
    }
    
    private static final String PHRED33 = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
    
    public static Pattern getPhred33Pattern(int minQual, int minLength) {
        return Pattern.compile("[\\Q" + PHRED33.substring(minQual) + "\\E]{" + Integer.toString(minLength) + ",}");
    }
    
    public static Pattern getNucleotideCharsPattern(int minLength) {
        return Pattern.compile("[ACGTU]{" + Integer.toString(minLength) + ",}", Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getHomoPolymerPattern(int length) {
        return Pattern.compile("(?:A{" + length + ",})" +
                              "|(?:C{" + length + ",})" +
                              "|(?:G{" + length + ",})" +
                              "|(?:T{" + length + ",})" +
                              "|(?:U{" + length + ",})", Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyATailPattern(int minLength) {
        return Pattern.compile("^.*A{" + minLength + ",}$", Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyATailMatchingPattern(int minLength) {
        return Pattern.compile("A{" + minLength + ",}$", Pattern.CASE_INSENSITIVE);
    }

    public static Pattern getPolyTHeadPattern(int minLength) {
        return Pattern.compile("^T{" + minLength + ",}.*$", Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyTHeadOrPolyATailPattern(int minLength) {
        return Pattern.compile("^T{" + minLength + ",}.*$|^.*A{" + minLength + ",}$", Pattern.CASE_INSENSITIVE);
    }
    
    public static final String[] POLY_A_SIGNALS = {"AATAAA", "ATTAAA", "AGTAAA", "TATAAA",
                                                   "CATAAA", "GATAAA", "AATATA", "AATACA",
                                                   "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA",
                                                   "AATGAA", "TTTAAA", "AAAACA", "GGGGCT",
                                                   "AATAAT", "AACAAA", "ATTACA", "ATTATA",
                                                   "AACAAG", "AATAAG", "TTTTTT"}; // PMID: 27382025

    public static final String[] POLY_A_SIGNALS_REV_COMP = {"TTTATT", "TTTAAT", "TTTACT", "TTTATA",
                                                   "TTTATG", "TTTATC", "TATATT", "TGTATT",
                                                   "TCTATT", "CTTTTT", "TTTAGT", "TTTCTT",
                                                   "TTCATT", "TTTAAA", "TGTTTT", "AGCCCC",
                                                   "ATTATT", "TTTGTT", "TGTAAT", "TATAAT",
                                                   "CTTGTT", "CTTATT", "AAAAAA"};
    
    public static Pattern getPolyASignalPattern() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS), Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyASignalReverseComplementPattern() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS_REV_COMP), Pattern.CASE_INSENSITIVE);
    }
    
    public static ArrayDeque<Integer> getPolyASignalPositions(String seq, Pattern pasPattern, Pattern tailPattern) {
        ArrayDeque<Integer> positions = new ArrayDeque<>();
        
        Matcher tailMatcher = tailPattern.matcher(seq);
        if (tailMatcher.find()) {
            int tailStartPos = tailMatcher.start();
            
            if (tailStartPos > 10) {
                Matcher pasMatcher = pasPattern.matcher(seq);
                pasMatcher.region(Math.max(0, tailStartPos-60), tailStartPos-10);

                while (pasMatcher.find()) {
                    positions.add(pasMatcher.start());
                }
            }
        }
        
        return positions;
    }
    
    public static int[] getPolyATailRegion(String seq, int searchRange, int window, int minWindowMatch, int seedSize) {
//        int searchRange = 100;
//        int window = 17;
//        int minWindowMatch = 15;
//        int seedSize = 2;
        
        int start = -1;
        int end = -1;
        
        int seqLen = seq.length();
        if (seqLen < window || searchRange < window) {
            return null;
        }
        
        int numMatch = 0;
        for (int i=seqLen-1; i>=seqLen-window; --i) {
            if (seq.charAt(i) == 'A') {
                ++numMatch;
            }
        }
        
        if (numMatch >= minWindowMatch) {
            start = seqLen-window;
            end = seqLen-1;
        }
        
        if (numMatch > 0 && seq.charAt(seqLen-1) == 'A') {
            --numMatch;
        }
        
        int min = Math.max(0, seqLen- searchRange);
        for (int i=seqLen-window-1; i>=min; --i) {            
            if (seq.charAt(i) == 'A') {
                ++numMatch;
            }
            
            if (numMatch >= minWindowMatch) {
                if (start < 0) {
                    start = i;
                    end = i+window;
                }
                else if (i+window>start) {
                    start = i;
                }
                else {
                    break;
                }
            }
            
            if (seq.charAt(i+window) == 'A') {
                --numMatch;
            }
        }
        
        if (start<0 || end<0) {
            return null;
        }
        
        for (int i=start; i<end; ++i) {
            boolean seedFound = false;
            
            for (int j=0; j<seedSize; ++j) {
                seedFound = true;
                if (seq.charAt(i+j) != 'A') {
                    seedFound = false;
                    i += j;
                    break;
                }
            }
            
            if (seedFound) {
                start = i;
                break;
            }
        }
        
        for (int i=end; i>start; --i) {
            boolean seedFound = false;

            for (int j=0; j<seedSize; ++j) {
                seedFound = true;
                if (seq.charAt(i-j) != 'A') {
                    seedFound = false;
                    i -= j;
                    break;
                }
            }
            
            if (seedFound) {
                end = i;
                break;
            }
        }
        
        return new int[]{start, end};
    }
    
    public static int[] getPolyTHeadRegion(String seq, int searchRange, int window, int minWindowMatch, int seedSize) {
        int start = -1;
        int end = -1;
        
        int seqLen = seq.length();
        if (seqLen < window || searchRange < window) {
            return null;
        }
        
        int numMatch = 0;
        for (int i=0; i<window; ++i) {
            if (seq.charAt(i) == 'T') {
                ++numMatch;
            }
        }
        
        if (numMatch >= minWindowMatch) {
            start = 0;
            end = window;
        }
        
        if (numMatch > 0 && seq.charAt(0) == 'T') {
            --numMatch;
        }
        
        int max = Math.min(seqLen, searchRange);
        for (int i=window; i<max; ++i) {            
            if (seq.charAt(i) == 'T') {
                ++numMatch;
            }
            
            if (numMatch >= minWindowMatch) {
                if (start < 0) {
                    start = i-window;
                    end = i;
                }
                else if (i-window<end) {
                    end = i;
                }
                else {
                    break;
                }
            }
            
            if (seq.charAt(i-window) == 'T') {
                --numMatch;
            }
        }
        
        if (start<0 || end<0) {
            return null;
        }
        
        for (int i=start; i<end; ++i) {
            boolean seedFound = false;
            
            for (int j=0; j<seedSize; ++j) {
                seedFound = true;
                if (seq.charAt(i+j) != 'T') {
                    seedFound = false;
                    i += j;
                    break;
                }
            }
            
            if (seedFound) {
                start = i;
                break;
            }
        }
        
        for (int i=end; i>start; --i) {
            boolean seedFound = false;
            
            for (int j=0; j<seedSize; ++j) {
                seedFound = true;
                if (seq.charAt(i-j) != 'T') {
                    seedFound = false;
                    i -=j ;
                    break;
                }
            }
            
            if (seedFound) {
                end = i;
                break;
            }
        }
        
        return new int[]{start, end};
    }
    
    public static boolean isHomoPolymer(String seq) {
        char c = seq.charAt(0);
        
        int len = seq.length();
        for (int i=1; i<len; ++i) {
            if (seq.charAt(i) != c) {
                return false;
            }
        }
        
        return true;
    }
    
    public static String compressHomoPolymers(String seq) {
        if (seq.isEmpty()) {
            return seq;
        }
        
        int len = seq.length();
        
        StringBuilder sb = new StringBuilder(len);
        
        char lastN = seq.charAt(0);
        sb.append(lastN);
        
        char n;
        for (int i=1; i<len; ++i) {
            n = seq.charAt(i);
            if (n != lastN) {
                sb.append(n);
                lastN = n;
            }
        }
        
        return sb.toString();
    }

    public static String longestSeq(String seq, Pattern seqPattern) {
        // filter sequence by keeping only ACGT characters
        Matcher m = seqPattern.matcher(seq);
        
        int longestStartPos = -1;
        int longestEndPos = -1;
        int longestLen = -1;
        while (m.find()) {
            int startPos = m.start();
            int endPos = m.end();
            int len = endPos - startPos;

            if (len > longestLen) {
                longestLen = len;
                longestStartPos = startPos;
                longestEndPos = endPos;
            }
        }
        
        if (longestStartPos < 0 || longestEndPos < 0 || longestLen < 0) {
            return "";
        }
        
        return seq.substring(longestStartPos, longestEndPos).toUpperCase();
    }
    
    public static String longestSeq(String seq, String qual, Pattern seqPattern, Pattern qualPattern) {
        // filter sequence by keeping only ACGT characters
        Matcher m = qualPattern.matcher(qual);
        
        int longestStartPos = -1;
        int longestEndPos = -1;
        int longestLen = -1;
        while (m.find()) {
            int startPos = m.start();
            int endPos = m.end();
            int len = endPos - startPos;

            if (len > longestLen) {
                longestLen = len;
                longestStartPos = startPos;
                longestEndPos = endPos;
            }
        }
        
        if (longestStartPos < 0 || longestEndPos < 0 || longestLen < 0) {
            return "";
        }
        
        return longestSeq(seq.substring(longestStartPos, longestEndPos), seqPattern);
    }
    
    public static ArrayList<String> filterFasta(String seq, Pattern seqPattern) {

        // filter sequence by keeping only ACGT characters
        Matcher m = seqPattern.matcher(seq);
        ArrayList<String> result = new ArrayList<>();

        int startPos;
        int len;
        int endPos = 0;
        while (m.find()) {
            startPos = m.start();
            
            if (endPos > 0) {
                // ensure that first item is not 'N'
                
                len = startPos - endPos;
                
                if (len == 1) {
                    result.add("N");
                }
                else {
                    char[] gap = new char[len];
                    Arrays.fill(gap, GAP_CHAR);
                    result.add(new String(gap));
                }
            }
            
            endPos = m.end();
            result.add(seq.substring(startPos, endPos).toUpperCase());
        }
        
        return result;
    }
    
    public static ArrayList<String> filterFastq(FastqRecord fq, Pattern seqPattern, Pattern qualPattern) {
        return filterFastq(fq.seq, fq.qual, seqPattern, qualPattern);
    }
    
    public static ArrayList<String> filterFastq(String seq, String qual, Pattern seqPattern, Pattern qualPattern) {
        // filter sequence by quality
        Matcher m = qualPattern.matcher(qual);
        StringBuilder qualFilteredSeq = new StringBuilder(seq.length());
        
        int startPos;
        int len;
        int endPos = 0;
        while (m.find()) {
            startPos = m.start();
            
            if (endPos > 0) {
                // ensure that first item is not 'N'
                
                len = startPos - endPos;

                if (len == 1) {
                    qualFilteredSeq.append('N');
                }
                else {
                    char[] gap = new char[len];
                    Arrays.fill(gap, GAP_CHAR);
                    qualFilteredSeq.append(gap);
                }
            }
            
            endPos = m.end();
            qualFilteredSeq.append(seq.substring(startPos, endPos).toUpperCase());
        }
        
        // filter sequence by keeping only ACGT characters
        seq = qualFilteredSeq.toString();
        m = seqPattern.matcher(seq);
        ArrayList<String> result = new ArrayList<>();

        endPos = 0;
        while (m.find()) {
            startPos = m.start();
            
            if (endPos > 0) {
                // ensure that first item is not 'N'
                
                len = startPos - endPos;
                
                if (len == 1) {
                    result.add("N");
                }
                else {
                    char[] gap = new char[len];
                    Arrays.fill(gap, GAP_CHAR);
                    result.add(new String(gap));
                }
            }
            
            endPos = m.end();
            result.add(seq.substring(startPos, endPos).toUpperCase());
        }
        
        return result;
    }
        
    public static void main(String[] args) {
        //debug
        String seq = "";
        ArrayList<String> segments = trimLowComplexityRegions(seq, 500);
        System.out.println(segments.size());
        for (String seg : segments) {
            System.out.println(seg);
        }
        System.out.println(isLowComplexityLongWindowed(seq));
    }
}
