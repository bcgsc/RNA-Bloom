/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
 * @author kmnip
 */
public final class SeqUtils {
    public static final char GAP_CHAR = 'N';
    public final static char[] NUCLEOTIDES = new char[] {'A','C','G','T'};
    public final static char[] A_ALT_NUCLEOTIDES = new char[] {'C','G','T'};
    public final static char[] C_ALT_NUCLEOTIDES = new char[] {'A','G','T'};
    public final static char[] G_ALT_NUCLEOTIDES = new char[] {'A','C','T'};
    public final static char[] T_ALT_NUCLEOTIDES = new char[] {'A','C','G'};
    
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
            default:
                return NUCLEOTIDES;
        }
    }
    
    public static float getPercentIdentity(String a, String b) {
        int aLen = a.length();
        int bLen = b.length();

        int d = getDistance(a, b, aLen, bLen);
        
        if (aLen <= bLen) {
            return ((float) (bLen - d))/bLen;
        }
        
        return ((float) (aLen - d))/aLen;
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
//            for (int j=0; j<=tLen; ++j) {
//                v0[j] = v1[j];
//            }
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
        
    private final static int CHAR_A_INT = (int) 'A';
    private final static int CHAR_C_INT = (int) 'C';
    private final static int CHAR_G_INT = (int) 'G';
    private final static int CHAR_T_INT = (int) 'T';
        
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
                default:
                    rc[--i] = 'N';
            }
        }
        
        return new String(rc);
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
        
        if (leftLength < rightLength && right.contains(left)) {
            // overlap: leftLength
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
        return Pattern.compile("[ACTG]{" + Integer.toString(minLength) + ",}");
    }
    
    public static Pattern getHomoPolymerPattern(int length) {
        return Pattern.compile("(?:A{" + length + "})" +
                              "|(?:C{" + length + "})" +
                              "|(?:G{" + length + "})" +
                              "|(?:T{" + length + "})");
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
        
    public static ArrayList<String> filterFastq(FastqRecord fq, Pattern qualPattern) {
        Matcher m = qualPattern.matcher(fq.qual);
        String seq = fq.seq;
        
        ArrayList<String> result = new ArrayList<>();
        
        int startPos;
        int len;
        int endPos = 0;
        while (m.find()) {
            startPos = m.start();
            
            if (endPos > 0) {
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
            result.add(seq.substring(startPos, endPos));
        }
        
        return result;
    }
       
//    public static void main(String[] args) {
//        String seq1 = "AATAA";
//        String seq2 = "AAAA";
//
//        System.out.println(getDistance(seq1, seq2));
//    }
}
