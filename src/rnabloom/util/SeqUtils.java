/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.PrimitiveIterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.io.FastqRecord;

/**
 *
 * @author kmnip
 */
public final class SeqUtils {
    
    public static final String[] kmerize(String seq, int k) {
        final int numKmers = seq.length() - k +1;
        final String[] kmers = new String[numKmers];
        
        for (int i=0; i<numKmers; ++i) {
            kmers[i] = seq.substring(i, i+k);
        }
        
        return kmers;
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
    
    public static Pattern getACTGPattern(int minLength) {
        return Pattern.compile("[ACTG]{" + Integer.toString(minLength) + ",}");
    }
        
    public static ArrayList<String> filterFastq(FastqRecord fq, Pattern qualPattern) {
        Matcher m = qualPattern.matcher(fq.qual);
        String seq = fq.seq;
        
        ArrayList<String> result = new ArrayList<>();
        
        while (m.find()) {
            result.add(seq.substring(m.start(), m.end()));
        }
        
        return result;
    }
        
}
