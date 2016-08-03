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
public final class SequenceOperations {
    
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
                    rc[--i] = (char) c;
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
    
    private static final String PHRED33 = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ";
    
    public static Pattern getPhred33Pattern(int minQual, int minLength) {
        return Pattern.compile("[" + PHRED33.substring(minQual) + "]{" + Integer.toString(minLength) + ",}");
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
