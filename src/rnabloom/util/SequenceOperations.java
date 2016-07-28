/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.PrimitiveIterator;

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
    
    public static final String reverseComplement(String seq) {
        int seqLen = seq.length();
        char[] rc = new char[seqLen];
        
        PrimitiveIterator.OfInt itr = seq.chars().iterator();
        int i = seqLen;
        
        while (itr.hasNext()) {
            char c = (char) itr.nextInt();
            switch(c) {
                case 'A':
                    rc[--i] = 'T';
                    break;
                case 'C':
                    rc[--i] = 'G';
                    break;
                case 'G':
                    rc[--i] = 'C';
                    break;
                case 'T':
                    rc[--i] = 'A';
                    break;
                default:
                    rc[--i] = c;
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
    
}
