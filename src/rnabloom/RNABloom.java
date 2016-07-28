/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import rnabloom.bloom.CountingBloomFilter;
import static rnabloom.util.SequenceOperations.*;
import static java.lang.Math.pow;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    public final static long NUM_BYTES_1GB = (long) pow(1024, 3);
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        String kmer = "AAAAATTTTTCCCCCGGGGG";
        System.out.println(smallestStrand(kmer));
    }
    
}
