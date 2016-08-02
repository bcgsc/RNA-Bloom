/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import rnabloom.bloom.CountingBloomFilter;
import static rnabloom.util.SequenceOperations.*;
import static java.lang.Math.pow;
import rnabloom.bloom.BloomFilter;

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
        String kmer2 = "AACAATTTTTCCCCCGGGGG";
        
        BloomFilter f1 = new BloomFilter(NUM_BITS_1GB, 3, 689, kmer.length());
        f1.add(kmer);
        System.out.println(f1.lookup(kmer));
        System.out.println(f1.lookup(kmer2));
        
        CountingBloomFilter f2 = new CountingBloomFilter((int) NUM_BYTES_1GB, 3, 689, kmer.length());
        
        for (int i=0; i < 20; ++i) {
            f2.increment(kmer);
        }
        
        System.out.println(f2.getCount(kmer));
    }
    
}
