/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom;

import rnabloom.bloom.BloomFilter;
import static java.lang.Math.pow;

/**
 *
 * @author kmnip
 */
public class RNABloom {
    
    
    public final static long NUM_BITS_1GB = (long) pow(1024, 3) * 8;
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        long size = NUM_BITS_1GB * 4;
        int num_hash = 4;
        int seed = 689;
        int key_length = 16;
        
        BloomFilter b = new BloomFilter(size, num_hash, seed, key_length);
        
        String key = "1234567890123456";
        
        b.add(key);
        
        System.out.println(b.lookup(key));
        System.out.println(b.lookup("0000000000123456"));
    }
    
}
