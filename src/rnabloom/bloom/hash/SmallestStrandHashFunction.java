/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.util.SequenceOperations.smallestStrand;

/**
 *
 * @author kmnip
 */
public class SmallestStrandHashFunction extends HashFunction {
    
    public SmallestStrandHashFunction(int numHash, int seed, int k) {
        super(numHash, seed, k);
    }

    public long[] getHashValues(final String kmer) {
        return super.getHashValues(smallestStrand(kmer));
    }    
}
