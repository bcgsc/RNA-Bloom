/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;

/**
 *
 * @author gengar
 */
public class HashFunction2 {
    private final int k;
    
    public HashFunction2(int k) {
        this.k = k;
    }
    
    public void getHashValues(final CharSequence kmer,
                              final int numHash,
                              final long[] out) {
        NTM64(kmer, k, numHash, out);
    }
    
    public NTHashIterator getHashIterator(final CharSequence seq,
                                          final int numHash) {
        return new NTHashIterator(k, numHash);
    }
}
