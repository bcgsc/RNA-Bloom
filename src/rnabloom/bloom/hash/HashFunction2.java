/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.bloom.hash.NTHash.NTM64B;
import static rnabloom.util.SeqUtils.NUCLEOTIDES;


/**
 *
 * @author gengar
 */
public class HashFunction2 {
    protected final int k;
    protected final int kMod64;
    protected final int kMinus1Mod64;
    
    public HashFunction2(int k) {
        this.k = k;
        this.kMod64 = k%64;
        this.kMinus1Mod64 = (k-1)%64;
    }
    
    public void getHashValues(final CharSequence kmer,
                              final int numHash,
                              final long[] out) {
        NTM64(kmer, k, numHash, out);
    }
    
    public long[][] getSuccessorsHashValues(final int numHash, final long[] hVals, final char leftMostNucleotide) {
        return NTM64(leftMostNucleotide, NUCLEOTIDES, k, numHash, hVals, kMod64);
    }

    public long[][] getPredecessorsHashValues(final int numHash, final long[] hVals, final char rightMostNucleotide) {
        return NTM64B(rightMostNucleotide, NUCLEOTIDES, k, numHash, hVals, kMinus1Mod64);
    }    
    
    public NTHashIterator getHashIterator(final int numHash) {
        return new NTHashIterator(k, numHash);
    }
    
    public long[] getHashValues(final CharSequence kmer1,
                                final CharSequence kmer2,
                                final int numHash) {
        final long[] hashVals1 = new long[numHash];
        //murmurhash3_x64_128(kmer1.getBytes(), 0, k, seed, numHash, hashVals1);
        getHashValues(kmer1, numHash, hashVals1);
        
        final long[] hashVals2 = new long[numHash];
        //murmurhash3_x64_128(kmer2.getBytes(), 0, k, seed, numHash, hashVals2);
        getHashValues(kmer2, numHash, hashVals2);
        
        final long[] hashVal = new long[numHash];
        for (int i=0; i<numHash; ++i) {
            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
        }
        
        return hashVal;
    }
    
    public long[] getHashValues(final long[] hashVals1,
                                final long[] hashVals2,
                                final int numHash) {
        final long[] hashVal = new long[numHash];
        for (int i=0; i<numHash; ++i) {
            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
        }
        
        return hashVal;
    }
    
    public static long combineHashValues(long a, long b) {
        // See: http://stackoverflow.com/a/27952689
        return a ^ (b + 0x9e3779b9 + (a << 6) + (b >>> 2));
        
//        a ^= b + 0x9e3779b9 + (a << 6) + (b >> 2);
//        return a < 0 ? -a : a;
    }
}
