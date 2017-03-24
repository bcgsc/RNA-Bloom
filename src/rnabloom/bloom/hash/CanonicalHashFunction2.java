/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.NTHash.NTMC64;
import static rnabloom.util.SeqUtils.NUCLEOTIDES;
import static rnabloom.util.SeqUtils.smallestStrand;

/**
 *
 * @author gengar
 */
public class CanonicalHashFunction2 extends HashFunction2 {
    public CanonicalHashFunction2(int k) {
        super(k);
    }
    
    @Override
    public void getHashValues(final CharSequence kmer,
                              final int numHash,
                              final long[] out) {
        NTMC64(kmer, k, numHash, out);
    }
    
    @Override
    public NTHashIterator getHashIterator(final int numHash) {
        return new CanonicalNTHashIterator(k, numHash);
    }
        
//    @Override
//    public long[] getHashValues(final CharSequence kmer1, final CharSequence kmer2, int numHash) {
//        String[] reorientedKmers = smallestStrand(kmer1.toString(), kmer2.toString());
//        
//        final long[] hashVals1 = new long[numHash];
//        super.getHashValues(reorientedKmers[0], numHash, hashVals1);
//        
//        final long[] hashVals2 = new long[numHash];
//        super.getHashValues(reorientedKmers[1], numHash, hashVals2);
//        
//        final long[] hashVal = new long[numHash];
//        for (int i=0; i<numHash; ++i) {
//            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
//        }
//        
//        return hashVal;
//    }
}
