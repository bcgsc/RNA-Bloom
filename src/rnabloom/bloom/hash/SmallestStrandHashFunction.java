/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.bloom.hash.MurmurHash3.murmurhash3_x64_128;
import static rnabloom.util.SeqUtils.smallestStrand;
import static rnabloom.bloom.hash.MurmurHash3.murmurhash3_x64_128;
import static rnabloom.util.SeqUtils.smallestStrand;

/**
 *
 * @author kmnip
 */
public class SmallestStrandHashFunction extends HashFunction {
    
    public SmallestStrandHashFunction(int seed, int k, int maxNumVals) {
        super(seed, k, maxNumVals);
    }
    
    @Override
    public long[] getHashValues(final CharSequence kmer, int numHash) {
        return super.getHashValues(smallestStrand(kmer.toString()), numHash);
    }
    
    @Override
    public void getHashValues(final CharSequence kmer, final int numHash, long[] out) {
        super.getHashValues(smallestStrand(kmer.toString()), numHash, out);
    }
    
    @Override
    public long[] getHashValues(final CharSequence kmer1, final CharSequence kmer2, int numHash) {
        String[] reorientedKmers = smallestStrand(kmer1.toString(), kmer2.toString());
        
        final long[] hashVals1 = new long[numHash];
        murmurhash3_x64_128(reorientedKmers[0].getBytes(), 0, k, seed, numHash, hashVals1);
        
        final long[] hashVals2 = new long[numHash];
        murmurhash3_x64_128(reorientedKmers[1].getBytes(), 0, k, seed, numHash, hashVals2);
        
        final long[] hashVal = new long[numHash];
        for (int i=0; i<numHash; ++i) {
            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
        }
        
        return hashVal;
    }
    
    /**@TODO long[][] getHashValues(final CharSequence cs)*/
    
}
