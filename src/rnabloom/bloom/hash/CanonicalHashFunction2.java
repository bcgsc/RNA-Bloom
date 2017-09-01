/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import java.util.ArrayList;
import java.util.Arrays;
import static rnabloom.bloom.hash.NTHash.NTMC64;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.CanonicalKmer;
import rnabloom.graph.Kmer2;
import static rnabloom.util.SeqUtils.stringToBytes;

/**
 *
 * @author Ka Ming Nip, Genome Sciences Centre, BC Cancer Agency
 */
public class CanonicalHashFunction2 extends HashFunction2 {
    public CanonicalHashFunction2(int k) {
        super(k);
    }
    
    @Override
    public Kmer2 getKmer(final String kmer, final int numHash, BloomFilterDeBruijnGraph graph) {
        long[] frhval = new long[2];
        long[] hVals = new long[numHash];
        NTMC64(kmer, k, numHash, frhval, hVals);
        return new CanonicalKmer(kmer, k, graph.getCount(hVals), frhval[0], frhval[1]);
    }
    
    @Override
    public ArrayList<Kmer2> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayList<Kmer2> result = new ArrayList<>();
        
        byte[] bytes = stringToBytes(seq, seq.length());
        
        CanonicalNTHashIterator itr = new CanonicalNTHashIterator(k, numHash);
        itr.start(seq);
        long[] hVals = itr.hVals;
        long[] frhval = itr.frhval;
        int i;
        while (itr.hasNext()) {
            itr.next();
            i = itr.getPos();
            result.add(new CanonicalKmer(Arrays.copyOfRange(bytes, i, i+k), graph.getCount(hVals), frhval[0], frhval[1]));
        }
        
        return result;
    }
    
    @Override
    public void getHashValues(final String kmer,
                              final int numHash,
                              final long[] out) {
        NTMC64(kmer, k, numHash, out);
    }
        
    @Override
    public NTHashIterator getHashIterator(final int numHash) {
        return new CanonicalNTHashIterator(k, numHash);
    }
    
    @Override
    public PairedNTHashIterator getPairedHashIterator(final int numHash, final int distance) {
        return new CanonicalPairedNTHashIterator(k, numHash, distance);
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
