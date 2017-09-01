/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import java.util.ArrayList;
import java.util.Arrays;
import static rnabloom.bloom.hash.NTHash.NTM64;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer2;
import static rnabloom.util.SeqUtils.stringToBytes;


/**
 *
 * @author Ka Ming Nip, Genome Sciences Centre, BC Cancer Agency
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
    
    public Kmer2 getKmer(final String kmer, final int numHash, BloomFilterDeBruijnGraph graph) {
        long[] hVals = new long[numHash];
        NTM64(kmer, k, numHash, hVals);
        return new Kmer2(kmer, k, graph.getCount(hVals), hVals[0]);
    }
    
    public ArrayList<Kmer2> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayList<Kmer2> result = new ArrayList<>();
        
        byte[] bytes = stringToBytes(seq, seq.length());
        
        NTHashIterator itr = new NTHashIterator(k, numHash);
        itr.start(seq);
        long[] hVals = itr.hVals;
        int i;
        while (itr.hasNext()) {
            itr.next();
            i = itr.getPos();
            result.add(new Kmer2(Arrays.copyOfRange(bytes, i, i+k), graph.getCount(hVals), hVals[0]));
        }
        
        return result;
    }
    
    public void getHashValues(final String kmer,
                              final int numHash,
                              final long[] out) {
        NTM64(kmer, k, numHash, out);
    }
    
    /**
     * Generate multiple hash values using the base hash value
     * @param hVal - the base hash value
     * @param numHash - the number of hash values in total, including the base value
     * @return - array of hash values
     */
    public long[] getHashValues(final long hVal, final int numHash) {
        long[] hVals = new long[numHash];
        NTM64(hVal, hVals, k, numHash);
        return hVals;
    }
            
//    public long[][] getSuccessorsHashValues(final int numHash, final long[] hVals, final char leftMostNucleotide) {
//        return NTM64(leftMostNucleotide, NUCLEOTIDES, k, numHash, hVals, kMod64);
//    }
//
//    public long[][] getPredecessorsHashValues(final int numHash, final long[] hVals, final char rightMostNucleotide) {
//        return NTM64B(rightMostNucleotide, NUCLEOTIDES, k, numHash, hVals, kMinus1Mod64);
//    }    
    
    public NTHashIterator getHashIterator(final int numHash) {
        return new NTHashIterator(k, numHash);
    }
    
    public PairedNTHashIterator getPairedHashIterator(final int numHash, final int distance) {
        return new PairedNTHashIterator(k, numHash, distance);
    }
    
//    public SuccessorsNTHashIterator getSuccessorsHashIterator(final int numHash) {
//        return new SuccessorsNTHashIterator(k, numHash);
//    }
//    
//    public PredecessorsNTHashIterator getPredecessorsNTHashIterator(final int numHash) {
//        return new PredecessorsNTHashIterator(k, numHash);
//    }
//    
//    public LeftVariantsNTHashIterator getLeftVariantsNTHashIterator(final int numHash) {
//        return new LeftVariantsNTHashIterator(k, numHash);
//    }
//    
//    public RightVariantsNTHashIterator getRightVariantsNTHashIterator(final int numHash) {
//        return new RightVariantsNTHashIterator(k, numHash);
//    }
        
    public long[] getHashValues(final String kmer1,
                                final String kmer2,
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
