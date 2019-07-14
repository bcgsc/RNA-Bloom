/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom.bloom.hash;

import java.util.ArrayList;
import java.util.Arrays;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.bloom.hash.NTHash.NTM64;


/**
 *
 * @author Ka Ming Nip
 */
public class HashFunction {
    protected int k;
    protected int kMod64;
    protected int kMinus1Mod64;
    
    public HashFunction(int k) {
        this.k = k;
        this.kMod64 = k%64;
        this.kMinus1Mod64 = (k-1)%64;
    }
    
    public void setK(int k) {
        this.k = k;
        this.kMod64 = k%64;
        this.kMinus1Mod64 = (k-1)%64;
    }
    
    public Kmer getKmer(final String kmer, final int numHash, BloomFilterDeBruijnGraph graph) {
        long[] hVals = new long[numHash];
        NTM64(kmer, k, numHash, hVals);
        return new Kmer(kmer, k, graph.getCount(hVals), hVals[0]);
    }
    
    public ArrayList<Kmer> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayList<Kmer> result = new ArrayList<>();
        
        int seqLength = seq.length();
        
        if (seqLength >= k) {
            byte[] bytes = stringToBytes(seq, seqLength);

            NTHashIterator itr = new NTHashIterator(k, numHash);
            itr.start(seq);
            long[] hVals = itr.hVals;
            int i;
            float c;
            while (itr.hasNext()) {                
                itr.next();
                c = graph.getCount(hVals);
                
//                if (c <= 0) {
//                    // not a valid sequence
//                    return new ArrayList<>();
//                }
                
                i = itr.getPos();
                result.add(new Kmer(Arrays.copyOfRange(bytes, i, i+k), c, hVals[0]));
            }
        }
        
        return result;
    }

    public ArrayList<Kmer> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph, float minCoverage) {       
        int seqLength = seq.length();
        
        ArrayList<Kmer> currentSegment = new ArrayList<>();
        ArrayList<Kmer> longestSegment = currentSegment;
                
        if (seqLength >= k) {
            float currentMinC = Float.MAX_VALUE;
            float longestMinC = Float.MAX_VALUE;
            int longestLen = 0;
            
            byte[] bytes = stringToBytes(seq, seqLength);

            NTHashIterator itr = new NTHashIterator(k, numHash);
            itr.start(seq);
            long[] hVals = itr.hVals;
            int i;
            float c;
            while (itr.hasNext()) {                
                itr.next();
                c = graph.getCount(hVals);
                
                if (c >= minCoverage) {
                    i = itr.getPos();
                    currentSegment.add(new Kmer(Arrays.copyOfRange(bytes, i, i+k), c, hVals[0]));
                    currentMinC = Math.min(currentMinC, c);
                }
                else if (!currentSegment.isEmpty()) {
                    if (longestSegment != currentSegment) {
                        int len = currentSegment.size();

                        if (len > longestLen || (len == longestLen && currentMinC > longestMinC)) {
                            longestSegment = currentSegment;
                            longestMinC = currentMinC;
                            longestLen = len;
                        }
                    }
                    
                    currentSegment = new ArrayList<>();
                    currentMinC = Float.MAX_VALUE;
                }
            }
        }
        
        return longestSegment;
    }
        
    public ArrayList<Kmer> getKmers(final String seq, final int start, final int end, final int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayList<Kmer> result = new ArrayList<>();
        
        int seqLength = seq.length();
        
        if (seqLength >= k) {
            byte[] bytes = stringToBytes(seq, seqLength);

            NTHashIterator itr = new NTHashIterator(k, numHash);
            itr.start(seq, start, end);
            long[] hVals = itr.hVals;
            int i;
            while (itr.hasNext()) {
                itr.next();
                i = itr.getPos();
                result.add(new Kmer(Arrays.copyOfRange(bytes, i, i+k), graph.getCount(hVals), hVals[0]));
            }
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
        
//    public long[] getHashValues(final String kmer1,
//                                final String kmer2,
//                                final int numHash) {
//        final long[] hashVals1 = new long[numHash];
//        //murmurhash3_x64_128(kmer1.getBytes(), 0, k, seed, numHash, hashVals1);
//        getHashValues(kmer1, numHash, hashVals1);
//        
//        final long[] hashVals2 = new long[numHash];
//        //murmurhash3_x64_128(kmer2.getBytes(), 0, k, seed, numHash, hashVals2);
//        getHashValues(kmer2, numHash, hashVals2);
//        
//        return getHashValues(hashVals1, hashVals2, numHash);
//    }
//    
//    public long[] getHashValues(final long[] hashVals1,
//                                final long[] hashVals2,
//                                final int numHash) {
//        
//        return getHashValues(combineHashValues(hashVals1[0], hashVals2[0]), numHash);
//    }
    
    public static long combineHashValues(long a, long b) {
        // See: http://stackoverflow.com/a/27952689
        return a ^ (b + 0x9e3779b9 + (a << 6) + (b >>> 2));
        
//        a ^= b + 0x9e3779b9 + (a << 6) + (b >> 2);
//        return a < 0 ? -a : a;
    }
}
