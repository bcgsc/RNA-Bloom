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
import rnabloom.graph.CanonicalKmer;
import rnabloom.graph.Kmer;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.bloom.hash.NTHash.NTMC64;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalHashFunction extends HashFunction {
    public CanonicalHashFunction(int k) {
        super(k);
    }
    
    @Override
    public Kmer getKmer(final String kmer, final int numHash, BloomFilterDeBruijnGraph graph) {
        long[] frhval = new long[2];
        long[] hVals = new long[numHash];
        NTMC64(kmer, k, numHash, frhval, hVals);
        return new CanonicalKmer(kmer, k, graph.getCount(hVals), frhval[0], frhval[1]);
    }
    
    @Override
    public ArrayList<Kmer> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayList<Kmer> result = new ArrayList<>();
        
        int seqLength = seq.length();
        
        if (seqLength >= k) {
            byte[] bytes = stringToBytes(seq, seqLength);

            CanonicalNTHashIterator itr = new CanonicalNTHashIterator(k, numHash);
            itr.start(seq);
            long[] hVals = itr.hVals;
            long[] frhval = itr.frhval;
            int i;
            float c;
            while (itr.hasNext()) {
                itr.next();
                c = graph.getCount(hVals);
                
                if (c <= 0) {
                    // not a valid sequence
                    return new ArrayList<>();
                }
                
                i = itr.getPos();
                result.add(new CanonicalKmer(Arrays.copyOfRange(bytes, i, i+k), c, frhval[0], frhval[1]));
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayList<Kmer> getKmers(final String seq, final int numHash, BloomFilterDeBruijnGraph graph, float minCoverage) {
        int seqLength = seq.length();
        
        ArrayList<Kmer> currentSegment = new ArrayList<>();
        ArrayList<Kmer> longestSegment = currentSegment;
        
        if (seqLength >= k) {
            float currentMinC = Float.MAX_VALUE;
            float longestMinC = Float.MAX_VALUE;
            int longestLen = 0;
            
            byte[] bytes = stringToBytes(seq, seqLength);

            CanonicalNTHashIterator itr = new CanonicalNTHashIterator(k, numHash);
            itr.start(seq);
            long[] hVals = itr.hVals;
            long[] frhval = itr.frhval;
            int i;
            float c;
            while (itr.hasNext()) {
                itr.next();
                c = graph.getCount(hVals);
                
                if (c >= minCoverage) {
                    i = itr.getPos();
                    currentSegment.add(new CanonicalKmer(Arrays.copyOfRange(bytes, i, i+k), c, frhval[0], frhval[1]));
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
