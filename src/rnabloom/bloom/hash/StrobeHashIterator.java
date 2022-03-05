/* 
 * Copyright (C) 2021-present BC Cancer Genome Sciences Centre
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

import java.util.Arrays;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;

/**
 *
 * @author Ka Ming Nip
 */
public class StrobeHashIterator implements StrobeHashIteratorInterface {
    private final int wMin;
    private final int wMax;
    private final int k;
    private final NTHashIterator itr;
    private int pos = -1;
    private int max = -2;
    private int numKmers = 0;
    private long[] kmerHashVals = null;
    private final int n;
    private int[] strobes = null;
    
    public StrobeHashIterator(int n, int k, int wMin, int wMax) {
        this.n = n;
        this.strobes = new int[n-1];
        this.wMin = wMin;
        this.wMax = wMax;
        this.k = k;
        this.itr = new NTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        Arrays.fill(strobes, 0);
        max = -2;
        kmerHashVals = null;
        if (itr.start(seq)) {
            numKmers = seq.length() - k + 1;
            if (numKmers > wMax * (n-1)) {
                kmerHashVals = new long[numKmers];
                max = numKmers - wMax * (n-2) - wMin -1;
                for (int i=0; i<numKmers; ++i) {
                    itr.next();
                    kmerHashVals[i] = itr.hVals[0];
                }
                
                return true;
            }
        }
        return false;
    }
    
    @Override
    public boolean hasNext() {
        return pos < max;
    }
    
    @Override
    public long next() {
        ++pos;        
        long strobemerHash = kmerHashVals[pos];
        
        for (int s=0; s<n-1; ++s) {
            int pos2 = pos + s*wMax + wMin;
            long h = combineHashValues(strobemerHash, kmerHashVals[pos2]);
            
            int end = Math.min(pos + s*wMax + wMax, numKmers);
            for (int i=pos2 +1; i<end; ++i) {
                long h2 = combineHashValues(strobemerHash, kmerHashVals[i]);
                if (Long.compareUnsigned(h, h2) >= 0) {
                    pos2 = i;
                    h = h2;
                }
            }
            
            strobemerHash = h;
            strobes[s] = pos2;
        }
        
        return strobemerHash;
    }
    
    @Override
    public HashedPositions get(int p) {
        long strobemerHash = kmerHashVals[p];
        int[] positions = new int[n];
        positions[0] = p;
        
        for (int s=0; s<n-1; ++s) {
            int pos2 = p + s*wMax + wMin;
            long pos2Kmer = kmerHashVals[pos2];
            long h = combineHashValues(strobemerHash, pos2Kmer);
            
            int end = Math.min(p + s*wMax + wMax, numKmers);
            for (int i=pos2 +1; i<end; ++i) {
                long altKmer = kmerHashVals[i];
                if (altKmer == pos2Kmer) {
                    pos2 = i;
                }
                else {
                    long h2 = combineHashValues(strobemerHash, altKmer);
                    if (Long.compareUnsigned(h, h2) >= 0) {
                        pos2 = i;
                        pos2Kmer = altKmer;
                        h = h2;
                    }
                }
            }
            
            strobemerHash = h;
            positions[s+1] = pos2;
        }
                
        return new HashedPositions(strobemerHash, positions);
    }
    
    @Override
    public int getPos() {
        return pos;
    }
    
    @Override
    public int[] getStrobes() {
        return strobes;
    }
    
    @Override
    public int getMin() {
        return 0;
    }
    
    @Override
    public int getMax() {
        return max;
    }
    
    @Override
    public int getNumStrobemers() {
        return max + 1;
    }
    
    public static void main(String[] args) {
        //debug
        
        StrobeHashIterator itr = new StrobeHashIterator(3, 11, 12, 50);
        String seq = "TCGAATCCGTCTGATGCCTGACTGTAGCTGCGACTGATCGTAGCTAGCGACGAGCAGTCGCCCCATCGTACGTAGTCATGCATGCATGCATGCAGTACTATCTGCACACATGATGCATGCAATCTATATATTTTTATAT";
        if (itr.start(seq)) {
            int[] strobes = itr.getStrobes();
            while(itr.hasNext()) {
                itr.next();
                int pos = itr.getPos();
                System.out.println(pos + "\t" + strobes[0] + "\t" + strobes[1]);
                
//                HashedPositions hp = itr.get(pos);
//                System.out.println(hp.pos[0] + "\t" + hp.pos[1] + "\t" + hp.pos[2]);
            }
        }
    }
}
