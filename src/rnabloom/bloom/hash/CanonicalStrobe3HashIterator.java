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
import java.util.HashSet;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalStrobe3HashIterator implements StrobeHashIteratorInterface {

    private final int wMin;
    private final int wMax;
    private final int k;
    private final CanonicalNTHashIterator itr;
    private int pos = -1;
    private int min = 0;
    private int max = -2;
    private long[] fHashVals = null;
    private long[] rHashVals = null;
    private final int[] strobes;
    
    public CanonicalStrobe3HashIterator(int k, int wMin, int wMax) {
        this.strobes = new int[2];
        this.k = k;
        this.wMin = wMin;
        this.wMax = wMax;
        this.itr = new CanonicalNTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        min = 0;
        max = -2;
        Arrays.fill(strobes, 0);
        fHashVals = null;
        rHashVals = null;
        if (itr.start(seq)) {
            int numKmers = seq.length() - k + 1;
            if (numKmers > wMax * 2) {
                fHashVals = new long[numKmers];
                rHashVals = new long[numKmers];
                min = wMax - 1;
                pos = min - 1;
                max = numKmers - wMax;
                for (int i=0; i<numKmers; ++i) {
                    itr.next();
                    fHashVals[i] = itr.frhval[0];
                    rHashVals[i] = itr.frhval[1];
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
        long fKmerHash = fHashVals[pos];
        long rKmerHash = rHashVals[pos];
        
        // fetch downstream strobe
        int pos3 = pos + wMin;
        long h3 = Math.min(combineHashValues(fKmerHash, fHashVals[pos3]), combineHashValues(rHashVals[pos3], rKmerHash));
        int end = pos + wMax;
        for (int i=pos3+1; i<end; ++i) {
            long h = Math.min(combineHashValues(fKmerHash, fHashVals[i]), combineHashValues(rHashVals[i], rKmerHash));
            if (Long.compareUnsigned(h3, h) > 0) {
                pos3 = i;
                h3 = h;
            }
        }
        
        // fetch upstream strobe
        int pos1 = pos - wMax + 1;
        long h1 = Math.min(combineHashValues(fHashVals[pos1], fKmerHash), combineHashValues(rKmerHash, rHashVals[pos1]));
        end = pos - wMin + 1;
        for (int i=pos1+1; i<end; ++i) {
            long h = Math.min(combineHashValues(fHashVals[i], fKmerHash), combineHashValues(rKmerHash, rHashVals[i]));
            if (Long.compareUnsigned(h1, h) > 0) {
                pos1 = i;
                h1 = h;
            }
        }
        
        strobes[0] = pos1;
        strobes[1] = pos3;
        
        long fhash = combineHashValues(fHashVals[pos1], fKmerHash, fHashVals[pos3]);
        long rhash = combineHashValues(rHashVals[pos3], rKmerHash, rHashVals[pos1]);
        
        return Math.min(fhash, rhash);
    }

    @Override
    public HashedPositions get(int p) {
        if (p < min || p > max) {
            return null;
        }
        
        long fKmerHash = fHashVals[p];
        long rKmerHash = rHashVals[p];
        
        // fetch downstream strobe
        int pos3 = p + wMin;
        long h3 = Math.min(combineHashValues(fKmerHash, fHashVals[pos3]), combineHashValues(rHashVals[pos3], rKmerHash));
        int end = p + wMax;
        for (int i=pos3+1; i<end; ++i) {
            long h = Math.min(combineHashValues(fKmerHash, fHashVals[i]), combineHashValues(rHashVals[i], rKmerHash));
            if (Long.compareUnsigned(h3, h) > 0) {
                pos3 = i;
                h3 = h;
            }
        }
        
        // fetch upstream strobe
        int pos1 = p - wMax + 1;
        long h1 = Math.min(combineHashValues(fHashVals[pos1], fKmerHash), combineHashValues(rKmerHash, rHashVals[pos1]));
        end = p - wMin + 1;
        for (int i=pos1+1; i<end; ++i) {
            long h = Math.min(combineHashValues(fHashVals[i], fKmerHash), combineHashValues(rKmerHash, rHashVals[i]));
            if (Long.compareUnsigned(h1, h) > 0) {
                pos1 = i;
                h1 = h;
            }
        }
        
        long fhash = combineHashValues(fHashVals[pos1], fKmerHash, fHashVals[pos3]);
        long rhash = combineHashValues(rHashVals[pos3], rKmerHash, rHashVals[pos1]);
        
        return new HashedPositions(Math.min(fhash, rhash), new int[]{pos1, p, pos3});
    }

    @Override
    public int getPos() {
        return pos;
    }

    @Override
    public int[] getStrobes() {
        return strobes;
    }

    public int getMin() {
        return min;
    }    
    
    @Override
    public int getMax() {
        return max;
    }
    
    public static void main(String[] args) {
        //debug
        
        CanonicalStrobe3HashIterator itr = new CanonicalStrobe3HashIterator(11, 12, 50);
        String seq   = "TCGAATCCGTCTGATGCCTGACTGTAGCTGCGACTGATCGTAGCTAGCGACGAGCAGTCGCCCCATCGTACGTAGTCATGCATGCATGCATGCAGTACTATCTGCACACATGATGCATGCAATCTATATATTTTTATAT";
        String seqRC = "ATATAAAAATATATAGATTGCATGCATCATGTGTGCAGATAGTACTGCATGCATGCATGCATGACTACGTACGATGGGGCGACTGCTCGTCGCTAGCTACGATCAGTCGCAGCTACAGTCAGGCATCAGACGGATTCGA";
        HashSet<Long> hashVals = new HashSet<>();
        if (itr.start(seq)) {
            while(itr.hasNext()) {
                hashVals.add(itr.next());
            }
        }
        
        if (itr.start(seqRC)) {
            while(itr.hasNext()) {
                if (!hashVals.contains(itr.next())) {
                    System.out.println(itr.getPos());
                }
            }
        }
        
        if (itr.start(seq)) {
            int min = itr.getMin();
            int max = itr.getMax();
            for (int i=min; i<=max; ++i) {
                if (!hashVals.contains(itr.get(i).hash)) {
                    System.out.println(i);
                }
            }
        }
    }
}
