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
        
        // forward strand
        long fKmerHash = fHashVals[pos];

        // fetch upstream strobe
        int fPos1 = pos - wMax + 1;
        long fh1 = combineHashValues(fHashVals[fPos1], fKmerHash);
        int fEnd = pos - wMin + 1;
        for (int i=fPos1+1; i<fEnd; ++i) {
            long h = combineHashValues(fHashVals[i], fKmerHash);
            if (Long.compareUnsigned(fh1, h) > 0) {
                fPos1 = i;
                fh1 = h;
            }
        }
        
        // fetch downstream strobe
        int fPos3 = pos + wMin;
        long fh3 = combineHashValues(fh1, fHashVals[fPos3]);
        fEnd = pos + wMax;
        for (int i=fPos3+1; i<fEnd; ++i) {
            long h = combineHashValues(fh1, fHashVals[i]);
            if (Long.compareUnsigned(fh3, h) > 0) {
                fPos3 = i;
                fh3 = h;
            }
        }

        // reverse strand
        long rKmerHash = rHashVals[pos];
        
        // fetch downstream strobe
        int rPos3 = pos + wMin;
        long rh3 = combineHashValues(rHashVals[rPos3], rKmerHash);
        int rEnd = pos + wMax;
        for (int i=rPos3+1; i<rEnd; ++i) {
            long h = combineHashValues(rHashVals[i], rKmerHash);
            if (Long.compareUnsigned(rh3, h) > 0) {
                rPos3 = i;
                rh3 = h;
            }
        }
        
        // fetch upstream strobe
        int rPos1 = pos - wMax + 1;
        long rh1 = combineHashValues(rh3, rHashVals[rPos1]);
        rEnd = pos - wMin + 1;
        for (int i=rPos1+1; i<rEnd; ++i) {
            long h = combineHashValues(rh3, rHashVals[i]);
            if (Long.compareUnsigned(rh1, h) > 0) {
                rPos1 = i;
                rh1 = h;
            }
        }
        
        if (Long.compareUnsigned(fh3, rh1) > 0) {
            strobes[0] = rPos1;
            strobes[1] = rPos3;
            return rh1;
        }

        strobes[0] = fPos1;
        strobes[1] = fPos3;
        return fh3;
    }

    @Override
    public HashedPositions get(int p) {
        if (p < min || p > max) {
            return null;
        }
        
        // forward strand
        long fKmerHash = fHashVals[p];

        // fetch upstream strobe
        int fPos1 = p - wMax + 1;
        long fh1 = combineHashValues(fHashVals[fPos1], fKmerHash);
        int fEnd = p - wMin + 1;
        for (int i=fPos1+1; i<fEnd; ++i) {
            long h = combineHashValues(fHashVals[i], fKmerHash);
            if (Long.compareUnsigned(fh1, h) > 0) {
                fPos1 = i;
                fh1 = h;
            }
        }
        
        // fetch downstream strobe
        int fPos3 = p + wMin;
        long fh3 = combineHashValues(fh1, fHashVals[fPos3]);
        fEnd = p + wMax;
        for (int i=fPos3+1; i<fEnd; ++i) {
            long h = combineHashValues(fh1, fHashVals[i]);
            if (Long.compareUnsigned(fh3, h) > 0) {
                fPos3 = i;
                fh3 = h;
            }
        }

        // reverse strand
        long rKmerHash = rHashVals[p];
        
        // fetch downstream strobe
        int rPos3 = p + wMin;
        long rh3 = combineHashValues(rHashVals[rPos3], rKmerHash);
        int rEnd = p + wMax;
        for (int i=rPos3+1; i<rEnd; ++i) {
            long h = combineHashValues(rHashVals[i], rKmerHash);
            if (Long.compareUnsigned(rh3, h) > 0) {
                rPos3 = i;
                rh3 = h;
            }
        }
        
        // fetch upstream strobe
        int rPos1 = p - wMax + 1;
        long rh1 = combineHashValues(rh3, rHashVals[rPos1]);
        rEnd = p - wMin + 1;
        for (int i=rPos1+1; i<rEnd; ++i) {
            long h = combineHashValues(rh3, rHashVals[i]);
            if (Long.compareUnsigned(rh1, h) > 0) {
                rPos1 = i;
                rh1 = h;
            }
        }
        
        if (Long.compareUnsigned(fh3, rh1) > 0) {
            return new HashedPositions(rh1, new int[]{rPos1, p, rPos3});
        }
        
        return new HashedPositions(fh3, new int[]{fPos1, p, fPos3});
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
                System.out.println(itr.strobes[0] + " " + itr.pos + " " + itr.strobes[1]);
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
