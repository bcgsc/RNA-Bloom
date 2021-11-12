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
public class Strobe3HashIterator implements StrobeHashIteratorInterface {

    private final int wMin;
    private final int wMax;
    private final int k;
    private final NTHashIterator itr;
    private int pos = -1;
    private int min = 0;
    private int max = -2;
    private long[] fHashVals = null;
    private final int[] strobes;
    
    public Strobe3HashIterator(int k, int wMin, int wMax) {
        this.strobes = new int[2];
        this.k = k;
        this.wMin = wMin;
        this.wMax = wMax;
        this.itr = new NTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        min = 0;
        max = -2;
        Arrays.fill(strobes, 0);
        fHashVals = null;
        if (itr.start(seq)) {
            int numKmers = seq.length() - k + 1;
            if (numKmers > wMax * 2) {
                fHashVals = new long[numKmers];
                min = wMax - 1;
                pos = min - 1;
                max = numKmers - wMax;
                for (int i=0; i<numKmers; ++i) {
                    itr.next();
                    fHashVals[i] = itr.hVals[0];
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
        
        // fetch downstream strobe
        int pos3 = pos + wMin;
        long h3 = combineHashValues(fKmerHash, fHashVals[pos3]);
        int end = pos + wMax;
        for (int i=pos3+1; i<end; ++i) {
            long h = combineHashValues(fKmerHash, fHashVals[i]);
            if (Long.compareUnsigned(h3, h) > 0) {
                pos3 = i;
                h3 = h;
            }
        }
        
        // fetch upstream strobe
        int pos1 = pos - wMax + 1;
        long h1 = combineHashValues(fHashVals[pos1], fKmerHash);
        end = pos - wMin + 1;
        for (int i=pos1+1; i<end; ++i) {
            long h = combineHashValues(fHashVals[i], fKmerHash);
            if (Long.compareUnsigned(h1, h) > 0) {
                pos1 = i;
                h1 = h;
            }
        }
        
        strobes[0] = pos1;
        strobes[1] = pos3;
        
        return combineHashValues(fHashVals[pos1], fKmerHash, fHashVals[pos3]);
    }
    
    @Override
    public HashedPositions get(int p) {
        if (p < min || p > max) {
            return null;
        }
        
        long fKmerHash = fHashVals[p];
        
        // fetch downstream strobe
        int pos3 = p + wMin;
        long h3 = combineHashValues(fKmerHash, fHashVals[pos3]);
        int end = p + wMax;
        for (int i=pos3+1; i<end; ++i) {
            long h = combineHashValues(fKmerHash, fHashVals[i]);
            if (Long.compareUnsigned(h3, h) > 0) {
                pos3 = i;
                h3 = h;
            }
        }
        
        // fetch upstream strobe
        int pos1 = p - wMax + 1;
        long h1 = combineHashValues(fHashVals[pos1], fKmerHash);
        end = p - wMin + 1;
        for (int i=pos1+1; i<end; ++i) {
            long h = combineHashValues(fHashVals[i], fKmerHash);
            if (Long.compareUnsigned(h1, h) > 0) {
                pos1 = i;
                h1 = h;
            }
        }
        
        strobes[0] = pos1;
        strobes[1] = pos3;
        
        return new HashedPositions(combineHashValues(fHashVals[pos1], fKmerHash, fHashVals[pos3]), new int[]{pos1, p, pos3});
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
        
        Strobe3HashIterator itr = new Strobe3HashIterator(11, 12, 50);
        String seq   = "TCGAATCCGTCTGATGCCTGACTGTAGCTGCGACTGATCGTAGCTAGCGACGAGCAGTCGCCCCATCGTACGTAGTCATGCATGCATGCATGCAGTACTATCTGCACACATGATGCATGCAATCTATATATTTTTATAT";
        HashSet<Long> hashVals = new HashSet<>();
        if (itr.start(seq)) {
            while(itr.hasNext()) {
                hashVals.add(itr.next());
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
