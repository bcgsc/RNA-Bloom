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
public class CanonicalStrobeHashIterator implements StrobeHashIteratorInterface {
    private final int wMin;
    private final int wMax;
    private final int k;
    private final CanonicalNTHashIterator itr;
    private int pos = -1;
    private int max = -2;
    private int numKmers = 0;
    private long[] fHashVals = null;
    private long[] rHashVals = null;
    private int n;
    private int[] strobes;
    
    public CanonicalStrobeHashIterator(int n, int k, int wMin, int wMax) {
        this.n = n;
        this.strobes = new int[n-1];
        this.wMin = wMin;
        this.wMax = wMax;
        this.k = k;
        this.itr = new CanonicalNTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        Arrays.fill(strobes, 0);
        max = -2;
        fHashVals = null;
        rHashVals = null;
        if (itr.start(seq)) {
            numKmers = seq.length() - k + 1;
            if (numKmers > wMax * (n-1)) {
                fHashVals = new long[numKmers];
                rHashVals = new long[numKmers];
                max = numKmers - wMax * (n-2) - wMin -1;
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
        long strobemerHash = fHashVals[pos];
        
        for (int s=0; s<n-1; ++s) {
            int pos2 = pos + s*wMax + wMin;
            long h = combineHashValues(strobemerHash, fHashVals[pos2]);
            
            int end = Math.min(pos + s*wMax + wMax, numKmers);
            for (int i=pos2 +1; i<end; ++i) {
                long h2 = combineHashValues(strobemerHash, fHashVals[i]);
                if (Long.compareUnsigned(h, h2) >= 0) {
                    pos2 = i;
                    h = h2;
                }
            }
            
            strobemerHash = h;
            strobes[s] = pos2;
        }
        
        long rStrobemerHash = rHashVals[strobes[n-2]];
        for (int s=n-3; s>=0; --s) {
            rStrobemerHash = combineHashValues(rHashVals[strobes[s]], rStrobemerHash);
        }
        rStrobemerHash = combineHashValues(rHashVals[pos], rStrobemerHash);
        
        return Math.min(strobemerHash, rStrobemerHash);
    }
    
    @Override
    public HashedPositions get(int p) {
        if (p > max) {
            return null;
        }
        
        long strobemerHash = fHashVals[p];
        int[] positions = new int[n];
        positions[0] = p;
        
        for (int s=0; s<n-1; ++s) {
            int pos2 = p + s*wMax + wMin;
            long h = combineHashValues(strobemerHash, fHashVals[pos2]);
            
            int end = Math.min(p + s*wMax + wMax, numKmers);
            for (int i=pos2 +1; i<end; ++i) {
                long h2 = combineHashValues(strobemerHash, fHashVals[i]);
                if (Long.compareUnsigned(h, h2) >= 0) {
                    pos2 = i;
                    h = h2;
                }
            }
            
            strobemerHash = h;
            positions[s+1] = pos2;
        }
        
        long rStrobemerHash = rHashVals[positions[n-1]];
        for (int s=n-2; s>=0; --s) {
            rStrobemerHash = combineHashValues(rHashVals[positions[s]], rStrobemerHash);
        }
        
        return new HashedPositions(Math.min(strobemerHash, rStrobemerHash), positions);
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
        
        CanonicalStrobeHashIterator itr = new CanonicalStrobeHashIterator(3, 15, 20, 50);
        String seq = "TCGAATCCGTCTGATGCCTGACTGTAGCTGCGACTGATCGTAGCTAGCGACGAGCAGTCGCCCCATCGTACGTAGTCATGCATGCATGCATGCAGTACTATCTGCACACATGATGCATGCAATCTATATATTTTTATAT";
        String seqRC = "ATATAAAAATATATAGATTGCATGCATCATGTGTGCAGATAGTACTGCATGCATGCATGCATGACTACGTACGATGGGGCGACTGCTCGTCGCTAGCTACGATCAGTCGCAGCTACAGTCAGGCATCAGACGGATTCGA";
        HashSet<Long> hashVals = new HashSet<>();
        if (itr.start(seq)) {
            int[] strobes = itr.getStrobes();
            while(itr.hasNext()) {
                hashVals.add(itr.next());
                //int pos = itr.getPos();
                //System.out.println(pos + "\t" + strobes[0] + "\t" + strobes[1]);
            }
        }
        
        HashSet<Long> hashValsRC = new HashSet<>();
        if (itr.start(seqRC)) {
            int[] strobes = itr.getStrobes();
            while(itr.hasNext()) {
                hashValsRC.add(itr.next());
                //int pos = itr.getPos();
                //System.out.println(pos + "\t" + strobes[0] + "\t" + strobes[1]);
            }
        }
        
        int numF = hashVals.size();
        int numR = hashValsRC.size();
        boolean changed = hashVals.retainAll(hashValsRC);
        int numI = hashVals.size();
        System.out.println(numF + "\n" + numR + "\n" + numI);
    }
}
