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

import static rnabloom.bloom.hash.HashFunction.combineHashValues;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalStrobeHashIterator implements StrobeHashIteratorInterface {
    private int wMin = -1;
    private int wMax = -1;
    private int k = -1;
    private CanonicalNTHashIterator itr = null;
    private int pos = -1;
    private int pos2 = -1;
    private int max = -2;
    private long[] fHashVals = null;
    private long[] rHashVals = null;
    
    public CanonicalStrobeHashIterator(int k, int wMin, int wMax) {
        this.wMin = wMin;
        this.wMax = wMax;
        this.k = k;
        this.itr = new CanonicalNTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        max = -2;
        fHashVals = null;
        rHashVals = null;
        if (itr.start(seq)) {
            int numKmers = seq.length() - k + 1;
            if (numKmers > wMax) {
                fHashVals = new long[numKmers];
                rHashVals = new long[numKmers];
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
    
    public long next() {
        ++pos;
        long s1f = fHashVals[pos];
        long s1r = rHashVals[pos];
        
        pos2 = pos+wMin;
        long strobe = Math.min(combineHashValues(s1f, fHashVals[pos2]), combineHashValues(rHashVals[pos2], s1r));
        
        for (int i=pos+wMin+1; i<pos+wMax; ++i) {
            long strobe2 = Math.min(combineHashValues(s1f, fHashVals[i]), combineHashValues(rHashVals[i], s1r));
            if (strobe2 < strobe) {
                pos2 = i;
                strobe = strobe2;
            }
        }
        return strobe;
    }
    
    @Override
    public int getPos() {
        return pos;
    }
    
    @Override
    public int getPos2() {
        return pos2;
    }
    
    @Override
    public int getMax() {
        return max;
    }
    
    public static void main(String[] args) {
        //debug
        
        StrobeHashIterator itr = new StrobeHashIterator(15, 20, 70);
        String seq = "";
        itr.start(seq);
        while(itr.hasNext()) {
            itr.next();
            System.out.println(itr.getPos() + "\t" + itr.getPos2());
        }
    }
}
