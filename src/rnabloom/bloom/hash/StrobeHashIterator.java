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
import rnabloom.olc.HashedInterval;

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
    private int pos2 = -1;
    private int max = -2;
    private long[] hashVals = null;
    
    public StrobeHashIterator(int k, int wMin, int wMax) {
        this.wMin = wMin;
        this.wMax = wMax;
        this.k = k;
        this.itr = new NTHashIterator(k, 1);
    }
    
    @Override
    public boolean start(String seq) {
        pos = -1;
        pos2 = -1;
        max = -2;
        hashVals = null;
        if (itr.start(seq)) {
            int numKmers = seq.length() - k + 1;
            if (numKmers > wMax) {
                hashVals = new long[numKmers];
                max = numKmers - wMax;
                for (int i=0; i<numKmers; ++i) {
                    itr.next();
                    hashVals[i] = itr.hVals[0];
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
        long s1 = hashVals[pos];
        
        pos2 = pos + wMin;
        long strobe = combineHashValues(s1, hashVals[pos2]);
        
        int end = pos + wMax;
        for (int i=pos+wMin+1; i<end; ++i) {
            long strobe2 = combineHashValues(s1, hashVals[i]);
            if (strobe2 < strobe) {
                pos2 = i;
                strobe = strobe2;
            }
        }
        return strobe;
    }
    
    @Override
    public HashedInterval get(int p1) {
        long s1 = hashVals[p1];
        
        int p2 = p1 + wMin;
        long strobe = combineHashValues(s1, hashVals[p2]);
        
        int end = p1 + wMax;
        for (int i=p1+wMin+1; i<end; ++i) {
            long strobe2 = combineHashValues(s1, hashVals[i]);
            if (strobe2 < strobe) {
                p2 = i;
                strobe = strobe2;
            }
        }
        
        return new HashedInterval(p1, p2, strobe);
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
        if (itr.start(seq)) {
            while(itr.hasNext()) {
                itr.next();
                System.out.println(itr.getPos() + "\t" + itr.getPos2());
            }
        }
    }
}
