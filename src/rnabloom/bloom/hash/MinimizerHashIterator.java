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

import static rnabloom.bloom.hash.NTHash.NTM64;
import rnabloom.util.LongRollingWindow;

/**
 *
 * @author Ka Ming Nip
 */
public class MinimizerHashIterator {
    protected int pos;
    protected int max;
    protected long prev = 0;
    protected final int k;
    protected final int w;
    private LongRollingWindow window;
    protected final NTHashIterator itr;
    
    public MinimizerHashIterator(int k, int w, NTHashIterator itr) {
        this.k = k;
        this.w = w;
        this.itr = itr;
    }
    
    public boolean start(CharSequence seq) {
        boolean ok = itr.start(seq);
        if (!ok) {
            return false;
        }
        
        pos = -1;
        int numKmers = seq.length() - k + 1;
        max = numKmers - w + 1;
        window = null;
        
        if (max > 0) {
            long[] firstWindow = new long[w];
            
            for (int i=0; i<w-1; ++i) {
                itr.next();
                firstWindow[i] = itr.hVals[0];
            }
            
            window = new LongRollingWindow(firstWindow);
            window.setIndex(w-2, w-2);
            
            // the first `next` call will set the last value in the window
            return true;
        }
        
        return false;
    }
    
    public int getPos() {
        return pos;
    }
    
    public boolean hasNext() {
        return pos < max;
    }
    
    public long next() {
        if (++pos < max) {
            itr.next();
            window.roll(itr.hVals[0]);
        }
        prev = window.getMin();
        return prev;
    }
    
    public int getMinimizerPos() {
        return window.getMinPos();
    }
    
    public void getMultipleHashValues(long bVal, long[] hVals) {
        NTM64(bVal, hVals, k, hVals.length);
    }
    
    public static void main(String[] args) {
        //debug
    }
}
